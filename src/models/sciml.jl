module SciML

using ...GeneRegulatorySystems: GeneRegulatorySystems
import ..Models: Models, Model, FlatState
import Catalyst
import JumpProcesses
import ModelingToolkit
import OrdinaryDiffEqRosenbrock


using Random
using Logging: LogLevel, @logmsg

Progress = LogLevel(-2)

normalize_name(s) = Symbol(
	replace(String(ModelingToolkit.getname(s)), '₊' => '.')
)

@kwdef mutable struct TriggerProgress
	i::Int = 0
end

function (trigger::TriggerProgress)(_u, _t, _integrator)
	trigger.i += 1
	trigger.i % 10000 == 0
end

@kwdef mutable struct EmitProgress
	t0::Float64 = 0.0
end

(progress::EmitProgress)(integrator) = @logmsg(
	Progress,
	:stepping,
	at = "JumpModel",
	done = integrator.t - progress.t0,
)

"""
	JumpState

Contains the prepared `JumpProcesses.JumpProblem` `problem` and an associated
`JumpProcesses.SSAIntegrator` `integrator` to be advanced by a `JumpModel`.

Since in JumpProcesses.jl the problem and integrator objects are tightly coupled
(because they alias various components), `JumpState` holds both of them. The
current time and state as well as the (potentially) recorded trajectory are
contained in the integrator, while the random number generator is part of the
problem object.

To be compatible with a `JumpModel` `f!`, the `JumpState`'s `problem` (and
`integrator`) must have been constructed by remaking `f!.problem`. To check
whether this is actually the case, `JumpState` also contains a reference to the
corresponding `f!`.
"""
@kwdef struct JumpState
	f!::Model{JumpState}
	problem::JumpProcesses.JumpProblem
	integrator::JumpProcesses.SSAIntegrator = ModelingToolkit.init(
		problem,
		JumpProcesses.SSAStepper(),
		save_start = false,
        seed = GeneRegulatorySystems.seed(problem.rng),
        alias_jump = false,
		callback = JumpProcesses.DiscreteCallback(
			TriggerProgress(),
			EmitProgress(),
			save_positions = (false, false),
		),
	)
end

Models.t(x::JumpState) = x.integrator.t
Models.randomness(x::JumpState) = x.problem.rng

FlatState(x::JumpState) = FlatState(
	counts = Dict(
		normalize_name(s) => x.integrator[s]
		for s in ModelingToolkit.SymbolicIndexingInterface.variable_symbols(
			x.integrator
		)
	),
	randomness = x.problem.rng;
	x.integrator.t,
)

"""
	JumpModel <: Model{JumpState}

Represents the stochastic dynamics of a compiled `JumpProcesses.JumpProblem`
`problem`. The problem is built once from a `ModelingToolkit.System` and a `JumpProcesses.AbstractAggregatorAlgorithm` method.

Gene regulation models in this package ultimately get compiled to `JumpModel`s.

# Specification

In JSON, `JumpModel`s can only be defined indirectly, such as via
[`Models.V1`](@ref).

# Invocation

	(f!::JumpModel)(x::JumpState, Δt::Float64; dense = false, _...)

Advance the simulation by applying the stochastic dynamics `f!` to `x` for `Δt`
time units, realizing a segment of the state trajectory.

`x` must be compatible with `f!`, that is, the `JumpProcesses.JumpProblem` (and
corresponding integrator) in `x` must have been produced by remaking `f!.problem`.
This is conservatively checked as `f! === x.f!`. If necessary, users can call
`adapt!(x, f!)` to convert `x` appropriately.

If `record === true`, `x.integrator` will record all jumps, otherwise the
trajectory will not be retained (and only the final state will be available).
Either way, the recorded trajectory will be initially cleared, so it needs to be
extracted before the next invocation of `f!`.

Unfortunately, JumpProcesses.jl always uses a dense trajectory encoding, so that
the recorded trajectory information is highly redundant and needs to be filtered
by `each_event` for output in sparse long format.
"""
@kwdef struct JumpModel <: Model{JumpState}
	problem::JumpProcesses.JumpProblem
end

# SymbolicUtils hash consing is thread unsafe
# https://github.com/SciML/ModelingToolkit.jl/issues/3315
const JUMP_PROBLEM_LOCK = ReentrantLock()

JumpModel(system::ModelingToolkit.System, method::JumpProcesses.AbstractAggregatorAlgorithm) =
    lock(JUMP_PROBLEM_LOCK) do
        JumpModel(problem = ModelingToolkit.JumpProblem(
            system,
            [s => 0 for s in ModelingToolkit.unknowns(system)],
            (0.0, Inf),
            aggregator = method,
            u0_eltype = Int,
        ))
    end

system(f!::JumpModel) = f!.problem.prob.f.sys
method(f!::JumpModel) = f!.problem.aggregator

variable_symbols(f!::JumpModel) = ModelingToolkit.SymbolicIndexingInterface.variable_symbols(f!.problem)
parameter_symbols(f!::JumpModel) = ModelingToolkit.SymbolicIndexingInterface.parameter_symbols(f!.problem)

Base.getindex(f!::JumpModel, s) = f!.problem.ps[s]

Models.parameters(f!::JumpModel) = Dict(
	normalize_name(s) => f![s] for s in parameter_symbols(f!)
)

Models.describe(::SciML.JumpModel) = Models.Label("SciML JumpSystem")

Models.adapt!(x::JumpState, f!::JumpModel, ::Val{Copy}) where {Copy} =
	if x.f! === f! && !Copy
		x
	else
		# Since SciML problems and integrators are tightly coupled we need to
		# remake the problem and then reinitialize the integrator if we want
		# a Model copy. Remaking JumpProblem only allows changing a limited
		# subset of the properties, and I am unsure which ones are aliased in
		# the process. To avoid trouble, we choose to simply extract the
		# current state to a FlatState and then proceed as if this were a new
		# model. Presumably this is slower than calling remake, yet safer, and
		# anyway could only be avoided when we are branching the simulation
		# without changing models.
		Models.adapt!(FlatState(x), f!)
	end

Models.adapt!(x::FlatState, f!::JumpModel, _copy) = JumpState(
    # because we are only modifying u0 here and not p, f!.problem is not mutated.
    problem = with_rng(JumpProcesses.remake(f!.problem;
        u0 = [
            s => get(x.counts, normalize_name(s), 0)
 			for s in variable_symbols(f!)
        ],
        tspan = (x.t, Inf));
        rng = x.randomness
    );
	f!,
)

# HACK: JumpProcesses.remake currently doesn't support rng as kwarg
function with_rng(jp::JumpProcesses.JumpProblem; rng)
    T = JumpProcesses.remaker_of(jp)
    T(jp.prob, jp.aggregator, jp.discrete_jump_aggregation, jp.jump_callback,
      jp.constant_jumps, jp.variable_jumps, jp.regular_jump,
      jp.massaction_jump, rng, jp.kwargs)
end

# HACK: JumpProcesses.remake mutate the original problem.
# see: https://github.com/SciML/JumpProcesses.jl/issues/416
# and https://github.com/SciML/JumpProcesses.jl/issues/554
function remake_p(jp::JumpProcesses.JumpProblem; p)
    new_inner = JumpProcesses.remake(jp.prob; p = p)
    new_maj = deepcopy(jp.massaction_jump)
    JumpProcesses.update_parameters!(new_maj, new_inner.p)
    T = JumpProcesses.remaker_of(jp)
    T(new_inner, jp.aggregator, jp.discrete_jump_aggregation, jp.jump_callback,
        jp.constant_jumps, jp.variable_jumps, jp.regular_jump,
        new_maj, jp.rng, jp.kwargs)
end

Models.remake(f!::JumpModel, parameters::AbstractDict{Symbol, <:Real}) = JumpModel(
    problem = remake_p(f!.problem; p = [
        s => get(parameters, normalize_name(s), f![s])
        for s in parameter_symbols(f!)
    ])
)

function Models.each_event(callback::Function, x::JumpState)
	solution = x.integrator.sol

	names = normalize_name.(
		ModelingToolkit.SymbolicIndexingInterface.variable_symbols(solution)
	)
	# ^ We assume that this access is safe and the order agrees with the values
	# in x.integrator.sol.u because this is how SciMLBase constructs the Table
	# reinterpretation in Tables.rows(::AbstractTimeseriesSolution).

	isempty(solution.u) && return nothing
	(t, previous), rest = Iterators.peel(zip(solution.t, solution.u))

	# We generate events for all variables at the beginning of the segment...
	for i in LinearIndices(previous)
		callback(t, names[i], previous[i])
	end

	# ...and only for changes at later timepoints.
	for (t, current) in rest
		for i in LinearIndices(current)
			if current[i] != previous[i]
				callback(t, names[i], current[i])
			end
		end
		previous = current
	end
end

function (f!::JumpModel)(x::JumpState, Δt::Float64; dense = false, _...)
	f! === x.f! || error("incompatible JumpState, must call adapt!(x, f!)")
	isfinite(Δt) || error("cannot do this forever")

	empty!(x.integrator.sol.u)
	empty!(x.integrator.sol.t)
	x.integrator.opts.callback.discrete_callbacks[1].affect!.t0 = Models.t(x)

	@logmsg Progress :stepping at = "JumpModel" todo = Δt
	if dense
		# Save every jump event into sol.u (for trajectory consumers like ExperimentTool).
		# Under active SSA dynamics this can allocate gigabytes per segment.
		x.integrator.save_everystep = true
		ModelingToolkit.savevalues!(x.integrator, true)
		ModelingToolkit.step!(x.integrator, Δt, true)
	else
		# Only retain the segment-end state. Sufficient when consumers read just the
		# current state (e.g. via FlatState) rather than iterating sol.u.
		x.integrator.save_everystep = false
		ModelingToolkit.step!(x.integrator, Δt, true)
		ModelingToolkit.savevalues!(x.integrator, true)
	end
	@logmsg Progress :done at = "JumpModel"

	x
end

"""
	ODEState

Holds the `ModelingToolkit.ODEProblem` `problem` and its `ODEIntegrator`
`integrator`, advanced by an [`ODEModel`](@ref). Mirrors [`JumpState`](@ref).

The model is deterministic, but the `Model`/`State` interface requires a random
number generator to be present and aliased through `adapt!`, so `ODEState`
carries the inbound `randomness` unchanged and never consumes it.
"""
@kwdef struct ODEState
	f!::Model{ODEState}
	problem::ModelingToolkit.ODEProblem
	randomness::AbstractRNG = Random.Xoshiro()
	integrator = ModelingToolkit.init(
		problem,
		f!.solver;
		save_everystep = false,
		save_start = false,
		# Species are concentrations: reject any adaptive trial step that probes
		# negative values. Without this the stiff solver evaluates Hill rate laws
		# (and their Jacobian) at negative concentrations — fractional powers there
		# return NaN and abort the solve at t=0.
		isoutofdomain = (u, _p, _t) -> any(<(0), u),
	)
end

"""
	ODEModel <: Model{ODEState}

Represents the deterministic **mean-field** dynamics of a regulation network
compiled to a `ModelingToolkit.ODESystem` (via `convert(ODESystem, …)` on the
`Catalyst.ReactionSystem`). It is the continuous sibling of [`JumpModel`](@ref):
the same `Model`/`State` interface and schedule integration, advancing an
`ODEIntegrator` by `Δt` instead of stepping an SSA.

An ODE solve yields the mean trajectory only — there is no cell-to-cell variance.
A recorded `ODEState` is therefore a population of one (mean = the value,
variance = 0). This is intended: `ODEModel` is a fast, differentiable
approximation of the network's mean behaviour, not a stochastic simulator.

Selected through [`build`](@ref)'s `method` flag as `method = :ode`.
"""
@kwdef struct ODEModel <: Model{ODEState}
	problem::ModelingToolkit.ODEProblem
	solver = OrdinaryDiffEqRosenbrock.Rosenbrock23(autodiff = false)
end

# `system` is the `Catalyst.ReactionSystem` (or any system `ODEProblem` accepts);
# `ODEProblem` builds the mean-field reaction-rate-equation ODE from it directly.
# `build_initializeprob = false`: the reaction-rate ODE has no algebraic constraints, so MTK's
# initialization system is spurious — and rebuilding it on each u0/p `remake` both wastes work
# and trips an MTK codegen bug. We always supply explicit initial conditions, so skip it.
ODEModel(system; solver = OrdinaryDiffEqRosenbrock.Rosenbrock23(autodiff = false), tspan = (0.0, Inf)) =
	ODEModel(
		problem = ModelingToolkit.ODEProblem(
			system,
			[s => 0.0 for s in ModelingToolkit.unknowns(system)],
			tspan;
			build_initializeprob = false,
		),
		solver = solver,
	)

Models.t(x::ODEState) = x.integrator.t
Models.randomness(x::ODEState) = x.randomness

FlatState(x::ODEState) = FlatState(
	counts = Dict{Symbol, Real}(
		normalize_name(s) => x.integrator[s]
		for s in ModelingToolkit.SymbolicIndexingInterface.variable_symbols(
			x.integrator
		)
	),
	randomness = x.randomness;
	x.integrator.t,
)

system(f!::ODEModel) = f!.problem.f.sys
method(f!::ODEModel) = f!.solver

variable_symbols(f!::ODEModel) = ModelingToolkit.SymbolicIndexingInterface.variable_symbols(f!.problem)
parameter_symbols(f!::ODEModel) = ModelingToolkit.SymbolicIndexingInterface.parameter_symbols(f!.problem)

Base.getindex(f!::ODEModel, s) = f!.problem.ps[s]

# An ODE problem's parameter set includes MTK's auto-generated `Initial(x)` initialization
# parameters alongside the genuine scalar rate constants. The former are composite call
# expressions (`Initial` applied to a variable); the latter are leaf symbols. Only the leaves
# are tunable model parameters, so restrict to those (a JumpProblem has no `Initial` params).
tunable_parameters(f!::ODEModel) = Iterators.filter(
	s -> !ModelingToolkit.iscall(s), parameter_symbols(f!)
)

Models.parameters(f!::ODEModel) = Dict(
	normalize_name(s) => f![s] for s in tunable_parameters(f!)
)

Models.describe(::ODEModel) = Models.Label("SciML ODESystem")

Models.adapt!(x::ODEState, f!::ODEModel, ::Val{Copy}) where {Copy} =
	if x.f! === f! && !Copy
		x
	else
		Models.adapt!(FlatState(x), f!)
	end

Models.adapt!(x::FlatState, f!::ODEModel, _copy) = ODEState(
	# ODEProblem.remake has none of JumpProblem's aliasing landmines, so we can
	# remake u0/tspan directly without the with_rng/remake_p workarounds.
	problem = ModelingToolkit.remake(f!.problem;
		u0 = [
			s => Float64(get(x.counts, normalize_name(s), 0))
			for s in variable_symbols(f!)
		],
		tspan = (x.t, Inf),
		build_initializeprob = false,
	),
	randomness = x.randomness;
	f!,
)

Models.remake(f!::ODEModel, parameters::AbstractDict{Symbol, <:Real}) = ODEModel(
	problem = ModelingToolkit.remake(f!.problem; p = [
		s => get(parameters, normalize_name(s), f![s])
		for s in tunable_parameters(f!)
	], build_initializeprob = false),
	solver = f!.solver,
)

# Mean-field has no per-event history; emit the current (segment-end) state.
function Models.each_event(callback::Function, x::ODEState)
	t = x.integrator.t
	for s in ModelingToolkit.SymbolicIndexingInterface.variable_symbols(x.integrator)
		callback(t, normalize_name(s), x.integrator[s])
	end
end

function (f!::ODEModel)(x::ODEState, Δt::Float64; _...)
	f! === x.f! || error("incompatible ODEState, must call adapt!(x, f!)")
	isfinite(Δt) || error("cannot do this forever")

	@logmsg Progress :stepping at = "ODEModel" todo = Δt
	ModelingToolkit.step!(x.integrator, Δt, true)
	@logmsg Progress :done at = "ODEModel"

	x
end

end
