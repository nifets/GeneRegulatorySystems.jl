module SciML

using ...GeneRegulatorySystems: GeneRegulatorySystems
import ..Models: Models, Model, FlatState
import Catalyst
import JumpProcesses
import ModelingToolkit


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
`integrator`) must have been constructed by remaking `f!.template`. To check
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
`template`. The template is built once from a `ModelingToolkit.System` and a `JumpProcesses.AbstractAggregatorAlgorithm` method.

Gene regulation models in this package ultimately get compiled to `JumpModel`s.

# Specification

In JSON, `JumpModel`s can only be defined indirectly, such as via
[`Models.V1`](@ref).

# Invocation

	(f!::JumpModel)(x::JumpState, Δt::Float64; dense = false, _...)

Advance the simulation by applying the stochastic dynamics `f!` to `x` for `Δt`
time units, realizing a segment of the state trajectory.

`x` must be compatible with `f!`, that is, the `JumpProcesses.JumpProblem` (and
corresponding integrator) in `x` must have been produced by remaking `f!.template`.
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
	template::JumpProcesses.JumpProblem
end

# SymbolicUtils hash consing is thread unsafe
# https://github.com/SciML/ModelingToolkit.jl/issues/3315
const JUMP_PROBLEM_LOCK = ReentrantLock()

JumpModel(system::ModelingToolkit.System, method::JumpProcesses.AbstractAggregatorAlgorithm) =
    lock(JUMP_PROBLEM_LOCK) do
        JumpModel(template = ModelingToolkit.JumpProblem(
            system,
            [s => 0 for s in ModelingToolkit.unknowns(system)],
            (0.0, Inf),
            aggregator = method,
            u0_eltype = Int,
        ))
    end

system(f!::JumpModel) = f!.template.prob.f.sys
method(f!::JumpModel) = f!.template.aggregator

variable_symbols(f!::JumpModel) = ModelingToolkit.SymbolicIndexingInterface.variable_symbols(f!.template)
parameter_symbols(f!::JumpModel) = ModelingToolkit.SymbolicIndexingInterface.parameter_symbols(f!.template)

Base.getindex(f!::JumpModel, s) = f!.template.ps[s]

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
    # because we are only modifying u0 here and not p, f!.template is not mutated.
    problem = with_rng(JumpProcesses.remake(f!.template;
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
    template = remake_p(f!.template; p = [
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

end
