module SciML

import ..Models: Models, Model, FlatState
import Catalyst
import JumpProcesses
using ModelingToolkit: ModelingToolkit, SymbolicIndexingInterface

import Random
using Logging: LogLevel, @logmsg

Progress = LogLevel(-2)

# Integrator reuse kill-switch. When enabled (the default), a single-owner
# (copy = false) continuation across an Instant/Model{FlatState} intervention
# resets the existing integrator in place instead of rebuilding one
# (clone_problem + init). Only ever engaged for single-owner continuations;
# branches (copy = true) and traces never carry the reuse handle, so they still
# build fresh independent integrators. Set to false to fall back to always
# rebuilding, which is behaviourally identical but allocates more.
const USE_REUSE = Ref(true)

# Share the immutable compiled inner problem (prob.f.sys and friends) via
# remake, and deepcopy ONLY the mutable jump state. `discrete_jump_aggregation`
# and `jump_callback.discrete_callbacks[1].condition` are the SAME object, so
# they are deepcopied in a single call to preserve that aliasing.
function lightclone(jp::JumpProcesses.JumpProblem)
    new_inner = JumpProcesses.remake(jp.prob; u0 = copy(jp.prob.u0))
    agg, cb, maj, rng = deepcopy((
        jp.discrete_jump_aggregation,
        jp.jump_callback,
        jp.massaction_jump,
        jp.rng,
    ))
    T = JumpProcesses.remaker_of(jp)
    T(new_inner, jp.aggregator, agg, cb,
        jp.constant_jumps, jp.variable_jumps, jp.regular_jump,
        maj, rng, jp.kwargs)
end

clone_problem(f!) = lightclone(f!.problem)

normalize_name(s) =
    Symbol(replace(String(ModelingToolkit.getname(s)), '₊' => '.'))

@kwdef mutable struct TriggerProgress
    i::Int = 0
end

function (trigger::TriggerProgress)(_u, _t, _integrator)
    trigger.i += 1
    trigger.i % 100000 == 0
end

# Since SciMLBase.DiscreteCallback is immutable, but we want to adjust progress
# reporting after construction, we wrap the progress reporter (emit directly or
# call back) in a mutable affect! of type ReportProgress.

@kwdef mutable struct EmitProgress
    t0::Float64 = 0.0
end

(progress::EmitProgress)(integrator) = @logmsg(
    Progress,
    :advancing,
    at = "JumpModel",
    done = integrator.t - progress.t0,
)

@kwdef struct CallbackProgress
    callback
end

function (progress::CallbackProgress)(integrator)
    progress.callback(:advancing, done = integrator.t)
    JumpProcesses.u_modified!(integrator, false)
end

@kwdef mutable struct ReportProgress
    reporter::Union{EmitProgress, CallbackProgress} = EmitProgress()
end

function (progress::ReportProgress)(integrator)
    progress.reporter(integrator)
    JumpProcesses.u_modified!(integrator, false)
end

"""
    JumpState

Contains the prepared `JumpProcesses.SSAIntegrator` `integrator` to be advanced
by a `JumpModel`.

In JumpProcesses.jl, the problem and integrator objects are tightly coupled. The
current time, values, RNG instance and (potentially) the recorded trajectories
are all contained in the integrator, but many mutable properties are usually
aliased from the problem, so we would typically consider both to be state.

However, since `JumpProblem` construction is expensive, we instead instruct
the `integrator` to perform a deepcopy of the whole `JumpProblem` before using
any of its resources. Except for the `JumpProblem`'s RNG, which we will re-seed
before `integrator` initialization, we can thus consider the `JumpProblem`s
immutable. This is still faster than rebuilding the whole `JumpProblem`.

To be compatible with a `JumpModel` `f!`, the `integrator` must have been
constructed for the same `f!.problem`. To check whether this is actually the
case, `JumpState` also contains a reference to the corresponding `f!`.
"""
@kwdef struct JumpState
    f!::Model{JumpState}
    integrator::JumpProcesses.SSAIntegrator
end

Models.t(x::JumpState) = x.integrator.t
Models.randomness(x::JumpState) = x.integrator.cb.affect!.rng

function Models.empty_trajectory!(x::JumpState)
    empty!(x.integrator.sol.u)
    empty!(x.integrator.sol.t)
end

FlatState(x::JumpState) = FlatState(
    t = Models.t(x),
    counts = Dict(
        normalize_name(s) => x.integrator[s]
        for s in SymbolicIndexingInterface.variable_symbols(x.integrator)
    ),
    randomness = copy(Models.randomness(x)),
)

"""
    JumpModel <: Model{JumpState}

Represents the stochastic dynamics of a compiled `JumpProcesses.JumpProblem`
`problem`. The problem is built once from a `ModelingToolkit.System` and a
`JumpProcesses.AbstractAggregatorAlgorithm` method.

Gene regulation models in this package ultimately get compiled to `JumpModel`s.

Parameter values can be changed cheaply via [`Models.remake`](@ref), which
patches the compiled `problem` in place instead of rebuilding it from the
symbolic system (avoiding recompilation).

# Specification

In JSON, `JumpModel`s can only be defined indirectly, such as via
[`Models.V1`](@ref).

# Invocation

    (f!::JumpModel)(x::JumpState, Δt::Float64; record = false, _...)

Advance the simulation by applying the stochastic dynamics `f!` to `x` for `Δt`
time units, realizing a segment of the state trajectory.

`x` must be compatible with `f!`, that is, the `JumpProcesses.JumpProblem` (and
corresponding integrator) in `x` must have been produced from `f!.problem`. This
is conservatively checked as `f! === x.f!`. If necessary, users can call
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

JumpModel(
    system::ModelingToolkit.System,
    method::JumpProcesses.AbstractAggregatorAlgorithm;
    # Forwarded to the aggregation constructor; e.g. `bracket_data` for the
    # RSSA family (see `V1.promoter_bracket_data`).
    problem...,
) =
    lock(JUMP_PROBLEM_LOCK) do
        JumpModel(problem = ModelingToolkit.JumpProblem(
            system,
            [s => 0 for s in ModelingToolkit.unknowns(system)],
            (0.0, Inf);
            aggregator = method,
            problem...,
            # Independent, deterministically seeded RNG (instead of the default
            # TaskLocalRNG): `JumpState` reseeds this before initializing any
            # integrator and relies on the deepcopy producing independent
            # instances, so we can treat the `JumpProblem` as effectively
            # immutable. See `JumpState` and SciML issue #554.
            rng = Random.Xoshiro(),
            u0_eltype = Float64,
        ))
    end

system(f!::JumpModel) = f!.problem.prob.f.sys
method(f!::JumpModel) = f!.problem.aggregator

variable_symbols(f!::JumpModel) = SymbolicIndexingInterface.variable_symbols(f!.problem)
parameter_symbols(f!::JumpModel) = SymbolicIndexingInterface.parameter_symbols(f!.problem)

Base.getindex(f!::JumpModel, s) = f!.problem.ps[s]

Models.parameters(f!::JumpModel) = Dict(
    normalize_name(s) => f![s] for s in parameter_symbols(f!)
)

Models.describe(::SciML.JumpModel) = Models.Label("SciML JumpSystem")

function JumpState(x::FlatState; f!::JumpModel)
    problem = clone_problem(f!)
    # ^ Given that ModelingToolkit.init below with alias_jump = false would just
    # deepcopy the whole JumpProblem anyway, we will do it ourselves here so we
    # can rest easy in the knowledge that f!.problem is never mutated.

    Random.setstate!(
        problem.jump_callback.discrete_callbacks[1].condition.rng,
        Random.getstate(Models.randomness(x)),
    )
    # ^ The integrator below will alias this RNG instance, which will then
    # become the returned JumpState's authoritative randomness. We need to
    # control it here already because integrator initialization below will
    # consume entropy! (Part of it will even taint problem, but we will then
    # make no further reference to it.) We know that the two RNGs have the same
    # type because we ensured that at construction. The Models.randomness(x)
    # instance will no longer be carried forward.

    integrator = ModelingToolkit.init(
        problem,
        JumpProcesses.SSAStepper(),
        save_start = false,
        callback = JumpProcesses.DiscreteCallback(
            TriggerProgress(),
            ReportProgress(),
            save_positions = (false, false),
        ),
        alias_jump = true,
        # ^ Might as well allow aliasing since we cloned the JumpProblem anyway.
    )
    integrator.t = Models.t(x)
    for s in SymbolicIndexingInterface.variable_symbols(integrator)
        integrator[s] = get(x.counts, normalize_name(s), 0)
    end
    # Re-set the integrator's authoritative RNG to the FlatState position AFTER
    # init has consumed its entropy, so the authoritative draws below start from
    # exactly the same position as the reuse fast path (reset_integrator!),
    # which never runs init. This makes fresh-build and reuse produce identical
    # trajectories. `integrator.cb.affect!.rng` is what Models.randomness reads.
    Random.setstate!(integrator.cb.affect!.rng, Random.getstate(x.randomness))
    JumpProcesses.reset_aggregated_jumps!(integrator)

    JumpState(; f!, integrator)
end

# Reset an existing integrator to represent `x`, reusing all of its backing
# storage instead of rebuilding a fresh one. Mirrors the tail of
# `JumpState(::FlatState)` exactly (RNG state, time, counts, aggregator reset)
# so a reused integrator is behaviourally identical to a freshly built one.
function reset_integrator!(src::JumpState, x::FlatState)
    integrator = src.integrator
    Random.setstate!(Models.randomness(src), Random.getstate(Models.randomness(x)))
    Models.empty_trajectory!(src)   # match a fresh integrator's empty trajectory
    # Reset the SSAStepper bookkeeping that accumulates across a segment's
    # step!(…, Δt, true) so a reused integrator starts a segment exactly like a
    # freshly `init`-ed one (which has empty tstops, indices at 1, etc.).
    empty!(integrator.tstops)
    integrator.tstops_idx = 1
    integrator.i = 1
    integrator.cur_saveat = 1
    integrator.u_modified = false
    integrator.keep_stepping = true
    integrator.t = Models.t(x)
    for s in SymbolicIndexingInterface.variable_symbols(integrator)
        integrator[s] = get(x.counts, normalize_name(s), 0)
    end
    reset_search_order!(integrator.cb.condition)
    JumpProcesses.reset_aggregated_jumps!(integrator)
    JumpState(; src.f!, integrator)
end

# Some SSA aggregators (notably SortingDirect) keep a self-organizing reaction
# search order that reorders as jumps fire, so a stepped integrator's order is
# path dependent. `reset_aggregated_jumps!` refreshes rates but not this order,
# so restore it to the canonical identity a freshly built aggregator starts
# with, making a reused integrator's next_jump selection identical to a fresh
# one's. Correctness does not depend on the order (any order samples reaction j
# with probability rate_j / sum_rate); this is purely to keep reuse an exact
# match for the from-scratch build. No-op for aggregators without such state.
function reset_search_order!(agg)
    if hasproperty(agg, :jump_search_order)
        jso = agg.jump_search_order
        @inbounds for i in eachindex(jso)
            jso[i] = i
        end
        agg.jump_search_idx = 0
    end
    nothing
end

function Models.adapt!(x::FlatState, f!::JumpModel, ::Val{false})
    src = x.source
    if USE_REUSE[] && src isa JumpState && src.f! === f!
        x.source = nothing   # consume the handle so it can never be reused twice
        return reset_integrator!(src, x)
    end
    JumpState(x; f!)
end

Models.adapt!(x::FlatState, f!::JumpModel, ::Val{true}) = JumpState(x; f!)

Models.adapt!(x::JumpState, f!::JumpModel, ::Val{Copy}) where {Copy} =
    x.f! === f! && !Copy ? x : JumpState(FlatState(x); f!)

# Flattening a live JumpState to feed a Model{FlatState} intervention (division,
# add/set, resampling, ...). When this is a single-owner continuation
# (copy = false), remember the live JumpState on the resulting FlatState so the
# following dynamic segment can reset its integrator in place rather than clone.
function Models.adapt!(x::JumpState, f!::Model{FlatState}, ::Val{false})
    fs = FlatState(x)
    USE_REUSE[] && (fs.source = x)
    fs
end

# HACK: JumpProcesses.remake mutates the original problem.
# see: https://github.com/SciML/JumpProcesses.jl/issues/416
# and https://github.com/SciML/JumpProcesses.jl/issues/554
function remake_p(jp::JumpProcesses.JumpProblem; p)
    new_inner = JumpProcesses.remake(jp.prob; p = p)
    # `update_parameters!` mutates only `scaled_rates`; the stoichiometry
    # vectors and `param_mapper` are read-only, so we share them (like the
    # light clone shares the compiled system) and only give the new jump a
    # fresh mutable `scaled_rates`. A full deepcopy here would drag the whole
    # compiled system via `param_mapper` (~3 GB / ~16% of a run measured).
    maj = jp.massaction_jump
    new_maj = JumpProcesses.MassActionJump(
        copy(maj.scaled_rates),
        maj.reactant_stoch,
        maj.net_stoch,
        maj.param_mapper;
        scale_rates = false,
        useiszero = false,
        nocopy = true,
        rescale_rates_on_update = maj.rescale_rates_on_update,
    )
    JumpProcesses.update_parameters!(new_maj, new_inner.p)
    T = JumpProcesses.remaker_of(jp)
    T(new_inner, jp.aggregator, jp.discrete_jump_aggregation, jp.jump_callback,
        jp.constant_jumps, jp.variable_jumps, jp.regular_jump,
        new_maj, jp.rng, jp.kwargs)
end

# Cheaply patch parameter values into the compiled `problem` (no recompilation
# of the symbolic system), returning a new `JumpModel`. Used e.g. by soft
# knockouts and other parameter interventions.
Models.remake(f!::JumpModel, parameters::AbstractDict{Symbol, <:Real}) = JumpModel(
    problem = remake_p(f!.problem; p = [
        s => get(parameters, normalize_name(s), f![s])
        for s in parameter_symbols(f!)
    ])
)

function Models.each_event(callback::Function, x::JumpState)
    solution = x.integrator.sol

    names = normalize_name.(
        SymbolicIndexingInterface.variable_symbols(solution)
    )
    # ^ We assume that this access is safe and the order agrees with the values
    # in x.integrator.sol.u because this is how SciMLBase constructs the Table
    # reinterpretation in Tables.rows(::AbstractTimeseriesSolution).
    # ^ TODO: Check if this is still the case for newer versions of SciML!

    isempty(solution.u) && return
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

function (f!::JumpModel)(
    x::JumpState,
    Δt::Float64;
    record = false,
    consolidated_progress = nothing,
    verbose = consolidated_progress === nothing,
    _...
)
    f! === x.f! || error("incompatible JumpState, must call adapt!(x, f!)")
    isfinite(Δt) || error("cannot do this forever")

    progress = x.integrator.opts.callback.discrete_callbacks[1].affect!
    if verbose
        progress.reporter.t0 = Models.t(x)
    else
        progress.reporter = CallbackProgress(consolidated_progress)
    end

    verbose && @logmsg Progress :advancing at = "JumpModel" todo = Δt
    x.integrator.save_everystep = record
    if record
        ModelingToolkit.savevalues!(x.integrator, true)
        ModelingToolkit.step!(x.integrator, Δt, true)
    else
        ModelingToolkit.step!(x.integrator, Δt, true)
        ModelingToolkit.savevalues!(x.integrator, true)
    end
    verbose && @logmsg Progress :done at = "JumpModel"

    x
end

end
