module SciML

import ..Models: Models, Model, FlatState
import Catalyst
import JumpProcesses
using ModelingToolkit: ModelingToolkit, SymbolicIndexingInterface
const MTKB = ModelingToolkit.ModelingToolkitBase
import JumpProcesses.DiffEqBase
const ArrayPartition = MTKB.ArrayPartition
import JumpProcesses.SciMLBase

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

# Replace only the last '.' in the name with the scope separator '₊' because
# gene names may contain dots while the kind names certainly do not.
symbolic_name(name::Symbol) =
    Symbol(replace(String(name), r"\.(?=[^.]*$)" => '₊'))

function Models.Reaction(reaction::Catalyst.Reaction)
    function reagents(species, stoichiometries)
        Models.Reagents(Dict(
            normalize_name(species) => Int(stoichiometry)
            for (species, stoichiometry) in zip(species, stoichiometries)
        ))
    end
    function reaction_name(properties)
        haskey(properties, :name) &&
            return Symbol("reaction.", properties[:name])

        kind = properties[:kind]
        owner = properties[:owner]
        kind === :proteolysis && return Symbol(owner, ".", kind, ".", properties[:from])

        Symbol(owner, ".", kind)
    end
    properties = Dict{Symbol, Any}(reaction.metadata)
    properties[:only_use_rate] = reaction.only_use_rate
    Models.Reaction(
        name = reaction_name(properties),
        from = reagents(reaction.substrates, reaction.substoich),
        to = reagents(reaction.products, reaction.prodstoich),
        k⁺ = reaction.rate,
        k⁻ = zero(reaction.rate),
        properties = properties
    )
end

function Models.describe(system::Catalyst.ReactionSystem)
    function merge_reactions(reactions)
        rxs = Models.Reaction[]
        groups = Dict{Symbol, Dict{Symbol, Models.Reaction}}()
        for r in reactions
            dir = get(r.properties, :direction, nothing)
            if isnothing(r.name) || isnothing(dir)
                push!(rxs, r)
            else
                directions = get!(groups, r.name) do
                    Dict{Symbol, Models.Reaction}()
                end
                directions[dir] = r
            end
        end
        for (name, directions) in groups
            forward = get(directions, :forward, nothing)
            reverse = get(directions, :reverse, nothing)
            base = something(forward, reverse)
            properties = copy(base.properties)
            delete!(properties, :direction)
            properties[:parameters] = merge(
                isnothing(forward) ? Dict{Symbol,Symbol}() :
                    forward.properties[:parameters],
                isnothing(reverse) ? Dict{Symbol,Symbol}() :
                    reverse.properties[:parameters],
            )
            push!(rxs, Models.Reaction(
                name=name,
                from=isnothing(forward) ? reverse.to : forward.from,
                to = isnothing(forward) ? reverse.from : forward.to,
                k⁺ = isnothing(forward) ? zero(reverse.k⁺) : forward.k⁺,
                k⁻ = isnothing(reverse) ? zero(forward.k⁺) : reverse.k⁺,
                properties=properties
            ))
        end
        rxs
    end
    Models.ReactionNetwork(reactions = merge_reactions(Models.Reaction.(Catalyst.reactions(system))))
end


@kwdef mutable struct TriggerProgress
    i::Int = 0
    last::Float64 = 0.0
end

function (trigger::TriggerProgress)(_u, _t, _integrator)
    trigger.i += 1
    trigger.i % 1000 == 0 || return false
    now = time()
    now - trigger.last < 0.1 && return false
    trigger.last = now
    true
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
    path = nothing
end

function (progress::CallbackProgress)(integrator)
    progress.callback(:advancing, done = integrator.t, path = progress.path)
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

function discrete_problem(system, op, tspan; kwargs...)
    inner, u0, p = MTKB.process_SciMLProblem(
        MTKB.EmptySciMLFunction{true}, system, op;
        t = tspan[1], check_length = false, build_initializeprob = false, kwargs...
    )
    SciMLBase.DiscreteProblem(
        SciMLBase.DiscreteFunction{true, true}(
            DiffEqBase.DISCRETE_INPLACE_DEFAULT; # no dynamics between jumps
            sys = system,
            observed = MTKB.ObservedFunctionCache(system),
            initialization_data = get(inner.kwargs, :initialization_data, nothing)
        ),
        u0, tspan, p; kwargs...
    )
end

function jump_problem(system, aggregator, source;
    op = [s => 0 for s in ModelingToolkit.unknowns(system)],
    tspan = (0.0, Inf),
    # `JumpProcesses` defaults to `(false, true)` for a `DiscreteProblem`;
    # `ModelingToolkit.JumpProblem` passed `(true, true)`, which `each_event`
    # and the trajectory sinks rely on.
    save_positions = (true, true),
    options...)
    problem = discrete_problem(system, op, tspan; u0_eltype=Float64)
    ids = Dict(ModelingToolkit.value(u) => i
        for (i, u) in enumerate(ModelingToolkit.unknowns(system)))
    jumps = MTKB.jumps(system)
    majs = Vector{JumpProcesses.MassActionJump}(
        filter(j -> j isa JumpProcesses.MassActionJump, jumps))
    crjs = Vector{JumpProcesses.ConstantRateJump}(
        filter(j -> j isa JumpProcesses.ConstantRateJump, jumps))
    graphs = dependency_graphs(system, aggregator, ArrayPartition(majs, crjs), ids)
    stoichiometry = net_stoichiometry(aggregator, ids, crjs)

    JumpProcesses.JumpProblem(
        problem, aggregator,
        JumpProcesses.JumpSet(
            massaction_jumps = isempty(majs) ? nothing : MTKB.assemble_maj(
                majs, ids,
                MTKB.JumpSysMajParamMapper(system, problem.p;
                    jseqs = jumps, rateconsttype = Float64)
            ),
            constant_jumps = constant_rate_jumps(source, problem, system, ids, crjs)
        );
        graphs..., stoichiometry...,
        callback = MTKB.process_events(system;
            op = MTKB.operating_point_preprocess(system, op),
            reset_jumps = true, tspan
        ),
        scale_rates = false, nocopy = true, save_positions,
        # Independent, deterministically seeded RNG (instead of the default
        # TaskLocalRNG): `JumpState` reseeds this before initializing any
        # integrator and relies on the deepcopy producing independent
        # instances, so we can treat the `JumpProblem` as effectively
        # immutable. See `JumpState` and SciML issue #554.
        rng = Random.Xoshiro(),
        options...
    )
end

function dependency_graphs(system, aggregator, jumps, ids)
    JumpProcesses.needs_vartojumps_map(aggregator) ||
        JumpProcesses.needs_depgraph(aggregator) ||
        aggregator isa JumpProcesses.NullAggregator ||
        return (;)

    forward = MTKB.asgraph([reads(jump, ids) for jump in jumps], ids)
    backward = variable_dependencies([writes(jump, ids) for jump in jumps], ids)
    (;
        vartojumps_map = forward.badjlist,
        jumptovars_map = backward.badjlist,
        dep_graph = JumpProcesses.needs_depgraph(aggregator) ?
            MTKB.eqeq_dependencies(forward, backward).fadjlist : nothing,
    )
end

function reads(jump, ids)
    buffer = Set()
    ModelingToolkit.Symbolics.get_variables!(buffer, jump.rate)
    [variable for variable in buffer if haskey(ids, ModelingToolkit.value(variable))]
end

function reads(jump::JumpProcesses.MassActionJump, ids)
    buffer = Set()
    rates = ModelingToolkit.value(jump.scaled_rates)
    rates isa Number || ModelingToolkit.Symbolics.get_variables!(buffer, rates)
    for (species, _) in jump.reactant_stoch
        push!(buffer, species)
    end
    [variable for variable in buffer if haskey(ids, ModelingToolkit.value(variable))]
end

writes(jump::JumpProcesses.MassActionJump, ids) =
    [species for (species, _) in jump.net_stoch
     if haskey(ids, ModelingToolkit.value(species))]
writes(jump, ids) =
    [affect.lhs for affect in jump.affect! # ::AbstractVector{Symbolics.Equation}
    # here MTK makes assumptions about the jump.affect! type that are not really documented anywhere
     if haskey(ids, ModelingToolkit.value(affect.lhs))]

function variable_dependencies(modified, ids)
    badjlist = [unique!(sort!([ids[ModelingToolkit.value(v)] for v in vars])) for vars in modified]
    fadjlist = [Vector{Int}() for _ in 1:length(ids)]
    edges = 0
    for (jump, variables) in enumerate(badjlist)
        foreach(variable -> push!(fadjlist[variable], jump), variables)
        edges += length(variables)
    end
    MTKB.BipartiteGraph(edges, fadjlist, badjlist)
end

net_stoichiometry(::JumpProcesses.AbstractAggregatorAlgorithm, ids, jumps) = (;)

net_stoichiometry(::JumpProcesses.HybridTau, ids, jumps) =
    (; crj_stoich = crj_stoichiometry(ids, jumps))

net_stoichiometry(::JumpProcesses.TauSplitting, ids, jumps) = (;
    jumptostoich_map = crj_stoichiometry(ids, jumps))

crj_stoichiometry(ids, jumps) = [
        Pair{Int, Int}[
            ids[ModelingToolkit.value(affect.lhs)] => stoichiometry(affect)
            for affect in jump.affect!
        ]
        for jump in jumps
    ]

function stoichiometry(affect::ModelingToolkit.Equation)
    change = ModelingToolkit.value(ModelingToolkit.Symbolics.expand(
        ModelingToolkit.value(affect.rhs) - MTKB.Pre(ModelingToolkit.value(affect.lhs))
    ))
    change isa Number || error(
        "affect `$(affect)` does not change its species by a constant amount"
    )
    Int(change)
end

const ParameterIndex = MTKB.ParameterIndex{MTKB.SciMLStructures.Tunable, Int}

struct Affect{N}
    unknowns::NTuple{N, Int}
    change::NTuple{N, Int8}
end

function (a::Affect{N})(integrator) where N
    @inbounds for i in 1:N
        integrator.u[a.unknowns[i]] += a.change[i]
    end
    nothing
end

struct Indices
    unknowns::Dict{Symbol, Int}
    parameters::Dict{Symbol, ParameterIndex}
end

Indices(system, ids::AbstractDict) = Indices(
    Dict{Symbol, Int}(
        normalize_name(unknown) => index for (unknown, index) in ids),
    Dict{Symbol, ParameterIndex}(
        normalize_name(parameter) =>
            SymbolicIndexingInterface.parameter_index(
                system, ModelingToolkit.getname(parameter))
        for parameter in ModelingToolkit.parameters(system)),
)

function lrate end
function urate end

bounded_jump(rate, affect) = JumpProcesses.ConstantRateJump(rate, affect;
    bounds = (ulow, uhigh, u, p, _) -> JumpProcesses.RateBounds(
        lrate = lrate(rate, ulow, uhigh, p), urate = urate(rate, ulow, uhigh, p)),
    lrate = (ulow, uhigh, u, p, _) -> JumpProcesses.RateBounds(
        lrate = lrate(rate, ulow, uhigh, p)),
    urate = (ulow, uhigh, u, p, _) -> JumpProcesses.RateBounds(
        urate = urate(rate, ulow, uhigh, p)))

constant_rate_jumps(::Nothing, problem, system, ids, jumps) =
    [MTKB.assemble_crj(system, j, ids) for j in jumps]

function constant_rate_jumps(directions::AbstractVector, problem, system, ids, jumps)
    unbounded = constant_rate_jumps(nothing, problem, system, ids, jumps)
    length(directions) == length(unbounded) || error(
        "got $(length(directions)) direction vectors for $(length(unbounded)) ConstantRateJumps"
    )
    [
        JumpProcesses.ConstantRateJump(c.rate, c.affect!;
            bounds = bracket_bound(c.rate, dirs, problem.u0),
            lrate = lower_bound(c.rate, dirs, problem.u0),
            urate = upper_bound(c.rate, dirs, problem.u0)
        ) for (c, dirs) in zip(unbounded, directions)
    ]
end

observed_activity(rate, u, p) = nothing

# SymbolicUtils hash consing is thread unsafe
# https://github.com/SciML/ModelingToolkit.jl/issues/3315
const JUMP_PROBLEM_LOCK = ReentrantLock()

JumpModel(
    system::ModelingToolkit.System,
    method::JumpProcesses.AbstractAggregatorAlgorithm,
    # Dispatches `constant_rate_jumps`; `nothing` builds unbounded jumps from
    # the symbolic system, a vector of propensity directions adds rate bounds.
    source = nothing;
    # Forwarded to the aggregation constructor; e.g. `bracket_data` for the
    # RSSA family (see `V1.promoter_bracket_data`).
    problem...,
) = lock(JUMP_PROBLEM_LOCK) do
    JumpModel(problem = jump_problem(system, method, source; problem...))
end

function corner_states!(ulow_corner, uhigh_corner, ulow, uhigh, directions)
    @inbounds for (species, direction) in directions
        if direction > 0
            ulow_corner[species] = ulow[species]
            uhigh_corner[species] = uhigh[species]
        else
            ulow_corner[species] = uhigh[species]
            uhigh_corner[species] = ulow[species]
        end
    end
    nothing
end

function bracket_bound(rate, directions, u)
    ulow_corner = copy(u)
    uhigh_corner = copy(u)
    function (ulow, uhigh, u, p, t)
        corner_states!(ulow_corner, uhigh_corner, ulow, uhigh, directions)
        JumpProcesses.RateBounds(
            lrate = rate(ulow_corner, p, t),
            urate = rate(uhigh_corner, p, t)
        )
    end
end

function lower_bound(rate, directions, u)
    ulow_corner = copy(u)
    uhigh_corner = copy(u)
    function (ulow, uhigh, u, p, t)
        corner_states!(ulow_corner, uhigh_corner, ulow, uhigh, directions)
        JumpProcesses.RateBounds(
            lrate = rate(ulow_corner, p, t)
        )
    end
end

function upper_bound(rate, directions, u)
    ulow_corner = copy(u)
    uhigh_corner = copy(u)
    function (ulow, uhigh, u, p, t)
        corner_states!(ulow_corner, uhigh_corner, ulow, uhigh, directions)
        JumpProcesses.RateBounds(
            urate = rate(uhigh_corner, p, t)
        )
    end
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

    integrator = JumpProcesses.init(
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

    p = x.integrator.p
    model = Models.unwrap(x.f!)
    covered = Set{Symbol}()
    for jump in model.problem.constant_jumps
        rate = jump.rate
        first_observed = observed_activity(rate, first(solution.u), p)
        isnothing(first_observed) && continue
        name, observed = first_observed
        name in covered && continue
        push!(covered, name)
        callback(first(solution.t), name, observed)

        for (t, u) in Iterators.drop(zip(solution.t, solution.u), 1)
            current = last(observed_activity(rate, u, p))
            current != observed && callback(t, name, current)
            observed = current
        end
    end

    sys = system(model)
    uncovered = filter(
        observable -> normalize_name(observable) ∉ covered,
        ModelingToolkit.observables(sys),
    )
    if !isempty(uncovered)
        evaluate = ModelingToolkit.build_explicit_observed_function(sys, uncovered)
        uncovered_names = normalize_name.(uncovered)

        (t, u), rest = Iterators.peel(zip(solution.t, solution.u))
        prev = evaluate(u, p, t)
        for i in eachindex(prev)
            callback(t, uncovered_names[i], prev[i])
        end

        for (t, u) in rest
            curr = evaluate(u, p, t)
            for i in eachindex(curr)
                curr[i] != prev[i] && callback(t, uncovered_names[i], curr[i])
            end
            prev = curr
        end
    end
end

function (f!::JumpModel)(
    x::JumpState,
    Δt::Float64;
    record = false,
    consolidated_progress = nothing,
    verbose = consolidated_progress === nothing,
    path = nothing,
    _...
)
    f! === x.f! || error("incompatible JumpState, must call adapt!(x, f!)")
    isfinite(Δt) || error("cannot do this forever")

    progress = x.integrator.opts.callback.discrete_callbacks[1].affect!
    if verbose
        progress.reporter.t0 = Models.t(x)
    else
        progress.reporter = CallbackProgress(consolidated_progress, path)
    end

    verbose && @logmsg Progress :advancing at = "JumpModel" todo = Δt
    x.integrator.save_everystep = record
    if record
        JumpProcesses.savevalues!(x.integrator, true)
        JumpProcesses.step!(x.integrator, Δt, true)
    else
        JumpProcesses.step!(x.integrator, Δt, true)
        JumpProcesses.savevalues!(x.integrator, true)
    end
    verbose && @logmsg Progress :done at = "JumpModel"

    x
end

end
