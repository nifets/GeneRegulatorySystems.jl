module PEtabExt

using GeneRegulatorySystems: Models
using GeneRegulatorySystems.Models: Scheduling, Plumbing, SciML
using Catalyst: ode_model, get_iv, species
using ModelingToolkit: ModelingToolkit, getname, substitute, complete
using PEtab: PEtab, PEtabCondition, PEtabEvent

intervention(primitive!::Scheduling.Primitive) = intervention(primitive!.f!)
intervention(f!::Models.Wrapped) = intervention(f!.model)
intervention(f!::Plumbing.Adjust) = f!
intervention(::Plumbing.Pass) = nothing
intervention(f!::Models.Instant) = error("cannot convert $(nameof(typeof(f!))) to PEtab")
intervention(_) = nothing

petab_name(name) = Symbol(replace(String(name), '.' => ModelingToolkit.NAMESPACE_SEPARATOR))

function assignments(segments, instants, symbolic)
    states = Dict{Float64, Dict{Symbol, Any}}()
    for segment in segments
        f! = instants[segment.model]
        isnothing(f!) && continue
        close = filter(t -> isapprox(t, segment.from; atol = 1e-9), collect(keys(states)))
        time = isempty(close) ? segment.from : only(close)
        state = get!(Dict{Symbol, Any}, states, time)
        for (name, value) in f!.adjustment
            haskey(symbolic, name) || continue
            initial = iszero(time) ? 0.0 : symbolic[name]
            state[name] = f!.adjust(get(state, name, initial), value)
        end
    end
    Dict(
        time => [symbolic[name] => value for (name, value) in state]
        for (time, state) in states if !isempty(state)
    )
end

function only_reaction_system(schedule!::Scheduling.Schedule, models = Scheduling.models(schedule!))
    dynamic_models = filter(((_, f!),) -> !(Models.unwrap(f!) isa Models.Instant), collect(models))
    length(dynamic_models) == 1 ||
        error("PEtab needs exactly one time-advancing model, schedule has $(length(dynamic_models))")
    path, f! = only(dynamic_models)
    system = SciML.reaction_system(f!)
    isnothing(system) && error("$(nameof(typeof(Models.unwrap(f!)))) is not a reaction system")
    invalid = filter(!Base.isidentifier, string.(getname.(species(system))))
    isempty(invalid) || error("species names not supported by PEtab: $(join(invalid, ", "))")
    foreach(intervention, values(models))
    (; path, system)
end

function petab_system(system)
    values = ModelingToolkit.initial_conditions(system)
    renames = Dict(
        p => ModelingToolkit.toparam(ModelingToolkit.Symbolics.variable(petab_name(getname(p))))
        for p in ModelingToolkit.parameters(system)
    )
    parametermap = [
        ModelingToolkit.Num(new) => ModelingToolkit.value(values[old]) for (old, new) in renames
    ]
    complete(substitute(ode_model(system), renames)), parametermap
end

function condition(index, label)
    positions = findall(s -> s.label == label, index)
    isempty(positions) && error("unknown label '$label'")
    branches = unique(Scheduling.branch(index[i].path) for i in positions)
    length(branches) == 1 || error("label '$label' lies in several branches: $branches")
    branch = only(branches)
    stretch = index[first(positions):last(positions)]
    any(s -> s.label != label && Scheduling.branch(s.path) == branch, stretch) &&
        error("label '$label' marks more than one stretch in branch $branch")
    (; branch, start = index[first(positions)].from)
end

"""
    PEtabParameter(schedule!::Schedule, id::Symbol; kwargs...)

Build a `PEtabParameter` for the model parameter `id` whose `value` is taken from `schedule!`.
"""
function PEtab.PEtabParameter(schedule!::Scheduling.Schedule, id::Symbol; kwargs...)
    (; system) = only_reaction_system(schedule!)
    parameters = filter(p -> petab_name(getname(p)) == id, ModelingToolkit.parameters(system))
    isempty(parameters) && error("schedule has no model parameter '$id'")
    value = ModelingToolkit.value(ModelingToolkit.initial_conditions(system)[only(parameters)])
    PEtab.PEtabParameter(id; value, kwargs...)
end

"""
    PEtabModel(schedule!::Schedule, observables, measurements, parameters;
        shared_parameters = Dict{Symbol, Symbol}(), kwargs...)

Build a `PEtabModel` that simulates `schedule!` as an ODE and compares it to `measurements`.

`schedule!` must advance through exactly one regulation model. Its `Instant`s must be `Adjust`s, which become initial values or events.

Each stage label used in `measurements` becomes a simulation condition that runs the label's branch of the schedule. The same label should not be used across parallel branches.

# Arguments
- `schedule!`: the `Schedule` to simulate.
- `observables`: `PEtabObservable`s, e.g. `[PEtabObservable("a", "scale * a₊mrnas", 0.1)]`.
- `measurements`: a `DataFrame` in PEtab measurement format whose `simulation_id` is a stage label and whose `time` is relative to the start of that label.
- `parameters`: the `PEtabParameter`s to estimate. This includes schedule model parameters to refit, parameters used in observables and any shared parameters, e.g. `PEtabParameter.([:a₊activation₊b₊at, :decay, :scale])`

# Keyword arguments
- `shared_parameters`: ties schedule model parameters to shared PEtab parameters, e.g. `Dict(:a₊mrna_decay => :decay, :b₊mrna_decay => :decay)`.
- `kwargs...`: passed on to `PEtabModel`.
"""
function PEtab.PEtabModel(
    schedule!::Scheduling.Schedule, observables, measurements, parameters;
    shared_parameters = Dict{Symbol, Symbol}(), kwargs...
)
    index = Scheduling.segments(schedule!)
    models = Scheduling.models(schedule!, index)
    (; system) = only_reaction_system(schedule!, models)
    instants = Dict(path => intervention(f!) for (path, f!) in models)
    symbolic = Dict(SciML.normalize_name(s) => ModelingToolkit.Num(s) for s in species(system))
    shared_assignments = [petab_name(name) => string(target) for (name, target) in shared_parameters]

    t = ModelingToolkit.Num(get_iv(system))

    labels = String.(string.(measurements.simulation_id))
    measurements = copy(measurements)
    measurements.time = Float64.(measurements.time)
    initials = Dict{String, Vector{Pair}}()
    events = PEtabEvent[]
    for label in unique(labels)
        (; branch, start) = condition(index, label)
        segments = filter(s -> Scheduling.ispathprefix(Scheduling.branch(s.path), branch), index)
        in_condition = labels .== label
        measurements.time[in_condition] .+= start
        stop = maximum(s.to for s in segments)
        all(time -> 0.0 <= time <= stop, measurements.time[in_condition]) ||
            error("times for label '$label' fall outside its branch (-$start to $(stop - start))")

        changes = assignments(segments, instants, symbolic)
        initials[label] = get(changes, 0.0, Pair[])
        for (time, changed) in changes
            time == 0.0 && continue
            push!(events, PEtabEvent(t == time, changed...; condition_ids = [label]))
            # PEtab.jl breaks on measurements exactly at an event time
            at_event = in_condition .& isapprox.(measurements.time, time; atol = 1e-9)
            measurements.time[at_event] .-= 1e-3
        end
    end

    targets = unique(first(pair) for pairs in values(initials) for pair in pairs)
    conditions = [
        PEtabCondition(id, (s => get(Dict(initial), s, 0.0) for s in targets)..., shared_assignments...)
        for (id, initial) in initials
    ]

    speciemap = [ModelingToolkit.Num(s) => 0.0 for s in species(system)]
    ode, parametermap = petab_system(system)
    PEtab.PEtabModel(ode, observables, measurements, parameters;
        simulation_conditions = conditions, events, speciemap, parametermap, kwargs...)
end

"""
    remake(schedule!::Schedule, result::EstimationResult, prob::PEtabODEProblem;
        condition = first(prob.model_info.simulation_info.conditionids[:experiment]))

Return a copy of `schedule!` that uses parameters from `result` in simulation `condition` of `prob`. `prob` must have been built from `schedule!` by `PEtabModel`.
"""
function Models.remake(
    schedule!::Scheduling.Schedule, result::PEtab.EstimationResult, prob::PEtab.PEtabODEProblem;
    condition = first(prob.model_info.simulation_info.conditionids[:experiment])
)
    (; path) = only_reaction_system(schedule!)
    values = Dict(SciML.normalize_name(p) => value for (p, value) in PEtab.get_ps(result, prob; condition))

    Models.remake(Scheduling.prebuild(schedule!, path), path, values)
end

end
