@kwdef mutable struct Sink
    lock::ReentrantLock = ReentrantLock()
    i::Int = 0
    records = []
    dimensions::Dict{Symbol, Dimension} = Dict{Symbol, Dimension}()
end

(sink::Sink)(; context...) = nothing


function (sink::Sink)(
    state;
    path,
    from,
    primitive!,
    into=nothing,
    context...,
)
    lock(sink.lock) do
        sink.i += 1
        to = Models.t(state)
        model = primitive!.path
        f! = Models.unwrap(primitive!.f!)
        instant = f! isa Models.Instant
        description = instant ? Models.describe(f!) : nothing
        bindings = primitive!.bindings
        label = get(bindings, :label, "")
        color = get(bindings, :color, get(bindings, :colour, ""))
        count = 0
        trajectories = Dict{Dimension, Series}()
        if !isnothing(into)
            Models.each_event(state) do t, name, value
                dimension = get!(() -> Dimension(Symbol(name)), sink.dimensions, Symbol(name))
                series = get!(() -> seriestype(dimension)(), trajectories, dimension)
                push!(series.ts, Float64(t))
                push!(series.ys, Float64(value))
                count += 1
            end
        end
        push!(
            sink.records,
            (
                i=sink.i,
                path=String(path),
                from=Float64(from),
                to=Float64(to),
                model=string(model),
                label=string(label),
                color=string(color),
                count,
                into=isnothing(into) ? nothing : String(into),
                instant,
                description=description isa Models.Label ? description.label : "",
                trajectories,
            ))
    end
    nothing
end

struct Trace{C}
    records::Vector
    catenations::C
end

describe(segment) = describe(segment.path, segment.label)

function catenation_color(records, catenation::Catenation)
    color = ""
    for position in catenation.segments
        isempty(records[position].color) || (color = String(records[position].color))
    end
    color
end

lineage(records, catenation::Catenation) =
    branch(records[first(catenation.segments)].path)

function describe(records, catenation::Catenation)
    label = ""
    for position in catenation.segments
        isempty(records[position].label) || (label = String(records[position].label))
    end
    describe(lineage(records, catenation), label)
end

function describer(records, catenation::Catenation)
    starts = [records[i].from for i in catenation.segments]
    t -> let record = records[catenation.segments[max(1, searchsortedlast(starts, t))]]
        describe(branch(record.path), record.label)
    end
end

function assemble(records, catenation::Catenation)
    length(catenation.segments) == 1 && return Catenation(
        catenation.segments,
        records[only(catenation.segments)].trajectories,
    )

    trajectories = Dict{Dimension, Series}()
    for position in catenation.segments
        for (dimension, series) in records[position].trajectories
            target = get!(() -> seriestype(dimension)(), trajectories, dimension)
            append!(target.ts, series.ts)
            append!(target.ys, series.ys)
        end
    end
    Catenation(catenation.segments, trajectories)
end

channels(sink::Sink) = sort!(unique(
    segment.into for segment in sink.records if segment.count > 0 && !isnothing(segment.into)
))

function instants(sink::Sink; path="", channel=nothing)
    marks = Dict{Float64, Vector{String}}()
    for segment in sink.records
        segment.instant || continue
        isnothing(channel) || segment.into == channel || continue
        Scheduling.ispathprefix(branch(segment.path), path) || continue
        entry = isempty(segment.label) ? segment.description : "$(segment.description) ($(segment.label))"
        push!(get!(Vector{String}, marks, segment.from), entry)
    end
    sort!([(; t, entries) for (t, entries) in marks]; by = mark -> mark.t)
end

function catenate(sink::Sink; path="", channel=nothing)
    records = filter(sink.records) do segment
        Scheduling.ispathprefix(path, segment.path) &&
            (isnothing(channel) || segment.into == channel)
    end
    Trace(records, [assemble(records, catenation) for catenation in cut(records)])
end

function within(trace::Trace, catenation::Catenation, path)
    segment = trace.records[first(catenation.segments)].path
    Scheduling.ispathprefix(path, segment) ||
        Scheduling.ispathprefix(branch(segment), branch(path))
end

select(trace::Trace{<:AbstractVector}, path) = isempty(path) ? trace : Trace(
    trace.records,
    filter(catenation -> within(trace, catenation, path), trace.catenations),
)

select(trace::Trace{<:AbstractDict}, path) = isempty(path) ? trace : Trace(
    trace.records,
    Dict(
        kind => filter(
            pair -> within(trace, last(pair), path),
            group,
        )
        for (kind, group) in trace.catenations
    ),
)
