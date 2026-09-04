@kwdef mutable struct Sink
    lock::ReentrantLock = ReentrantLock()
    i::Int = 0
    index = []

    is::Vector{Int} = Int[]
    ts::Vector{Float64} = Float64[]
    names::Vector{Symbol} = Symbol[]
    values::Vector{Float64} = Float64[]
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
        bindings = primitive!.bindings
        label = get(bindings, :label, "")
        color = get(bindings, :color, get(bindings, :colour, ""))
        count = 0
        if !isnothing(into)
            Models.each_event(state) do t, name, value
                push!(sink.is, sink.i)
                push!(sink.ts, Float64(t))
                push!(sink.names, Symbol(name))
                push!(sink.values, Float64(value))
                count += 1
            end
        end
        push!(
            sink.index,
            (
                i=sink.i,
                path=String(path),
                from=Float64(from),
                to=Float64(to),
                model=string(model),
                label=string(label),
                color=string(color),
                count,
                into=isnothing(into) ? "" : String(into)
            ))
    end
    nothing
end

struct Trace{C}
    index::Vector
    catenations::C
end

describe(segment) = describe(segment.path, segment.label)

function catenation_color(index, catenation::Catenation)
    color = ""
    for position in catenation.segments
        isempty(index[position].color) || (color = String(index[position].color))
    end
    color
end

lineage(index, catenation::Catenation) =
    branch(index[first(catenation.segments)].path)

function describe(index, catenation::Catenation)
    label = ""
    for position in catenation.segments
        isempty(index[position].label) || (label = String(index[position].label))
    end
    describe(lineage(index, catenation), label)
end

function catenate(sink::Sink; path="")
    index = filter(sink.index) do segment
        Scheduling.ispathprefix(path, segment.path)
    end
    catenations = cut(index)
    by_segment = Dict(
        index[position].i => catenation
        for catenation in catenations
        for position in catenation.segments
    )
    for j in eachindex(sink.is)
        segment = sink.is[j]
        haskey(by_segment, segment) || continue
        place!(
            by_segment[segment],
            sink.ts[j],
            sink.names[j],
            sink.values[j]
        )
    end
    Trace(index, catenations)
end
