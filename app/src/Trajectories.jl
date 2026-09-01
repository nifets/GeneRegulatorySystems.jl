module Trajectories

using GeneRegulatorySystems: Models
using GeneRegulatorySystems.Visualisation: Catenation, cut, place!

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
        label = get(primitive!.bindings, :label, "")
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
                count,
                into=isnothing(into) ? "" : String(into)
            ))
    end
    nothing
end

function prepare(sink::Sink)
    catenations = cut(sink.index)
    by_segment = Dict(
        segment => catenation
        for catenation in catenations
        for segment in catenation.segments
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
    by_kind = Dict{Symbol, Dict{Int, Catenation}}()
    for catenation in catenations
        split = Dict{Symbol, Catenation}()
        for (dimension, series) in catenation.series
            selected = get!(split, dimension.kind) do
                Catenation(segments=catenation.segments)
            end
            selected.series[dimension] = series
        end
        for (kind, selected) in split
            get!(
                Dict{Int, Catenation},
                by_kind,
                kind
            )[last(selected.segments)] = selected
        end
    end
    (; index=copy(sink.index), catenations=by_kind)
end

end
