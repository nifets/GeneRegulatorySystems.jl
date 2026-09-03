const DIMENSION_PATTERN = r"(?<group>.+)\.(?<kind>.+)"

struct Dimension
    kind::Symbol
    group::String
end

function Dimension(name::Symbol)
    matched = match(DIMENSION_PATTERN, String(name))
    isnothing(matched) && return Dimension(name, "")
    Dimension(Symbol(matched[:kind]), String(matched[:group]))
end

abstract type Series end

@kwdef struct CountSeries <: Series
    ts::Vector{Float64} = Float64[]
    ys::Vector{Float64} = Float64[]
end

@kwdef struct FractionSeries <: Series
    ts::Vector{Float64} = Float64[]
    ys::Vector{Float64} = Float64[]
end


seriestype(dimension::Dimension) = seriestype(dimension.kind)
seriestype(kind::Symbol) = seriestype(Val(kind))
seriestype(::Val) = CountSeries
seriestype(::Val{:activity}) = FractionSeries

struct Catenation{T}
    segments::Vector{Int}
    trajectories::Dict{Dimension, T}
end

Catenation(; segments, trajectories=Dict{Dimension, Series}()) = Catenation(segments, trajectories)

paths(index) = sort!(unique(
        prefix
        for segment in index
        for prefix in [ancestors(segment.path); string(segment.path)]
        if !isempty(prefix)
    ); by=pathorder)

describe(path, label) = isempty(label) ? String(path) : "$path ($label)"

function path_labels(index)
    labels = Dict{String, String}()
    ambiguous = Set{String}()

    for segment in index
        isempty(segment.label) && continue
        label = string(segment.label)
        for prefix in [ancestors(segment.path); string(segment.path)]
            isempty(prefix) && continue
            prefix in ambiguous && continue
            if get(labels, prefix, label) == label
                labels[prefix] = label
            else
                push!(ambiguous, prefix)
                delete!(labels, prefix)
            end
        end
    end

    labels
end

function cut(index)
    segments = collect(index)
    branches = Dict{String, Vector{Int}}()

    for (i, segment) in enumerate(segments)
        segment.count > 0 || continue
        push!(
            get!(Vector{Int}, branches, branch(segment.path)),
            i
        )
    end

    result = Catenation[]

    for segment_ids in values(branches)
        pending = Int[]

        for i in segment_ids
            segment = segments[i]

            if segment.from < segment.to
                if !isempty(pending)
                    push!(result, Catenation(segments=copy(pending)))
                    empty!(pending)
                end
                push!(result, Catenation(segments=[i]))
            else
                push!(pending, i)
            end
        end
        isempty(pending) || push!(result, Catenation(segments=copy(pending)))
    end
    sort!(result, by=catenation -> first(catenation.segments))
end

function place!(catenation::Catenation, t::Real, name::Symbol, value::Real)
    dimension = Dimension(name)
    series = get!(catenation.trajectories, dimension) do
        seriestype(dimension)()
    end
    push!(series.ts, Float64(t))
    push!(series.ys, Float64(value))
    nothing
end


place!(catenation::Catenation, event) =
    place!(catenation, event.t, event.name, event.value)

(series::Series)(t) =
    (i = searchsortedlast(series.ts, t); iszero(i) ? NaN : series.ys[i])
