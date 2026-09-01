const BRANCH_PATTERN = r"^(?:.*/\d+)?"
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

function FractionSeries(these::CountSeries, others::CountSeries)
    ts = Float64[]
    ys = Float64[]
    xs1 = zip(these.ts, these.ys)
    xs2 = zip(others.ts, others.ys)
    t1 = t2 = 0.0
    x1 = x2 = x1′ = x2′ = 0.0
    next1 = next2 = (1, 1)
    sentinel = ((Inf, 0.0), (0, 0))
    while true
        t1′ = t1
        if t1 ≤ t2
            ((t1, x1′), next1) = something(iterate(xs1, next1), sentinel)
        end
        if t2 ≤ t1′
            ((t2, x2′), next2) = something(iterate(xs2, next2), sentinel)
        end
        isfinite(t1) || isfinite(t2) || break
        t1 ≤ t2 && (x1 = x1′)
        t2 ≤ t1 && (x2 = x2′)
        push!(ts, min(t1, t2))
        push!(ys, iszero(x1) ? 0.0 : x1 / (x1 + x2))
    end

    FractionSeries(; ts, ys)
end

seriestype(dimension::Dimension) = seriestype(dimension.kind)
seriestype(kind::Symbol) = seriestype(Val(kind))
seriestype(::Val) = CountSeries
seriestype(::Val{:activity}) = FractionSeries

@kwdef struct Catenation
    segments::Vector{Int}
    series::Dict{Dimension, Series} = Dict{Dimension, Series}()
end

branch(path::AbstractString) = string(match(BRANCH_PATTERN, path).match)

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
    series = get!(catenation.series, dimension) do
        CountSeries()
    end
    push!(series.ts, Float64(t))
    push!(series.ys, Float64(value))
    nothing
end

function augment!(catenation::Catenation)
    for (dimension, active) in collect(catenation.series)
        dimension.kind == :active || continue
        inactive_dimension = Dimension(:inactive, dimension.group)
        catenation.series[Dimension(:activity, dimension.group)] =
            if haskey(catenation.series, inactive_dimension)
                FractionSeries(
                    active,
                    catenation.series[inactive_dimension],
                )
            else
                FractionSeries(
                    ts=active.ts,
                    ys=copy(active.ys),
                )
            end
    end

    nothing
end

place!(catenation::Catenation, event) =
    place!(catenation, event.t, event.name, event.value)
