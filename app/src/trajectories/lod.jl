const PIXELS_PER_BIN = 4

abstract type TrajectoryLOD end

struct CountSample
    position::Int
    t::Float64
    value::Float64
end

struct CountBin # face
    index::Int
    first::CountSample
    minimum::CountSample
    maximum::CountSample
    last::CountSample
end

struct CountPyramid <: TrajectoryLOD
    series::CountSeries
    origin::Float64
    resolution::Float64
    levels::Vector{Vector{CountBin}}
end

struct FractionIntegral <: TrajectoryLOD
    series::FractionSeries
    area::Vector{Float64}
end

function bound_series(series::Series, to)
    isempty(series.ts) && return series
    last(series.ts) >= to && return series

    typeof(series)(
        ts=[series.ts; to],
        ys=[series.ys; last(series.ys)],
    )
end

function lod(series::CountSeries; from, to)
    series = bound_series(series, to)
    isempty(series.ts) && return CountPyramid(
        series,
        from,
        1.0,
        Vector{CountBin}[],
    )

    span = max(to - from, eps(Float64))
    levels = max(1, ceil(Int, log2(max(length(series.ts), 2))))
    resolution = span / 2.0^levels
    result = Vector{Vector{CountBin}}(undef, levels)
    result[1] = build_bins(series, from, resolution)
    for level in 2:levels
        result[level] = coarsen_bins(result[level - 1])
    end
    CountPyramid(series, from, resolution, result)
end

function lod(series::FractionSeries; from, to)
    series = bound_series(series, to)
    area = zeros(length(series.ts))
    for i in 2:length(series.ts)
        duration = series.ts[i] - series.ts[i - 1]
        area[i] = area[i - 1] + duration * series.ys[i - 1]
    end
    FractionIntegral(series, area)
end

function build_bins(series::CountSeries, origin, resolution)
    result = CountBin[]
    first_position = 1

    while first_position <= length(series.ts)
        index = floor(
            Int,
            (series.ts[first_position] - origin) / resolution,
        )
        last_position = first_position

        while last_position < length(series.ts)
            next_index = floor(
                Int,
                (series.ts[last_position + 1] - origin) / resolution,
            )
            next_index == index || break
            last_position += 1
        end

        first_sample = CountSample(
            first_position,
            series.ts[first_position],
            series.ys[first_position],
        )
        minimum_sample = first_sample
        maximum_sample = first_sample

        for position in (first_position + 1):last_position
            sample = CountSample(
                position,
                series.ts[position],
                series.ys[position],
            )
            sample.value < minimum_sample.value && (minimum_sample = sample)
            sample.value > maximum_sample.value && (maximum_sample = sample)
        end

        push!(
            result,
            CountBin(
                index,
                first_sample,
                minimum_sample,
                maximum_sample,
                CountSample(
                    last_position,
                    series.ts[last_position],
                    series.ys[last_position],
                ),
            ),
        )
        first_position = last_position + 1
    end

    result
end

function combine_bins(left::CountBin, right::CountBin, index)
    CountBin(
        index,
        left.first,
        left.minimum.value <= right.minimum.value ?
            left.minimum : right.minimum,
        left.maximum.value >= right.maximum.value ?
            left.maximum : right.maximum,
        right.last,
    )
end

function coarsen_bins(level)
    result = CountBin[]
    sizehint!(result, cld(length(level), 2))

    for bin in level
        index = bin.index >> 1
        if !isempty(result) && last(result).index == index
            result[end] = combine_bins(last(result), bin, index)
        else
            push!(result, CountBin(index, bin.first, bin.minimum, bin.maximum, bin.last))
        end
    end

    result
end

function series_slice(series::Series, from, to)
    isempty(series.ts) && return typeof(series)()
    (to < first(series.ts) || from > last(series.ts)) && return typeof(series)()

    first_position = max(searchsortedlast(series.ts, from), 1)
    last_position = min(
        searchsortedfirst(series.ts, to),
        length(series.ts),
    )

    typeof(series)(
        ts=series.ts[first_position:last_position],
        ys=series.ys[first_position:last_position],
    )
end

function emit_bin!(ts, ys, bin::CountBin)
    samples = [bin.first, bin.minimum, bin.maximum, bin.last]
    sort!(samples; by=sample -> sample.position)

    previous = 0
    for sample in samples
        sample.position == previous && continue
        push!(ts, sample.t)
        push!(ys, sample.value)
        previous = sample.position
    end

    nothing
end

bin_index(bin::CountBin) = bin.index
bin_index(index::Integer) = index

function resample(detail::CountPyramid; from, to, width)
    series = detail.series
    isempty(series.ts) && return CountSeries()
    (!isfinite(from) || !isfinite(to) ||
        to < first(series.ts) || from > last(series.ts)) && return CountSeries()
    from = max(from, first(series.ts))
    to = min(to, last(series.ts))

    width = max(round(Int, width), 1)
    target = (to - from) / width
    if isempty(detail.levels) || target <= detail.resolution
        return series_slice(series, from, to)
    end

    level = clamp(
        ceil(
            Int,
            log2(PIXELS_PER_BIN * target / detail.resolution),
        ) + 1,
        1,
        length(detail.levels),
    )
    bins = detail.levels[level]
    resolution = detail.resolution * 2.0^(level - 1)
    first_index = floor(Int, (from - detail.origin) / resolution)
    last_index = floor(Int, (to - detail.origin) / resolution)
    first_bin = searchsortedfirst(bins, first_index; by=bin_index)
    first_bin > 1 && (first_bin -= 1)

    ts = Float64[]
    ys = Float64[]
    for position in first_bin:length(bins)
        bin = bins[position]
        bin.index > last_index + 1 && break
        emit_bin!(ts, ys, bin)
    end

    CountSeries(; ts, ys)
end

function fraction_integral(detail::FractionIntegral, t)
    series = detail.series
    isempty(series.ts) && return (0.0, 0.0)

    bounded = clamp(t, first(series.ts), last(series.ts))
    position = searchsortedlast(series.ts, bounded)
    area = detail.area[position]
    area += (bounded - series.ts[position]) * series.ys[position]
    area, bounded - first(series.ts)
end

function resample(detail::FractionIntegral; from, to, width)
    series = detail.series
    isempty(series.ts) && return FractionSeries()
    (!isfinite(from) || !isfinite(to) ||
        to < first(series.ts) || from > last(series.ts)) && return FractionSeries()
    from = max(from, first(series.ts))
    to = min(to, last(series.ts))

    width = max(round(Int, width), 1)
    bin_count = max(1, cld(width, PIXELS_PER_BIN))
    exact = series_slice(series, from, to)
    length(exact.ts) <= bin_count && return exact

    edges = collect(range(from, to; length=bin_count + 1))
    values = Float64[]
    sizehint!(values, bin_count)

    for (left, right) in zip(edges, @view(edges[2:end]))
        left_area, left_covered = fraction_integral(detail, left)
        right_area, right_covered = fraction_integral(detail, right)
        covered = right_covered - left_covered
        push!(
            values,
            covered > 0 ?
                (right_area - left_area) / covered :
                NaN,
        )
    end

    FractionSeries(
        ts=edges,
        ys=[values; last(values)],
    )
end

function lod(catenation::Catenation{<:Series}; from, to)
    trajectories = Dict{Dimension, TrajectoryLOD}(
        dimension => lod(series; from, to)
        for (dimension, series) in catenation.trajectories
    )
    Catenation(catenation.segments, trajectories)
end
