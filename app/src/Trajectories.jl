module Trajectories

using WGLMakie
using GeneRegulatorySystems: Models, Scheduling
using GeneRegulatorySystems.Visualisation:
    Catenation,
    Series,
    CountSeries,
    FractionSeries,
    Dimension,
    cut, place!, seriestype

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

function prepare(sink::Sink; path="")
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
    by_kind = Dict{Symbol, Dict{Int, Catenation}}()
    for catenation in catenations
        split = Dict{Symbol, Catenation}()
        for (dimension, series) in catenation.trajectories
            selected = get!(split, dimension.kind) do
                Catenation(segments=catenation.segments)
            end
            selected.trajectories[dimension] = series
        end
        for (kind, selected) in split
            from = index[first(selected.segments)].from
            to = index[last(selected.segments)].to
            get!(
                Dict{Int, Catenation},
                by_kind,
                kind
            )[last(selected.segments)] = lod(selected; from, to)
        end
    end
    (; index, catenations=by_kind)
end

include("TrajectoryLOD.jl")

function sampling_window(axis)
    window = map(
        axis.finallimits,
        axis.scene.viewport;
        ignore_equal_values=true,
    ) do limits, viewport
        lower, upper = extrema(limits)
        from, to = lower[1], upper[1]
        width = max(widths(viewport)[1], 1)
        (; from, to, width)
    end
    WGLMakie.Observables.async_latest(window)
end

function bounded_xlimits(limits, from, to)
    lower, upper = extrema(limits)
    span = min(upper[1] - lower[1], to - from)
    left = clamp(lower[1], from, to - span)
    Rect2d((left, lower[2]), (span, upper[2] - lower[2]))
end

function bounded_interaction(interaction, from, to; stop_at_extent=false)
    function (event, axis)
        if stop_at_extent && event isa ScrollEvent && event.y < 0 &&
                widths(axis.targetlimits[])[1] >= to - from
            return Consume(true)
        end

        consumed = WGLMakie.Makie.process_interaction(interaction, event, axis)
        limits = axis.targetlimits[]
        bounded = bounded_xlimits(limits, from, to)
        limits == bounded || (axis.targetlimits[] = bounded)
        consumed
    end
end

function constrain_x!(axis, from, to)
    from < to || return

    for name in (:scrollzoom, :dragpan)
        interaction = deregister_interaction!(axis, name)
        register_interaction!(
            axis,
            name,
            bounded_interaction(
                interaction,
                from,
                to;
                stop_at_extent=name === :scrollzoom,
            ),
        )
    end

    nothing
end

function fit_y!(axis, values)
    minimum_y = Inf
    maximum_y = -Inf

    for value in values
        isfinite(value) || continue
        minimum_y = min(minimum_y, value)
        maximum_y = max(maximum_y, value)
    end

    isfinite(minimum_y) || return
    padding = minimum_y == maximum_y ?
        max(abs(minimum_y) * 0.05, 0.5) :
        (maximum_y - minimum_y) * 0.05
    ylims!(axis, minimum_y - padding, maximum_y + padding)
    nothing
end

function viewport_values(points, from, to)
    values = Float64[]
    left_value = nothing

    for point in points
        if point[1] <= from
            left_value = Float64(point[2])
        elseif point[1] <= to
            push!(values, point[2])
        else
            break
        end
    end

    isnothing(left_value) || pushfirst!(values, left_value)
    values
end

function aggregate(::Type{FractionSeries}, catenations; index, from, to, width)
    rows = index isa AbstractVector ? index : eachrow(index)
    grouped = Dict()

    for catenation in values(catenations)
        for (dimension, detail) in catenation.trajectories
            detail isa FractionIntegral || continue
            push!(get!(grouped, dimension, []), (catenation, detail))
        end
    end

    bin_count = max(1, cld(max(round(Int, width), 1), PIXELS_PER_BIN))
    edges = collect(range(from, to; length=bin_count + 1))
    result = Dict()
    for (dimension, members) in grouped
        values = Float64[]
        nlive = Int[]
        paths = Vector{String}[]

        for (left, right) in zip(edges, @view(edges[2:end]))
            active = Float64[]
            live_paths = String[]

            for (catenation, detail) in members
                branch_from = rows[first(catenation.segments)].from
                branch_to = rows[last(catenation.segments)].to
                covered_from = max(left, branch_from)
                covered_to = min(right, branch_to)
                covered_from < covered_to || continue

                left_area, left_covered = fraction_integral(detail, covered_from)
                right_area, right_covered = fraction_integral(detail, covered_to)
                duration = right_covered - left_covered
                duration > 0 || continue

                value = (right_area - left_area) / duration
                isfinite(value) || continue
                push!(active, value)

                for position in catenation.segments
                    segment = rows[position]
                    if segment.from < covered_to && segment.to > covered_from
                        push!(live_paths, String(segment.path))
                    end
                end
            end

            push!(values, isempty(active) ? NaN : sum(active) / length(active))
            push!(nlive, length(active))
            push!(paths, unique!(live_paths))
        end

        result[dimension] = (; ts=copy(edges), values, nlive, paths)
    end
    result
end

function linear_quantile(values, q)
    position = 1 + (length(values) - 1) * q
    lower = floor(Int, position)
    upper = ceil(Int, position)
    lower == upper && return values[lower]

    weight = position - lower
    values[lower] * (1 - weight) + values[upper] * weight
end

function aggregate(::Type{CountSeries}, catenations; index, from, to, width)
    rows = index isa AbstractVector ? index : eachrow(index)
    grouped = Dict()

    for catenation in values(catenations)
        for (dimension, detail) in catenation.trajectories
            detail isa CountPyramid || continue
            push!(get!(grouped, dimension, []), (catenation, detail))
        end
    end

    bin_count = max(1, cld(max(round(Int, width), 1), PIXELS_PER_BIN))
    step = (to - from) / bin_count
    ts = [from + (i - 0.5) * step for i in 1:bin_count]
    result = Dict()

    for (dimension, members) in grouped
        median = Float64[]
        q25 = Float64[]
        q75 = Float64[]
        nlive = Int[]
        paths = String[]

        for (catenation, _) in members
            append!(
                paths,
                String(rows[position].path)
                for position in catenation.segments
            )
        end
        unique!(paths)

        for t in ts
            values = Float64[]

            for (catenation, detail) in members
                branch_from = rows[first(catenation.segments)].from
                branch_to = rows[last(catenation.segments)].to
                branch_from <= t <= branch_to || continue

                position = searchsortedlast(detail.series.ts, t)
                iszero(position) && continue
                push!(values, detail.series.ys[position])
            end

            sort!(values)
            push!(nlive, length(values))
            if isempty(values)
                push!(median, NaN)
                push!(q25, NaN)
                push!(q75, NaN)
            else
                push!(median, linear_quantile(values, 0.5))
                push!(q25, linear_quantile(values, 0.25))
                push!(q75, linear_quantile(values, 0.75))
            end
        end

        result[dimension] = (; ts, median, q25, q75, nlive, paths)
    end

    result
end

function render!(
    axis,
    ::Type{FractionSeries},
    catenations;
    index,
    group_colors,
    visible,
    window,
)
    rows = map(window) do view
        aggregate(FractionSeries, catenations; index, view...)
    end
    dimensions = sort!(collect(keys(rows[])); by=dimension -> dimension.group)
    visibilities = [visible(dimension.group) for dimension in dimensions]
    layout = if isempty(dimensions)
        Observable(Dict{Dimension, Int}())
    else
        map(visibilities...) do shown...
            visible_dimensions = (
                dimension
                for (dimension, isshown) in zip(dimensions, shown)
                if isshown
            )
            Dict(
                dimension => row
                for (row, dimension) in enumerate(visible_dimensions)
            )
        end
    end

    visual = map(rows, layout) do all_rows, row_numbers
        shown = sort!(collect(keys(row_numbers)); by=dimension -> row_numbers[dimension])
        row_count = max(length(shown), 1)
        edges = isempty(all_rows) ? [0.0, 1.0] : first(values(all_rows)).ts
        colors = fill(RGBAf(0, 0, 0, 0), length(edges) - 1, row_count)
        tooltips = fill("", size(colors))

        for dimension in shown
            data = all_rows[dimension]
            gene = dimension.group
            row = row_numbers[dimension]
            color = to_color(get(group_colors, gene, :gray))

            for column in eachindex(data.values)
                value = data.values[column]
                isnan(value) && continue

                colors[column, row] = RGBAf(
                    color,
                    clamp(value, 0.0, 1.0),
                )
                tooltips[column, row] = join([
                    "gene: $gene",
                    "activity: $(round(value; digits=2))",
                    "paths: $(join(data.paths[column], ", "))",
                    "time: $(round(data.ts[column]; digits=2))–$(round(data.ts[column + 1]; digits=2))",
                ], "\n")
            end
        end

        yedges = collect(range(0.5; step=1.0, length=row_count + 1))
        (; edges, yedges, colors, tooltips, shown)
    end

    current = Ref(visual[])
    plot = heatmap!(
        axis,
        current[].edges,
        current[].yedges,
        current[].colors;
        interpolate=false,
        inspector_label=(_, (column, row), _) ->
            current[].tooltips[column, row],
    )

    on(visual; update=true) do updated
        current[] = updated
        update!(plot, updated.edges, updated.yedges, updated.colors)
        axis.yticks[] = (
            eachindex(updated.shown),
            getproperty.(updated.shown, :group),
        )
        ylims!(axis, 0.5, max(length(updated.shown), 1) + 0.5)
    end
    axis
end

function render!(
    axis,
    ::Type{CountSeries},
    ::Val{:raw},
    catenations;
    index,
    group_colors,
    visible,
    window
)
    point_sets = Observable[]
    visibilities = Observable[]

    for catenation in values(catenations)
        paths = unique(
            index[position].path
            for position in catenation.segments
        )

        catenation_to = index[last(catenation.segments)].to

        for (dimension, detail) in catenation.trajectories
            detail isa CountPyramid || continue
            points = map(window) do view
                series = resample(detail; view...)
                isempty(series.ts) && return Point2f[]

                ts = copy(series.ts)
                ys = copy(series.ys)

                display_to = min(catenation_to, view.to)
                if last(ts) < display_to
                    push!(ts, display_to)
                    push!(ys, last(ys))
                end
                Point2f.(ts, ys)
            end

            gene = dimension.group
            visibility = visible(gene)
            tooltip = join([
                "gene: $gene",
                "path: $(join(paths, " → "))",
            ], "\n")

            push!(point_sets, points)
            push!(visibilities, visibility)

            stairs!(
                axis,
                points;
                step=:post,
                color=get(group_colors, gene, :gray),
                visible=visibility,
                inspector_label=(_, _, _) -> tooltip,
            )
        end
    end

    update_y = function (_...)
        view = window[]
        visible_points = (
            viewport_values(points[], view.from, view.to)
            for (points, visibility) in zip(point_sets, visibilities)
            if visibility[]
        )
        fit_y!(
            axis,
            (value for values in visible_points for value in values),
        )
    end
    onany(update_y, window, visibilities...; update=true)

    axis
end

function render!(
    axis,
    ::Type{CountSeries},
    ::Val{:aggregate},
    catenations;
    index,
    group_colors,
    visible,
    window,
)
    rows = map(window) do view
        aggregate(CountSeries, catenations; index, view...)
    end
    dimensions = sort!(collect(keys(rows[])); by=dimension -> dimension.group)
    visibilities = Observable[]

    for dimension in dimensions
        gene = dimension.group
        color = get(group_colors, gene, :gray)
        visibility = visible(gene)
        data = Ref(rows[][dimension])

        lower = Point2f.(data[].ts, data[].q25)
        upper = Point2f.(data[].ts, data[].q75)
        center = Point2f.(data[].ts, data[].median)

        band_plot = band!(
            axis,
            lower,
            upper;
            color=(color, 0.18),
            visible=visibility,
            inspectable=false,
        )
        line_plot = lines!(
            axis,
            center;
            color,
            linewidth=1.5,
            visible=visibility,
            inspector_label=(_, i, _) -> join([
                "gene: $gene",
                "median: $(round(data[].median[i]; digits=2))",
                "q25–q75: $(round(data[].q25[i]; digits=2))–$(round(data[].q75[i]; digits=2))",
                "branches: $(data[].nlive[i])",
            ], "\n"),
        )
        push!(visibilities, visibility)

        on(rows) do all_rows
            data[] = all_rows[dimension]
            update!(
                band_plot,
                Point2f.(data[].ts, data[].q25),
                Point2f.(data[].ts, data[].q75),
            )
            update!(line_plot, Point2f.(data[].ts, data[].median))
        end
    end

    update_y = function (_...)
        visible_data = (
            rows[][dimension]
            for (dimension, visibility) in zip(dimensions, visibilities)
            if visibility[]
        )
        fit_y!(
            axis,
            (
                value
                for data in visible_data
                for boundary in (data.q25, data.q75)
                for value in boundary
            ),
        )
    end
    onany(update_y, rows, visibilities...; update=true)

    axis
end

function render(
    trajectories;
    tracks,
    selected_genes,
    default_genes,
    group_colors,
    aggregate_mode=:raw,
)
    aggregate_mode = Symbol(aggregate_mode)
    aggregate_mode in (:raw, :aggregate) || throw(
        ArgumentError("aggregate_mode must be :raw or :aggregate"),
    )
    tracks = collect(tracks)
    figure = Figure(size=(1200, 180 * max(1, length(tracks))))
    axes = Axis[]
    from = isempty(trajectories.index) ? 0.0 :
        minimum(segment.from for segment in trajectories.index)
    to = isempty(trajectories.index) ? 1.0 :
        maximum(segment.to for segment in trajectories.index)

    visible(gene) = lift(selected_genes) do selected
        if selected isa AbstractString
            gene == selected
        elseif isempty(selected)
            gene in default_genes
        else
            gene in selected
        end
    end

    for (row, kind) in enumerate(tracks)
        axis = Axis(
            figure[row, 1];
            ylabel=kind,
            xlabel=row == length(tracks) ? "time" : "",
            backgroundcolor=:transparent,
            xzoomlock=false,
            yzoomlock=true,
            xpanlock=false,
            ypanlock=true,
        )
        push!(axes, axis)
        constrain_x!(axis, from, to)

        kind = Symbol(kind)
        catenations = get(
            trajectories.catenations,
            kind,
            Dict{Int, Catenation}(),
        )
        type = seriestype(kind)
        window = sampling_window(axis)
        if type === CountSeries
            render!(
                axis,
                type,
                Val(aggregate_mode),
                catenations;
                index=trajectories.index,
                group_colors,
                visible,
                window,
            )
        else
            render!(
                axis,
                type,
                catenations;
                index=trajectories.index,
                group_colors,
                visible,
                window,
            )
        end
    end

    length(axes) > 1 && linkxaxes!(axes...)
    if !isempty(trajectories.index) && !isempty(axes)
        from < to && xlims!(first(axes), from, to)
    end

    DataInspector(figure; fontsize=12, show_bbox_indicators=false)
    figure
end

end
