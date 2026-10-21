function sampling_window(axis; interval=0.05)
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
    WGLMakie.Makie.Observables.throttle(interval, window)
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
    register_interaction!(axis, :reset_view) do event, axis
        if event isa MouseEvent &&
                event.type == WGLMakie.Makie.MouseEventTypes.leftdoubleclick
            lower, upper = extrema(axis.targetlimits[])
            axis.targetlimits[] = Rect2d(
                (from, lower[2]),
                (to - from, upper[2] - lower[2]),
            )
            return Consume(true)
        end
        Consume(false)
    end

    nothing
end

function set_ylimits!(axis, lower, upper)
    limits = axis.targetlimits[]
    axis.targetlimits[] = Rect2d(
        (limits.origin[1], lower),
        (limits.widths[1], upper - lower),
    )
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
    set_ylimits!(axis, minimum_y - padding, maximum_y + padding)
    nothing
end

function branches_by_dimension(::Type{T}, catenations; records, selected_genes, max_branches) where {T}
    rows = records isa AbstractVector ? records : eachrow(records)
    result = Dict{
        Dimension,
        Vector{@NamedTuple{
            detail::T,
            from::Float64,
            to::Float64,
            weight::Float64,
            sampled::Bool,
        }},
    }()

    for catenation in values(catenations)
        detail_from = rows[first(catenation.segments)].from
        detail_to = rows[last(catenation.segments)].to
        sampled = all(i -> rows[i].from == rows[i].to, catenation.segments)

        for (dimension, detail) in catenation.trajectories
            detail isa T || continue
            dimension.group in selected_genes || continue
            push!(
                get!(valtype(result), result, dimension),
                (;
                    detail,
                    from=Float64(detail_from),
                    to=Float64(detail_to),
                    weight=1.0,
                    sampled,
                ),
            )
        end
    end

    for (dimension, branches) in result
        length(branches) <= max_branches && continue
        sampled = branches[begin:cld(length(branches), max_branches):end]
        weight = length(branches) / length(sampled)
        result[dimension] = [(; branch..., weight) for branch in sampled]
    end

    result
end

function quantiles(samples)
    n = length(samples)
    n == 1 && return (samples[1], samples[1], samples[1])
    at(q) = begin
        position = 1 + (n - 1) * q
        lower = floor(Int, position)
        upper = ceil(Int, position)
        lower == upper && return partialsort!(samples, lower)
        weight = position - lower
        partialsort!(samples, lower) * (1 - weight) +
            partialsort!(samples, upper) * weight
    end
    (at(0.5), at(0.25), at(0.75))
end

function aggregate(::Type{CountSeries}, grouped; from, to, width)
    bin_count = max(1, cld(max(round(Int, width), 1), PIXELS_PER_BIN))
    step = (to - from) / bin_count
    ts = [from + (i - 0.5) * step for i in 1:bin_count]
    result = Dict{Dimension, @NamedTuple{
        ts::Vector{Float64},
        median::Vector{Float64},
        q25::Vector{Float64},
        q75::Vector{Float64},
        nlive::Vector{Int},
    }}()

    for (dimension, branches) in grouped
        windows = [
            (branch, resample(branch.detail; from, to, width))
            for branch in branches
        ]
        filter!(pair -> !isempty(last(pair).ts), windows)

        cursors = ones(Int, length(windows))
        samples = Vector{Float64}(undef, length(windows))
        median = Float64[]
        q25 = Float64[]
        q75 = Float64[]
        nlive = Int[]

        for t in ts
            count = 0
            live = 0.0
            right = t + step / 2

            for (position, (branch, series)) in enumerate(windows)
                cursor = cursors[position]
                if branch.sampled
                    while cursor < length(series.ts) && series.ts[cursor + 1] < right
                        cursor += 1
                    end
                    cursors[position] = cursor
                    right - step <= series.ts[cursor] < right || continue
                else
                    branch.from <= t <= branch.to || continue
                    while cursor < length(series.ts) && series.ts[cursor + 1] <= t
                        cursor += 1
                    end
                    cursors[position] = cursor
                    series.ts[cursor] <= t || continue
                end
                count += 1
                live += branch.weight
                samples[count] = series.ys[cursor]
            end

            push!(nlive, round(Int, live))
            if iszero(count)
                push!(median, NaN)
                push!(q25, NaN)
                push!(q75, NaN)
            else
                middle, lower, upper = quantiles(@view samples[1:count])
                push!(median, middle)
                push!(q25, lower)
                push!(q75, upper)
            end
        end

        keep = all(branch -> branch.sampled, branches) ?
            findall(!isnan, median) : eachindex(ts)
        result[dimension] = (;
            ts=ts[keep],
            median=median[keep],
            q25=q25[keep],
            q75=q75[keep],
            nlive=nlive[keep],
        )
    end

    result
end

function integral_at(detail::FractionIntegral, t, cursor)
    series = detail.series
    bounded = clamp(t, first(series.ts), last(series.ts))
    while cursor < length(series.ts) && series.ts[cursor + 1] <= bounded
        cursor += 1
    end
    area = detail.area[cursor] + (bounded - series.ts[cursor]) * series.ys[cursor]
    (area, bounded - first(series.ts), cursor)
end

function aggregate(::Type{FractionSeries}, grouped; from, to, width)
    bin_count = max(1, cld(max(round(Int, width), 1), PIXELS_PER_BIN))
    edges = collect(range(from, to; length=bin_count + 1))
    result = Dict{Dimension, @NamedTuple{
        ts::Vector{Float64},
        values::Vector{Float64},
        nlive::Vector{Int},
    }}()

    for (dimension, branches) in grouped
        visible = filter(
            branch -> branch.to > from && branch.from < to && !isempty(branch.detail.series.ts),
            branches,
        )
        cursors = ones(Int, length(visible))
        values = Float64[]
        nlive = Int[]

        for (left, right) in zip(edges, @view(edges[2:end]))
            total = 0.0
            count = 0
            live = 0.0

            for (position, branch) in enumerate(visible)
                covered_from = max(left, branch.from)
                covered_to = min(right, branch.to)
                covered_from < covered_to || continue

                left_area, left_covered, cursor =
                    integral_at(branch.detail, covered_from, cursors[position])
                right_area, right_covered, cursor =
                    integral_at(branch.detail, covered_to, cursor)
                cursors[position] = cursor

                duration = right_covered - left_covered
                duration > 0 || continue

                value = (right_area - left_area) / duration
                isfinite(value) || continue
                total += value
                count += 1
                live += branch.weight
            end

            push!(values, iszero(count) ? NaN : total / count)
            push!(nlive, round(Int, live))
        end

        result[dimension] = (; ts=copy(edges), values, nlive)
    end

    result
end

function render!(
    axis,
    ::Type{FractionSeries},
    catenations;
    records,
    group_colors,
    selected_genes,
    window,
    max_branches,
)
    grouped = branches_by_dimension(
        FractionIntegral,
        catenations;
        records,
        selected_genes,
        max_branches,
    )
    rows = map(window) do view
        aggregate(FractionSeries, grouped; view...)
    end
    dimensions = sort!(collect(keys(rows[])); by=dimension -> dimension.group)
    row_numbers = Dict(
        dimension => row
        for (row, dimension) in enumerate(dimensions)
    )

    visual = map(rows) do all_rows
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
                    "branches: $(data.nlive[column])",
                    "time: $(round((data.ts[column] + data.ts[column + 1]) / 2; digits=2))",
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
        set_ylimits!(axis, 0.5, max(length(updated.shown), 1) + 0.5)
    end
    axis
end

function render!(
    axis,
    ::Type{CountSeries},
    ::Val{:raw},
    catenations;
    records,
    group_colors,
    selected_genes,
    window,
    max_branches,
    dimming,
    connect_samples=false,
)
    strands = Dict{Tuple{String, Bool}, Vector{NamedTuple}}()

    for catenation in values(catenations)
        describe_at = describer(records, catenation)
        from = records[first(catenation.segments)].from
        to = records[last(catenation.segments)].to
        sampled = all(i -> records[i].from == records[i].to, catenation.segments)

        for (dimension, detail) in catenation.trajectories
            detail isa CountPyramid || continue
            dimension.group in selected_genes || continue
            push!(
                get!(Vector{NamedTuple}, strands, (dimension.group, sampled)),
                (; detail, from, to, describe_at),
            )
        end
    end

    buffers = Ref[]

    for ((gene, sampled), entries) in strands
        data = map(window) do view
            points = Point2f[]
            owners = String[]

            for entry in entries
                (entry.to < view.from || entry.from > view.to) && continue
                series = resample(entry.detail; view...)
                isempty(series.ts) && continue

                if !isempty(points)
                    push!(points, Point2f(NaN, NaN))
                    push!(owners, "")
                end
                for (t, y) in zip(series.ts, series.ys)
                    push!(points, Point2f(t, y))
                    push!(owners, entry.describe_at(t))
                end

                display_to = min(entry.to, view.to)
                if !sampled && last(series.ts) < display_to
                    push!(points, Point2f(display_to, last(series.ys)))
                    push!(owners, entry.describe_at(display_to))
                end
            end

            (; points, owners)
        end

        current = Ref(data[])
        push!(buffers, current)

        inspector_label = (_, i, _) -> begin
            buffer = current[]
            checkbounds(Bool, buffer.points, i) && !isempty(buffer.owners[i]) || return ""
            point = buffer.points[i]
            join([
                "gene: $gene",
                "path: $(buffer.owners[i])",
                "time: $(round(point[1]; digits=2))",
                "value: $(round(point[2]; digits=2))",
            ], "\n")
        end
        color = get(group_colors, gene, :gray)
        plot = sampled && !connect_samples ?
            scatter!(
                axis,
                current[].points;
                color,
                alpha=dimming,
                markersize=3,
                inspector_label,
            ) :
            sampled ?
            scatterlines!(
                axis,
                current[].points;
                color,
                alpha=dimming,
                linewidth=0.75,
                markersize=3,
                inspector_label,
            ) :
            lines!(
                axis,
                current[].points;
                color,
                alpha=dimming,
                inspector_label,
            )

        on(data) do updated
            current[] = updated
            update!(plot; arg1=updated.points)
        end
    end

    update_y = function (_...)
        view = window[]
        fit_y!(
            axis,
            (
                point[2]
                for buffer in buffers
                for point in buffer[].points
                if isfinite(point[2]) && view.from <= point[1] <= view.to
            ),
        )
    end
    on(update_y, window; update=true)

    axis
end

function render!(
    axis,
    ::Type{CountSeries},
    ::Val{:aggregate},
    catenations;
    records,
    group_colors,
    selected_genes,
    window,
    max_branches,
    dimming,
    connect_samples=false,
)
    grouped = branches_by_dimension(
        CountPyramid,
        catenations;
        records,
        selected_genes,
        max_branches,
    )
    rows = map(window) do view
        aggregate(CountSeries, grouped; view...)
    end
    dimensions = sort!(collect(keys(rows[])); by=dimension -> dimension.group)

    for dimension in dimensions
        gene = dimension.group
        color = get(group_colors, gene, :gray)
        data = Ref(rows[][dimension])

        lower = Point2f.(data[].ts, data[].q25)
        upper = Point2f.(data[].ts, data[].q75)
        center = Point2f.(data[].ts, data[].median)

        band_plot = band!(
            axis,
            lower,
            upper;
            color=(color, 0.18),
            inspectable=false,
        )
        line_plot = lines!(
            axis,
            center;
            color,
            linewidth=1.5,
            inspector_label=(_, i, _) -> join([
                "gene: $gene",
                "median: $(round(data[].median[i]; digits=2))",
                "q25–q75: $(round(data[].q25[i]; digits=2))–$(round(data[].q75[i]; digits=2))",
                "branches: $(data[].nlive[i])",
                "time: $(round(data[].ts[i]; digits=2))",
            ], "\n"),
        )
        on(rows) do all_rows
            data[] = all_rows[dimension]
            update!(
                band_plot,
                Point2f.(data[].ts, data[].q25),
                Point2f.(data[].ts, data[].q75),
            )
            update!(line_plot; arg1=Point2f.(data[].ts, data[].median))
        end
    end

    update_y = function (_...)
        visible_data = (
            rows[][dimension]
            for dimension in dimensions
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
    on(update_y, rows; update=true)

    axis
end

function extent(trajectories::Trace{<:AbstractDict})
    segments = unique(
        segment
        for group in values(trajectories.catenations)
        for catenation in values(group)
        for segment in catenation.segments
    )
    isempty(segments) && return 0.0, 1.0
    records = trajectories.records[segments]
    minimum(record.from for record in records), maximum(record.to for record in records)
end

function render(
    trajectories::Trace{<:AbstractDict};
    tracks,
    selected_genes,
    group_colors,
    aggregate_mode=:raw,
    gene_limit=10,
    max_branches=400,
    dimming=0.7,
    connect_samples=false,
    instants=(),
)
    aggregate_mode = Symbol(aggregate_mode)
    aggregate_mode in (:raw, :aggregate) || throw(
        ArgumentError("aggregate_mode must be :raw or :aggregate"),
    )
    tracks = Symbol.(collect(tracks))
    selected_genes = selected_genes isa Observable ? selected_genes[] : selected_genes
    selected_genes = selected_genes isa AbstractString ?
        [String(selected_genes)] : unique(string.(collect(selected_genes)))
    selected_genes = Set(last(selected_genes, gene_limit))
    figure = Figure(size=(1200, 180 * max(1, length(tracks))))
    axes = Axis[]
    from, to = extent(trajectories)
    margin = 0.005 * (to - from)
    from, to = from - margin, to + margin

    for (row, kind) in enumerate(tracks)
        type = seriestype(kind)
        axis = Axis(
            figure[row, 1];
            ylabel=string(kind),
            xlabel=row == length(tracks) ? "time" : "",
            backgroundcolor=:transparent,
            limits=(from < to ? (from, to) : nothing, nothing),
            xzoomlock=false,
            yzoomlock=true,
            xpanlock=false,
            ypanlock=true,
            panbutton=Mouse.left,
            (type === CountSeries ? (ytickformat = values -> [
                v == round(v) ? string(round(Int, v)) : string(v)
                for v in values
            ],) : NamedTuple())...
        )
        push!(axes, axis)
        constrain_x!(axis, from, to)

        catenations = get(
            trajectories.catenations,
            kind,
            Dict{Int, Catenation}(),
        )
        window = sampling_window(axis)
        if type === CountSeries
            render!(
                axis,
                type,
                Val(aggregate_mode),
                catenations;
                records=trajectories.records,
                group_colors,
                selected_genes,
                window,
                max_branches,
                dimming,
                connect_samples,
            )
        else
            render!(
                axis,
                type,
                catenations;
                records=trajectories.records,
                group_colors,
                selected_genes,
                window,
                max_branches,
            )
        end
    end

    for axis in axes, mark in instants
        text = join(["time: $(round(mark.t; digits=2))"; mark.entries], "\n")
        vlines!(
            axis,
            [mark.t];
            color=(:gray, 0.6),
            linestyle=:dash,
            linewidth=1,
            inspector_label=(_, _, _) -> text,
        )
    end

    length(axes) > 1 && linkxaxes!(axes...)

    DataInspector(figure; fontsize=12, show_bbox_indicators=false)
    figure
end
