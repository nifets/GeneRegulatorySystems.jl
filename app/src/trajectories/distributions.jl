struct Samples
    X::Matrix{Float64}
    ids::Vector{Int}
    colors::Dict{Int, String}
    dimensions::Vector{Dimension}
end

function counts(
    trace::Trace{<:AbstractVector}, genes, track; resolution=200, budget=200_000,
)
    track = Symbol(track)
    dimensions = [Dimension(track, string(gene)) for gene in genes]
    (; records, catenations) = trace

    resolution = max(2, min(
        resolution,
        fld(budget, max(1, length(dimensions)) * max(1, length(catenations))),
    ))

    values = Float64[]
    ids = Int[]
    colors = Dict{Int, String}()
    column = Vector{Float64}(undef, length(dimensions))

    for (id, catenation) in enumerate(catenations)
        from = records[first(catenation.segments)].from
        to = records[last(catenation.segments)].to
        from < to || continue
        colors[id] = catenation_color(records, catenation)

        for t in range(from, to; length=resolution)
            for (i, dimension) in enumerate(dimensions)
                series = get(catenation.trajectories, dimension, nothing)
                column[i] = isnothing(series) ? NaN : series(t)
            end
            all(isfinite, column) || continue
            append!(values, column)
            push!(ids, id)
        end
    end

    Samples(reshape(values, length(dimensions), :), ids, colors, dimensions)
end

function path_color(samples::Samples)
    for id in sort!(collect(keys(samples.colors)))
        color = samples.colors[id]
        isempty(color) || return color
    end
    "#9ca3af"
end

gene_color(dimension::Dimension, group_colors) =
    to_color(get(group_colors, dimension.group, :gray))

moments(samples::Samples) =
    (vec(mean(samples.X; dims=2)), vec(var(samples.X; dims=2)))

function moment_axis(figure)
    reset_on_doubleclick!(Axis(
        figure[1, 1];
        xlabel="mean",
        ylabel="variance",
        xscale=log10,
        yscale=log10,
        backgroundcolor=:transparent,
    ))
end

function poisson!(axis, means)
    isempty(means) && return
    low, high = extrema(means)
    high > low || return
    lines!(axis, [low, high], [low, high]; color=:gray, linestyle=:dash)
end

function reset_on_doubleclick!(axis)
    register_interaction!(axis, :reset_view) do event, axis
        if event isa MouseEvent &&
                event.type == WGLMakie.Makie.MouseEventTypes.leftdoubleclick
            autolimits!(axis)
            return Consume(true)
        end
        Consume(false)
    end
    axis
end

function joint_axis(position, samples::Samples; limits=(nothing, nothing))
    reset_on_doubleclick!(Axis(
        position;
        xlabel=samples.dimensions[1].group,
        ylabel=samples.dimensions[2].group,
        backgroundcolor=:transparent,
        limits,
    ))
end

cellsize(((xlow, xhigh), (ylow, yhigh)), bins) =
    ((xhigh - xlow) / bins, (yhigh - ylow) / bins)

flat(values) = isempty(values) || (extrema(values) |> ((low, high),) -> high <= low)

function edges(values...; maxbins=40)
    low, high = extrema(Iterators.flatten(values))
    integral = all(v -> v == round(v), Iterators.flatten(values))
    integral && high - low <= maxbins ?
        range(low - 0.5, high + 0.5; length=Int(high - low) + 2) :
        range(low, high; length=maxbins + 1)
end

degenerate(samples::Samples) =
    any(row -> flat(view(samples.X, row, :)), 1:size(samples.X, 1))

function scatter_fallback!(axis, samples::Samples)
    scatter!(
        axis,
        view(samples.X, 1, :),
        view(samples.X, 2, :);
        color=path_color(samples),
        markersize=8,
    )
end

function joint_limits(a::Samples, b::Samples)
    x = extrema(vcat(view(a.X, 1, :), view(b.X, 1, :)))
    y = extrema(vcat(view(a.X, 2, :), view(b.X, 2, :)))
    (x, y)
end

render(samples::Samples; group_colors=Dict()) =
    render(Val(length(samples.dimensions)), samples; group_colors)

render(primary::Samples, secondary::Samples; group_colors=Dict()) =
    render(Val(length(primary.dimensions)), primary, secondary; group_colors)

render(::Val{0}, ::Samples; group_colors=nothing) = Figure()
render(::Val{0}, ::Samples, ::Samples; group_colors=nothing) = Figure()

function render(::Val{1}, samples::Samples; group_colors)
    figure = Figure()
    axis = reset_on_doubleclick!(Axis(
        figure[1, 1];
        xlabel=samples.dimensions[1].group,
        ylabel="density",
        backgroundcolor=:transparent,
    ))
    flat(vec(samples.X)) || hist!(
        axis,
        vec(samples.X);
        bins=edges(vec(samples.X)),
        normalization=:pdf,
        color=(gene_color(samples.dimensions[1], group_colors), 0.55),
        strokecolor=path_color(samples),
        strokewidth=1,
    )
    figure
end

function render(::Val{1}, primary::Samples, secondary::Samples; group_colors)
    figure = Figure()
    axis = reset_on_doubleclick!(Axis(
        figure[1, 1];
        xlabel=primary.dimensions[1].group,
        ylabel="density",
        backgroundcolor=:transparent,
    ))
    shared = edges(vec(primary.X), vec(secondary.X))
    flat(vec(primary.X)) || hist!(
        axis,
        vec(primary.X);
        bins=shared,
        normalization=:pdf,
        color=(gene_color(primary.dimensions[1], group_colors), 0.35),
    )
    for samples in (primary, secondary)
        flat(vec(samples.X)) && continue
        stephist!(
            axis,
            vec(samples.X);
            bins=shared,
            normalization=:pdf,
            color=path_color(samples),
            linewidth=2,
        )
    end
    figure
end

function render(::Val{2}, samples::Samples; group_colors=nothing)
    figure = Figure()
    axis = joint_axis(figure[1, 1], samples)
    if degenerate(samples)
        scatter_fallback!(axis, samples)
    else
        plot = hexbin!(axis, view(samples.X, 1, :), view(samples.X, 2, :); bins=40)
        Colorbar(figure[1, 2], plot)
    end
    figure
end

function render(::Val{2}, primary::Samples, secondary::Samples; group_colors=nothing)
    figure = Figure()
    limits = joint_limits(primary, secondary)
    size = cellsize(limits, 40)
    for (column, samples) in enumerate((primary, secondary))
        axis = joint_axis(figure[1, column], samples; limits)
        axis.title = column == 1 ? "primary" : "secondary"
        axis.titlecolor = path_color(samples)
        if degenerate(samples)
            scatter_fallback!(axis, samples)
        else
            hexbin!(
                axis,
                view(samples.X, 1, :),
                view(samples.X, 2, :);
                cellsize=size,
            )
        end
    end
    figure
end

label(dimension::Dimension, mean, variance) = join([
    dimension.group,
    "mean $(round(mean; digits=2))",
    "variance $(round(variance; digits=2))",
    "fano $(round(variance / mean; digits=2))",
], "\n")

function render(::Val, samples::Samples; group_colors)
    means, variances = moments(samples)
    keep = means .> 0 .&& variances .> 0
    shown = samples.dimensions[keep]
    m = means[keep]
    v = variances[keep]
    figure = Figure()
    axis = moment_axis(figure)
    scatter!(
        axis,
        m,
        v;
        color=[gene_color(d, group_colors) for d in shown],
        strokecolor=path_color(samples),
        strokewidth=1,
        markersize=8,
        inspector_label=(_, i, _) -> label(shown[i], m[i], v[i]),
    )
    poisson!(axis, m)
    DataInspector(figure; fontsize=12, show_bbox_indicators=false)
    figure
end

function render(::Val, primary::Samples, secondary::Samples; group_colors)
    primary_means, primary_variances = moments(primary)
    secondary_means, secondary_variances = moments(secondary)
    keep = primary_means .> 0 .&& primary_variances .> 0 .&&
        secondary_means .> 0 .&& secondary_variances .> 0

    segments = Point2f[]
    for i in findall(keep)
        push!(segments, Point2f(primary_means[i], primary_variances[i]))
        push!(segments, Point2f(secondary_means[i], secondary_variances[i]))
    end

    shown = primary.dimensions[keep]
    pm = primary_means[keep]
    pv = primary_variances[keep]
    sm = secondary_means[keep]
    sv = secondary_variances[keep]
    colors = [gene_color(d, group_colors) for d in shown]

    figure = Figure()
    axis = moment_axis(figure)
    linesegments!(axis, segments; color=:gray, linewidth=0.75)
    scatter!(
        axis,
        pm,
        pv;
        color=colors,
        marker=:circle,
        strokecolor=path_color(primary),
        strokewidth=1,
        markersize=8,
        inspector_label=(_, i, _) ->
            label(shown[i], pm[i], pv[i]),
    )
    scatter!(
        axis,
        sm,
        sv;
        color=colors,
        marker=:rect,
        strokecolor=path_color(secondary),
        strokewidth=1,
        markersize=8,
        inspector_label=(_, i, _) ->
            label(shown[i], sm[i], sv[i]),
    )
    poisson!(axis, vcat(pm, sm))
    DataInspector(figure; fontsize=12, show_bbox_indicators=false)
    figure
end
