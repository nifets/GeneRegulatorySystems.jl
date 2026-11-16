transform(::CountSeries, y) = log1p(y)
transform(::FractionSeries, y) = y

function sample(catenation::Catenation, dimension::Dimension, t)
    series = get(catenation, dimension)
    transform(series, series(t))
end

function snapshots(trace::Trace{<:AbstractVector}, genes, track; resolution=200)
    track = Symbol(track)
    dimensions = [Dimension(track, string(gene)) for gene in genes]
    (; records, catenations) = trace

    values = Float64[]
    column = Vector{Float64}(undef, length(dimensions))
    states = @NamedTuple{catenation::Int, t::Float64}[]
    descriptions = String[]
    declared = String[]
    lineages = Dict{Int, String}()
    sampled = Dict{Int, Bool}()

    for (id, catenation) in enumerate(catenations)
        from = records[first(catenation.segments)].from
        to = records[last(catenation.segments)].to
        snapshot = all(i -> records[i].from == records[i].to, catenation.segments)
        snapshot || from < to || continue
        describe_at = describer(records, catenation)
        at = record_at(records, catenation)
        lineages[id] = lineage(records, catenation)
        sampled[id] = snapshot

        times = snapshot ?
            unique(records[i].from for i in catenation.segments) :
            range(from, to; length=resolution)
        for t in times
            map!(dimension -> sample(catenation, dimension, t), column, dimensions)
            all(isfinite, column) || continue
            append!(values, column)
            push!(states, (; catenation=id, t))
            push!(descriptions, describe_at(t))
            push!(declared, at(t).color)
        end
    end

    X = reshape(values, length(dimensions), :)
    (; X, states, descriptions, declared, lineages, sampled, dimensions)
end

struct Projection
    coords::Matrix{Float64}
    states::Vector{@NamedTuple{catenation::Int, t::Float64}}
    colors::Vector{RGBf}
    descriptions::Vector{String}
    lineages::Dict{Int, String}
    sampled::Dict{Int, Bool}
    method::Symbol
    labels::Vector{String}
end

unval(::Val{S}) where {S} = S

method(X, components) =
    size(X, 1) <= components ? Val(:direct) :
    size(X, 1) <= 50         ? Val(:pca) : Val(:umap)

function standardise(X)
    Z = copy(X)
    for row in eachrow(Z)
        μ, σ = mean(row), std(row)
        row .= σ > 0 ? (row .- μ) ./ σ : row .- μ
    end
    Z
end

project(::Val{:direct}, X, components) = X

function project(::Val{:pca}, X, components)
    Z = standardise(X)
    predict(fit(PCA, Z; maxoutdim=components, method=:svd), Z)
end

project(::Val{:umap}, X, components) =
    UMAP.fit(
        project(Val(:pca), X, min(50, size(X, 1))),
        components;
        n_neighbors=min(15, size(X, 2) - 1),
    ).embedding

labels(::Val{:direct}, dimensions, _) =
    [dimension.group for dimension in dimensions]
labels(::Val{:pca}, _, n) = ["PC $i" for i in 1:n]
labels(::Val{:umap}, _, n) = ["UMAP $i" for i in 1:n]

function blend(X, dimensions, group_colors; temperature=nothing)
    colors = [
        Colors.RGB(to_color(get(group_colors, dimension.group, :gray)))
        for dimension in dimensions
    ]
    scale = isnothing(temperature) ? std(X) : temperature
    isfinite(scale) && scale > 0 || (scale = 1.0)
    weights = similar(X, size(X, 1))

    function mix(column)
        weights .= exp.((column .- maximum(column)) ./ scale)
        weights ./= sum(weights)
        RGBf(
            sum(w * Colors.red(c) for (w, c) in zip(weights, colors)),
            sum(w * Colors.green(c) for (w, c) in zip(weights, colors)),
            sum(w * Colors.blue(c) for (w, c) in zip(weights, colors)),
        )
    end

    [mix(column) for column in eachcol(X)]
end

palette(::Val{:genes}, snapshot, group_colors, temperature) =
    blend(snapshot.X, snapshot.dimensions, group_colors; temperature)

function palette(::Val{:time}, snapshot, _group_colors, _temperature)
    scheme = [
        RGBf(convert(Colors.RGB, Colors.HSL(h, 0.65, 0.58)))
        for h in range(140, 300; length=256)
    ]
    ts = [state.t for state in snapshot.states]
    low, high = extrema(ts)
    span = high > low ? high - low : 1.0
    [
        RGBf(scheme[clamp(
            round(Int, (t - low) / span * (length(scheme) - 1)) + 1,
            1,
            length(scheme),
        )])
        for t in ts
    ]
end

palette(::Val{:schedule}, snapshot, _group_colors, _temperature) = [
    isempty(color) ? RGBf(0.6, 0.6, 0.6) : RGBf(Colors.RGB(to_color(color)))
    for color in snapshot.declared
]

Projection(::Nothing, components) = Projection(
    zeros(Float64, max(components, 2), 0),
    @NamedTuple{catenation::Int, t::Float64}[],
    RGBf[],
    String[],
    Dict{Int, String}(),
    Dict{Int, Bool}(),
    Symbol(""),
    String[],
)

function project(
    snapshot;
    components=2,
    group_colors=Dict(),
    temperature=nothing,
    coloring=:genes,
)
    size(snapshot.X, 2) < 3 && return Projection(nothing, components)
    coloring = Symbol(coloring)
    colors = palette(Val(coloring), snapshot, group_colors, temperature)
    chosen = method(snapshot.X, components)
    coords = try
        project(chosen, snapshot.X, components)
    catch error
        @warn "phase projection failed" error method=unval(chosen) size=size(snapshot.X)
        return Projection(nothing, components)
    end
    Projection(
        coords,
        snapshot.states,
        colors,
        snapshot.descriptions,
        snapshot.lineages,
        snapshot.sampled,
        unval(chosen),
        labels(chosen, snapshot.dimensions, size(coords, 1)),
    )
end

axislabel(projection, i) = get(projection.labels, i, "")

function phase_axis(figure, projection, n)
    n == 3 && return Axis3(
        figure[1, 1];
        xlabel=axislabel(projection, 1),
        ylabel=axislabel(projection, 2),
        zlabel=axislabel(projection, 3),
        backgroundcolor=:transparent,
        viewmode=:free,
    )
    Axis(
        figure[1, 1];
        xlabel=axislabel(projection, 1),
        ylabel=axislabel(projection, 2),
        backgroundcolor=:transparent,
    )
end

related(a, b) =
    Scheduling.ispathprefix(a, b) || Scheduling.ispathprefix(b, a)

function polyline(projection::Projection, n; budget, sampled)
    grouped = Dict{Int, Vector{Int}}()
    for (i, state) in enumerate(projection.states)
        get(projection.sampled, state.catenation, false) == sampled || continue
        push!(get!(Vector{Int}, grouped, state.catenation), i)
    end
    stride = max(1, cld(length(projection.states), budget))
    coords = projection.coords
    at(j) = n == 3 ?
        Point3f(coords[1, j], coords[2, j], coords[3, j]) :
        Point2f(coords[1, j], coords[2, j])
    blank = n == 3 ? Point3f(NaN, NaN, NaN) : Point2f(NaN, NaN)

    points = typeof(blank)[]
    colors = RGBf[]
    sources = Int[]
    spans = Tuple{String, UnitRange{Int}}[]

    for id in sort!(collect(keys(grouped)))
        mask = grouped[id][begin:stride:end]
        isempty(mask) && continue
        if !isempty(points)
            push!(points, blank)
            push!(colors, RGBf(0, 0, 0))
            push!(sources, 0)
        end

        start = length(points) + 1
        for j in mask
            push!(points, at(j))
            push!(colors, projection.colors[j])
            push!(sources, j)
        end
        push!(spans, (get(projection.lineages, id, ""), start:length(points)))
    end

    (; points, colors, sources, spans, blank)
end

function render(
    projection::Projection;
    emphasis=3,
    thickness=2,
    budget=20_000,
    connect_samples=false,
    selection=Observable(""),
    hovered=Observable{Union{Nothing, String}}(nothing),
)
    n = min(size(projection.coords, 1), 3)
    figure = Figure()
    axis = phase_axis(figure, projection, max(n, 2))

    register_interaction!(axis, :reset_view) do event, axis
        if event isa MouseEvent &&
                event.type == WGLMakie.Makie.MouseEventTypes.leftdoubleclick
            autolimits!(axis)
            return Consume(true)
        end
        Consume(false)
    end

    n >= 2 && !isempty(projection.states) || return figure

    continuous = polyline(projection, n; budget, sampled=false)
    discrete = polyline(projection, n; budget, sampled=true)
    blank = continuous.blank

    highlight_points = Observable(typeof(blank)[])
    highlight_colors = Observable(RGBf[])
    highlighted = (; sources=Int[])

    function emphasise()
        current = something(hovered[], selection[])
        highlight = typeof(blank)[]
        shades = RGBf[]
        empty!(highlighted.sources)

        if !isnothing(hovered[]) || !isempty(current)
            for polyline in (continuous, discrete), (key, range) in polyline.spans
                related(current, key) || continue
                if !isempty(highlight)
                    push!(highlight, blank)
                    push!(shades, RGBf(0, 0, 0))
                    push!(highlighted.sources, 0)
                end
                append!(highlight, @view polyline.points[range])
                append!(shades, @view polyline.colors[range])
                append!(highlighted.sources, @view polyline.sources[range])
            end
        end

        highlight_points[] = highlight
        highlight_colors[] = shades
    end

    on(_ -> emphasise(), axis.scene, selection)
    on(_ -> emphasise(), axis.scene, hovered)

    selectable!(axis, hovered, selection)

    halo = contrast(WGLMakie.Makie.theme(:backgroundcolor)[])
    dimming = opacity(length(unique(values(projection.lineages))))

    hover_with(source; follow=true) = (inspector, plot, i) -> begin
        j = get(source.sources, i, 0)
        iszero(j) && return false
        current = get(projection.lineages, projection.states[j].catenation, "")
        isempty(selection[]) || related(selection[], current) || return false
        deeper = Scheduling.ispathprefix(something(hovered[], selection[]), current)
        (follow || deeper) && hovered[] != current && (hovered[] = current)
        WGLMakie.Makie.show_data(inspector, plot, i)
    end
    tooltip_for(source) = (_, i, _) -> begin
        j = source.sources[i]
        "branch: $(projection.descriptions[j])\ntime: $(round(projection.states[j].t; digits=2))"
    end
    inspection(source; follow=true) =
        (; inspector_hover=hover_with(source; follow), inspector_label=tooltip_for(source))

    lines!(
        axis,
        continuous.points;
        color=continuous.colors,
        alpha=dimming,
        linewidth=thickness,
        inspection(continuous)...,
    )
    connect_samples ?
        lines!(
            axis,
            discrete.points;
            color=discrete.colors,
            alpha=dimming,
            linewidth=thickness,
            inspection(discrete)...,
        ) :
        scatter!(
            axis,
            discrete.points;
            color=discrete.colors,
            alpha=min(2dimming, 0.9),
            markersize=2thickness,
            inspection(discrete)...,
        )
    lines!(
        axis,
        highlight_points;
        color=halo,
        linewidth=emphasis + 3,
        inspection(highlighted; follow=false)...,
    )
    lines!(
        axis,
        highlight_points;
        color=highlight_colors,
        linewidth=emphasis,
        inspection(highlighted; follow=false)...,
    )

    emphasise()
    inspector = DataInspector(figure; fontsize=12)
    on(inspector.plot.visible) do visible
        visible || isnothing(hovered[]) || (hovered[] = nothing)
    end
    figure
end
