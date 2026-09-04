transform(::CountSeries, y) = log1p(y)
transform(::FractionSeries, y) = y

function sample(catenation::Catenation, dimension::Dimension, t)
    series = get(catenation.trajectories, dimension, nothing)
    isnothing(series) ? NaN : transform(series, series(t))
end

function snapshots(trace::Trace{<:AbstractVector}, genes, track; resolution=200)
    track = Symbol(track)
    dimensions = [Dimension(track, string(gene)) for gene in genes]
    (; records, catenations) = trace

    values = Float64[]
    column = Vector{Float64}(undef, length(dimensions))
    states = @NamedTuple{catenation::Int, t::Float64}[]
    branches = Dict{Int, String}()
    lineages = Dict{Int, String}()
    colors = Dict{Int, String}()

    for (id, catenation) in enumerate(catenations)
        from = records[first(catenation.segments)].from
        to = records[last(catenation.segments)].to
        from < to || continue
        branches[id] = describe(records, catenation)
        lineages[id] = lineage(records, catenation)
        colors[id] = catenation_color(records, catenation)

        for t in range(from, to; length=resolution)
            map!(dimension -> sample(catenation, dimension, t), column, dimensions)
            all(isfinite, column) || continue
            append!(values, column)
            push!(states, (; catenation=id, t))
        end
    end

    X = reshape(values, length(dimensions), :)
    (; X, states, branches, lineages, colors, dimensions)
end

struct Projection
    coords::Matrix{Float64}
    states::Vector{@NamedTuple{catenation::Int, t::Float64}}
    colors::Vector{RGBf}
    branches::Dict{Int, String}
    lineages::Dict{Int, String}
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
    predict(fit(PCA, Z; maxoutdim=components), Z)
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
    scheme = to_colormap(:viridis)
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

function project(
    snapshot;
    components=2,
    group_colors=Dict(),
    temperature=nothing,
    coloring=:genes,
)
    coloring = Symbol(coloring)
    colors = palette(Val(coloring), snapshot, group_colors, temperature)
    swatches = coloring === :time ? Dict{Int, Nothing}() : Dict(
        id => try
            isempty(declared) ? nothing : RGBf(Colors.RGB(to_color(declared)))
        catch
            nothing
        end
        for (id, declared) in snapshot.colors
    )
    for (i, state) in enumerate(snapshot.states)
        declared = get(swatches, state.catenation, nothing)
        isnothing(declared) || (colors[i] = declared)
    end
    chosen = method(snapshot.X, components)
    coords = project(chosen, snapshot.X, components)
    Projection(
        coords,
        snapshot.states,
        colors,
        snapshot.branches,
        snapshot.lineages,
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

function polyline(projection::Projection, n; budget)
    grouped = Dict{Int, Vector{Int}}()
    for (i, state) in enumerate(projection.states)
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

function render(projection::Projection; emphasis=3, dimming=0.4, budget=20_000)
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

    (; points, colors, sources, spans, blank) = polyline(projection, n; budget)

    line = lines!(
        axis,
        points;
        color=colors,
        alpha=dimming,
    )

    highlight_points = Observable(typeof(blank)[])
    highlight_colors = Observable(RGBf[])
    lines!(
        axis,
        highlight_points;
        color=highlight_colors,
        linewidth=emphasis,
        depth_shift=-0.05f0,
    )

    overlay = Scene(figure.scene; clear=false)
    campixel!(overlay)
    tip_position = Observable(Point2f(0, 0))
    tip_text = Observable("")
    tip_visible = Observable(false)
    tooltip!(
        overlay,
        tip_position,
        tip_text;
        visible=tip_visible,
        fontsize=12,
        overdraw=true,
        depth_shift=-0.1f0,
    )

    hovered = Ref{Union{Nothing, String}}(nothing)

    on(events(axis.scene).mouseposition) do _
        picked, i = pick(axis.scene)
        hit = picked === line &&
            checkbounds(Bool, points, i) && !iszero(sources[i])
        state = hit ? projection.states[sources[i]] : nothing

        if hit
            tip_position[] = Point2f(events(axis.scene).mouseposition[])
            tip_text[] = join([
                "branch: $(get(projection.branches, state.catenation, ""))",
                "time: $(round(state.t; digits=2))",
            ], "\n")
        end
        tip_visible[] == hit || (tip_visible[] = hit)

        current = hit ? get(projection.lineages, state.catenation, "") : nothing
        hovered[] == current && return
        hovered[] = current

        highlight = typeof(blank)[]
        shades = RGBf[]

        if !isnothing(current)
            for (key, range) in spans
                related(current, key) || continue
                if !isempty(highlight)
                    push!(highlight, blank)
                    push!(shades, RGBf(0, 0, 0))
                end
                append!(highlight, @view points[range])
                append!(shades, @view colors[range])
            end
        end

        highlight_points[] = highlight
        highlight_colors[] = shades
    end

    figure
end
