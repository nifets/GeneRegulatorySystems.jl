transform(::CountSeries, y) = log1p(y)
transform(::FractionSeries, y) = y

function sample(catenation::Catenation, dimension::Dimension, t)
    series = get(catenation.trajectories, dimension, nothing)
    isnothing(series) ? NaN : transform(series, series(t))
end

function snapshots(sink::Sink, genes, track; path="", resolution=200)
    track = Symbol(track)
    dimensions = [Dimension(track, string(gene)) for gene in genes]
    (; index, catenations) = catenate(sink; path)

    columns = Vector{Float64}[]
    states = @NamedTuple{catenation::Int, t::Float64}[]
    branches = Dict{Int, String}()

    for (id, catenation) in enumerate(catenations)
        from = index[first(catenation.segments)].from
        to = index[last(catenation.segments)].to
        from < to || continue
        branches[id] = describe(index, catenation)

        for t in range(from, to; length=resolution)
            column = [sample(catenation, dimension, t) for dimension in dimensions]
            all(isfinite, column) || continue
            push!(columns, column)
            push!(states, (; catenation=id, t))
        end
    end

    X = reduce(hcat, columns; init=zeros(length(dimensions), 0))
    (; X, states, branches, dimensions)
end

struct Projection
    coords::Matrix{Float64}
    states::Vector{@NamedTuple{catenation::Int, t::Float64}}
    colors::Vector{RGBf}
    branches::Dict{Int, String}
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

function mix(values, colors, temperature)
    weights = exp.((values .- maximum(values)) ./ temperature)
    weights ./= sum(weights)
    RGBf(
        sum(w * Colors.red(c) for (w, c) in zip(weights, colors)),
        sum(w * Colors.green(c) for (w, c) in zip(weights, colors)),
        sum(w * Colors.blue(c) for (w, c) in zip(weights, colors)),
    )
end

function blend(X, dimensions, group_colors; temperature=nothing)
    colors = [
        Colors.RGB(to_color(get(group_colors, dimension.group, :gray)))
        for dimension in dimensions
    ]
    scale = isnothing(temperature) ? std(X) : temperature
    isfinite(scale) && scale > 0 || (scale = 1.0)
    [mix(column, colors, scale) for column in eachcol(X)]
end

function project(snapshot; components=2, group_colors=Dict(), temperature=nothing)
    chosen = method(snapshot.X, components)
    coords = project(chosen, snapshot.X, components)
    Projection(
        coords,
        snapshot.states,
        blend(snapshot.X, snapshot.dimensions, group_colors; temperature),
        snapshot.branches,
        unval(chosen),
        labels(chosen, snapshot.dimensions, size(coords, 1)),
    )
end

function phase_axis(figure, projection, n)
    n == 3 && return Axis3(
        figure[1, 1];
        xlabel=projection.labels[1],
        ylabel=projection.labels[2],
        zlabel=projection.labels[3],
        backgroundcolor=:transparent,
    )
    Axis(
        figure[1, 1];
        xlabel=projection.labels[1],
        ylabel=projection.labels[2],
        backgroundcolor=:transparent,
    )
end

function render(projection::Projection; dimmed=0.12)
    n = min(size(projection.coords, 1), 3)
    figure = Figure()
    n >= 2 && !isempty(projection.states) || return figure

    axis = phase_axis(figure, projection, n)
    point = n == 3 ? Point3f : Point2f
    hovered = Observable{Union{Nothing, Int}}(nothing)
    owners = IdDict{Any, Int}()

    for id in sort!(unique(state.catenation for state in projection.states))
        mask = findall(state -> state.catenation == id, projection.states)
        positions = [point(view(projection.coords, 1:n, j)...) for j in mask]
        colors = projection.colors[mask]
        shaded = map(hovered) do current
            alpha = isnothing(current) || current == id ? 1.0 : dimmed
            [RGBAf(Colors.red(c), Colors.green(c), Colors.blue(c), alpha) for c in colors]
        end
        label = get(projection.branches, id, "")

        line = lines!(axis, positions; color=shaded, inspectable=false)
        dot = scatter!(
            axis,
            positions;
            color=shaded,
            markersize=4,
            inspector_label=(_, i, _) -> join([
                "branch: $label",
                "time: $(round(projection.states[mask[i]].t; digits=2))",
            ], "\n"),
        )
        owners[line] = id
        owners[dot] = id
    end

    on(events(axis.scene).mouseposition) do _
        picked, _ = pick(axis.scene)
        current = get(owners, picked, nothing)
        hovered[] == current || (hovered[] = current)
    end

    DataInspector(figure; fontsize=12, show_bbox_indicators=false)
    figure
end
