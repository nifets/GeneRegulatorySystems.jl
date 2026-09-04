transform(::CountSeries, y) = log1p(y)
transform(::FractionSeries, y) = y

function sample(catenation::Catenation, dimension::Dimension, t)
    series = get(catenation.trajectories, dimension, nothing)
    isnothing(series) ? NaN : transform(series, series(t))
end

function snapshots(trace::Trace{<:AbstractVector}, genes, track; resolution=200)
    track = Symbol(track)
    dimensions = [Dimension(track, string(gene)) for gene in genes]
    (; index, catenations) = trace

    values = Float64[]
    column = Vector{Float64}(undef, length(dimensions))
    states = @NamedTuple{catenation::Int, t::Float64}[]
    branches = Dict{Int, String}()
    lineages = Dict{Int, String}()
    colors = Dict{Int, String}()

    for (id, catenation) in enumerate(catenations)
        from = index[first(catenation.segments)].from
        to = index[last(catenation.segments)].to
        from < to || continue
        branches[id] = describe(index, catenation)
        lineages[id] = lineage(index, catenation)
        colors[id] = catenation_color(index, catenation)

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

function project(snapshot; components=2, group_colors=Dict(), temperature=nothing)
    colors = blend(snapshot.X, snapshot.dimensions, group_colors; temperature)
    swatches = Dict(
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

function render(projection::Projection; emphasis=3, budget=20_000)
    n = min(size(projection.coords, 1), 3)
    figure = Figure()
    axis = phase_axis(figure, projection, max(n, 2))
    n >= 2 && !isempty(projection.states) || return figure

    grouped = Dict{Int, Vector{Int}}()
    for (i, state) in enumerate(projection.states)
        push!(get!(Vector{Int}, grouped, state.catenation), i)
    end
    stride = max(1, cld(length(projection.states), budget))
    coords = projection.coords
    at(j) = n == 3 ?
        Point3f(coords[1, j], coords[2, j], coords[3, j]) :
        Point2f(coords[1, j], coords[2, j])

    owners = IdDict{Any, String}()
    strands = Dict{String, Vector{Any}}()

    for id in sort!(collect(keys(grouped)))
        mask = grouped[id][begin:stride:end]
        points = [at(j) for j in mask]
        colors = projection.colors[mask]
        label = get(projection.branches, id, "")
        key = get(projection.lineages, id, "")

        annotate(_, i, _) = join([
            "branch: $label",
            "time: $(round(projection.states[mask[i]].t; digits=2))",
        ], "\n")

        line = lines!(axis, points; color=colors, inspector_label=annotate)
        dot = scatter!(
            axis,
            points;
            color=colors,
            markersize=4,
            inspector_label=annotate,
        )
        owners[line] = key
        owners[dot] = key
        push!(get!(Vector{Any}, strands, key), line)
    end

    hovered = Ref{Union{Nothing, String}}(nothing)
    active = Any[]

    on(events(axis.scene).mouseposition) do _
        picked, _ = pick(axis.scene)
        current = get(owners, picked, nothing)
        hovered[] == current && return
        hovered[] = current

        for line in active
            line.linewidth[] = 1.5
        end
        empty!(active)
        isnothing(current) && return

        for (key, lines) in strands
            related(current, key) || continue
            for line in lines
                line.linewidth[] = emphasis
                push!(active, line)
            end
        end
    end

    DataInspector(figure; fontsize=12, show_bbox_indicators=false)
    figure
end
