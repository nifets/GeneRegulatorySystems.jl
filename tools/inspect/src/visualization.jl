module Visualization

using ..Common: Dimension

import Colors: Colors, Color, @colorant_str
using DataFrames
using Makie
using GeneRegulatorySystems: Models, Scheduling
import Graphs
import GraphMakie

using Printf

struct GroupColors
    colors::Dict{String, Color}
end
GroupColors(::Nothing; _...) = GroupColors(Dict{String, Color}())
GroupColors(
    groups::AbstractVector{String};
    reserved = [colorant"white", colorant"black", colorant"crimson"],
    fixed = [
        # Light group: In Oklch color space, starting with colorant"crimson",
        # split the hue circle into 12 equidistant colors, increase luminosity
        # (by 0.2 on the logit scale), then pick 8 colors manually.
        colorant"#4196FF"
        colorant"#3CBC0B"
        colorant"#E38600"
        colorant"#D864F1"
        colorant"#A07CFF"
        colorant"#00C1D3"
        colorant"#A9A600"
        colorant"#FB53B1"

        # Dark group: Take the same 8 colors, shift their hue by -15° and
        # decrease luminosity (by 0.5 on the logit scale).
        colorant"#004D9D"
        colorant"#335D00"
        colorant"#8E2400"
        colorant"#6C00AB"
        colorant"#3D30AC"
        colorant"#006A67"
        colorant"#6B4900"
        colorant"#8D007C"
    ],
    seed = vcat(reserved, fixed),
    drop = length(reserved),
) = GroupColors(Dict(zip(
    groups,
    Colors.distinguishable_colors(length(groups) + drop, seed)[(drop + 1):end]
)))

Base.getindex(colors::GroupColors, group::Symbol) = colors[string(group)]
Base.getindex(colors::GroupColors, group::String) =
    get(colors.colors, group, colorant"gray")

kindname(kind::Symbol) = kindname(Val(kind))
kindname(::Val{Kind}) where {Kind} = String(Kind)
kindname(::Val{:activity}) = "promoter"
kindname(::Val{:mrnas}) = "mRNAs"
kindname(::Val{:premrnas}) = "pre-mRNAs"

abstract type Series end

@kwdef struct CountSeries <: Series
    ts::Vector{Float64} = Float64[]
    ys::Vector{Int} = Int[]
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
    x1 = x2 = x1′ = x2′ = 0
    next1 = next2 = (1, 1)
    sentinel = ((Inf, 0), (0, 0))
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
    front::Int
    back::Int
    series::Dict{Dimension, Series} = Dict{Dimension, Series}()
end

function attach_trajectory_label!(figure; kind, yscale)
    label = Label(
        figure,
        kindname(kind),
        rotation = π / 2,
        tellheight = false,
    )

    mouseevents = addmouseevents!(
        label.blockscene,
        label.layoutobservables.computedbbox,
    )
    onmouseleftdown(mouseevents) do _
        yscale[] = yscale[] == log10 ? identity : log10
    end

    label
end

function attach_trajectory_components!(
    figure,
    ::Type{CountSeries};
    index,
    catenations,
    group_colors,
    yscale,
)
    axis = Axis(
    	figure,
    	xticksvisible = false,
    	xticklabelsvisible = false;
    	yscale,
    )

    top = 0.0
    right = 1.0
    for catenation in values(catenations)
        to = index[catenation.back, :to]
        right = max(right, to)
        for (dimension, series) in catenation.series
            top = max(top, maximum(series.ys))
            color = group_colors[dimension.group]

            previous_i = index[catenation.front, :previous]
            if previous_i > 0 && haskey(catenations, previous_i)
                previous_t = index[previous_i, :to]
                previous_series = catenations[previous_i].series
                previous_y =
                if haskey(previous_series, dimension)
                    previous_y = last(previous_series[dimension].ys)
                    scatterlines!(
                        axis,
                        [previous_t, first(series.ts)],
                        [previous_y, first(series.ys)],
                        markersize = 2,
                        linewidth = 0.5,
                        linestyle = :dash;
                        color,
                    )
                end
            end

            if catenation.front == catenation.back
                stairs!(
                    axis,
                    series.ts,
                    series.ys,
                    step = :post,
                    linewidth = 1;
                    color,
                )

                if last(series.ts) < to
                    stairs!(
                        axis,
                        [last(series.ts), to],
                        [last(series.ys), last(series.ys)],
                        step = :post,
                        linewidth = 1;
                        color,
                    )
                end
            else
                scatterlines!(
                    axis,
                    series.ts,
                    series.ys,
                    markersize = 2,
                    linewidth = 0.5,
                    linestyle = :dash;
                    color,
                )
            end
        end
    end

    limits!(axis, 0.0, right, 0.5, top + 1.0)

    axis
end

function attach_trajectory_components!(
    figure,
    ::Type{FractionSeries};
    index,
    catenations,
    group_colors,
    yscale,
)
    axis = Axis(
        figure;
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        yreversed = true,
        tellheight = false,
    )

    right = 1.0
    for catenation in values(catenations)
        catenation.front == catenation.back || continue
        segment = index[catenation.back, :]
        segment.from < segment.to || continue
        right = max(right, segment.to)
        s = 1 / length(catenation.series)
        sortedseries = sort(collect(catenation.series), by = x -> x.first.group)
        for (j, (dimension, series)) in enumerate(sortedseries)
            if segment.previous > 0
                previous_t = index[segment.previous, :to]
                previous_y = index[segment.previous, :track]
                joint = scatterlines!(
                    axis,
                    [previous_t, segment.from],
                    [previous_y, segment.track + 1.0],
                    [-1.0, -1.0],
                    markersize = 2,
                    linewidth = 0.5,
                    linestyle = :dash,
                    color = colorant"black",
                )
                translate!(joint, 0, 0, 1)
            end
            y = segment.track + j * s
            ts = [repeat(series.ts, inner = 2)[2 : end]; segment.to]
            ys = repeat(series.ys, inner = 2)
            band!(
                axis,
                ts,
                y - 0.5s .- (0.5s .* ys),
                y - 0.5s .+ (0.5s .* ys),
                color = group_colors[dimension.group],
            )
        end
    end

    xlims!(axis, 0.0, right)

    axis
end

attach_trajectory_components!(figure; events, kind, rest...) =
    if haskey(events, kind)
        attach_trajectory_components!(
            figure,
            seriestype(kind);
            catenations = events[kind],
            rest...,
        )
    else
        Label(figure, "(no data)", tellheight = false, tellwidth = false)
    end

function attach_trajectory!(figure; index, events, kinds, group_colors)
    grid = GridLayout(tellheight = false)

    transform_pattern = r"""
        (?<transform>[[:word:]]+)? (?(transform)\(|)
            (?<kind>.+)
        (?(transform)\)|)
    """x
    for (i, kind) in enumerate(kinds)
        m = match(transform_pattern, String(kind))
        kind = Symbol(m[:kind])
        yscale = Observable{Function}(
            isnothing(m[:transform]) ? identity : log10
        )
        # ^ for now, every non-empty transform is interpreted to mean `log10`
        # TODO: clean this up

        grid[i, 1] = attach_trajectory_label!(figure; kind, yscale)
        grid[i, 2] = attach_trajectory_components!(
            figure;
            index,
            events,
            kind,
            group_colors,
            yscale,
        )
    end

    axes = [x for x in contents(grid[:, 2]) if x isa Axis]
    if !isempty(axes)
        bottom = last(axes)
        bottom.xlabel = L"model time $t$"
        bottom.xticksvisible = true
        bottom.xticklabelsvisible = true
        linkxaxes!(axes...)
    end

    grid
end

provenance_chain(description::Models.Description) = [description]
provenance_chain(provenance::Models.Provenance) = [
    provenance.description
    provenance_chain(provenance.source)
]

attach_model!(figure, ::Models.Description; tellheight = true, _...) =
    Label(figure, "(no description)", tellwidth = false; tellheight)

attach_model!(figure, label::Models.Label; tellheight = true, _...) =
    Label(figure, label.label, tellwidth = false; tellheight)

function attach_model!(
    figure,
    provenance::Models.Provenance;
    group_colors,
    _...,
)
    grid = GridLayout()

    actual, sources = Iterators.peel(provenance_chain(provenance))
    i = 1
    grid[i, 1] = attach_model!(figure, actual; group_colors)
    for source in sources
        i += 1
        if source isa Models.Label
            grid[i, 1] = Label(
                figure,
                "constructed from: $(source.label)",
                fontsize = 10,
                tellwidth = false,
                tellheight = true,
            )
        else
            grid[i, 1] = Label(
                figure,
                "constructed from:",
                fontsize = 10,
                tellwidth = false,
                tellheight = true,
            )
            i += 1
            grid[i, 1] = attach_model!(figure, source; group_colors)
        end
    end

    grid
end

function attach_model!(
    figure,
    composite::Models.Descriptions;
    group_colors,
    _...,
)
    components = Dict(typeof(d) => d for d in composite.descriptions)
    attach_model!(
        figure,
        components[Models.Network],
        get(Models.Label, components, Models.Label),
        get(Models.ReactionNetwork, components, Models.ReactionNetwork);
        group_colors,
    )
end

function attach_model!(
    figure,
    network::Models.Network,
    label::Models.Label,
    reactions::Models.ReactionNetwork;
    group_colors,
    _...,
)
    axis = Axis(figure, autolimitaspect = 1, title = label.label)

    node_style = (color = :gray, strokecolor = :gray)
    edge_styles = Dict(
        :activation => (color = :black, linestyle = :solid),
        :repression => (color = :red, linestyle = :solid),
        :proteolysis => (color = :red, linestyle = :dash),
        :multiple => (color = :darkred, linestyle = :dash),
        :massaction => (color = :gray, linestyle = :solid),
    )

    ungrouped = Symbol[]
    for reaction in reactions.reactions
        reagents = Iterators.flatten((reaction.from.counts, reaction.to.counts))
        for (species, _count) in reagents
            union!(ungrouped, (get(network.aliases, species, species),))
        end
    end
    setdiff!(ungrouped, network.species_groups)

    nodes = [
        map(network.species_groups) do g
            properties = (;
                node_style...,
                label = "",
                color = group_colors[g],
                strokecolor = :black,
                size = 16,
            )
            (:species, g, properties)
        end

        map(ungrouped) do s
            properties = (;
                node_style...,
                label = "",
                color = :white,
                size = 8,
            )
            (:species, s, properties)
        end

        [
            (
                :forward_reaction,
                i,
                (; node_style..., label = @sprintf("%.2g", r.k₊), size = 3)
            )
            for (i, r) in enumerate(reactions.reactions)
            if r.k₊ > 0.0
        ]

        [
            (
                :backward_reaction,
                i,
                (; node_style..., label = @sprintf("%.2g", r.k₋), size = 3)
            )
            for (i, r) in enumerate(reactions.reactions)
            if r.k₋ > 0.0
        ]
    ]

    nodes_index = Dict(
        (kind, node) => i
        for (i, (kind, node, _properties)) in enumerate(nodes)
    )
    for (alias, species) in network.aliases
        nodes_index[(:species, alias)] = nodes_index[(:species, species)]
    end

    links_by_edge = reduce(network.links, init = Dict()) do by_edge, link
        source = get(nodes_index, (:species, link.from), nothing)
        target = get(nodes_index, (:species, link.to), nothing)
        if source !== nothing && target !== nothing
            push!(get!(Vector, by_edge, source => target), link)
        end
        by_edge
    end

    edges = Dict()
    for (edge, links) in links_by_edge
        edges[edge] =
            if length(links) == 1
                link = only(links)
                value = link.properties[link.kind == :proteolysis ? :k : :at]
                (; label = @sprintf("%.2g", value), edge_styles[link.kind]...)
            else
                values = ["", "", ""]
                for link in links
                    if link.kind == :activation
                        values[1] = @sprintf("%.2g", link.properties[:at])
                    elseif link.kind == :repression
                        values[2] = @sprintf("%.2g", link.properties[:at])
                    elseif link.kind == :proteolysis
                        values[3] = @sprintf("%.2g", link.properties[:k])
                    end
                end
                (; label = join(values, "/"), edge_styles[:multiple]...)
            end
    end

    for (i, reaction) in enumerate(reactions.reactions)
        forward_node = get(nodes_index, (:forward_reaction, i), nothing)
        if forward_node !== nothing
            for (source, count) in reaction.from.counts
                source_node = nodes_index[(:species, source)]
                edges[source_node => forward_node] =
                    (; label = string(count), edge_styles[:massaction]...)
            end
            for (target, count) in reaction.to.counts
                target_node = nodes_index[(:species, target)]
                edges[forward_node => target_node] =
                    (; label = string(count), edge_styles[:massaction]...)
            end
        end

        backward_node = get(nodes_index, (:backward_reaction, i), nothing)
        if backward_node !== nothing
            for (source, count) in reaction.to.counts
                source_node = nodes_index[(:species, source)]
                edges[source_node => backward_node] =
                    (; label = string(count), edge_styles[:massaction]...)
            end
            for (target, count) in reaction.from.counts
                target_node = nodes_index[(:species, target)]
                edges[backward_node => target_node] =
                    (; label = string(count), edge_styles[:massaction]...)
            end
        end
    end

    graph = Graphs.DiGraph(length(nodes))
    Graphs.add_edge!.(Ref(graph), keys(edges))

    node_properties(property) = map(nodes) do (_, _, node)
        getproperty(node, property)
    end

    edge_properties(property) = map(Graphs.edges(graph)) do edge
        getproperty(edges[Graphs.src(edge) => Graphs.dst(edge)], property)
    end

    edge_attributes = Graphs.ne(graph) == 0 ? (;) : (
        elabels = edge_properties(:label),
        elabels_color = edge_properties(:color),
        elabels_distance = 12,
        elabels_fontsize = 8,
        edge_plottype = :beziersegments,
        edge_attr = (
            linestyle = edge_properties(:linestyle),
            color = edge_properties(:color),
        ),
        arrow_attr = (
            markersize = 16,
            color = edge_properties(:color),
        ),
    )

    GraphMakie.graphplot!(
        axis,
        graph,
        nlabels = node_properties(:label),
        nlabels_fontsize = 8,
        nlabels_distance = 3,
        nlabels_color = :gray,
        node_size = node_properties(:size),
        node_strokewidth = 2,
        node_attr = (; strokecolor = node_properties(:strokecolor),),
        node_color = node_properties(:color);
        edge_attributes...
    )

    hidespines!(axis)
    hidedecorations!(axis)
    deregister_interaction!(axis, :rectanglezoom)
    register_interaction!(
        axis,
        :nodedrag,
        GraphMakie.NodeDrag(only(axis.scene.plots))
    )

    axis
end

end
