module Visualization

using ..Common: Dimension

import Colors: Colors, Color, @colorant_str
using DataFrames
using Makie
using GeneRegulatorySystems: Models, Scheduling, NetworkRepresentation
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
    network=NetworkRepresentation.Network(provenance),
    _...,
)
    grid = GridLayout()

    actual = network
    sources = provenance_chain(provenance.source)
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
    descriptions::Models.Descriptions;
    group_colors,
    network=NetworkRepresentation.Network(descriptions),
    _...,
)
    label = findfirst(x -> x isa(Models.Label), descriptions.descriptions)
    title = label === nothing ? "" : descriptions.descriptions[label].label
    attach_model!(figure, network; group_colors, title)
end

function attach_model!(
    figure,
    description::Union{
        Models.RegulatoryNetwork,
        Models.ReactionNetwork
    };
    group_colors,
    network=NetworkRepresentation.Network(description),
    _...,
)
    attach_model!(figure, network; group_colors)
end

function attach_model!(
    figure,
    network::NetworkRepresentation.Network;
    group_colors,
    title="",
    _...,
)
    gene_axis = Axis(figure; autolimitaspect=1, title)
    species_axis = Axis(figure; autolimitaspect=1, title)
    panel = GridLayout()
    panel[1,1] = gene_axis
    panel[1,1] = species_axis
    axes = Dict(:gene_view => gene_axis, :species_view => species_axis)

    networks = Dict(
        :gene_view => NetworkRepresentation.gene_view(network),
        :species_view => NetworkRepresentation.species_view(network),
    )

    for (view, axis) in axes
        plot = attach_network!(axis, networks[view], group_colors; full_network=network)
        hidespines!(axis)
        hidedecorations!(axis)
        deregister_interaction!(axis, :rectanglezoom)
        register_interaction!(axis, :nodedrag, GraphMakie.NodeDrag(plot))
    end

    selector = Menu(
        figure;
        options=[
            ("gene view", :gene_view),
            ("species view", :species_view),
        ],
        default="gene view",
        width=120,
    )

    controls = GridLayout(panel[1, 1]; tellwidth=false, tellheight=false, halign=:right, valign=:top)
    controls[1, 1] = selector

    function show_view(view)
        view === nothing && return
        for (name, axis) in axes
            name === view ? Makie.unhide!(axis) : Makie.hide!(axis)
        end
        autolimits!(axes[view])
    end

    on(show_view, selector.selection)
    show_view(:gene_view)

    panel
end

function attach_network!(
    axis,
    network::NetworkRepresentation.Network,
    group_colors;
    full_network=network
)
    nodes = network.nodes
    node_index = Dict(node.name => i for (i, node) in enumerate(nodes))
    length(node_index) == length(nodes) || error("network node names must be unique")

    for link in network.links
        haskey(node_index, link.from) || error("missing network node $(link.from)")
        haskey(node_index, link.to) || error("missing network node $(link.to)")
    end


    format_property(value, ::Val) = string(value)
    format_property(value::Real, ::Val) = @sprintf("%.2g", value)
    format_property(value::Real, ::Val{:stoichiometry}) = string(Int(value))
    format_property(value, key::Symbol) = format_property(value, Val(key))

    function properties_label(properties)
        entries = [key => value for (key, value) in properties if key !== :parameters]
        isempty(entries) && return ""
        if length(entries) == 1
            key, value = only(entries)
            return format_property(value, key)
        end
        join(("$key=$(format_property(value, key))" for (key, value) in entries), " ")
    end

    node_label(node) = node_label(Val(node.kind), node)
    node_label(::Val, node) = string(node.name)
    node_label(::Val{:reaction}, node) = string(get(node.properties, :kind, node.name))

    parent_color(::Nothing, _colors, fallback) = fallback
    parent_color(parent::Symbol, colors, _fallback) = colors[parent]

    node_color(node) = node_color(Val(node.kind), node)
    node_color(::Val, _node) = :white
    node_color(::Val{:gene}, node) = group_colors[node.name]
    node_color(::Val{:species}, node) = parent_color(node.parent, group_colors, :white)
    node_color(::Val{:reaction}, _node) = :black

    node_styles = Dict(
        :gene     => (; marker=:circle,  size=0.70, fontsize=0.4),
        :species  => (; marker=:circle,  size=0.25, fontsize=0.25),
        :reaction => (; marker=:diamond, size=0.05, fontsize=0.2),
    )
    node_style(node) = get(node_styles, node.kind, node_styles[:species])
    styles = node_style.(nodes)

    function reaction_side(network, reaction, kind)
        links = filter(network.links) do link
            link.kind === kind &&
                (kind === :substrate ? link.to : link.from) === reaction.name
        end
        isempty(links) && return "∅"
        join((
            begin
                species = kind === :substrate ? link.from : link.to
                stoichiometry = get(link.properties, :stoichiometry, 1)
                stoichiometry == 1 ? string(species) : "$stoichiometry $species"
            end
            for link in links
        ), " + ")
    end

    function parameter_lines(item, network)
        associations = get(item.properties, :parameters, Dict())
        paths = isempty(item.present_in) ? keys(network.parameters) : item.present_in
        lines = String[]
        for (label, parameter) in associations
            values = unique(
                network.parameters[path][parameter]
                for path in paths
                if haskey(network.parameters, path) &&
                    haskey(network.parameters[path], parameter)
            )
            isempty(values) && continue
            push!(lines, "$label = $(join(format_property.(values, Ref(label)), ", "))")
        end
        lines
    end

    node_tooltip(node, network) = node_tooltip(Val(node.kind), node, network)
    node_tooltip(::Val, node, _network) = string(node.name)
    node_tooltip(::Val{:gene}, node, _network) = "gene $(node.name)"
    function node_tooltip(::Val{:reaction}, node, network)
        kind = get(node.properties, :kind, :reaction)
        heading = node.parent === nothing ? string(kind) : "$kind on $(node.parent)"
        arrow = iszero(get(node.properties, :k⁻, 0)) ? "->" : "⇌"
        equation = join((
            reaction_side(network, node, :substrate),
            reaction_side(network, node, :product),
        ), " $arrow ")
        join((heading, equation, parameter_lines(node, network)...), "\n")
    end

    links_by_edge = Dict{Pair{Int, Int}, Vector{NetworkRepresentation.Link}}()
    for link in network.links
        edge = node_index[link.from] => node_index[link.to]
        push!(get!(links_by_edge, edge, NetworkRepresentation.Link[]), link)
    end

    edge_styles = Dict(
        :activation  => (; color=:darkgreen,  linestyle=:solid, linewidth=3.0, marker=:rtriangle, fontsize=0.15),
        :repression  => (; color=:darkred,    linestyle=:solid, linewidth=3.0, marker=:rect, fontsize=0.15),
        :proteolysis => (; color=:darkred,    linestyle=:solid, linewidth=2.5, marker=:diamond, fontsize=0.15),
        :promotes    => (; color=:green4,     linestyle=:dash,  linewidth=1.3, marker=:rtriangle, fontsize=0.1),
        :inhibits    => (; color=:darkred,    linestyle=:dash,  linewidth=1.3, marker=:rect, fontsize=0.1),
        :catalyses   => (; color=:darkorange, linestyle=:dash,  linewidth=1.3, marker=:circle, fontsize=0.1),
        :substrate   => (; color=:gray,       linestyle=:solid, linewidth=1.3, marker=:rtriangle, fontsize=0.1),
        :product     => (; color=:gray,       linestyle=:solid, linewidth=1.3, marker=:rtriangle, fontsize=0.1),
        :multiple    => (; color=:black,      linestyle=:solid, linewidth=5.0, marker=:rtriangle, fontsize=0.3),
    )


    edge_properties = Dict(
        edge => if length(links) == 1
            link = only(links)
            (;
                label = properties_label(link.properties),
                get(edge_styles, link.kind, edge_styles[:substrate])...
            )
        else
            (;
                label=join(("$(link.kind): $(properties_label(link.properties))" for link in links), "\n"),
                edge_styles[:multiple]...
            )
        end
        for (edge, links) in links_by_edge
    )

    edge_tooltip(links) = join(unique(string(link.kind) for link in links), "\n")



    graph = Graphs.DiGraph(length(nodes))
    Graphs.add_edge!.(Ref(graph), keys(edge_properties))
    edges = collect(Graphs.edges(graph))
    edge_property(edge) = edge_properties[Graphs.src(edge) => Graphs.dst(edge)]
    edge_colors = [edge_property(edge).color for edge in edges]
    edge_linestyles = [edge_property(edge).linestyle for edge in edges]
    edge_widths = [edge_property(edge).linewidth for edge in edges]
    edge_labels = [edge_property(edge).label for edge in edges]
    edge_markers = [edge_property(edge).marker for edge in edges]
    edge_markersizes = [
        marker === :rect ? Vec2f(0.08, 0.2) : Vec2f(0.2, 0.2)
        for marker in edge_markers
    ]
    edge_fontsizes = [edge_property(edge).fontsize for edge in edges]

    edge_attributes = isempty(edges) ? (;) : (
        elabels=edge_labels,
        elabels_color=edge_colors,
        elabels_distance=0.2,
        elabels_fontsize=edge_fontsizes,
        elabels_attr=(markerspace=:data,),
        edge_plottype=:linesegments,
        edge_attr=(color=edge_colors, linestyle=edge_linestyles, linewidth=edge_widths),
        arrow_attr=(
            markerspace=:data,
            color=edge_colors,
            marker=edge_markers,
            markersize=edge_markersizes,
        )
    )

    function arrow_shifts(positions)
        map(edges) do edge
            source = Graphs.src(edge)
            target = Graphs.dst(edge)
            marker = edge_property(edge).marker
            source == target && return marker === :rtriangle ? 0.94 : 0.9
            distance = sqrt(sum(abs2, positions[target] - positions[source]))
            radius = styles[target].size / 2
            extent = marker === :rtriangle ? 0.1 : marker === :rect ? 0.02 : 0.15
            clamp(1 - (radius + extent) / distance, 0.0, 1.0)
        end
    end

    plot = GraphMakie.graphplot!(
        axis, graph;
        arrow_shift=0.92,
        nlabels=node_label.(nodes),
        nlabels_align=(:center, :bottom),
        nlabels_distance=0.2,
        nlabels_fontsize=getproperty.(styles, :fontsize),
        nlabels_color=:black,
        nlabels_attr=(markerspace=:data,),
        node_size=getproperty.(styles, :size),
        node_color=node_color.(nodes), node_strokewidth=2,
        node_attr=(
            markerspace=:data,
            marker=getproperty.(styles, :marker),
            strokecolor=:black
        ),
        edge_attributes...
    )
    update_arrows(positions) = plot[:arrow_shift][] = arrow_shifts(positions)
    on(update_arrows, plot[:node_pos])
    update_arrows(plot[:node_pos][])
    node_tip = Makie.tooltip!(axis, Observable(Point2f(0)); text=Observable(""), visible=false, fontsize=10, backgroundcolor=:white, transparency=false, overdraw=true, textpadding=(3,3,2,2), outline_linewidth=1)
    edge_tip = Makie.tooltip!(axis, Observable(Point2f(0)); text=Observable(""), visible=false, fontsize=10, backgroundcolor=:white, transparency=false, overdraw=true, textpadding=(3,3,2,2), outline_linewidth=1)
    for tip in (node_tip, edge_tip)
        tip.plots[1].draw_on_top[] = true
        tip.plots[1].alpha[] = 1.0
    end
    register_interaction!(axis, :node_tooltip,
        GraphMakie.NodeHoverHandler() do state, index, _event, _axis
            if state
                node_tip[1][] = plot[:node_pos][][index]
                node_tip.text[] = node_tooltip(nodes[index], full_network)
            end
            node_tip.visible[] = state
        end
    )
    register_interaction!(axis, :edge_tooltip,
        GraphMakie.EdgeHoverHandler() do state, index, event, _axis
            if state
                edge = edges[index]
                links = links_by_edge[
                    Graphs.src(edge) => Graphs.dst(edge)
                ]
                edge_tip[1][] = event.data
                edge_tip.text[] = edge_tooltip(links)
            end
            edge_tip.visible[] = state
        end
    )
    plot
end

end
