
@kwdef struct Node
    name::Symbol
    kind::Symbol
    parent::Union{Symbol, Nothing} = nothing
    properties::Dict{Symbol, Any} = Dict{Symbol, Any}()
    present_in::Set{String} = Set{String}()
end

@kwdef struct Link
    kind::Symbol
    from::Symbol
    to::Symbol
    scope::Symbol = :all
    properties::Dict{Symbol, Any} = Dict{Symbol, Any}()
    present_in::Set{String} = Set{String}()
end

@kwdef struct Network
    nodes::Vector{Node} = Node[]
    links::Vector{Link} = Link[]
    groups::Vector{Symbol} = Symbol[]
    parameters::Dict{String, Dict{Symbol, Float64}} = Dict{String, Dict{Symbol, Float64}}()
end

paths(network::Network) = sort!(collect(keys(network.parameters)))

Network(::Models.Description) = Network()
Network(descriptions::Models.Descriptions) =
    merge_networks(Network.(descriptions.descriptions)...)

function Network(provenance::Models.Provenance)
    function latest(type, description)
        description isa type && return Models.Description[description]
        if description isa Models.Descriptions
            return mapreduce(
                item -> latest(type, item),
                append!,
                description.descriptions;
                init=Models.Description[],
            )
        end
        description isa Models.Provenance || return Models.Description[]
        result = latest(type, description.description)
        isempty(result) ? latest(type, description.source) : result
    end

    merge_networks(
        Network.(latest(Models.RegulatoryNetwork, provenance))...,
        Network.(latest(Models.ReactionNetwork, provenance))...,
    )
end

function Network(network::Models.RegulatoryNetwork)
    nodes = vcat(
        [Node(name=g, kind=:gene) for g in network.species_groups],
        [Node(name=species,kind=:species, properties=Dict{Symbol,Any}(:shared=>true)) for species in network.shared_species]
    )
    gene_links = [Link(kind=link.kind, from=link.from, to=link.to, scope=:gene, properties=link.properties) for link in network.links]
    modulation_links = [Link(
        kind=link.modulation.kind,
        from=link.modulation.from,
        to=link.modulation.to,
        scope=:species,
        properties=merge(link.properties, Dict(:regulation => link.kind)),
    ) for link in network.links if hasproperty(link, :modulation)]
    Network(
        nodes=nodes,
        links=vcat(gene_links, modulation_links),
        groups = copy(network.species_groups)
    )
end

function Network(network::Models.ReactionNetwork)
    nodes = Node[Node(name=species, kind=:species) for species in Models.species(network)]
    links = Link[]
    for reaction in network.reactions
        parent = get(reaction.properties, :owner, nothing)
        scope = parent === nothing ? :all : :species
        properties = Dict{Symbol, Any}(
            :kind => get(reaction.properties, :kind, :reaction),
            :name => get(
                reaction.properties,
                :name,
                get(reaction.properties, :kind, reaction.name),
            ),
            :k⁺ => reaction.k⁺,
            :k⁻ => reaction.k⁻
        )
        for key in (:parameters, :gene_link)
            haskey(reaction.properties, key) && (properties[key] = reaction.properties[key])
        end
        push!(nodes, Node(
            name=reaction.name,
            kind=:reaction,
            parent=parent,
            properties=properties
        ))
        for (species, stoichiometry) in reaction.from.counts
            push!(links, Link(
                kind=:substrate,
                from=species,
                to=reaction.name,
                scope=scope,
                properties=Dict(:stoichiometry => stoichiometry)
            ))
        end
        for (species, stoichiometry) in reaction.to.counts
            push!(links, Link(
                kind=:product,
                from=reaction.name,
                to=species,
                scope=scope,
                properties=Dict(:stoichiometry => stoichiometry)
            ))
        end
    end
    Network(; nodes, links)
end

function Network(path::String, description::Models.Description;parameters=Dict{Symbol, Float64}())
    network = Network(description)
    Network(
        nodes = [Node(; name=node.name, kind=node.kind, parent=node.parent, properties=node.properties, present_in=Set([path])) for node in network.nodes],
        links = [Link(; kind=link.kind, from=link.from, to=link.to, scope=link.scope, properties=link.properties, present_in=Set([path])) for link in network.links],
        groups=network.groups,
        parameters=Dict(path => parameters)
    )
end

function Network(schedule!::Models.Scheduling.Schedule)
    networks = Dict{String, Network}()
    function collect_network!(primitive!, x, Δt; kwargs...)
        isfinite(Δt) && Δt > 0 || return nothing
        path = primitive!.path
        haskey(networks, path) && return nothing
        networks[path] = Network(
            path,
            Models.describe(primitive!);
            parameters=Dict{Symbol, Float64}(Models.parameters(primitive!))
        )
        nothing
    end
    schedule!(
        Models.FlatState();
        dryrun=collect_network!,
        parallel=false
    )
    merge_networks(values(networks)...)
end

function gene_of(species::Symbol, groups)
    index = findfirst(gene -> startswith(String(species), string(gene, ".")), groups)
    isnothing(index) ? nothing : groups[index]
end

function infer_parents(network::Network)
    parents = Dict{Symbol, Symbol}()
    for node in network.nodes
        node.kind === :species || continue
        parent = isnothing(node.parent) ?
            gene_of(node.name, network.groups) : node.parent
        isnothing(parent) || (parents[node.name] = parent)
    end

    adopted = Dict{Symbol, Set{Union{Symbol, Nothing}}}()
    for link in network.links
        link.kind in (:substrate, :product) || continue
        reaction, species = link.kind === :substrate ?
            (link.to, link.from) : (link.from, link.to)
        push!(
            get!(Set{Union{Symbol, Nothing}}, adopted, reaction),
            get(parents, species, nothing),
        )
    end

    nodes = map(network.nodes) do node
        parent = if node.kind === :species
            get(parents, node.name, nothing)
        elseif node.kind === :reaction && isnothing(node.parent)
            owners = get(adopted, node.name, ())
            length(owners) == 1 ? only(owners) : nothing
        else
            node.parent
        end
        Node(node.name, node.kind, parent, node.properties, node.present_in)
    end

    Network(
        nodes=nodes,
        links=network.links,
        groups=network.groups,
        parameters=network.parameters,
    )
end

function merge_networks(networks::Network...)
    isempty(networks) && return Network()
    nodes = Dict{Tuple{Symbol, Symbol}, Node}()
    links = Dict{Tuple{Symbol, Symbol, Symbol}, Link}()

    groups = unique(Symbol[group for n in networks for group in n.groups])

    for network in networks
        for node in network.nodes
            key = (node.kind, node.name)
            if haskey(nodes, key)
                prev = nodes[key]
                nodes[key] = Node(
                    name=node.name,
                    kind=node.kind,
                    parent=prev.parent === nothing ? node.parent : prev.parent,
                    properties = merge(prev.properties, node.properties),
                    present_in = union(prev.present_in, node.present_in)
                )
            else
                nodes[key] = node
            end
        end
        for link in network.links
            key = (link.kind, link.from, link.to)
            if haskey(links, key)
                prev = links[key]
                links[key] = Link(
                    kind=link.kind, from=link.from, to=link.to,
                    scope=prev.scope === link.scope ? link.scope : :all,
                    properties=merge(prev.properties, link.properties),
                    present_in=union(prev.present_in, link.present_in)
                )
            else
                links[key] = link
            end
        end
    end
    infer_parents(Network(
        nodes=collect(values(nodes)),
        links=collect(values(links)),
        groups=groups,
        parameters=merge((network.parameters for network in networks)...)
    ))
end

models(index, path::AbstractString="") = sort!(unique(
    segment.model
    for segment in index
    if Models.Scheduling.ispathprefix(path, segment.path)
))

models(network::Network, index, path::AbstractString="") =
    filter(in(keys(network.parameters)), models(index, path))

present_at(item, path) = isempty(item.present_in) || path in item.present_in

function path_view(network::Network, path::String)
    parameters = haskey(network.parameters, path) ?
        Dict(path => network.parameters[path]) :
        Dict{String, Dict{Symbol, Float64}}()

    Network(
        nodes=filter(node -> present_at(node, path), network.nodes),
        links=filter(link -> present_at(link, path), network.links),
        groups=network.groups,
        parameters=parameters,
    )
end

function species_view(network::Network; include_shared=false)
    nodes = [n for n in network.nodes if n.kind !== :gene && (include_shared || !get(n.properties, :shared, false))]
    visible = Set(node.name for node in nodes)
    Network(
        nodes=nodes,
        links=[link for link in network.links
            if link.scope !== :gene && link.from in visible && link.to in visible],
        groups=network.groups,
        parameters=network.parameters
    )
end

function gene_view(network::Network; include_shared=false)
    lookup = Dict(node.name => node for node in network.nodes)
    representative(name) = let node = lookup[name]
        something(node.parent, node.name)
    end
    nodes = [node for node in network.nodes if (node.kind === :gene || node.parent === nothing) && (include_shared || !get(node.properties, :shared, false))]
    visible = Set(node.name for node in nodes)
    links = Link[
        if link.scope === :gene
            link
        else
            let
                from = representative(link.from)
                to = representative(link.to)
                Link(kind=link.kind, from=from, to=to, scope=link.scope, present_in=link.present_in,
                properties = from === link.from && to === link.to ? link.properties : Dict{Symbol, Any}())
            end
        end
        for link in network.links if link.scope !== :species
    ]
    substrates = Dict{Symbol, Dict{Symbol, Int}}()
    products = Dict{Symbol, Dict{Symbol, Int}}()
    for link in network.links
        link.kind === :substrate &&
            (get!(substrates, link.to, Dict{Symbol, Int}())[link.from] = link.properties[:stoichiometry])
        link.kind === :product &&
            (get!(products, link.from, Dict{Symbol, Int}())[link.to] = link.properties[:stoichiometry])
    end

    for reaction in network.nodes
        reaction.kind === :reaction || continue
        reaction.parent === nothing && continue
        haskey(reaction.properties, :gene_link) && continue
        from = get(substrates, reaction.name, Dict{Symbol, Int}())
        to = get(products, reaction.name, Dict{Symbol, Int}())
        catalysts = Set(
            species for (species, stoichiometry) in from
            if get(to, species, 0) == stoichiometry
        )
        for catalyst in catalysts
            from = representative(catalyst)
            from === reaction.parent && continue
            push!(links, Link(
                kind=:affects,
                from=from,
                to=reaction.parent,
                scope=:gene,
                present_in=reaction.present_in,
            ))
        end
    end
    links = merge_networks(Network(links=links)).links

    links = [link for link in links if link.from in visible && link.to in visible]
    Network(nodes=nodes, links=links, groups=network.groups, parameters=network.parameters)
end


format_property(value, ::Val) = string(value)
format_property(value::Real, ::Val) = string(round(value; sigdigits=2))
format_property(value::Integer, ::Val{:stoichiometry}) = string(value)
format_property(value, key::Symbol) = format_property(value, Val(key))

function properties_label(properties)
    entries = [
        key => value
        for (key, value) in properties
        if key ∉ (:parameters, :gene_link, :regulation)
    ]

    isempty(entries) && return ""

    if length(entries) == 1
        key, value = only(entries)
        return format_property(value, key)
    end

    join(
        ("$key=$(format_property(value, key))" for (key, value) in entries),
        " ",
    )
end

node_label(node::Node) = node_label(Val(node.kind), node)
node_label(::Val, node) = string(node.name)
node_label(::Val{:reaction}, node) = string(node.properties[:name])

link_label(link::Link) = link_label(Val(link.kind), link)
link_label(::Val, link) = properties_label(link.properties)
link_label(::Val{:substrate}, link) = nothing
link_label(::Val{:product}, link) = nothing



function parameter_lines(item, network)
    associations = get(item.properties, :parameters, Dict())
    paths = isempty(item.present_in) ?
        keys(network.parameters) :
        item.present_in

    lines = String[]

    for (label, parameter) in associations
        values = unique(
            network.parameters[path][parameter]
            for path in paths
            if haskey(network.parameters, path) &&
               haskey(network.parameters[path], parameter)
        )

        isempty(values) && continue

        push!(
            lines,
            "$label = $(join(format_property.(values, Ref(label)), ", "))",
        )
    end

    lines
end

node_tooltip(node::Node, network::Network) =
    node_tooltip(Val(node.kind), node, network)

node_tooltip(::Val, node, network) = string(node.name)
node_tooltip(::Val{:gene}, node, network) = "gene $(node.name)"

function node_tooltip(::Val{:reaction}, node, network)
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
                stoichiometry == 1 ? string(species) : "[$stoichiometry]$species"
            end
            for link in links
        ), " + ")
    end
    label = node_label(node)
    heading = node.parent === nothing ?
        label :
        "$label on $(node.parent)"

    k⁻ = get(node.properties, :k⁻, 0)
    arrow = isequal(k⁻, zero(k⁻)) ? "→" : "⇌"

    equation = join((
        reaction_side(network, node, :substrate),
        reaction_side(network, node, :product),
    ), " $arrow ")

    join((heading, equation, parameter_lines(node, network)...), "\n")
end

link_tooltip(link::Link, network::Network) =
    link_tooltip(Val(link.kind), link, network)

link_tooltip(::Val{:substrate}, link, network) = nothing
link_tooltip(::Val{:product}, link, network) = nothing

function link_tooltip(::Val, link, network)
    species_level = link.kind in (:promotes, :inhibits)
    from = species_level ? link.from : something(gene_of(link.from, network.groups), link.from)
    to = species_level ? link.to : something(gene_of(link.to, network.groups), link.to)
    heading = "$(link.kind): $from → $to"
    parameters = ("  $line" for line in parameter_lines(link, network))
    join((heading, parameters...), "\n")
end

function node_variants(node::Node, network::Network)
    Dict(path => begin
        view = path_view(network, path)
        (;
            label=node_label(node),
            tooltip=node_tooltip(node, view),
            parameters=parameter_lines(node, view),
        )
    end for path in keys(network.parameters) if present_at(node, path))
end

function link_variants(link::Link, network::Network)
    Dict(path => begin
        view = path_view(network, path)
        parameters = parameter_lines(link, view)
        (;
            label=isempty(parameters) ? link_label(link) : join(parameters, " "),
            tooltip=link_tooltip(link, view),
        )
    end for path in keys(network.parameters) if present_at(link, path))
end
