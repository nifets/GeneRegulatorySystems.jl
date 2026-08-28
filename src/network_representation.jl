module NetworkRepresentation

using ..Models

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

Network(::Models.Description) = Network()
Network(descriptions::Models.Descriptions) =
    merge_networks(Network.(descriptions.descriptions)...)
Network(provenance::Models.Provenance) =
    merge_networks(Network(provenance.source), Network(provenance.description))

function Network(network::Models.RegulatoryNetwork)
    nodes = vcat(
        [Node(name=g, kind=:gene) for g in network.species_groups],
        [Node(name=species,kind=:species, properties=Dict{Symbol,Any}(:shared=>true)) for species in network.shared_species]
    )
    gene_links = [Link(kind=link.kind, from=link.from, to=link.to, scope=:gene, properties=link.properties) for link in network.links]
    modulation_links = [Link(kind=link.modulation.kind, from=link.modulation.from, to=link.modulation.to, scope=:species, properties=link.properties) for link in network.links if hasproperty(link, :modulation)]
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
            :kind => get(reaction.properties, :kind, :auxiliary),
            :k⁺ => reaction.k⁺,
            :k⁻ => reaction.k⁻
        )
        haskey(reaction.properties, :parameters) && (properties[:parameters] = reaction.properties[:parameters])
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

function merge_networks(networks::Network...)
    isempty(networks) && return Network()
    nodes = Dict{Tuple{Symbol, Symbol}, Node}()
    links = Dict{Tuple{Symbol, Symbol, Symbol}, Link}()

    groups = unique(Symbol[group for n in networks for group in n.groups])

    function gene_of(species::Symbol, groups)
        index = findfirst(gene -> startswith(String(species), string(gene, ".")), groups)
        isnothing(index) ? nothing : groups[index]
    end

    for network in networks
        for node in network.nodes
            parent = node.kind === :species && node.parent === nothing ?
                gene_of(node.name, groups) : node.parent
            node = Node(node.name, node.kind, parent, node.properties, node.present_in)
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
    Network(
        nodes=collect(values(nodes)),
        links=collect(values(links)),
        groups=groups,
        parameters=merge((network.parameters for network in networks)...)
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
        from = get(substrates, reaction.name, Dict{Symbol, Int}())
        to = get(products, reaction.name, Dict{Symbol, Int}())
        catalysts = Set(
            species for (species, stoichiometry) in from
            if get(to, species, 0) == stoichiometry
        )
        for catalyst in catalysts
            from = representative(catalyst)
            from === reaction.parent && continue
            push!(links, Link(; kind=:catalyses, from, to=reaction.parent, scope=:gene))
        end
    end
    links = merge_networks(Network(links=links)).links

    links = [link for link in links if link.from in visible && link.to in visible]
    Network(nodes=nodes, links=links, groups=network.groups, parameters=network.parameters)
end

end
