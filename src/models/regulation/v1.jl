"""
Contains components to build `SciML.JumpModel`s for gene regulation.

These models can be constructed by [`build`](@ref), either directly (from a
[`Definition`](@ref) or from the corresponding JSON specification), or
indirectly from a model template such as `KroneckerNetworks.Definition` or
`Differentiation.Definition` via their respective `build` functions.

The regulation models require setting up the initial state to include required
molecular species (polymerases, ribosomes and proteasomes) before they can be
invoked. If a model is specified as part of a
[`Schedule`](@ref Models.Scheduling.Schedule) and therefore the defaults are
available, this can be achieved by prefixing the model specification with an
`Instant` adjustment step such as:
```
{"{add}": {"\$": ["defaults", "bootstrap"]}}
```
See `examples/specification/simple.schedule.json` for an example.
"""
module V1

using ..Models: Models, SciML
import ...Specifications: constructor, cast, representation

using Base: @invoke

using Catalyst
using JumpProcesses
using ModelingToolkit
using StatsBase

"""
    ProkaryoteBaseRates

Defines a prokaryotic [`Gene`](@ref) reaction cascade's base rates.

# Specification

They are specified in JSON as a JSON object
```
{
    "activation": <...>,
    "deactivation": <...>,
    "trigger": <...>,
    "transcription": <...>,
    "translation": <...>,
    "abortion": <...>,
    "mrna_decay": <...>,
    "protein_decay": <...>
}
```
where `<...>` are JSON numbers setting the corresponding reaction rate
constants.
"""
@kwdef struct ProkaryoteBaseRates
    activation::Float64
    deactivation::Float64
    trigger::Float64
    transcription::Float64
    translation::Float64
    abortion::Float64
    mrna_decay::Float64
    protein_decay::Float64
end

"""
    EukaryoteBaseRates

Defines a eukaryotic [`Gene`](@ref) reaction cascade's base rates.

# Specification

They are specified in JSON as a JSON object
```
{
    "activation": <...>,
    "deactivation": <...>,
    "trigger": <...>,
    "transcription": <...>,
    "processing": <...>,
    "translation": <...>,
    "abortion": <...>,
    "premrna_decay": <...>,
    "mrna_decay": <...>,
    "protein_decay": <...>
}
```
where `<...>` are JSON numbers setting the corresponding reaction rate
constants.
"""
@kwdef struct EukaryoteBaseRates
    activation::Float64
    deactivation::Float64
    trigger::Float64
    transcription::Float64
    processing::Float64
    translation::Float64
    abortion::Float64
    premrna_decay::Float64
    mrna_decay::Float64
    protein_decay::Float64
end

@kwdef struct DirectRegulator
    from::Symbol
    k::Float64
end

@kwdef struct HillRegulator
    from::Symbol
    at::Float64
    k::Float64 = -1.0
    w::Float64 = 1.0
end

abstract type Regulation end

"""
    Activation <: Regulation

Defines how a `Gene` should be transcriptionally regulated by tempering its
`deactivation` rate.

For details, see [`Gene`](@ref).

# Specification

In JSON, a V1 `Activation` is specified by either one of:
- A JSON array of *slots* `[<slot>...]`, each `<slot>` a JSON object
  `{"from": <from>, "at": <at>}` specifying a single inbound transcriptional
  regulation link. `<from>` is a JSON string referring to the regulating
  chemical species (potenially containing `.` separators), or if it names a
  gene, that gene's proteins. `<at>` is a JSON number specifying the number of
  `<from>` molecules where the link is at half-saturation (assuming quasi-steady
  state). If more than one slot is specified, default aggregation will apply
  (see below).
- A JSON object `{"slots": [<slot>...], "aggregate": <aggregation>}` where
  `<aggregation>` is a JSON string specifying how to aggregate inbound
  regulatory links if more than one of them is present by type -- one of the
  following:
  - `"neutral"`
  - `"minimum"` (default)
  - `"maximum"`
  - `"mean"`
  - `"geometric_mean"`
  - `"harmonic_mean"`
  - `"generalized_mean"`
  If `"aggregate"` is set to `"generalized_mean"`, `"p"` may additionally set
  to a numeric value to control the generalized mean parameter. It is set to
  `0.0` by default, corresponding to a geometric mean. The defaults for
  `"aggregate"` and `"p"` may be overridden by specifying them in
  `"activation"` or `"repression"` JSON objects at the `Definition` level; see
  also `examples/specification/aggregations.schedule.json`.
"""
@kwdef struct Activation <: Regulation
    slots::Vector{HillRegulator} = []
    aggregate::Function = minimum
end
(activation::Activation)(xs; T) =
    isempty(activation.slots) ? one(T) : activation.aggregate(xs)

"""
    Repression <: Regulation

Defines how a `Gene` should be transcriptionally regulated by tempering its
`activation` rate.

For details, see [`Gene`](@ref).

# Specification

Equivalently to [`Activation`](@ref).
"""
@kwdef struct Repression <: Regulation
    slots::Vector{HillRegulator} = []
    aggregate::Function = minimum
end
(repression::Repression)(xs; T) =
    isempty(repression.slots) ? one(T) : repression.aggregate(xs)

"""
    Proteolysis <: Regulation

Defines how a `Gene` should be proteolytically regulated.

For each inbound link, `build` synthesizes a reaction of the form
`@reaction k, proteases + proteins --> proteases`.

(This is in addition to the degradation
`@reaction protein_decay, proteasomes + proteins --> proteasomes` that is
implicitly defined for each gene.)

# Specification

In JSON, a V1 `Proteolysis` is specified by a JSON array of *slots*
`[<slot>...]`, each `<slot>` a JSON object `{"from": <from>, "k": <k>}`
specifying a single inbound proteolytical repression link. `<from>` is a JSON
string referring to the protease chemical species (potenially containing `.`
separators), or if it names a gene, that gene's proteins. `<k>` is a JSON
number specifying the decay reaction propensity.
"""
@kwdef struct Proteolysis <: Regulation
    slots::Vector{DirectRegulator} = []
end

"""
    Gene{BaseRates}

Defines a single gene within a V1 `Definition`, including the reaction rate
constants of its reaction cascade, and optionally inbound regulation and other
information.

# Specification

In JSON, a V1 `Gene` is specified as a JSON object
```
{
    "name": <name>,
    "base_rates": <base_rates>,
    "unique": <unique>,
    "activation": <activation>,
    "repression": <repression>,
    "proteolysis": <proteolysis>
}
```
where only `<base_rates>` is required and the other mappings are optional.

If present, `<name>` must be a JSON string, otherwise it will be set
automatically to the gene's index. (All-digits names are therefore reserved and
should not be specified.) The gene can then be referred to by name in other
`Gene` specifications as a transcription factor (see below) or in additional
mass-action reactions.

`<base_rates>` specifies either [`ProkaryoteBaseRates`](@ref) or
[`EukaryoteBaseRates`](@ref), depending on which [`build`](@ref) will
instantiate a corresponding reaction cascade for this `Gene`. The cascade will
include the following reactions:
```
@reaction_network begin
    trigger, active + \$polymerases --> active + elongations
    transcription, elongations --> premrnas + \$polymerases
    processing, premrnas --> mrnas
    translation, mrnas + \$ribosomes --> mrnas + proteins + \$ribosomes
    abortion, elongations --> \$polymerases
    premrna_decay, premrnas --> 0
    mrna_decay, mrnas --> 0
    protein_decay, proteins + \$proteasomes --> \$proteasomes
end
```
Here, `polymerases`, `ribosomes` and `proteasomes` are chemical species that
will be shared by all genes. `build` will then add reactions for the inter-gene
regulation network, in which the `activation` and `deactivation` base rates will
be respectively tempered by repression and activation links as defined below.
For prokaryotic genes, the `processing` and `premrna_decay` reactions will be
omitted, and `transcription` will directly produce `mrnas` instead.

If present, `<unique>` must be a JSON boolean, otherwise it defaults to `true`.
Setting it to `false` makes the promoter's copy number a dynamic quantity
(starting at 0) instead of the constant 1 -- see also
`examples/specification/copies.schedule.json`.

In the constructed model, genes' promoters dynamically change their
configuration (between `active` and `inactive`) depending on the base
`activation` and `deactivation` rates as well as potentially on the presence of
transcription factors (although only `active` is actually tracked if `<unique>`
is `true`).

If present, `<activation>` must specify inbound transcriptional regulation by
[`Activation`](@ref).

If present, `<repression>` must specify inbound transcriptional regulation by
[`Repression`](@ref).

Transcriptionally regulating links correspond to abstract promoter binding sites
for the regulated gene. They are defined by the choice of transcription factor
(`from`, typically another gene's protein) and a binding affinity parameter
`at`. The binding and unbinding is assumed to occur on a fast time scale and
therefore to always be in quasi-steady state, and `at` specifies the amount of
`from` at which the promoter is bound half of the time. [`build`](@ref)
synthesizes `Catalyst.Reaction`s for promoter activation and deactivation, and
their respective reaction rates are tempered from the specified `activation` and
`deactivation` base rates by the fractional binding of the transcription factor.

If multiple transcriptionally regulating links of the same type are specified
for this `Gene`, in V1, they are simply aggregated as chosen in `<activation>`
and `<repression>`.

If present, `<proteolysis>` must specify inbound regulation by
[`Proteolysis`](@ref).
"""
@kwdef struct Gene{BaseRates}
    name::Symbol
    base_rates::BaseRates
    unique::Bool = true
    activation::Activation = Activation()
    repression::Repression = Repression()
    proteolysis::Proteolysis = Proteolysis()
end

Gene(gene::Gene{BaseRates}; name) where {BaseRates} =
    Gene{BaseRates}(;
        name,
        gene.base_rates,
        gene.activation,
        gene.repression,
        gene.proteolysis,
    )

"""
    Definition

Defines how to construct a `JumpModel` with transcriptionally and
proteolytically inter-regulated genes, and optionally additional reactions with
independent rates constants.

A `Definition` can be constructed directly or parsed from a JSON specification,
for example as part of a `Schedule` execution. It contains instructions for
[`build`](@ref) to assemble a concrete `Catalyst.ReactionSystem` and embed it
into a resulting `JumpModel`, which will interpret the contained reactions as
constituting a mass-action jump process.

# Specification

In JSON, a V1 `Definition` is specified as a JSON object

```
{
    "genes": [<gene>...],
    "reactions": [<Models.Reaction>...],
    "polymerases": <polymerases>,
    "ribosomes": <ribosomes>,
    "proteasomes": <proteasomes>
}}
```
where `[<gene>...]` is a JSON array of [`Gene`](@ref) specifications,
`[<reaction>...]` is a JSON array of [`Models.Reaction`](@ref)s, and
`<polymerases>`, `<ribosomes>` and `<proteasomes>` are JSON strings specifying
the names the species taking that function in the resulting model. All of
these mappings are optional, with the following defaults:
- `[<gene>...]`: `[]`
- `[<reaction>...]`: `[]`
- `<polymerases>`: `"polymerases"`
- `<ribosomes>`: `"ribosomes"`
- `<proteasomes>`: `"proteasomes"`
However, at least one gene or at least one reaction must be specified so that
the system is not empty.
"""
@kwdef struct Definition
    polymerases::Symbol = :polymerases
    ribosomes::Symbol = :ribosomes
    proteasomes::Symbol = :proteasomes
    genes::Vector{Gene} = Gene[]
    reactions::Vector{Models.Reaction} = Models.Reaction[]
    promoter_model::Symbol = :switching
end

cast(::Type{Vector{Gene}}, xs::AbstractVector; context) = [
    cast(
        Gene,
        merge(Dict(:name => lpad(i, ndigits(length(xs)), '0')), x);
        context,
    )
    for (i, x) in enumerate(xs)
]

function cast(::Type{Gene}, x::AbstractDict{Symbol}; context)
    Rates =
        if haskey(x[:base_rates], :processing)
            EukaryoteBaseRates
        else
            ProkaryoteBaseRates
        end

    @invoke cast(
        Gene{Rates}::Type,
        # Ensure we descend on these, even if they are not in x, because we will
        # look up model-wide defaults further down:
        merge(
            Dict(:activation => empty(x), :repression => empty(x)),
            x
        )::AbstractDict{Symbol};
        context,
    )
end

cast(T::Type{<:Regulation}, xs::AbstractVector; context) =
    cast(T, Dict(:slots => xs); context)

cast(::Type{Activation}, x::AbstractDict{Symbol}; context) = @invoke cast(
    Activation::Type,
    merge(get(context, :activation, empty(x)), x)::AbstractDict{Symbol};
    context
)

cast(::Type{Repression}, x::AbstractDict{Symbol}; context) = @invoke cast(
    Repression::Type,
    merge(get(context, :repression, empty(x)), x)::AbstractDict{Symbol};
    context
)

cast(::Type{<:Regulation}, x::AbstractDict{Symbol}, ::Val{:aggregate}; _...) =
    aggregation(Val(Symbol(x[:aggregate])), x)

function cast(::Type{HillRegulator}, x::AbstractDict{Symbol}; _...)
    if haskey(x, Symbol("-k"))
        haskey(x, :k) && error("ambiguous HillRegulator definition")
        x = merge(x, Dict(:k => -x[Symbol("-k")]))
    end

    @invoke cast(HillRegulator::Type, x::AbstractDict{Symbol})
end

aggregation(::Val{:neutral}, _) = one ∘ typeof ∘ first
aggregation(::Val{:minimum}, _) = minimum
aggregation(::Val{:maximum}, _) = maximum
aggregation(::Val{:mean}, _) = mean
aggregation(::Val{:geometric_mean}, _) = geomean
aggregation(::Val{:harmonic_mean}, _) = harmmean
aggregation(k::Val{:generalized_mean}, x) =
    aggregation(k, cast(Float64, get(x, :p, 0.0)))
aggregation(::Val{:generalized_mean}, p::Float64) =
    p == -Inf ? minimum :
    p == -1.0 ? harmmean :
    p == 0.0 ? geomean :
    p == 1.0 ? mean :
    p == Inf ? maximum :
    Base.Fix2(genmean, p)

aggregation_name(::typeof(one ∘ typeof ∘ first)) = "neutral"
aggregation_name(::typeof(geomean)) = "geometric_mean"
aggregation_name(::typeof(harmmean)) = "harmonic_mean"
aggregation_name(f::Function) = nameof(f)

regulation_representation(slots, ::typeof(minimum)) = representation(slots)
regulation_representation(slots, aggregation::Function) = Dict(
    :slots => representation(slots),
    :aggregation => aggregation_name(aggregation),
)
regulation_representation(slots, ::Base.Fix2{typeof(genmean), P}) where {P} =
    Dict(
        :slots => representation(slots),
        :aggregation => "generalized_mean",
        :p => P,
    )

representation(x::ProkaryoteBaseRates) = representation(x, simple = true)
representation(x::EukaryoteBaseRates) = representation(x, simple = true)
representation(x::DirectRegulator) = representation(x, simple = true)
representation(x::HillRegulator) =
    representation(x, simple = true, omit_defaults = [:k => -1.0, :w => 1.0])
representation(x::Activation) = regulation_representation(x.slots, x.aggregate)
representation(x::Repression) = regulation_representation(x.slots, x.aggregate)
representation(x::Proteolysis) = representation(x.slots)
representation(x::Gene) = representation(
    x,
    simple = true,
    omit_defaults = [
        :name => "",
        :activation => [],
        :repression => [],
        :proteolysis => [],
        :unique => true,
    ],
)
representation(x::Definition) = Dict{Symbol, Any}(
    Symbol("{regulation/v1}") => representation(
        x,
        simple = true,
        omit_defaults = [
            :polymerases => "polymerases",
            :ribosomes => "ribosomes",
            :proteasomes => "proteasomes",
            :genes => [],
            :reactions => [],
        ],
    )
)

function Models.describe(definition::Definition)
    genes = Set(gene.name for gene in definition.genes)
    regulator(name) = name in genes ? Symbol("$(name).proteins") : name
    function modulation(gene, kind, from)
        if definition.promoter_model === :equilibrium
            edge_kind = kind === :activation ? :promotes : :inhibits
            reaction = Symbol("$(gene).trigger")
        else
            edge_kind = :inhibits
            reaction = Symbol("$(gene).$(kind === :activation ? "deactivation" : "activation")")
        end
        (; kind=edge_kind, from=regulator(from), to=reaction)
    end
    links = mapreduce(vcat, definition.genes; init=NamedTuple[]) do gene
        vcat(
            map(gene.activation.slots) do (; from, at, k)
                properties = Dict(:at => at, :k => k,
                    :parameters => Dict(
                        :at => Symbol("$(gene.name).activation.$from.at"),
                        :k => Symbol("$(gene.name).activation.$from.k")
                ))
                (; to = gene.name, from, kind = :activation, modulation=modulation(gene.name, :activation, from), properties)
            end,
            map(gene.repression.slots) do (; from, at, k)
                properties = Dict(:at => at, :k => k,
                    :parameters => Dict(
                        :at => Symbol("$(gene.name).repression.$from.at"),
                        :k => Symbol("$(gene.name).repression.$from.k")
                ))
                (; to = gene.name, from, kind = :repression, modulation=modulation(gene.name, :repression, from), properties)
            end,
            map(gene.proteolysis.slots) do (; from, k)
                properties = Dict(:k => k,
                    :parameters => Dict(
                        :k => Symbol("$(gene.name).proteolysis.$from.k")
                ))
                (; to = gene.name, from, kind = :proteolysis, properties)
            end,
        )
    end
    Models.Descriptions([
        Models.Label(
            "'regulation/v1' network with $(length(definition.genes)) genes"
        ),
        Models.RegulatoryNetwork(;
            species_groups = [gene.name for gene in definition.genes],
            links,
            shared_species=Set([definition.polymerases, definition.ribosomes, definition.proteasomes])
        ),
        Models.ReactionNetwork(;reactions=definition.reactions)
    ])
end

function annotate(reaction, kind::Symbol; owner=nothing, metadata...)
    Symbolics.setmetadata(reaction, :kind, kind)
    owner === nothing || Symbolics.setmetadata(reaction, :owner, owner)
    for (key, value) in pairs(metadata)
        Symbolics.setmetadata(reaction, key, value)
    end
    reaction
end

function cascade(definition::Gene{Rates}, ::Val{PromoterModel}; polymerases, ribosomes, proteasomes) where {Rates, PromoterModel}
    name = definition.name
    rxs = Reaction[]

    add(kind, reaction) = push!(rxs, annotate(
        reaction,
        kind;
        owner=name,
        parameters=Dict(:rate => Symbol("$(name).$(kind)"))
    ))

    if PromoterModel === :switching
        add(:trigger, @reaction trigger, active + $polymerases --> active + elongations)
    end

    if Rates === ProkaryoteBaseRates
        add(:transcription, @reaction transcription, elongations --> mrnas + $polymerases)
    elseif Rates === EukaryoteBaseRates
        add(:transcription, @reaction transcription, elongations --> premrnas + $polymerases)
        add(:processing, @reaction processing, premrnas --> mrnas)
        add(:premrna_decay, @reaction premrna_decay, premrnas --> 0)
    end

    add(:translation, @reaction translation, mrnas + $ribosomes --> mrnas + proteins + $ribosomes),
    add(:abortion, @reaction abortion, elongations --> $polymerases),
    add(:mrna_decay, @reaction mrna_decay, mrnas --> 0),
    add(:protein_decay, @reaction protein_decay, proteins + $proteasomes --> $proteasomes)

    ReactionSystem(rxs, default_t(); name)
end

function gene(definition::Gene; polymerases, ribosomes, proteasomes, t, promoter_model)
    result = cascade(definition, Val(promoter_model); polymerases, ribosomes, proteasomes)
    if promoter_model === :switching && !definition.unique
        result = extend(result, @network_component (@species inactive(t);))
    end
    result
end

function species_variable(name::Symbol; t)
    # Replace only the last '.' in the name with the scope separator '₊'
    # because the gene names may contain dots while the kind names certainly do
    # not.
    name = Symbol(replace(String(name), r"\.(?=[^.]*$)" => '₊'))
    only(@species $name(t))
end

species_reference(name::Symbol; t, genes) =
    haskey(genes, name) ? genes[name].proteins : species_variable(name; t)

# Catalyst.hill returns NaN if either both X and K are 0, or if X is 0 and
# n < 0. Since X may become (and often initially is) 0, this restricts what both
# K (which here is the activation/repression half-activation concentration) and
# n may be specified as by the user. To allow both of these violations in the
# limit case, we define a slightly more flexible hill function.
safe_div(K::Real, X::Real) = iszero(K) ? zero(K / X) : K / X
# Allow hill2 to also work with symbolic K
@register_symbolic safe_div(K, X)

hill2(X, v, K, n) = v / (1.0 + safe_div(K, X) ^ n)

make_parameter(name::Symbol) = ModelingToolkit.toparam(Symbolics.variable(name))
make_parameter(name::Symbol, default::Float64) = Symbolics.setmetadata(make_parameter(name), Symbolics.VariableDefaultValue, default)

function each_parameter(callback::Function, definition::Definition)
    for g in definition.genes
        for kind in fieldnames(typeof(g.base_rates))
            callback(Symbol("$(g.name).$(kind)"), getfield(g.base_rates, kind))
        end
        for slot in g.activation.slots
            callback(Symbol("$(g.name).activation.$(slot.from).at"), slot.at)
            callback(Symbol("$(g.name).activation.$(slot.from).k"),  slot.k)
        end
        for slot in g.repression.slots
            callback(Symbol("$(g.name).repression.$(slot.from).at"), slot.at)
            callback(Symbol("$(g.name).repression.$(slot.from).k"),  slot.k)
        end
        for slot in g.proteolysis.slots
            callback(Symbol("$(g.name).proteolysis.$(slot.from).k"), slot.k)
        end

        for (kind, reg) in ((:activation, g.activation), (:repression, g.repression))
            reg.aggregate isa Base.Fix2{typeof(genmean)} &&
                callback(Symbol("$(g.name).$(kind).p"), reg.aggregate.x)
        end
    end
    for rxn in definition.reactions
        callback(Symbol("reaction.$(rxn.name).k⁺"), rxn.k⁺)
        callback(Symbol("reaction.$(rxn.name).k⁻"), rxn.k⁻)
    end
end

Models.parameters(definition::Definition) = let result = Dict{Symbol, Float64}()
    each_parameter(definition) do name, default
        result[name] = default
    end
    result
end

Models.remake(definition::Definition, parameters::AbstractDict{Symbol, <:Real}) = Definition(;
    definition.polymerases,
    definition.ribosomes,
    definition.proteasomes,
    genes=[Models.remake(g, parameters) for g in definition.genes],
    reactions=[
        Models.Reaction(;
            rxn.name, rxn.from, rxn.to,
            k⁺ = get(parameters, Symbol("reaction.$(rxn.name).k⁺"), rxn.k⁺),
            k⁻ = get(parameters, Symbol("reaction.$(rxn.name).k⁻"), rxn.k⁻),
        )
        for rxn in definition.reactions
    ]
)


function Models.remake(gene::Gene, parameters::AbstractDict{Symbol, <:Real})
    T = typeof(gene.base_rates)
    Gene(;
        gene.name, gene.unique,
        base_rates = T(; (
            f => get(parameters, Symbol("$(gene.name).$(f)"), getfield(gene.base_rates, f))
            for f in fieldnames(T)
        )...),
        activation = Activation(;
            gene.activation.aggregate,
            slots = [
                HillRegulator(; slot.from,
                    at = get(parameters, Symbol("$(gene.name).activation.$(slot.from).at"), slot.at),
                    k = get(parameters, Symbol("$(gene.name).activation.$(slot.from).k"), slot.k)
                )
                for slot in gene.activation.slots
            ]
        ),
        repression = Repression(;
            gene.repression.aggregate,
            slots = [
                HillRegulator(; slot.from,
                    at = get(parameters, Symbol("$(gene.name).repression.$(slot.from).at"), slot.at),
                    k = get(parameters, Symbol("$(gene.name).repression.$(slot.from).k"), slot.k)
                )
                for slot in gene.repression.slots
            ]
        ),
        proteolysis = Proteolysis(;
            slots = [
                DirectRegulator(; slot.from,
                    k = get(parameters, Symbol("$(gene.name).proteolysis.$(slot.from).k"), slot.k),
                )
                for slot in gene.proteolysis.slots
            ]
        ),
    )
end

function regulation(
    genes::Dict{Symbol,<:ModelingToolkit.AbstractSystem};
    definition::Definition,
    t::Num,
)

    inactive(target::Gene) =
        if target.unique
            1 - genes[target.name].active
        else
            genes[target.name].inactive
        end

    aggregate(reg::Regulation, xs, name::Symbol, kind::String) =
        isempty(reg.slots) ? one(Num) :
        reg.aggregate isa Base.Fix2{typeof(genmean)} ?
            let hs = collect(xs),
                ws = [make_parameter(Symbol("$(name).$(kind).$(s.from).w"), s.w) for s in reg.slots],
                p = make_parameter(Symbol("$(name).$(kind).p"), reg.aggregate.x),
                W = sum(ws),
                num = sum(i -> ws[i] * hs[i]^p, eachindex(hs)),
                logmean = sum(i -> ws[i] * log(hs[i]), eachindex(hs)) / W
                ifelse(W > 0.5, ifelse(abs(p) < 1e-6, exp(logmean), (num / W)^inv(p)), one(Num))
            end :
            reg.aggregate(collect(xs))

    k_on(target::Gene) = (
        make_parameter(Symbol("$(target.name).activation"), target.base_rates.activation)
        * aggregate(target.repression,
            (  # ^ arguments and value go towards 0 as repression increases
                hill2(species_reference(from; t, genes), 1.0,
                    make_parameter(Symbol("$(target.name).repression.$(from).at"), at),
                    make_parameter(Symbol("$(target.name).repression.$(from).k"), k))
                for (; from, k, at) in target.repression.slots
            ),
            target.name, "repression"
        )
    )

    k_off(target::Gene) = (
        make_parameter(Symbol("$(target.name).deactivation"), target.base_rates.deactivation)
        * aggregate(target.activation,
            (  # ^ arguments and value go towards 0 as activation increases
                hill2(species_reference(from; t, genes), 1.0,
                    make_parameter(Symbol("$(target.name).activation.$(from).at"), at),
                    make_parameter(Symbol("$(target.name).activation.$(from).k"), k))
                for (; from, k, at) in target.activation.slots
            ),
            target.name, "activation"
        )
    )

    function p_active(target::Gene)
        kon = k_on(target)
        koff = k_off(target)
        kon / (kon + koff)
    end

    activation_rate(target::Gene) = k_on(target) * inactive(target)
    deactivation_rate(target::Gene) = k_off(target) * genes[target.name].active

    trigger_rate(target::Gene) = make_parameter(Symbol("$(target.name).trigger"), target.base_rates.trigger)


    # Regulation for the whole network:
    [
        # For each gene...
        mapreduce(vcat, definition.genes, init = Reaction[]) do target::Gene
            vcat(
                if definition.promoter_model === :switching
                    [
                        # ...activation (by tempering promoter deactivation)
                        annotate(Reaction(
                            deactivation_rate(target),
                            [genes[target.name].active],
                            target.unique ? nothing : [genes[target.name].inactive],
                            only_use_rate = true
                        ), :deactivation;
                            owner=target.name,
                            parameters=Dict(:rate => Symbol("$(target.name).deactivation")))

                        # ...repression (by tempering promoter activation)
                        annotate(Reaction(
                            activation_rate(target),
                            target.unique ? nothing : [genes[target.name].inactive],
                            [genes[target.name].active],
                            only_use_rate = true
                        ), :activation;
                            owner=target.name,
                            parameters=Dict(:rate => Symbol("$(target.name).activation")))
                    ]
                elseif definition.promoter_model === :equilibrium
                    polymerases = species_variable(definition.polymerases; t)
                    directions = [
                        polymerases => Int8(1)
                        (species_reference(slot.from; t, genes) => Int8(1) for slot in target.activation.slots)...
                        (species_reference(slot.from; t, genes) => Int8(-1) for slot in target.repression.slots)...
                    ]
                    annotate(Reaction(
                        trigger_rate(target) * p_active(target),
                        [polymerases],
                        [genes[target.name].elongations];
                        metadata = [:propensity_directions => directions]
                    ), :trigger;
                        owner=target.name,
                        parameters=Dict(:rate => Symbol("$(target.name).trigger")))
                end,
                # ...repression (by proteolysis)
                map(target.proteolysis.slots) do (; from, k)
                    proteases = species_reference(from; t, genes)
                    proteins = genes[target.name].proteins
                    k_symbolic = make_parameter(Symbol("$(target.name).proteolysis.$(from).k"), k)
                    reaction =
                        if from == target.name
                            # This is a loop in the proteolysis repression network
                            # and means that the protein decays without another
                            # protease.
                            Reaction(k_symbolic, [proteins], [proteins], [2], [1])
                        else
                            Reaction(k_symbolic, [proteases, proteins], [proteases])
                        end
                    annotate(reaction, :proteolysis;
                        owner=target.name,
                        from,
                        parameters=Dict(:rate => Symbol("$(target.name).proteolysis.$(from).k")))
                end
            )
        end

        # Additionally, we add arbitrary mass-action reactions as specified.
        # Bidirectional pairs are broken up, and reactions are only included if
        # their rate is nonzero.
        [
            annotate(Reaction(
                # need to use :k⁺ instead of :k⁺ because ₊ is used as a scope separator in MTK
                make_parameter(Symbol("reaction.$(name).k⁺"), k⁺),
                species_reference.(keys(from.counts); t, genes),
                species_reference.(keys(to.counts); t, genes),
                collect(values(from.counts)),
                collect(values(to.counts)),
            ), :auxiliary;
                name,
                direction=:forward,
                parameters=Dict(
                    :k⁺ => Symbol("reaction.$(name).k⁺")
                ))
            for (; name, from, k⁺, to) in definition.reactions
            if k⁺ > 0.0
        ]

        [
            annotate(Reaction(
                make_parameter(Symbol("reaction.$(name).k⁻"), k⁻),
                species_reference.(keys(to.counts); t, genes),
                species_reference.(keys(from.counts); t, genes),
                collect(values(to.counts)),
                collect(values(from.counts)),
            ), :auxiliary;
                name,
                direction=:reverse,
                parameters=Dict(
                    :k⁻ => Symbol("reaction.$(name).k⁻")
                ))
            for (; name, from, k⁻, to) in definition.reactions
            if k⁻ > 0.0
        ]
    ]
end

const JUMP_PROCESSES_METHODS = Dict(
    :Direct => Direct,
    :SortingDirect => SortingDirect,
    :RSSA => RSSA,
    :RSSACR => RSSACR,
    :TauSplitting => TauSplitting
)

pick_method(system; method) = get(JUMP_PROCESSES_METHODS, method) do
    if numspecies(system) < 100 && numreactions(system) < 1000
        SortingDirect
    else
        RSSACR
    end
end


"""
Bracket data pinning the promoter state species to exact (zero-width) brackets.

Every regulated rate is linear increasing in its gene's promoter state species
while being decreasing in the regulator proteins, which makes it non-monotone
over a box that gives the promoter any width. Since the promoter species are
booleans, bracketing them at all is pointless anyway; pinning them to zero width
removes that coordinate from the box and restores monotonicity in the rest.

The cost is that a promoter toggle always falls outside its bracket and forces a
propensity recomputation. Promoter transitions are a small minority of events
relative to protein birth/death, so most of RSSA's saving is retained.
"""
function promoter_bracket_data(system)
    names = [String(ModelingToolkit.getname(s)) for s in ModelingToolkit.unknowns(system)]
    exact = [endswith(n, "active") for n in names]  # matches `active` and `inactive`
    # The vector-valued accessors in JumpProcesses dispatch on the field type
    # being exactly `AbstractVector`, so the type parameters are spelled out.
    BracketData{AbstractVector{Float64}, AbstractVector{Int}}(
        [e ? 0.0 : 0.1 for e in exact],
        [e ? 0 : 25 for e in exact],
        [e ? 0 : 4 for e in exact],
    )
end

function aggregator_options(algorithm, reaction_system, jump_system, definition)
    algorithm in (RSSA, RSSACR, TauSplitting) || return (;)
    all(definition.genes) do gene
        all(slot -> slot.k <= 0.0, gene.activation.slots) &&
            all(slot -> slot.k <= 0.0, gene.repression.slots)
    end || error("$(nameof(algorithm)) requires every activation/repression `k` to be non-positive so that propensities are monotone in the regulator counts; this definition has a positive `k`.")

    if definition.promoter_model === :equilibrium
        all(definition.genes) do gene
            activators = Set(slot.from for slot in gene.activation.slots)
            repressors = Set(slot.from for slot in gene.repression.slots)
            isdisjoint(activators, repressors)
        end || error("$(nameof(algorithm)) does not support a species both activating and repressing the same equilibrium promoter")
    end

    system = Catalyst.flatten(reaction_system)
    reactions = Catalyst.reactions(system)
    unknowns = ModelingToolkit.unknowns(jump_system)
    species_index = Dict(ModelingToolkit.value(species) => i
        for (i, species) in enumerate(unknowns)
    )
    constant_reactions = filter(reactions) do rx
        !Catalyst.ismassaction(rx, system)
    end

    function propensity_directions(rx, species_index)
        directions = get(Dict(rx.metadata), :propensity_directions, nothing)
        directions === nothing || return [
            species_index[ModelingToolkit.value(species)] => dir
            for (species, dir) in directions
        ]

        substrates = Set(ModelingToolkit.value.(rx.substrates))
        species = union(substrates, ModelingToolkit.value.(Catalyst.get_variables(rx.rate)))
        Pair{Int, Int8}[
            species_index[s] => Int8(s in substrates ? 1 : -1)
            for s in species
            if haskey(species_index, s)
        ]
    end
    propensity_bounds = DirectionalBounds{Float64}([propensity_directions(rx, species_index)
        for rx in constant_reactions], length(unknowns))

    if algorithm in (RSSA, RSSACR)
        (;
            bracket_data = promoter_bracket_data(jump_system),
            propensity_bounds
        )
    elseif algorithm === TauSplitting
        jumptostoich_map = [
            Pair{Int, Int}[
                species_index[species] => Int(stoch)
                for (raw_species, stoch) in rx.netstoich
                for species in (ModelingToolkit.value(raw_species),)
                if haskey(species_index, species)
            ]
            for rx in constant_reactions
        ]
        (;
            propensity_bounds,
            jumptostoich_map
        )
    else
        (;)
    end
end

"""
    build(specification::AbstractDict{Symbol})
    build(definition::Definition; method::Symbol = :default)

Construct a [`SciML.JumpModel`](@ref) from a [`Definition`](@ref).

When interpreting a JSON specification, this function (in its first form) is
called to construct a concrete regulation model on encountering a
`{"{regulation/v1}": {...}}` literal. It will first destructure the parsed JSON
into a `Definition` and then proceed from there.

This function is also called (directly in its second form) when lowering a
higher-level model template (such as [`Models.Differentiation`](@ref)).

The result is constructed by first assembling a `Catalyst.ReactionSystem` as
specified by `definition`, interpreting it as a `JumpProcesses.JumpSystem`,
packaging that up as a `JumpModel`, wrapping that up with the `ReactionSystem`
in a [`Models.Wrapped`](@ref), and finally further wrapping that in another
`Wrapped` with `definition`. This will result in the following stack of
abstractions:
- [`SciML.JumpModel`](@ref Models.SciML.JumpModel), specified by a
- `Catalyst.ReactionSystem`, specified by
- `definition::`[`V1.Definition`](@ref), potentially specified by
- `specification`

The `method` to use for the `JumpModel` can be set in `specification[:method]`
or passed as a keyword argument. By default, a method will automatically be
chosen based on the size of the `ReactionSystem`: `SortingDirect` for small
systems (having less than 100 species and less than 1000 reactions), and
`RSSACR` otherwise.

# Specification

V1 models are specified in JSON as `{"{regulation/v1}": <definition>}` where
`<definition>` specifies a [`Definition`](@ref) as described there.

For an example, see `examples/specification/simple.schedule.json`.
"""
function build end

build(specification::AbstractDict{Symbol}) = build(
    cast(Definition, specification),
    method = Symbol(get(specification, :method, "default"))
)

function build(definition::Definition; method::Symbol = :default)
    allequal(typeof.(definition.genes)) ||
        error("mixing eukaryotic and prokaryotic genes is forbidden")

    let names = [rxn.name for rxn in definition.reactions]
        duplicates = unique([n for n in names if count(==(n), names) > 1])
        isempty(duplicates) ||
            error("reaction names must be unique; duplicated: ", join(duplicates, ", "))
    end

    definition.promoter_model in (:switching, :equilibrium) ||
        error("promoter_model must be :switching or :equilibrium, got $(definition.promoter_model)")
    definition.promoter_model === :equilibrium && any(g -> !g.unique, definition.genes) &&
        error("equilibrium promoter model does not support `unique=false` genes")

    t = default_t()
    polymerases = species_variable(definition.polymerases; t)
    ribosomes = species_variable(definition.ribosomes; t)
    proteasomes = species_variable(definition.proteasomes; t)

    genes = Dict{Symbol, ReactionSystem}(
        g.name => gene(
            g,
            polymerases = ParentScope(polymerases),
            ribosomes = ParentScope(ribosomes),
            proteasomes = ParentScope(proteasomes);
            t,
            definition.promoter_model
        )
        for g in definition.genes
    )

    @named reaction_system = ReactionSystem(
        regulation(genes; definition, t),
        t,
        systems=collect(values(genes)),
        initial_conditions=Dict(
            getproperty(genes[g.name], kind) => getfield(g.base_rates, kind)
            for g in definition.genes
            for kind in fieldnames(typeof(g.base_rates))
            if kind ∉ (:activation, :deactivation) &&
                (definition.promoter_model === :switching || kind !== :trigger)
        )
    )
    reaction_system = complete(reaction_system)
    jump_system = complete(jump_model(reaction_system))
    algorithm = pick_method(reaction_system; method)

    Models.Wrapped(;
        definition,
        model = Models.Wrapped(
            definition = reaction_system,
            model = SciML.JumpModel(
                jump_system,
                algorithm();
                aggregator_options(algorithm, reaction_system, jump_system, definition)...,
            ),
        ),
    )
end

constructor(::Val{Symbol("regulation/v1")}) = build

"""
    knockout(model; genes, soft=false)

Knock out `genes` from a V1 model.

With `soft=false` (default, *hard* knockout), the genes are structurally removed
from the `Definition` and the model is fully recompiled.

With `soft=true` (*soft* knockout), the `trigger` reaction rate is zeroed for
each knocked-out gene, without recompiling the model.

# Specification

`{"{regulation/v1/knockout}": {"of": {"\$": "do"}, "genes": ["A"], "soft": false}}`
"""
function knockout end

knockout(specification::AbstractDict{Symbol}) = knockout(
    specification[:of];
    genes = Symbol.(specification[:genes]),
    soft = get(specification, :soft, false),
)

function knockout(model::Models.Wrapped; genes, soft=false)
    if soft
        Models.remake(model, Dict(Symbol("$(g).trigger") => 0.0 for g in genes))
    else
        knockout(model.definition, model.model; genes)
    end
end
knockout(::Any, model::Models.Model; genes, soft=false) = knockout(model; genes, soft)
knockout(definition::Definition, ::Models.Model; genes, soft=false) = build(knockout(definition; genes))

function knockout(definition::Definition; genes)
    knocked_out = Set(genes)
    remaining = filter(g -> g.name ∉ knocked_out, definition.genes)
    cleaned = map(remaining) do gene
        Gene(;
            gene.name, gene.base_rates, gene.unique,
            activation = Activation(;
                gene.activation.aggregate,
                slots = filter(s -> s.from ∉ knocked_out, gene.activation.slots),
            ),
            repression = Repression(;
                gene.repression.aggregate,
                slots = filter(s -> s.from ∉ knocked_out, gene.repression.slots),
            ),
            proteolysis = Proteolysis(;
                slots = filter(s -> s.from ∉ knocked_out, gene.proteolysis.slots),
            ),
        )
    end
    Definition(; definition.polymerases, definition.ribosomes, definition.proteasomes,
        genes = cleaned, definition.reactions)
end

constructor(::Val{Symbol("regulation/v1/knockout")}) = knockout

end
