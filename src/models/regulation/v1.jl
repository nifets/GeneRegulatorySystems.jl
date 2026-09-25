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
    BaseRates

Defines a [`Gene`](@ref) reaction cascade's base rates.

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
constants. `"processing"` and `"premrna_decay"` are optional; omitting them
drops `premrnas` from the cascade, so that `transcription` directly produces
`mrnas`.
"""
@kwdef struct BaseRates
    activation::Float64
    deactivation::Float64
    trigger::Float64
    transcription::Float64
    translation::Float64
    abortion::Float64
    mrna_decay::Float64
    protein_decay::Float64
    processing::Union{Nothing, Float64} = nothing
    premrna_decay::Union{Nothing, Float64} = nothing
end

struct Stage
    species::Symbol
    birth::Symbol
    death::Symbol
    persists::Bool
end

const STAGES = (
    Stage(:elongations, :trigger, :abortion, false),
    Stage(:premrnas, :transcription, :premrna_decay, false),
    Stage(:mrnas, :processing, :mrna_decay, true),
    Stage(:proteins, :translation, :protein_decay, true),
)

available(rates::BaseRates, species::Symbol) = species !== :premrnas ||
    (rates.processing !== nothing && rates.premrna_decay !== nothing)

default_species(rates::BaseRates) = [
    :active,
    (stage.species for stage in STAGES if available(rates, stage.species))...,
]

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
    Gene

Defines a single gene within a V1 `Definition`, including the reaction rate
constants of its reaction cascade, and optionally inbound regulation and other
information.

# Specification

In JSON, a V1 `Gene` is specified as a JSON object
```
{
    "name": <name>,
    "base_rates": <base_rates>,
    "species": [<species>...],
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

`<base_rates>` specifies [`BaseRates`](@ref), from which [`build`](@ref) will
instantiate a reaction cascade for this `Gene`. The cascade will include the
following reactions:
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

If present, `[<species>...]` must be a JSON array naming a subset of
`"active"`, `"elongations"`, `"premrnas"`, `"mrnas"` and `"proteins"`; it
defaults to all of the ones the base rates provide, and may also be given once
for the whole model alongside `"genes"`. Omitting a species removes it from the
cascade and merges the reactions that produced and consumed it into one, whose
rate constant is scaled to preserve the steady-state flux reaching the next
species: by `k_next / (k_next + k_loss)` where the omitted species is consumed
(`elongations`, `premrnas`), and by `k_next / k_loss` where it acts as a
catalyst (`mrnas`). This preserves steady-state means but not the dwell time the
omitted species contributed, so a reduced gene is burstier than the full one.

Omitting `"active"` puts the promoter in quasi-steady state instead, folding its
occupancy into the `trigger` rate as before, and is incompatible with
`"unique": false`. Other genes may regulate this one by naming it, which refers
to the last species it keeps, or by naming a species explicitly as
`"<gene>.<species>"`.

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
@kwdef struct Gene
    name::Symbol
    base_rates::BaseRates
    species::Union{Nothing, Vector{Symbol}} = nothing
    unique::Bool = true
    activation::Activation = Activation()
    repression::Repression = Repression()
    proteolysis::Proteolysis = Proteolysis()
end

Gene(gene::Gene; kwargs...) = Gene(;
    (field => getfield(gene, field) for field in fieldnames(Gene))...,
    kwargs...,
)

species_of(gene::Gene) = @something(gene.species, default_species(gene.base_rates))

keeps(gene::Gene, species::Symbol) = gene.species === nothing ?
    available(gene.base_rates, species) : species in gene.species

switching(gene::Gene) = keeps(gene, :active)

function first_transcript(gene::Gene)
    for stage in STAGES
        keeps(gene, stage.species) && return stage.species
    end
end

function last_species(gene::Gene)
    for stage in reverse(STAGES)
        keeps(gene, stage.species) && return stage.species
    end
    :active
end

regulator_name(gene::Gene) = Symbol("$(gene.name).$(last_species(gene))")

collapse_factor(dropped::Stage, next::Stage, rates::BaseRates) = collapse_factor(
    getfield(rates, next.birth), getfield(rates, dropped.death), dropped.persists)
collapse_factor(::Nothing, ::Nothing, ::Bool) = 1.0
collapse_factor(out::Float64, loss::Float64, persists::Bool) =
    persists ? out / loss : out / (out + loss)

function rate_values(gene::Gene)
    rates = gene.base_rates
    effective = Dict{Symbol, Float64}()
    previous = 0
    for (i, stage) in enumerate(STAGES)
        keeps(gene, stage.species) || continue
        birth = STAGES[previous + 1].birth
        effective[birth] = foldl(
            (rate, j) -> rate * collapse_factor(STAGES[j], STAGES[j + 1], rates),
            (previous + 1):(i - 1);
            init = getfield(rates, birth),
        )
        effective[stage.death] = getfield(rates, stage.death)
        previous = i
    end
    effective
end

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
end

Definition(base::Definition; kwargs...) = Definition(;
    (field => getfield(base, field) for field in fieldnames(Definition))...,
    kwargs...
)

cast(::Type{Vector{Gene}}, xs::AbstractVector; context) = [
    cast(
        Gene,
        merge(Dict(:name => lpad(i, ndigits(length(xs)), '0')), x);
        context,
    )
    for (i, x) in enumerate(xs)
]

cast(::Type{Gene}, x::AbstractDict{Symbol}; context) = @invoke cast(
    Gene::Type,
    # Ensure we descend on these, even if they are not in x, because we will
    # look up model-wide defaults further down:
    merge(
        Dict(:activation => empty(x), :repression => empty(x)),
        haskey(context, :species) ? Dict(:species => context[:species]) : empty(x),
        x
    )::AbstractDict{Symbol};
    context,
)

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

function thermodynamic end

aggregation(::Val{:neutral}, _) = one ∘ typeof ∘ first
aggregation(::Val{:thermodynamic}, _) = thermodynamic
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
aggregation_name(::typeof(thermodynamic)) = "thermodynamic"
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

representation(x::BaseRates) = representation(x, simple = true, omit_defaults = [:processing => nothing, :premrna_decay => nothing])
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
        :species => nothing,
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
    genes = Dict(gene.name => gene for gene in definition.genes)
    regulator(name) = haskey(genes, name) ? regulator_name(genes[name]) : name
    present(slots) = filter(slot -> !iszero(slot.w), slots)
    function modulation(gene, kind, from)
        if switching(genes[gene])
            edge_kind = :inhibits
            reaction = Symbol("$(gene).$(kind === :activation ? "deactivation" : "activation")")
        else
            edge_kind = kind === :activation ? :promotes : :inhibits
            reaction = Symbol("$(gene).trigger")
        end
        (; kind=edge_kind, from=regulator(from), to=reaction)
    end
    links = mapreduce(vcat, definition.genes; init=NamedTuple[]) do gene
        vcat(
            map(present(gene.activation.slots)) do (; from, at, k, w)
                properties = Dict(:at => at, :k => k, :w => w,
                    :parameters => Dict(
                        :at => Symbol("$(gene.name).activation.$from.at"),
                        :k => Symbol("$(gene.name).activation.$from.k")
                ))
                (; to = gene.name, from, kind = :activation, modulation=modulation(gene.name, :activation, from), properties)
            end,
            map(present(gene.repression.slots)) do (; from, at, k, w)
                properties = Dict(:at => at, :k => k, :w => w,
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

death_reaction(::Val{:elongations}, x, k; polymerases, _...) =
    Reaction(k, [x], [polymerases])
death_reaction(::Val{:proteins}, x, k; proteasomes, _...) =
    Reaction(k, [x, proteasomes], [proteasomes])
death_reaction(::Val, x, k; _...) = Reaction(k, [x], nothing)

function cascade(definition::Gene; polymerases, ribosomes, proteasomes, t)
    name = definition.name
    rxs = Reaction[]

    add(kind, reaction) = push!(rxs, annotate(
        reaction,
        kind;
        owner=name,
        parameters=Dict(:rate => Symbol("$(name).$(kind)"))
    ))

    previous = 0
    for (i, stage) in enumerate(STAGES)
        keeps(definition, stage.species) || continue
        target = species_variable(stage.species; t)
        catalysts = stage.species === :proteins ? (ribosomes,) : ()

        if previous == 0
            if switching(definition)
                active = species_variable(:active; t)
                held = stage.species === :elongations
                add(STAGES[1].birth, Reaction(make_parameter(STAGES[1].birth),
                    [active, polymerases, catalysts...],
                    [active, target, catalysts...,
                        (held ? () : (polymerases,))...]))
            end
        else
            kind = STAGES[previous + 1].birth
            source = species_variable(STAGES[previous].species; t)
            released = STAGES[previous].species === :elongations ? (polymerases,) : ()
            add(kind, Reaction(make_parameter(kind),
                [source, catalysts...],
                [target, catalysts..., released...,
                    (STAGES[previous].persists ? (source,) : ())...]))
        end

        add(stage.death, death_reaction(Val(stage.species), target,
            make_parameter(stage.death); polymerases, proteasomes))
        previous = i
    end

    ReactionSystem(rxs, t; name)
end

function gene(definition::Gene; polymerases, ribosomes, proteasomes, t)
    result = cascade(definition; polymerases, ribosomes, proteasomes, t)
    switching(definition) || return result
    first_transcript(definition) === nothing &&
        (result = extend(result, @network_component (@species active(t);)))
    definition.unique && return result
    extend(result, @network_component (@species inactive(t);))
end

function species_variable(name::Symbol; t)
    name = SciML.symbolic_name(name)
    only(@species $name(t))
end

function observed_variable(name::Symbol; t)
    name = SciML.symbolic_name(name)
    only(@variables $name(t))
end

function hill2(X, v, K, n)
    m = abs(n)
    a = X^m
    b = K^m
    total = a + b
    ifelse(total > 0,
        ifelse(n >= 0, v * a / total, v * b / total), ifelse(n >= 0, v, zero(v)))
end

make_parameter(name::Symbol) = ModelingToolkit.toparam(Symbolics.variable(name))
make_parameter(name::Symbol, default::Float64) = Symbolics.setmetadata(make_parameter(name), Symbolics.VariableDefaultValue, default)

function each_parameter(callback::Function, definition::Definition)
    for g in definition.genes
        for kind in fieldnames(BaseRates)
            rate = getfield(g.base_rates, kind)
            rate === nothing ||
                callback(Symbol("$(g.name).$(kind)"), rate)
        end
        for slot in g.activation.slots
            callback(Symbol("$(g.name).activation.$(slot.from).at"), slot.at)
            callback(Symbol("$(g.name).activation.$(slot.from).k"),  slot.k)
            callback(Symbol("$(g.name).activation.$(slot.from).w"),  slot.w)
        end
        for slot in g.repression.slots
            callback(Symbol("$(g.name).repression.$(slot.from).at"), slot.at)
            callback(Symbol("$(g.name).repression.$(slot.from).k"),  slot.k)
            callback(Symbol("$(g.name).repression.$(slot.from).w"),  slot.w)
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
        gene.name, gene.unique, gene.species,
        base_rates = T(; (
            f => get(parameters, Symbol("$(gene.name).$(f)"), getfield(gene.base_rates, f))
            for f in fieldnames(T)
        )...),
        activation = Activation(;
            gene.activation.aggregate,
            slots = [
                HillRegulator(; slot.from,
                    at = get(parameters, Symbol("$(gene.name).activation.$(slot.from).at"), slot.at),
                    k = get(parameters, Symbol("$(gene.name).activation.$(slot.from).k"), slot.k),
                    w = get(parameters, Symbol("$(gene.name).activation.$(slot.from).w"), slot.w)
                )
                for slot in gene.activation.slots
            ]
        ),
        repression = Repression(;
            gene.repression.aggregate,
            slots = [
                HillRegulator(; slot.from,
                    at = get(parameters, Symbol("$(gene.name).repression.$(slot.from).at"), slot.at),
                    k = get(parameters, Symbol("$(gene.name).repression.$(slot.from).k"), slot.k),
                    w = get(parameters, Symbol("$(gene.name).repression.$(slot.from).w"), slot.w)
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

struct Regulators
    source::Vector{Int}
    at::Vector{SciML.ParameterIndex}
    k::Vector{SciML.ParameterIndex}
    w::Vector{SciML.ParameterIndex}
end

weight(r::Regulators, j, p) = @inbounds p[r.w[j]]
weights(r::Regulators, p) = sum(j -> weight(r, j, p), eachindex(r.w); init = 0.0)

hill(r::Regulators, j, u, p) =
    @inbounds hill2(u[r.source[j]], 1.0, p[r.at[j]], p[r.k[j]])

aggregate(::typeof(one ∘ typeof ∘ first), r::Regulators, u, p) = 1.0
aggregate(::typeof(thermodynamic), r::Regulators, u, p) =
    error("thermodynamic promoters have no per-class aggregate")
included(r::Regulators, j, p) = weight(r, j, p) > 0.5

aggregate(::typeof(minimum), r::Regulators, u, p) =
    minimum(j -> included(r, j, p) ? hill(r, j, u, p) : 1.0, eachindex(r.source))
function aggregate(::typeof(maximum), r::Regulators, u, p)
    weights(r, p) > 0.5 || return 1.0
    maximum(j -> included(r, j, p) ? hill(r, j, u, p) : 0.0, eachindex(r.source))
end

function aggregate(::typeof(mean), r::Regulators, u, p)
    total = weights(r, p)
    total > 0.5 || return 1.0
    sum(j -> weight(r, j, p) * hill(r, j, u, p), eachindex(r.w)) / total
end

function aggregate(::typeof(geomean), r::Regulators, u, p)
    total = weights(r, p)
    total > 0.5 || return 1.0
    exp(sum(j -> weight(r, j, p) * log(hill(r, j, u, p)), eachindex(r.w)) / total)
end

function aggregate(::typeof(harmmean), r::Regulators, u, p)
    total = weights(r, p)
    total > 0.5 || return 1.0
    total / sum(j -> weight(r, j, p) / hill(r, j, u, p), eachindex(r.w))
end

function aggregate(f::Base.Fix2{typeof(genmean)}, r::Regulators, u, p)
    q = f.x
    total = weights(r, p)
    total > 0.5 || return 1.0
    abs(q) < 1e-6 &&
        return exp(sum(j -> weight(r, j, p) * log(hill(r, j, u, p)), eachindex(r.w)) / total)
    (sum(j -> weight(r, j, p) * hill(r, j, u, p) ^ q, eachindex(r.w)) / total) ^ inv(q)
end

function apply(a, r::Regulators, u, p)
    n = length(r.source)
    n == 0 && return 1.0
    n == 1 && return included(r, 1, p) ? hill(r, 1, u, p) : 1.0
    aggregate(a, r, u, p)
end
struct SwitchingRate{A, N}
    name::Symbol
    k::SciML.ParameterIndex
    regulators::Regulators
    aggregation::A
    site::Int
    active::Int
    inactive::Int
    scale::Float64
    offset::Float64
    affect::SciML.Affect{N}
end

occupancy(f::SwitchingRate, u) = @inbounds f.offset + f.scale * u[f.site]
corner(f::SwitchingRate, lo, hi) = max(occupancy(f, f.scale > 0 ? lo : hi), 0.0)

function active_fraction(f::SwitchingRate, u)
    f.inactive == 0 && return Float64(u[f.active])
    total = u[f.active] + u[f.inactive]
    total > 0 ? u[f.active] / total : 0.0
end

SciML.observed_activity(f::SwitchingRate, u, _) = f.name => active_fraction(f, u)

(f::SwitchingRate)(u, p, _) = p[f.k] * apply(f.aggregation, f.regulators, u, p) * occupancy(f, u)

SciML.lrate(f::SwitchingRate, ulow, uhigh, p) =
    p[f.k] * apply(f.aggregation, f.regulators, uhigh, p) * corner(f, ulow, uhigh)
SciML.urate(f::SwitchingRate, ulow, uhigh, p) =
    p[f.k] * apply(f.aggregation, f.regulators, ulow, p) * corner(f, uhigh, ulow)

struct EquilibriumRate{A, B, N}
    name::Symbol
    k::SciML.ParameterIndex
    kon::SciML.ParameterIndex
    koff::SciML.ParameterIndex
    repression::Regulators
    activation::Regulators
    aggregation_repression::A
    aggregation_activation::B
    polymerases::Int
    affect::SciML.Affect{N}
end

function active_fraction(f::EquilibriumRate, urepression, uactivation, p)
    kon = p[f.kon] * apply(f.aggregation_repression, f.repression, urepression, p)
    koff = p[f.koff] * apply(f.aggregation_activation, f.activation, uactivation, p)
    kon / (kon + koff)
end

SciML.observed_activity(f::EquilibriumRate, u, p) =
    f.name => active_fraction(f, u, u, p)

(f::EquilibriumRate)(u, p, _) =
    @inbounds p[f.k] * active_fraction(f, u, u, p) * u[f.polymerases]

SciML.lrate(f::EquilibriumRate, ulow, uhigh, p) =
    @inbounds p[f.k] * active_fraction(f, uhigh, ulow, p) * ulow[f.polymerases]
SciML.urate(f::EquilibriumRate, ulow, uhigh, p) =
    @inbounds p[f.k] * active_fraction(f, ulow, uhigh, p) * uhigh[f.polymerases]

struct ThermodynamicRate{N}
    name::Symbol
    k::SciML.ParameterIndex
    kon::SciML.ParameterIndex
    koff::SciML.ParameterIndex
    repression::Regulators
    activation::Regulators
    polymerases::Int
    affect::SciML.Affect{N}
end

function logsum(r::Regulators, u, p)
    total = 0.0
    @inbounds for j in eachindex(r.source)
        total += log1p((u[r.source[j]] / p[r.at[j]]) ^ abs(p[r.k[j]]))
    end
    total
end

function active_fraction(f::ThermodynamicRate, urepression, uactivation, p)
    @inbounds kon, koff = p[f.kon], p[f.koff]
    la = logsum(f.activation, uactivation, p)
    lr = logsum(f.repression, urepression, p)
    (expm1(la) + kon / (kon + koff)) * exp(-(la + lr))
end

SciML.observed_activity(f::ThermodynamicRate, u, p) =
    f.name => active_fraction(f, u, u, p)

(f::ThermodynamicRate)(u, p, _) =
    @inbounds p[f.k] * active_fraction(f, u, u, p) * u[f.polymerases]

SciML.lrate(f::ThermodynamicRate, ulow, uhigh, p) =
    @inbounds p[f.k] * active_fraction(f, uhigh, ulow, p) * ulow[f.polymerases]
SciML.urate(f::ThermodynamicRate, ulow, uhigh, p) =
    @inbounds p[f.k] * active_fraction(f, ulow, uhigh, p) * uhigh[f.polymerases]

weighted(f, hs, ws, name, kind) = f(hs)

weighted(::typeof(minimum), hs, ws, name, kind) =
    reduce(min, (ifelse(ws[i] > 0, hs[i], one(Num)) for i in eachindex(hs)))

function weighted(::typeof(maximum), hs, ws, name, kind)
    W = sum(ws)
    best = reduce(max, (ifelse(ws[i] > 0, hs[i], zero(Num)) for i in eachindex(hs)))
    ifelse(W > 0, best, one(Num))
end

function weighted(::typeof(mean), hs, ws, name, kind)
    W = sum(ws)
    ifelse(W > 0, sum(i -> ws[i] * hs[i], eachindex(hs)) / W, one(Num))
end

function weighted(::typeof(geomean), hs, ws, name, kind)
    W = sum(ws)
    ifelse(W > 0, exp(sum(i -> ws[i] * log(hs[i]), eachindex(hs)) / W), one(Num))
end

function weighted(::typeof(harmmean), hs, ws, name, kind)
    W = sum(ws)
    ifelse(W > 0, W / sum(i -> ws[i] / hs[i], eachindex(hs)), one(Num))
end

function weighted(f::Base.Fix2{typeof(genmean)}, hs, ws, name, kind)
    p = make_parameter(Symbol("$(name).$(kind).p"), f.x)
    W = sum(ws)
    num = sum(i -> ws[i] * hs[i]^p, eachindex(hs))
    logmean = sum(i -> ws[i] * log(hs[i]), eachindex(hs)) / W
    ifelse(W > 0, ifelse(abs(p) < 1e-6, exp(logmean), (num / W)^inv(p)), one(Num))
end

function regulators(indices::SciML.Indices, genes, gene::Gene, kind::String, slots, aggregate)
    aggregate isa typeof(one ∘ typeof ∘ first) && return Regulators(
        Int[], SciML.ParameterIndex[], SciML.ParameterIndex[], SciML.ParameterIndex[])
    regulator(from) = haskey(genes, from) ? regulator_name(genes[from]) : from
    Regulators(
        [indices.unknowns[regulator(slot.from)] for slot in slots],
        [indices.parameters[Symbol("$(gene.name).$(kind).$(slot.from).at")] for slot in slots],
        [indices.parameters[Symbol("$(gene.name).$(kind).$(slot.from).k")] for slot in slots],
        [indices.parameters[Symbol("$(gene.name).$(kind).$(slot.from).w")] for slot in slots],
    )
end

function promoter_rate(jump, genes, definition::Definition, indices::SciML.Indices)
    net_stoich = [
        (SciML.normalize_name(ModelingToolkit.value(affect.lhs)),
            SciML.stoichiometry(affect))
        for affect in jump.affect!
    ]
    affect = SciML.Affect(
        Tuple(indices.unknowns[species] for (species, _) in net_stoich),
        Tuple(Int8(change) for (_, change) in net_stoich),
    )

    @something(
        promoter_rate(Val(:switching), net_stoich, affect, genes, definition, indices),
        promoter_rate(Val(:thermodynamic), net_stoich, affect, genes, definition, indices),
        promoter_rate(Val(:equilibrium), net_stoich, affect, genes, definition, indices),
        Some(nothing),
    )
end

# the MTK gernerated jumps are not labelled so we need to reverse engineer which gene they belong to
function gene_of(match, net_stoich::Vector{Tuple{Symbol, Int}}, genes)
    for (species, change) in net_stoich
        text = String(species)
        separator = findlast('.', text)
        separator === nothing && continue
        name = Symbol(SubString(text, 1, separator - 1))
        haskey(genes, name) || continue
        match(genes[name], Symbol(SubString(text, separator + 1)), change) &&
            return name
    end
end

function promoter_rate(
    ::Val{:equilibrium}, net_stoich, affect, genes, definition, indices,
)
    name = gene_of(net_stoich, genes) do gene, kind, change
        change > 0 && !switching(gene) && kind === first_transcript(gene)
    end
    name === nothing && return nothing
    gene = genes[name]
    EquilibriumRate(
        Symbol("$(name).activity"),
        indices.parameters[Symbol("$(name).trigger")],
        indices.parameters[Symbol("$(name).activation")],
        indices.parameters[Symbol("$(name).deactivation")],
        regulators(indices, genes, gene, "repression", gene.repression.slots,
            gene.repression.aggregate),
        regulators(indices, genes, gene, "activation", gene.activation.slots,
            gene.activation.aggregate),
        gene.repression.aggregate,
        gene.activation.aggregate,
        indices.unknowns[definition.polymerases],
        affect,
    )
end

function promoter_rate(
    ::Val{:thermodynamic}, net_stoich, affect, genes, definition, indices,
)
    name = gene_of(net_stoich, genes) do gene, kind, change
        change > 0 && !switching(gene) && kind === first_transcript(gene)
    end
    name === nothing && return nothing
    gene = genes[name]
    gene.activation.aggregate === thermodynamic || return nothing
    ThermodynamicRate(
        Symbol("$(name).activity"),
        indices.parameters[Symbol("$(name).trigger")],
        indices.parameters[Symbol("$(name).activation")],
        indices.parameters[Symbol("$(name).deactivation")],
        regulators(indices, genes, gene, "repression", gene.repression.slots,
            gene.repression.aggregate),
        regulators(indices, genes, gene, "activation", gene.activation.slots,
            gene.activation.aggregate),
        indices.unknowns[definition.polymerases],
        affect,
    )
end

function promoter_rate(
    ::Val{:switching}, net_stoich, affect, genes, definition, indices,
)
    name = gene_of((_, kind, _) -> kind === :active, net_stoich, genes)
    name === nothing && return nothing
    gene = genes[name]
    active_name = Symbol("$(name).active")
    active = indices.unknowns[active_name]
    activating = first(change for (species, change) in net_stoich
        if species === active_name) > 0

    kind = activating ? "repression" : "activation"
    regulation = activating ? gene.repression : gene.activation
    aggregate = regulation.aggregate

    site, scale, offset = if !activating
        active, 1.0, 0.0
    elseif gene.unique
        active, -1.0, 1.0
    else
        indices.unknowns[Symbol("$(name).inactive")], 1.0, 0.0
    end

    SwitchingRate(
        Symbol("$(name).activity"),
        indices.parameters[Symbol("$(name).$(activating ? "activation" : "deactivation")")],
        regulators(indices, genes, gene, kind, regulation.slots, aggregate),
        aggregate,
        site,
        active,
        gene.unique ? 0 : indices.unknowns[Symbol("$(name).inactive")],
        scale, offset, affect,
    )
end

function SciML.constant_rate_jumps(definition::Definition, problem, system, ids, jumps)
    genes = Dict(g.name => g for g in definition.genes)
    indices = SciML.Indices(system, ids)
    map(jumps) do jump
        promoter = promoter_rate(jump, genes, definition, indices)
        promoter === nothing ?
        only(SciML.constant_rate_jumps(nothing, problem, system, ids, [jump])) :
        SciML.bounded_jump(promoter, promoter.affect)
    end
end

function regulation(
    genes::Dict{Symbol,<:ModelingToolkit.AbstractSystem};
    definition::Definition,
    t::Num,
)

    reference = Dict(g.name => regulator_name(g) for g in definition.genes)
    regulator_of(from) = species_variable(get(reference, from, from); t)
    effective = Dict(g.name => rate_values(g) for g in definition.genes)

    inactive(target::Gene) =
        if target.unique
            1 - genes[target.name].active
        else
            genes[target.name].inactive
        end

    aggregate(reg::Regulation, xs, name::Symbol, kind::String) =
        isempty(reg.slots) ? one(Num) :
        weighted(reg.aggregate, collect(xs),
            [make_parameter(Symbol("$(name).$(kind).$(s.from).w"), s.w) for s in reg.slots],
            name, kind)

    k_on(target::Gene) = (
        make_parameter(Symbol("$(target.name).activation"), target.base_rates.activation)
        * aggregate(target.repression,
            (  # ^ arguments and value go towards 0 as repression increases
                hill2(regulator_of(from), 1.0,
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
                hill2(regulator_of(from), 1.0,
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

    partition(target::Gene, kind::String, slots) = prod(
        (
            1 + (regulator_of(from) /
                make_parameter(Symbol("$(target.name).$(kind).$(from).at"), at)) ^
                abs(make_parameter(Symbol("$(target.name).$(kind).$(from).k"), k))
            for (; from, k, at) in slots
        );
        init = one(Num),
    )

    function p_thermodynamic(target::Gene)
        kon = make_parameter(Symbol("$(target.name).activation"),
            target.base_rates.activation)
        koff = make_parameter(Symbol("$(target.name).deactivation"),
            target.base_rates.deactivation)
        P = partition(target, "activation", target.activation.slots)
        Q = partition(target, "repression", target.repression.slots)
        (P - koff / (kon + koff)) / (P * Q)
    end

    p_promoter(target::Gene) = target.activation.aggregate === thermodynamic ?
        p_thermodynamic(target) : p_active(target)

    function promoter_activity(target)
        if !switching(target)
            p_promoter(target)
        elseif target.unique
            genes[target.name].active
        else
            active_counts = genes[target.name].active
            inactive_counts = genes[target.name].inactive
            total = active_counts + inactive_counts
            ifelse(total > 0, active_counts / total, 0.0)
        end
    end

    observed = [
        observed_variable(Symbol("$(target.name).activity"); t) ~
            promoter_activity(target)
        for target in definition.genes
    ]

    activation_rate(target::Gene) = k_on(target) * inactive(target)
    deactivation_rate(target::Gene) = k_off(target) * genes[target.name].active

    trigger_rate(target::Gene) =
        make_parameter(Symbol("$(target.name).trigger"), effective[target.name][:trigger])


    # Regulation for the whole network:
    reactions = [
        # For each gene...
        mapreduce(vcat, definition.genes, init = Reaction[]) do target::Gene
            vcat(
                if switching(target)
                    [
                        # ...activation (by tempering promoter deactivation)
                        annotate(Reaction(
                            deactivation_rate(target),
                            [genes[target.name].active],
                            target.unique ? nothing : [genes[target.name].inactive],
                            only_use_rate = true;
                            metadata = [:propensity_directions => [
                                (regulator_of(slot.from) => Int8(-1) for slot in target.activation.slots)...
                                genes[target.name].active => Int8(1)
                            ]]
                        ), :deactivation;
                            owner=target.name,
                            parameters=Dict(:rate => Symbol("$(target.name).deactivation")))

                        # ...repression (by tempering promoter activation)
                        annotate(Reaction(
                            activation_rate(target),
                            target.unique ? nothing : [genes[target.name].inactive],
                            [genes[target.name].active],
                            only_use_rate = true;
                            metadata = [:propensity_directions => [
                                (regulator_of(slot.from) => Int8(-1) for slot in target.repression.slots)...
                                (target.unique ? genes[target.name].active => Int8(-1) :
                                                 genes[target.name].inactive => Int8(1))
                            ]]
                        ), :activation;
                            owner=target.name,
                            parameters=Dict(:rate => Symbol("$(target.name).activation")))
                    ]
                else
                    transcript = first_transcript(target)
                    polymerases = species_variable(definition.polymerases; t)
                    catalysts = transcript === :proteins ?
                        [species_variable(definition.ribosomes; t)] : []
                    held = transcript === :elongations
                    directions = [
                        polymerases => Int8(1)
                        (catalyst => Int8(1) for catalyst in catalysts)...
                        (regulator_of(slot.from) => Int8(1) for slot in target.activation.slots)...
                        (regulator_of(slot.from) => Int8(-1) for slot in target.repression.slots)...
                    ]
                    annotate(Reaction(
                        trigger_rate(target) * p_promoter(target),
                        [polymerases, catalysts...],
                        [getproperty(genes[target.name], transcript), catalysts...,
                            (held ? () : (polymerases,))...];
                        metadata = [:propensity_directions => directions]
                    ), :trigger;
                        owner=target.name,
                        parameters=Dict(:rate => Symbol("$(target.name).trigger")))
                end,
                # ...repression (by proteolysis)
                map(target.proteolysis.slots) do (; from, k)
                    proteases = regulator_of(from)
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
                        gene_link=:proteolysis,
                        parameters=Dict(:rate => Symbol("$(target.name).proteolysis.$(from).k")))
                end
            )
        end

        # Additionally, we add arbitrary mass-action reactions as specified.
        # Bidirectional pairs are broken up, and reactions are only included if
        # their rate is nonzero.
        [
            annotate(Reaction(
                make_parameter(Symbol("reaction.$(name).k⁺"), k⁺),
                regulator_of.(keys(from.counts)),
                regulator_of.(keys(to.counts)),
                collect(values(from.counts)),
                collect(values(to.counts)),
            ), :reaction;
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
                regulator_of.(keys(to.counts)),
                regulator_of.(keys(from.counts)),
                collect(values(to.counts)),
                collect(values(from.counts)),
            ), :reaction;
                name,
                direction=:reverse,
                parameters=Dict(
                    :k⁻ => Symbol("reaction.$(name).k⁻")
                ))
            for (; name, from, k⁻, to) in definition.reactions
            if k⁻ > 0.0
        ]
    ]
    (; reactions, observed)
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
    BracketData{Vector{Float64}, Vector{Int}}(
        [e ? 0.0 : 0.1 for e in exact],
        [e ? 0 : 25 for e in exact],
        [e ? 0 : 4 for e in exact],
    )
end

hybrid(definition::Definition; exact = RSSACR(), dt = Inf, nc = nothing,
        policy = nothing, kwargs...) =
    JumpProcesses.HybridTau(exact,
        blending_policy(policy === nothing ? definition : Val(Symbol(policy)), nc),
        dt; kwargs...)

blending_policy(definition::Definition, nc) =
    nc === nothing ? blending_policy(definition) : JumpProcesses.CriticalBlend(nc)

blending_policy(::Val{:AlwaysLeap}, _) = JumpProcesses.AlwaysLeap()
blending_policy(::Val{:LinearBlend}, _) = JumpProcesses.LinearBlend()
blending_policy(::Val{:CriticalBlend}, nc) =
    JumpProcesses.CriticalBlend(something(nc, 10))

blending_policy(definition::Definition) =
    any(switching, definition.genes) ? JumpProcesses.CriticalBlend(2) :
    JumpProcesses.AlwaysLeap()

resolve_method(method::JumpProcesses.AbstractAggregatorAlgorithm, system, definition) =
    method

resolve_method(method::Symbol, system, definition) =
    method === :HybridTau ? hybrid(definition) : pick_method(system; method)

function resolve_method(method::AbstractDict{Symbol}, system, definition)
    name = Symbol(get(method, :name, "default"))
    name === :HybridTau ||
        error("only \"HybridTau\" takes options; got $(name)")
    options = (; (k => v for (k, v) in method if k !== :name && k !== :exact)...)
    hybrid(definition;
        exact = pick_method(system; method = Symbol(get(method, :exact, "RSSACR")))(),
        options...)
end

# `bounds = false` skips deriving the propensity directions, which only the
# symbolic path needs -- the fast path carries its own bounds.
function aggregator_options(algorithm::JumpProcesses.HybridTau, reaction_system,
        jump_system, definition; bounds = true)
    if any(switching, definition.genes) &&
       algorithm.policy isa JumpProcesses.CriticalBlend
        nc = algorithm.policy.nc
        if all(gene -> gene.unique, definition.genes)
            minimum(nc) <= 1 &&
                @warn "HybridTau with `active` genes and CriticalBlend(nc = $nc) leaps transcription with the promoter gate held fixed over the window. Pass CriticalBlend(2) or larger, or use `V1.hybrid(definition)`."
        else
            @warn "HybridTau with `unique = false` genes and CriticalBlend(nc = $nc): the promoter copy number is dynamic, so the gate is only kept exact while it stays below $nc. Raise nc above the largest promoter copy number you seed, or use `V1.hybrid(definition)`."
        end
    end
    aggregator_options(typeof(algorithm.exact), reaction_system, jump_system, definition;
        bounds)
end

function aggregator_options(algorithm, reaction_system, jump_system, definition;
    bounds = true)
    algorithm in (RSSA, RSSACR, TauSplitting) || return (nothing, (;))
    all(definition.genes) do gene
        all(slot -> slot.k <= 0.0, gene.activation.slots) &&
            all(slot -> slot.k <= 0.0, gene.repression.slots)
    end || error("$(nameof(algorithm)) requires every activation/repression `k` to be non-positive so that propensities are monotone in the regulator counts; this definition has a positive `k`.")

    bounds || return nothing, algorithm in (RSSA, RSSACR) ?
        (; bracket_data = promoter_bracket_data(jump_system)) : (;)

    if any(!switching, definition.genes)
        all(definition.genes) do gene
            activators = Set(slot.from for slot in gene.activation.slots)
            repressors = Set(slot.from for slot in gene.repression.slots)
            isdisjoint(activators, repressors)
        end || error("$(nameof(algorithm)) with compilation=:symbolic does not support a species both activating and repressing the same equilibrium promoter; use compilation=:fast")
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
    directions = [propensity_directions(rx, species_index) for rx in constant_reactions]

    algorithm in (RSSA, RSSACR) || return directions, (;)
    directions, (; bracket_data = promoter_bracket_data(jump_system))
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

`"HybridTau"` selects the hybrid exact/tau-leaping aggregator, with a blending
policy chosen from whether any gene keeps `active`. It may instead be given as a
JSON object
```
{"name": "HybridTau", "exact": <exact>, "epsilon": <epsilon>, "dt": <dt>,
 "policy": <policy>, "nc": <nc>}
```
where `<exact>` names the inner exact aggregator (`"RSSACR"` by default),
`<epsilon>` is the tau-selection tolerance, and `<dt>` an upper bound on the
leap window. Setting `<epsilon>` to `null` selects fixed steps of `<dt>` instead
of adaptive tau selection. `<policy>` overrides the blending policy with one of
`"AlwaysLeap"`, `"CriticalBlend"` or `"LinearBlend"`, and `<nc>` sets the
critical threshold.

# Specification

V1 models are specified in JSON as `{"{regulation/v1}": <definition>}` where
`<definition>` specifies a [`Definition`](@ref) as described there.

For an example, see `examples/specification/simple.schedule.json`.
"""
function build end

build(specification::AbstractDict{Symbol}) = build(
    cast(Definition, specification),
    method = method_specification(get(specification, :method, "default")),
    compilation = Symbol(get(specification, :compilation, "fast")),
)

method_specification(method) = Symbol(method)
method_specification(method::AbstractDict{Symbol}) = method

function build(definition::Definition;
    method::Union{Symbol, AbstractDict{Symbol},
        JumpProcesses.AbstractAggregatorAlgorithm} = :default,
    compilation::Symbol = :fast)

    let names = [rxn.name for rxn in definition.reactions]
        duplicates = unique([n for n in names if count(==(n), names) > 1])
        isempty(duplicates) ||
            error("reaction names must be unique; duplicated: ", join(duplicates, ", "))
    end

    for gene in definition.genes
        chosen = species_of(gene)
        isempty(chosen) &&
            error("gene $(gene.name) must keep at least one species")
        unknown = setdiff(chosen, default_species(gene.base_rates))
        isempty(unknown) ||
            error("gene $(gene.name) cannot keep $(join(unknown, ", ")); available: ",
                join(default_species(gene.base_rates), ", "))
        switching(gene) || gene.unique ||
            error("gene $(gene.name) without `active` does not support `unique=false`")
        (gene.base_rates.processing === nothing) ==
        (gene.base_rates.premrna_decay === nothing) ||
            error("gene $(gene.name) must set `processing` and `premrna_decay` together")
        isempty(gene.proteolysis.slots) || :proteins in chosen ||
            error("gene $(gene.name) cannot be regulated by proteolysis without `proteins`")
    end

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
        )
        for g in definition.genes
    )
    (; reactions, observed) = regulation(genes; definition, t)
    @named reaction_system = ReactionSystem(
        reactions,
        t;
        observed,
        systems=collect(values(genes)),
        initial_conditions=Dict(
            getproperty(genes[g.name], kind) => value
            for g in definition.genes
            for (kind, value) in rate_values(g)
            if switching(g) || kind !== :trigger
        )
    )
    reaction_system = complete(reaction_system)
    jump_system = complete(jump_model(reaction_system))
    algorithm = resolve_method(method, reaction_system, definition)

    compilation in (:symbolic, :fast) ||
        error("compilation must be :symbolic or :fast, got $(compilation)")

    directions, options = aggregator_options(
        algorithm, reaction_system, jump_system, definition; bounds = compilation !== :fast)
    source = compilation === :fast ? definition : directions

    Models.Wrapped(;
        definition,
        model = Models.Wrapped(
            definition = reaction_system,
            model = SciML.JumpModel(jump_system,
                algorithm isa JumpProcesses.AbstractAggregatorAlgorithm ? algorithm : algorithm(),
                source; options...),
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
            gene.name, gene.base_rates, gene.unique, gene.species,
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
