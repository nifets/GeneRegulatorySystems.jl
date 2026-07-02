"""
Contains components to build random `Models.V1`-based `SciML.JumpModel`s for
gene regulation with an SKG-linked network.

These models can be constructed by [`build`](@ref) from a [`Definition`](@ref)
or from the corresponding JSON specification. The link structure is sampled
from a stochastic Kronecker graph ("SKG") network model, and link types and
parameters are sampled independently for each link.
"""
module KroneckerNetworks

using ...GeneRegulatorySystems: σ, logit
using ..Models: Models, V1
using ..Sampling: Nonnegative, BaseRatesTemplate
import ..Specifications

using Random
using SparseArrays

using Distributions
using Kronecker
using Roots

function shadow_task_randomness!(callback, imposed::Random.Xoshiro)
    hidden = Random.getstate(Random.default_rng())
    Random.setstate!(Random.default_rng(), Random.getstate(imposed))
    result = callback(Random.default_rng())
    Random.setstate!(imposed, Random.getstate(Random.default_rng()))
    Random.setstate!(Random.default_rng(), hidden)
    result
end

"""
    NetworkTemplate

Abstract supertype of network templates.

These define how to sample a regulatory network of a specific link type, which
includes both link existence and parameters. Existence is sampled from a
`Bernoulli` distribution independently for each potential link according to a
stochastic `adjacency` matrix of unit-range values. Link parameters are
determined differently between `NetworkTemplate` subtypes, but the `adjacency`
specification is equivalent.

# Specification

In JSON, any `<:NetworkTemplate` is specified by a JSON object
`{"adjacency": <adjacency>, ...}` with subtype-specific other mappings.

The `<adjacency>` must be specified in one of the following forms:
- In the simplest case, the values are given directly. The specification must
  then be a JSON array of JSON arrays, each specifying a single row of the
  matrix. For example:
  ```
  [
      [0.7, 0.6],
      [0.4, 0.2]
  ]
  ```
- Otherwise, `<adjacency>` may be a JSON array of such matrices. They will then
  be Kronecker-multiplied to give an expanded adjacency matrix.
- Finally, `<adjacency>` may be defined as a Kronecker power. In this case, it
  is specified as a JSON object
  ```
  {
      "initiator": <initiator>,
      "power": <power>,
      "expected_links_per_gene": <links>
  }
  ```
  where `<initiator>` is a matrix specification like in the first case and
  `<power>` is a JSON (integer) number that specifies to which power the
  initiator should be raised by an iterated Kronecker product. This therefore
  specifies a stochastic Kronecker graph that shapes the regulatory network. The
  `"expected_links_per_gene"` specification is optional; if present, `<links>`
  must be a JSON number, and all entries of the final expanded adjacency matrix
  will then be shifted on the logit scale such that the expected total number of
  sampled links is equal to `<links>`. This can be used to control the density
  of the regulatory network without changing its nesting structure.
"""
abstract type NetworkTemplate end

"""
    ActivationNetworkTemplate <: NetworkTemplate

Defines how to sample a network of transcriptional activation.

# Specification

In JSON, an `ActivationNetworkTemplate` is specified by a JSON object
```
{
    "adjacency": <adjacency>,
    "at": <at>
}
```
where `<at>` specifies a [`Nonnegative{<:UnivariateDistribution}`](@ref) prior
for each link's parameter and `<adjacency>` is specified as described for
[`NetworkTemplate`](@ref).
"""
@kwdef struct ActivationNetworkTemplate <: NetworkTemplate
    adjacency::AbstractMatrix
    at::Nonnegative{UnivariateDistribution}
    k::UnivariateDistribution = Dirac(-1.0)
    # TODO: Should we support aggregation here?
end

"""
    RepressionNetworkTemplate <: NetworkTemplate

Defines how to sample a network of transcriptional repression.

# Specification

In JSON, an `RepressionNetworkTemplate` is specified by a JSON object
```
{
    "adjacency": <adjacency>,
    "at": <at>
}
```
where `<at>` specifies a [`Nonnegative{UnivariateDistribution}`](@ref) prior for
each link's parameter and `<adjacency>` is specified as described for
[`NetworkTemplate`](@ref).
"""
@kwdef struct RepressionNetworkTemplate <: NetworkTemplate
    adjacency::AbstractMatrix
    at::Nonnegative{UnivariateDistribution}
    k::UnivariateDistribution = Dirac(-1.0)
    # TODO: Should we support aggregation here?
end

"""
    ProteolysisNetworkTemplate <: NetworkTemplate

Defines how to sample a network of proteolytic repression.

# Specification

In JSON, an `ProteolysisNetworkTemplate` is specified by a JSON object
```
{
    "adjacency": <adjacency>,
    "k": <at>
}
```
where `<k>` specifies a [`Nonnegative{<:UnivariateDistribution}`](@ref) prior
for each link's parameter and `<adjacency>` is specified as described for
[`NetworkTemplate`](@ref).
"""
@kwdef struct ProteolysisNetworkTemplate <: NetworkTemplate
    adjacency::AbstractMatrix
    k::Nonnegative{UnivariateDistribution}
end

"""
    Template

Defines how to sample a [`V1.Definition`](@ref Models.V1.Definition) from given
base rate and regulation network templates.

[`build`](@ref) uses these definitions to generate genes independently
according to the specified base rate distribution and then adds links for the
`activation`, `repression` and `proteolysis` regulatory networks. The existence
of those links is determined by a stochastic Kronecker graph, and their
parameters are chosen independently according to the specified distributions.

# Specification

In JSON, a `KroneckerNetworks.Template` is specified as a JSON object
```
{
    "base_rates": <base_rates>,
    "activation": <activation>,
    "repression": <repression>,
    "proteolysis": <proteolysis>
}
```
where only `<base_rates>` is required and the other mappings are optional.

These define the priors that `build` uses to pick the generated genes' base
rates as well as the existence and parameters of inter-gene regulatory links.
The regulatory networks for `activation`, `repression` and `proteolysis` are
specified separately, and `build` will implicitly overlay them.

`<base_rates>` specifies a [`BaseRatesTemplate`](@ref).

If present, `<activation>` must specify a [`ActivationNetworkTemplate`](@ref).

If present, `<repression>` must specify a [`RepressionNetworkTemplate`](@ref).

If present, `<proteolysis>` must specify a [`ProteolysisNetworkTemplate`](@ref).

(Note that the JSON object specifying the `Template` will typically also
contain a `"seed"` mapping if it is specified as part of a
`KroneckerNetworks.Definition`.)
"""
@kwdef struct Template
    base_rates::BaseRatesTemplate
    activation::Union{Some{ActivationNetworkTemplate}, Nothing} = nothing
    repression::Union{Some{RepressionNetworkTemplate}, Nothing} = nothing
    proteolysis::Union{Some{ProteolysisNetworkTemplate}, Nothing} = nothing
    count::Int =
        first(size(@something(activation, repression, proteolysis).adjacency))
    prefix::String = ""
    approximate::Union{Some{Bool}, Nothing} = nothing
end

"""
    Definition

Defines how to construct a [`SciML.JumpModel`](@ref Models.SciML.JumpModel)
(via [`V1`](@ref Models.V1)) by sampling from a contained [`Template`](@ref)
with a fixed `seed`.

# Specification

In JSON, a `KroneckerNetworks.Definition` is specified as a JSON object
```
{
    "seed": <seed>,
    <template>...
}}
```
where `<seed>` is a JSON string that fixes a specific instance to sample from
the [`Template`](@ref) specified by the JSON object containing `<template>...`.
"""
@kwdef struct Definition
    seed::String
    template::Template
end

Models.describe(definition::Definition) = Models.Label("\
    'regulation/kronecker' definition with seed '$(definition.seed)' \
    and a $(definition.template.count)×$(definition.template.count) expansion\
")

gene_name(i; n, prefix) = Symbol("$prefix$(lpad(i, ndigits(n), '0'))")

regulator(
    template::Union{ActivationNetworkTemplate, RepressionNetworkTemplate};
    from::Symbol,
    randomness::AbstractRNG,
) = V1.HillRegulator(
    at = rand(randomness, template.at),
    k = rand(randomness, template.k);
    from
)

regulator(
    template::ProteolysisNetworkTemplate;
    from::Symbol,
    randomness::AbstractRNG,
) = V1.DirectRegulator(k = rand(randomness, template.k); from)

regulation(::ActivationNetworkTemplate; slots::Vector{V1.HillRegulator}) =
    V1.Activation(; slots)

regulation(::RepressionNetworkTemplate; slots::Vector{V1.HillRegulator}) =
    V1.Repression(; slots)

regulation(::ProteolysisNetworkTemplate; slots::Vector{V1.DirectRegulator}) =
    V1.Proteolysis(; slots)

function regulations(
    template::NetworkTemplate;
    n,
    approximate,
    prefix,
    randomness::AbstractRNG,
)
    size(template.adjacency) == (n, n) ||
        error("invalid expansion: Kronecker product must have size ($n, $n)")
    isprob(template.adjacency) ||
        error("invalid initiator: must be probabilities")

    if approximate
        realized =
            shadow_task_randomness!(randomness) do randomness
                fastsample(template.adjacency)
            end
        map(eachcol(realized)) do column
            regulation(template, slots = [
                regulator(template, from = gene_name(i; n, prefix); randomness)
                for i in first(findnz(column))
            ])
        end
    else
        map(eachcol(template.adjacency)) do column
            regulation(template, slots = [
                regulator(template; from = gene_name(i; n, prefix), randomness)
                for (i, cell) in enumerate(column)
                if rand(randomness) < cell
            ])
        end
    end
end

adjacency(initiator::AbstractMatrix; k::Int) =
    k == 1 ? initiator : kronecker(initiator, k)
adjacency(factors::AbstractVector) =
    length(factors) == 1 ? only(factors) : kronecker(factors...)

adjusted(initiator::AbstractMatrix, α::Float64) = σ.(α .+ logit.(initiator))
function adjusted(initiator::AbstractMatrix; k::Int, 𝔼links::Float64)
    n = size(initiator)[1]
    equivalent = 𝔼links ^ (1 / k)
    0 < equivalent < n || error("impossible expected link number target")
    target = n * equivalent

    residual(α) = sum(adjusted(initiator, α)) - target
    α₀ = find_zero(residual, 0.0)

    adjusted(initiator, α₀)
end

factor(xs::AbstractVector) = Float64.(mapreduce(adjoint, vcat, xs))

function Specifications.cast(::Type{AbstractMatrix}, xs::AbstractVector; _...)
    factors = factor.(xs)
    result = adjacency(factors)
    issquare(result) || error("expanded adjacency must be square")
    result
end

function Specifications.cast(
    T::Type{AbstractMatrix},
    x::AbstractDict{Symbol};
    _...,
)
    initiator = factor(x[:initiator])
    issquare(initiator) || error("initiator must be square")
    k = x[:power]
    if haskey(x, :expected_links_per_gene)
        𝔼links = x[:expected_links_per_gene]
        initiator = adjusted(initiator; k, 𝔼links)
    end
    adjacency(initiator; k)
end

function Base.rand(randomness::AbstractRNG, template::Template)
    n = template.count
    arguments = (;
        approximate = @something(template.approximate, n > 16),
        template.prefix,
        randomness,
    )

    activations =
        if isnothing(template.activation)
            [V1.Activation() for _ in 1:n]
        else
            regulations(something(template.activation); n, arguments...)
        end

    repressions =
        if isnothing(template.repression)
            [V1.Repression() for _ in 1:n]
        else
            regulations(something(template.repression); n, arguments...)
        end

    proteolyses =
        if isnothing(template.proteolysis)
            [V1.Proteolysis() for _ in 1:n]
        else
            regulations(something(template.proteolysis); n, arguments...)
        end

    V1.Definition(;
        genes = [
            V1.Gene(
                name = gene_name(i; n, template.prefix),
                base_rates = rand(randomness, template.base_rates);
                activation,
                repression,
                proteolysis,
            )
            for (i, activation, repression, proteolysis) in
                zip(1:n, activations, repressions, proteolyses)
        ],
    )
end

"""
    build(specification::AbstractDict{Symbol})
    build(definition::Definition; method::Symbol = :default)

Construct a `SciML.JumpModel` from a [`Definition`](@ref).

When interpreting a JSON specification, this function (in its first form) is
called to construct a concrete regulation model on encountering a
`{"{regulation/kronecker}": {...}}` literal. It will first destructure the
parsed JSON into a `Definition` and then proceed from there.

The result is constructed by first making a concrete [`V1.Definition`](@ref)
from [`definition.template`](@ref Template) to obtain the corresponding `Model`
and then wrapping that up in a [`Models.Wrapped`](@ref) with `definition`. This
will result in the following stack of abstractions:
- [`SciML.JumpModel`](@ref Models.SciML.JumpModel), specified by a
- `Catalyst.ReactionSystem`, specified by a
- [`V1.Definition`](@ref), specified by
- `definition::`[`KroneckerNetworks.Definition`](@ref), potentially specified by
- `specification`

# Specification

Kronecker-linked networks are specified in JSON as
`{"{regulation/kronecker}": <definition>}` where `<definition>` specifies a
[`Definition`](@ref) as described there.

For an example, see `examples/specification/kronecker.schedule.json`.
"""
function build end

build(specification::AbstractDict{Symbol}) = build(
    # Pick a specific model instance by affixing the randomness:
    Definition(
        seed = specification[:seed],
        template = Specifications.cast(Template, specification),
    ),
    method = Symbol(get(specification, :method, "default"))
)

build(definition::Definition; method::Symbol = :default) = Models.Wrapped(
    model = V1.build(
        # Deterministically fill in the template to create a concrete
        # V1.Definition from it:
        rand(Xoshiro(definition.seed), definition.template);
        method,
    );
    definition,
)

Specifications.constructor(::Val{Symbol("regulation/kronecker")}) = build

end
