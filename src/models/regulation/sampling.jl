module Sampling

import ...Specifications: cast
using ..Models: V1

using Distributions

using Random

"""
    struct Nonnegative{T <: UnivariateDistribution}

Wraps a Distributions.jl `UnivariateDistribution`, but requires it to have
nonnegative support.

# Specification

In JSON, a `Nonnegative{T <: UnivariateDistribution}` is specified as one of the
following:
- In the simple case, it is specified as a JSON number, which will result in
  a `Dirac` distribution parameterized by that number.
- Otherwise, it is specified as a JSON Array `[<distribution>, <parameters>...]`
  where `<distribution>` is a JSON string naming a nonnegative
  `UnivariateDistribution`, and `<parameters>...` will be passed to its
  constructor. For example: `["LogNormal", 2.0, 1.0]`.
"""
struct Nonnegative{T <: UnivariateDistribution}
    inner::T
end

function cast(::Type{Nonnegative{T}}, x; _...) where {T}
    result = cast(T, x)
    minimum(result) ≥ 0.0 ||
        error("distribution must have nonnegative support")
    Nonnegative{T}(result)
end

cast(::Type{<:UnivariateDistribution}, x::Real; _...) = Dirac(x)

function cast(::Type{<:UnivariateDistribution}, xs::AbstractVector; _...)
    T = getfield(Distributions, Symbol(first(xs)))
    T <: UnivariateDistribution || error("not a UnivariateDistribution")
    T((identity.(x) for x in xs[2:end])...)
end

cast(::Type{MultivariateDistribution}, xs::AbstractVector; _...) =
    if !isempty(xs) && first(xs) isa AbstractString
        T = getfield(Distributions, Symbol(first(xs)))
        T <: MultivariateDistribution || error("not a MultivariateDistribtuion")
        T((identity.(x) for x in xs[2:end])...)
    else
        Product(Dirac.(xs))
    end

Base.rand(randomness::AbstractRNG, d::Nonnegative{<:UnivariateDistribution}) =
    rand(randomness, d.inner)

"""
    BaseRatesTemplate

Defines how to sample [`V1.BaseRates`](@ref).

Base rates are sampled independently for each kind of rate.

# Specification

In JSON, a `BaseRatesTemplate` is specified as a JSON object
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
where each `<...>` specifies a [`Nonnegative{<:UnivariateDistribution}`](@ref)
from which that rate should be sampled from. `"processing"` and
`"premrna_decay"` are optional, as in [`V1.BaseRates`](@ref).
"""
@kwdef struct BaseRatesTemplate
    activation::Nonnegative{UnivariateDistribution}
    deactivation::Nonnegative{UnivariateDistribution}
    trigger::Nonnegative{UnivariateDistribution}
    transcription::Nonnegative{UnivariateDistribution}
    translation::Nonnegative{UnivariateDistribution}
    abortion::Nonnegative{UnivariateDistribution}
    mrna_decay::Nonnegative{UnivariateDistribution}
    protein_decay::Nonnegative{UnivariateDistribution}
    processing::Union{Nothing, Nonnegative{UnivariateDistribution}} = nothing
    premrna_decay::Union{Nothing, Nonnegative{UnivariateDistribution}} = nothing
end

sample(randomness::AbstractRNG, ::Nothing) = nothing
sample(randomness::AbstractRNG, d::Nonnegative) = rand(randomness, d)

Base.rand(randomness::AbstractRNG, template::BaseRatesTemplate) =
    V1.BaseRates(; (
        field => sample(randomness, getfield(template, field))
        for field in fieldnames(V1.BaseRates)
    )...)

end
