"""
Contains `Instant` models to be used in extraction schemes (observation models).
"""
module Resampling

using ..Models: Model, FlatState
import ..Specifications

using Distributions

"""
    ResampleEachBinomial <: Model{FlatState}

Retain each molecule independently with probability `p`, and drop the rest.

This replaces each species' count *n* with a value sampled from from a
binomial distribution with parameters *n* and `p`.

# Specification

In JSON, `ResampleEachBinomial` is specified as a JSON object
```
{"{resample-each-binomial}": <p>}
```
where `<p>` is a unit-range JSON number specifying the per-molecule retain
probability.
```
{"{resample-each-binomial}": {"p": <p>, "except": "<regex>"}}
```
"""
@kwdef struct ResampleEachBinomial <: Model{FlatState}
    p::Float64
    except::Union{Regex, Nothing} = nothing

    function ResampleEachBinomial(p, except = nothing)
        0.0 ≤ p ≤ 1.0 || error("invalid probability")
        new(p, except isa AbstractString ? Regex(except) : except)
    end
end

ResampleEachBinomial(spec::AbstractDict{Symbol}) =
    ResampleEachBinomial(spec[:p], get(spec, :except, nothing))

Specifications.constructor(::Val{Symbol("resample-each-binomial")}) =
    ResampleEachBinomial

function (f!::ResampleEachBinomial)(x::FlatState, _Δt::Float64; _...)
    for (name, count) in x.counts
        f!.except !== nothing && occursin(f!.except, String(name)) && continue
        x.counts[name] = rand(x.randomness, Binomial(floor(Int, count), f!.p))
    end
    x
end

"""
    ResampleEachAccumulate <: Model{FlatState}

Drop, retain or multiply each molecule independently with specified
probabilities `ps`.

The species are treated an exchangeable. For each molecule, the number of copies
it should be replaced by is sampled independently with probabilities given by
the `Vector` `ps`, where each `ps[i]` defines the probability of resulting in
`i` - 1 copies.

# Specification

In JSON, `ResampleEachAccumulate` is specified as a JSON object
```
{"{resample-each-accumulate}": <ps>}
```
where `<ps>` is a JSON array of unit-range JSON numbers that sum to 1 and
specify the per-molecule copy probabilities as defined above.
"""
@kwdef struct ResampleEachAccumulate <: Model{FlatState}
    ps::Vector{Float64}

    function ResampleEachAccumulate(ps::AbstractVector)
        ps = Specifications.cast(Vector{Float64}, ps)
        isprobvec(ps) ? new(ps) : error("invalid probability vector")
    end
end

Specifications.constructor(::Val{Symbol("resample-each-accumulate")}) =
    ResampleEachAccumulate

function (f!::ResampleEachAccumulate)(x::FlatState, _Δt::Float64; _...)
    map!(values(x.counts)) do count
        sum(
            enumerate(rand(x.randomness, Multinomial(floor(Int, count), f!.ps)))
        ) do (i, count′)
            (i - 1) * count′
        end
    end
    x
end

"""
    ResampleHypergeometric <: Model{FlatState}

Retain (sample without replacement) `n` molecules with equal probability, and
drop the rest.

This replaces the per-species counts with values sampled from a multivariate
hypergeometric distribution parametrized by `n` and the current counts. If `n`
exceeds the total count, the per-species counts are left unchanged.

# Specification

In JSON, `ResampleHypergeometric` is specified as a JSON object
```
{"{resample-hypergeometric}": <n>}
```
where `<n>` is a JSON number specifying the total count of molecules to retain.
"""
@kwdef struct ResampleHypergeometric <: Model{FlatState}
    n::Int

    ResampleHypergeometric(n) =
        n ≥ 0 ? new(n) : error("invalid sample size")
end

Specifications.constructor(::Val{Symbol("resample-hypergeometric")}) =
    ResampleHypergeometric

function (f!::ResampleHypergeometric)(x::FlatState, _Δt::Float64; _...)
    reservoir = sum(floor(Int, count) for count in values(x.counts))
    todo = f!.n
    reservoir ≤ todo && return x

    map!(values(x.counts)) do count
        count = floor(Int, count)
        reservoir -= count
        count = rand(x.randomness, Hypergeometric(count, reservoir, todo))
        todo -= count
        count
    end

    x
end

"""
    ResampleMultinomial <: Model{FlatState}

Sample `n` molecules (with replacement) with equal probability.

This replaces the per-species counts with values sampled from a multinomial
distribtution parametrized by `n` and the current counts. If the per-species
counts are all zero, they are left unchanged.

# Specification

In JSON, `ResampleMultinomial` is specified as a JSON object
```
{"{resample-multinomial}": <n>}
```
where `<n>` is a JSON number specifying the total count of molecules to sample.
"""
@kwdef struct ResampleMultinomial <: Model{FlatState}
    n::Int

    ResampleMultinomial(n) = n ≥ 0 ? new(n) : error("invalid sample size")
end

Specifications.constructor(::Val{Symbol("resample-multinomial")}) =
    ResampleMultinomial

function (f!::ResampleMultinomial)(x::FlatState, _Δt::Float64; _...)
    total = sum(values(x.counts))
    total > 0 || return x
    counts = rand(x.randomness, Multinomial(f!.n, values(x.counts) ./ total))
    for (kind, count) in zip(keys(x.counts), counts)
        x.counts[kind] = count
    end
    x
end

@kwdef struct ResampleDirichletMultinomial <: Model{FlatState}
    n::Int

    ResampleDirichletMultinomial(n) =
        n ≥ 0 ? new(n) : error("invalid sample size")
end

Specifications.constructor(::Val{Symbol("resample-Dirichlet-multinomial")}) =
    ResampleDirichletMultinomial

function (f!::ResampleDirichletMultinomial)(x::FlatState, _Δt::Float64; _...)
    counts = rand(
        x.randomness,
        DirichletMultinomial(f!.n, collect(values(x.counts))),
    )
    for (kind, count) in zip(keys(x.counts), counts)
        x.counts[kind] = count
    end
    x
end

"""
    ResampleTargetMeanEachBinomial <: Model{FlatState}

Retain each molecule independently with the same probability such that the
resulting expected total count is `target`.

This replaces each species ``i``'s count ``n_i`` with a value sampled from from
a binomial distribution ``B(n_i, p)`` with
``p = \\texttt{target} / \\sum_i{n_i}``.

# Specification

In JSON, `ResampleTargetMeanEachBinomial` is specified as a JSON object
```
{"{resample-target-mean-each-binomial}": <target>}
```
where `<target>` is a JSON number specifying the target expected total molecule
count.
"""
@kwdef struct ResampleTargetMeanEachBinomial <: Model{FlatState}
    μ::Float64

    ResampleTargetMeanEachBinomial(μ) =
        μ ≥ 0.0 ? new(μ) : error("invalid sample mean")
end

Specifications.constructor(
    ::Val{Symbol("resample-target-mean-each-binomial")}
) = ResampleTargetMeanEachBinomial

(f!::ResampleTargetMeanEachBinomial)(x::FlatState, _Δt::Float64; _...) =
    ResampleEachBinomial(p = min(1.0, f!.μ / sum(values(x.counts))))(x, 0.0)

end
