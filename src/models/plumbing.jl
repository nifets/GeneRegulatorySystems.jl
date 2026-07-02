"""
Contains miscellaneous `Model`s, including simple `Instant` interventions.

These are likely not useful on their own, but they are mostly intended for
connecting [regulation models](@ref "Regulation") within a `Schedule`.
"""
module Plumbing

import ....GeneRegulatorySystems
using ..Models: Models, Model, Instant, FlatState
import ..Specifications

"""
    Pass <: Instant{Any}

Do nothing for an instant.

# Specification

Specified in JSON as `{"{pass}": true}`.

# Invocation

    (::Pass)(x, _Δt::Float64; _...) = x
"""
struct Pass <: Instant{Any} end

Pass(_) = Pass()

Specifications.constructor(::Val{:pass}) = Pass

(::Pass)(x, _Δt::Float64; _...) = x

"""
    Wait <: Model{FlatState}

Do nothing for a while.

# Specification

Specified in JSON as `{"{wait}": true}`.

# Invocation

    (::Wait)(x, Δt::Float64; _...)

Advance simulation time by `Δt` but do not otherwise change state.
"""
struct Wait <: Model{FlatState} end

Wait(_) = Wait()

Specifications.constructor(::Val{:wait}) = Wait

function (::Wait)(x::FlatState, Δt::Float64; _...)
    x.t += Δt
    x
end

@doc raw"""
    Filter <: Instant{FlatState}

Instantly remove all species counts that have names not matching a pattern.

# Specification

Specified in JSON as `{"{filter}": "<regex>"}` where `<regex>` is any JSON
string that represents a valid Julia regular expression
([in PCRE2 syntax](https://www.pcre.org/current/doc/html/pcre2syntax.html)).
`\` characters must be escaped by `\\`. For example, `\\.(pre)?mrnas$` will only
retain species with names that end in either `.premrnas` or `.mrnas`.

# Invocation

    (f!::Filter)(x::FlatState, _Δt::Float64; _...)

Remove all mappings from `x.counts` whose key does not match `f!.kinds`.
"""
struct Filter <: Instant{FlatState}
    kinds::Regex
end
Filter(s::AbstractString) = Filter(Regex(s))

Specifications.constructor(::Val{:filter}) = Filter

function (f!::Filter)(x::FlatState, _Δt::Float64; _...)
    for key in keys(x.counts)
        if isnothing(match(f!.kinds, String(key)))
            delete!(x.counts, key)
        end
    end
    x
end

flatten(xs::AbstractDict{Symbol}; T = Any) =
    mapreduce(merge, xs) do (key, value)
        if value isa AbstractDict{Symbol}
            Dict{Symbol, T}(
                Symbol("$(key).$(key′)") => value′
                for (key′, value′) in flatten(value)
            )
        else
            Dict{Symbol, T}(key => value)
        end
    end

"""
    Adjust <: Instant{FlatState}

Instantly adjust counts.

# Specification

Specified in JSON as one of
- `{"{set}": <adjustment>}`, setting counts to constant values
- `{"{add}": <adjustment>}`, adding constant values to counts
- `{"{multiply}": <adjustment>}`, multiplying constant values with counts
where `<adjustment>` is a JSON object mapping species names to nonnegative
adjustment values. Nested JSON objects will be flattened, joining keys on `"."`
such that for example

    {"{multiply}": {
        "a": {
            "proteins": 0.5,
            "mnras": 0.5
        }
    }}

is equivalent to:

    {"{multiply}": {
        "a.proteins": 0.5,
        "a.mnras": 0.5
    }}

# Invocation

    (f!::Adjust)(x::FlatState, _Δt::Float64; _...)

Instantly apply flattened `f!.adjustment` to `x.counts` using `f!.adjust`.

Adjusted values will be rounded down to `Int`s (which is relevant if `f!.adjust`
is `*`).
"""
struct Adjust <: Instant{FlatState}
    adjust::Function
    adjustment::Dict{Symbol, <:Real}

    function Adjust(adjust, adjustment)
        all(values(adjustment) .≥ zero(valtype(adjustment))) ||
            error("adjustment must be nonnegative")
        new(adjust, adjustment)
    end
end

adder(counts::AbstractDict{Symbol}) = Adjust(+, flatten(counts, T = Int))
multiplier(counts::AbstractDict{Symbol}) = Adjust(*, flatten(counts, T = Real))
setter(counts::AbstractDict{Symbol}) =
    Adjust(last ∘ Pair, flatten(counts, T = Int))

Specifications.constructor(::Val{:set}) = setter
Specifications.constructor(::Val{:add}) = adder
Specifications.constructor(::Val{:multiply}) = multiplier

function (f!::Adjust)(x::FlatState, _Δt::Float64; _...)
    mergewith!(Base.Fix1(floor, Int) ∘ f!.adjust, x.counts, f!.adjustment)
    x
end

Models.parameters(f!::Adjust) = f!.adjustment
Models.remake(f!::Adjust, parameters::AbstractDict{Symbol}) = Adjust(
    f!.adjust,
    Dict(key => get(parameters, key, value) for (key, value) in f!.adjustment)
)

end
