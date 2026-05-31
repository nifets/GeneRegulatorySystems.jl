module Models

using ..Specifications: Specifications, Specification, Template

import JSON

using Random

const DEFAULTS = "$(@__DIR__)/defaults.specification.json"

"""
    FlatState

Represents species counts at a single time point.

Counts are stored as a flat `Dict{Symbol, Int}` mapping dimension names to
counts. For models that have dimensions associated with subsystems (such as
genes in [`Models.V1`](@ref GeneRegulatorySystems.Models.V1)), the corresponding
dimension names in `FlatState` are flattened by joining on `"."`.
"""
@kwdef mutable struct FlatState
    t::Float64 = 0.0
    counts::Dict{Symbol, Int} = Dict{Symbol, Int}()
    randomness::AbstractRNG = Random.Xoshiro()
end
FlatState(x::FlatState) = FlatState(
    counts = deepcopy(x.counts);
    x.t,
    x.randomness
)

struct Branched
    stem
    branches::Vector
end
Branched(x) = Branched(x, [])

"""
    t(x)

Access the current simulation time of state `x`.
"""
function t end
t(x::FlatState) = x.t
t(x::Branched) = t(x.stem)

"""
    randomness(x)

Access the random number generator instance of state `x`.
"""
function randomness end
randomness(x::FlatState) = x.randomness
randomness(x::Branched) = randomness(x.stem)

"""
    Model{State}

Abstract supertype of all models.

`Model`s are functors: subtypes `M <: Model{State}` must define a call overload
method `(::M)(x::State, Δt::Float64; arguments...)` that advances the
simulation. Subtypes indicate which `State` type they expect to receive, and
must accept at least such `x`.

If `State` is a newly defined type (that is, specific to `M` and e.g. not
[`FlatState`](@ref)), that also requires implementing at least the following
methods:
- [`each_event(callback::Function, x::State)`](@ref each_event) to extract the
  state (and, if contained, trajectory) in long format.
- [`t(x::State)`](@ref t) to access the current simulation time.
- [`randomness(x::State)`](@ref randomness) to access the contained random
  number generator.
- [`adapt!(x::FlatState, f!::M, copy::Val)`](@ref adapt!) to convert a
  `FlatState` to a `State` and return it. The result must alias `x.randomness`.
  If `copy` is `Val(true)`, the result must otherwise be an independent deep
  copy of `x`.
- [`FlatState(x::State)`](@ref FlatState) to allow `adapt!` for subsequent
  models to fall back on converting to `FlatState` and retrying if there is no
  more specific `adapt!` method defined.
However, implementing `adapt!` methods for more specific state-model-pairs may
allow for more efficient state conversion between model invocations, for example
because a copy isn't required or parts of another state type can be reused.
`adapt!` is used by
[`Schedule`](@ref GeneRegulatorySystems.Models.Scheduling.Schedule) to convert
between the various state representations required by its contained primitive
`Model`s, but it may also be called directly.

In addition, a `Specifications.constructor(::Val{<:Symbol})` method typically
needs to be defined to return a function that instantiates a new `M` instance
from its `JSON.parse`d specification. This allows `M` to be defined within JSON
schedules.

# Construction

    Model(definition; bindings = Dict{Symbol, Any}())

Construct a `Model` subtype from `definition`.

If `definition` is not a [`Specification`](@ref), it will first be converted to
one. Then, if `definition` is a `Template`, it will be expanded to (and asserted to
actually be) a `Model`, and otherwise it will be put into a `Schedule`.

The `bindings` will be used for template expansion (if it contains such
references) or bound to the `Schedule` on construction.

# Invocation

    (f!::Model)(x, Δt::Float64; arguments...) = error("unimplemented")

Advance the simulation state `x` by at most `Δt` time units according to the
dynamics defined by `f!`.

Subtypes `<: Model{State}` must define specialized methods for this function,
one for each type of `x` they support, but at least for `State`.

Methods must advance the simulation state and return it, and they may do so as
an arbitrary type. They may arbitrarily modify `x` in the process, and they may
return (i.e. alias) it. Methods may assume that `Δt ≥ 0.0` and may accept
additional keyword `arguments` to control their behavior.
"""
abstract type Model{State} end

Model(x; bindings = Dict{Symbol, Any}()) =
    Model(Specification(x; bound = Set(keys(bindings))); bindings)
Model(specification::Specification; bindings = Dict{Symbol, Any}()) =
    Scheduling.Schedule(; specification, bindings)
Model(template::Template; bindings = Dict{Symbol, Any}()) =
    Specifications.expand(template; bindings)::Model

(f!::Model)(_x, _Δt::Float64; _...) = error("unimplemented")

"""
    Instant{State} <: Model{State}

Abstract supertype of instant models that ignore their time budget parameter
(`Δt`) and never progress on the time axis.
"""
abstract type Instant{State} <: Model{State} end

"""
    Wrapped{State} <: Model{State}

Wraps a `Model` together with an arbitrary `definition` object.

Invoking a `Wrapped` model will just forward the call to the contained `model`.

This is used to annotate any `Model` with the `definition` that produced it,
which can either be a specification or model-specific definition object (such as
[`Models.V1.Definition`](@ref)) or a
[`Scheduling.Locator`](@ref GeneRegulatorySystems.Models.Scheduling.Locator)
indicating where the model was defined within a
[`Schedule`](@ref GeneRegulatorySystems.Models.Scheduling.Schedule). Because
`model` can be `Wrapped` itself, this mechanism can be used to annotate `Model`s with a
chain of provenance, and the intermediate definitions produced during model
construction can be recreated using [`Scheduling.reify`](@ref).
"""
@kwdef struct Wrapped{State} <: Model{State}
    definition
    model::Model{State}
end

(f!::Wrapped)(x, Δt::Float64; arguments...) = f!.model(x, Δt; arguments...)

unwrap(model) = model
unwrap(wrapped::Wrapped) = unwrap(wrapped.model)

"""
    parameters(model)

Returns an `AbstractDict{Symbol, Any}`, corresponding to the tunable parameters of the model.
"""
parameters(wrapped::Wrapped) = parameters(wrapped.model)

"""
    remake(model, parameters::AbstractDict{Symbol, Any})

Return a copy of `model` with updated parameters.

`parameters` stores named tunable values. Each model type defines
which parameter types it recognises and how it applies them.
"""
remake(model, parameters::AbstractDict) = model
remake(wrapped::Wrapped, parameters::AbstractDict) = Wrapped(
    remake(wrapped.definition, parameters),
    remake(wrapped.model, parameters)
)


"""
    adapt!(x, f!; copy = false)

Convert the simulation state `x` to a type accepted by the model `f!`.

If `copy` is set, `x` is not modified and the return value independent of it,
except for retaining the the same `randomness` instance. Otherwise, `adapt!` is
allowed to modify or alias parts of `x`, which should therefore no longer be
used.
"""
function adapt! end
adapt!(x, f!::Model; copy = false) = _adapt!(x, f!, Val(copy))
adapt!(x, f!::Model, _copy) = _adapt!(FlatState(x), f!, Val(false))
adapt!(x, f!::Wrapped, copy) = _adapt!(x, f!.model, copy)

_adapt!(x, f!::Model, copy::Val) = adapt!(x, f!, copy)
_adapt!(x::Branched, ::Model{Branched}, ::Val{false}) = x
_adapt!(x::Branched, f!::Model, copy::Val) = _adapt!(x.stem, f!, copy)
_adapt!(x::FlatState, ::Model{FlatState}, ::Val{false}) = x
_adapt!(x::FlatState, ::Model{Any}, ::Val{false}) = x
_adapt!(x::FlatState, f!::Model, ::Val{true}) = _adapt!(
    FlatState(counts = deepcopy(x.counts); x.t, x.randomness),
    f!,
    Val(false)
)

"""
    each_event(callback::Function, x)

Invoke `callback` once for each state change event captured in state `x`.

This will iterate the state in long format, with one event per dimension at the
earliest contained timepoint and then, if the trajectory is included in `x`,
every change event ordered by timepoint. The required signature is
`callback(t::Float64, key::Symbol, value::Int)`.

(I intend to reimplement this as a Tables.jl rows interface.)

# Examples

```jldoctest; setup = :(using GeneRegulatorySystems)
julia> x = Models.FlatState(counts = Dict(:a => 10));

julia> Models.each_event(println ∘ tuple, x)
(0.0, :a, 10)

julia> Models.each_event(x) do t, key, value
           println((; t, key, value))
       end
(t = 0.0, key = :a, value = 10)
```
"""
function each_event end

each_event(callback::Function, x::FlatState) =
    for (key, value) in x.counts
        callback(x.t, key, value)
    end

each_event(callback::Function, x::Branched) = each_event(callback, x.stem)

"""
    Reagents

Defines (integer) stoichiometries for a set of chemical species.

This is used to define reagents on either side of a [`Reaction`](@ref).

# Specification

In JSON, `Reagents` are specified as a JSON object mapping species names to
(integer JSON number) stoichiometric coefficients. Names may contain `.`
separators to refer to species of a subsystem (such as a gene).

As a shortcut, `Reagents` may alternatively be specified as a JSON Array of such
species names. Repeated entries will increment the corresponding stoichiometric
coefficient.
"""
@kwdef struct Reagents
    counts::Dict{Symbol, Int} = Dict{Symbol, Int}()
end

"""
    Reaction

Defines a constant-rate reaction pair between two sets of [`Reagents`](@ref)s as
part of a model.

# Specification

In JSON, a `Reaction` is specified as a JSON object
```
{
    "name": <name>,
    "from": <from>,
    "to": <to>,
    "rates": [<forward>, <reverse>]
}
```
where `<from>` and `<to>` each specify [`Reagents`](@ref) defining the (integer)
stoichiometries respectively for the reactants and products (of the forward
reaction), and `<forward>` and `<reverse>` must be JSON numbers defining the
corresponding rate constants.

If present, `<name>` must be a JSON string identifying the reaction, otherwise it
will be set automatically to `"rxn-<i>"` using the reaction's index
within its containing list. Names must be unique within a model.

For example,
```
{"from": ["A", "B.mrnas"], "rates": [0.00001, 0.00002], "to": ["AB"]}
```
will define the Catalyst.jl
```
@reaction (0.00001, 0.00002), A + B <--> AB
```
e.g. to be added to a `JumpModel`.

As a convenience, `"rates": [<forward>, <reverse>]` may alternatively by written
as `"rate": <forward>`, setting `<reverse>` to zero.
"""
@kwdef struct Reaction
    name::Symbol
    from::Reagents = Reagents()
    to::Reagents = Reagents()
    k₊::Float64 = 0.0
    k₋::Float64 = 0.0
end

Specifications.cast(::Type{Vector{Reaction}}, xs::AbstractVector; context) = [
    Specifications.cast(
        Reaction,
        merge(Dict(:name => "rxn-$i"), x);
        context,
    )
    for (i, x) in enumerate(xs)
]

Specifications.cast(::Type{Reaction}, x::AbstractDict{Symbol}; context) =
    @invoke Specifications.cast(
        Reaction::Type,
        if haskey(x, :rates)
            merge(x, Dict(zip((:k₊, :k₋), x[:rates])))
        elseif haskey(x, :rate)
            merge(x, Dict(:k₊ => x[:rate]))
        else
            error("missing rates in reaction specification")
        end::AbstractDict{Symbol};
        context,
    )

Specifications.cast(::Type{Reagents}, x::AbstractDict{Symbol}; _...) =
    Reagents(x)

function Specifications.cast(::Type{Reagents}, xs::AbstractVector; _...)
    result = Reagents()

    for x in xs
        reagent = Symbol(x)
        result.counts[reagent] = get(result.counts, reagent, 0) + 1
    end

    result
end

Specifications.representation(x::Reagents) = x.counts
Specifications.representation(x::Reaction) = Dict{Symbol, Any}(
    :name => Specifications.representation(x.name),
    :from => Specifications.representation(x.from),
    :to => Specifications.representation(x.to),
    :rates => [x.k₊, x.k₋]
)

abstract type Description end

struct EmptyDescription <: Description end

struct Descriptions <: Description
    descriptions::Vector{Description}
end

@kwdef struct Provenance <: Description
    source::Description
    description::Description
end

@kwdef struct Label <: Description
    label::String = ""
end

@kwdef struct Network <: Description
    species_groups::Vector{Symbol}
    links
    aliases::Dict{Symbol, Symbol} = Dict{Symbol, Symbol}()
end

@kwdef struct ReactionNetwork <: Description
    reactions::Vector{Reaction} = Reaction[]
end

describe(::Any) = EmptyDescription()
describe(wrapped::Wrapped) = Provenance(
    source = describe(wrapped.definition),
    description = describe(wrapped.model),
)

load_defaults() = JSON.parsefile(DEFAULTS, dicttype = Dict{Symbol, Any})

include("plumbing.jl")
include("scheduling.jl")
include("resampling.jl")
include("sciml.jl")
include("regulation/v1.jl")
include("regulation/differentiation.jl")
include("regulation/sampling.jl")
include("regulation/kronecker_networks.jl")
include("regulation/random_differentiation.jl")
include("extraction.jl")

"""
    parse(
        definition;
        into = "",
        channel = "",
        defaults = load_defaults(),
        others...,
    )

Parse a JSON `definition` to construct a [`Model`](@ref).

The `definition` can be a `String` or an `IO`, and keyword arguments will be
collected and passed as bindings for the model construction.

# Examples

```jldoctest; setup = :(using GeneRegulatorySystems)
julia> f! = Models.parse(\"""
           {"{add}": {"a": 10}}
       \""")
GeneRegulatorySystems.Models.Plumbing.Adjust(+, Dict(:a => 10))
```

This is equivalent to the following:

```jldoctest; setup = :(using GeneRegulatorySystems)
julia> import JSON

julia> source = \"""
           {"{add}": {"a": 10}}
       \""";

julia> definition = JSON.parse(source, dicttype = Dict{Symbol, Any})
Dict{Symbol, Any} with 1 entry:
  Symbol("{add}") => Dict{Symbol, Any}(:a=>10)

julia> f! = Model(definition)
GeneRegulatorySystems.Models.Plumbing.Adjust(+, Dict(:a => 10))
```
"""
parse(
    definition;
    into = "",
    channel = "",
    defaults = load_defaults(),
    others...,
) = Model(
    JSON.parse(definition, dicttype = Dict{Symbol, Any}),
    bindings = (; into, channel, defaults, others...) |> pairs |> Dict
)

"""
    load(path::AbstractString; bindings...)

Load a `Model` from a JSON file.

This is just a convenience wrapper around [`parse`](@ref); keyword arguments
will be forwarded.
"""
load(path::AbstractString; bindings...) = open(path) do file
    parse(file; bindings...)
end

end
