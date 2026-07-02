"""
Contains components to build `Differentiation`-based `SciML.JumpModel`s for
random differentiation trees.

These models can be constructed by [`build`](@ref) from a [`Definition`](@ref)
or from the corresponding JSON specification. The differentiation tree is
chosen as the Huffman tree of a randomly sampled unit vector of target ratios
for the terminal states, with other parameters sampled independently for each
inner node. The corresponding [`Differentiation`](@ref) core network may
optionally regulate a peripheral network that is constructed from a
[`KroneckerNetworks.Definition`](@ref).

Like for `Differentiation`, the initial state must be set up before regulation
can start, which can be achieved by [`Differentiation.bootstrap`].
"""
module RandomDifferentiation

using ..Models: Models, V1, Differentiation, KroneckerNetworks
using ..Sampling: Nonnegative, BaseRatesTemplate
import ..Specifications

using Random

using DataStructures: PriorityQueue, enqueue!, dequeue!
using Distributions

"""
    DifferentiationTemplate

Defines how to sample a [`Differentiation.Definition`](@ref) from the specified
priors.

# Specification

In JSON, a `DifferentiationTemplate` is specified by a JSON object
```
{
    "ratios": <ratios>,
    "differentiator_base_rates": <differentiator_base_rates>,
    "timer_base_rates": <timer_base_rates>,
    "duration": <duration>,
    "trigger": <trigger>,
    "trigger_deposit": <trigger_deposit>,
    "buffer": <buffer>,
    "timer_deposit": <timer_deposit>,
    "buffer_deposit": <buffer_deposit>
}
```
where:
- `<ratios>` specifies a multivariate distribution by a JSON array, either
  - of the form `[<distribution>, <parameters>...]` where `<distribution>` is a
    JSON string naming a Distributions.jl `MultivariateDistribution`, and
    `<parameters>...` will be passed to its constructor, for example
    `["Dirichlet", 6, 5.0]`, or otherwise
  - listing the JSON numbers to be jointly returned on sampling, that is,
    specifying a `Product` of `Dirac` distributions, for example
    `[1, 2, 3, 4, 5]`.
  The `Vector{Float64}` sampled from this distribution will then define the
  differentiation tree and split ratios as described in [`Template`](@ref).
- `<differentiator_base_rates>` specifies the differentiator gene base rates'
  prior by a [`BaseRatesTemplate`](@ref).
- `<timer_base_rates>` specifies the timer genes' prior by a
  [`BaseRatesTemplate`](@ref).
- `<duration>` specifies created [`Differentiation.Transient`](@ref)'s
  `duration` prior by a [`Nonnegative{UnivariateDistribution}`](@ref).
- `<trigger>` optionally specifies the name of the root `Transient`'s upstream
  `differentiator` by a JSON string, defaulting to `"trigger"`.
- The optional mappings for `<trigger_deposit>`, `<buffer>`, `<timer_deposit>`
  and `<buffer_deposit>` all specify (exact) values for the respective
  properties shared by all the created [`Differentiation.Transient`](@ref)s.
"""
@kwdef struct DifferentiationTemplate
    ratios::MultivariateDistribution
    differentiator_base_rates::BaseRatesTemplate
    timer_base_rates::BaseRatesTemplate
    duration::Nonnegative{UnivariateDistribution}
    trigger::Symbol = :trigger
    trigger_deposit::Int = 0
    buffer::Vector{Float64} = Float64[]
    timer_deposit::Dict{Symbol, Int} = Dict{Symbol, Int}()
    buffer_deposit::Int = 0
end

"""
    InterRegulationTemplate

Defines how to sample regulatory links from the core differentiation network
into the peripheral network.

The total number of links will be sampled from `count`, and their parameters
will be sampled independently from `at`.

# Specification

In JSON, an `InterRegulationTemplate` is specified by a JSON object
```
{
    "count": <count>,
    "at": <at>
}
```
where `<count>` specifies a
[`Nonnegative{DiscreteUnivariateDistribution}`](@ref) prior for the number of
inbound links independently for each peripheral gene, and `<at>` specifies a
[`Nonnegative{UnivariateDistribution}`](@ref) prior for each link's parameter.
"""
@kwdef struct InterRegulationTemplate
    count::Nonnegative{DiscreteUnivariateDistribution} =
        Nonnegative{DiscreteUnivariateDistribution}(Dirac(0))
    at::Nonnegative{UnivariateDistribution} =
        Nonnegative{UnivariateDistribution}(Dirac(Inf))
    k::UnivariateDistribution = Dirac(-1.0)
end

"""
    Template

Defines how to sample a randomly differenting model from the specified priors.

[`build`](@ref) uses these definitions to construct a differentiation tree as
well as a peripheral network and its regulation by the differentiator genes. The
construction is as follows:
1. A `Vector{Float64}` is sampled from `differentiation.ratios` and normalized.
   It defines the target ratios (and count) of the terminal states, each of
   which will be represented by a signalling "differentiator" gene.
2. The differentiators will be organized into a Huffman tree (resulting in
   mostly balanced split ratios in the upper level) with one
   [`Differentiation.Transient`](@ref) for each inner node. The other parameters
   are set or sampled independently according to the specifications in
   `differentiation`. This constitutes the "core" network facilitating the
   differentiation.
3. Optionally, a Kronecker-linked peripheral network is sampled from
   `peripheral`.
4. Optionally, regulation of the peripheral genes is added from the set of all
   differentiator genes, transient or terminal, according to the specifications
   in `activation` and `repression`. For each peripheral gene and both kinds,
   the link count is sampled from the specified prior, and then that many
   inbound links are placed, sampling the regulator uniformly from the set of
   differentiator genes (with replacement) and the regulation constants
   independently from the corresponding priors.
5. The deposit instructions for the (top-level) trigger species are attached to
   the resulting [`Differentiation.Definition`](@ref).

# Specification

In JSON, a `RandomDifferentiation.Template` is specified as a JSON object
```
{
    "differentiation": <differentiation>,
    "peripheral": <peripheral>,
    "activation": <activation>,
    "repression": <repression>
}
```
where `<differentiation>` specifies a [`DifferentiationTemplate`](@ref),
`<peripheral>` specifies a [`KroneckerNetworks.Template`](@ref) and
`<activation>` and `<repression>` each specify an
[`InterRegulationTemplate`](@ref) respectively for transcriptional activation
and repression. Only `<differentiation>` is mandatory, the other mappings are
optional.

(Note that the JSON object specifying the `Template` will typically also
contain a `"seed"` mapping if it is specified as part of a
`RandomDifferentiation.Definition`.)
"""
@kwdef struct Template
    differentiation::DifferentiationTemplate
    peripheral::Union{Some{KroneckerNetworks.Template}, Nothing} = nothing
    activation::InterRegulationTemplate = InterRegulationTemplate()
    repression::InterRegulationTemplate = InterRegulationTemplate()
end

"""
    Definition

Defines how to construct a [`SciML.JumpModel`](@ref Models.SciML.JumpModel)
(via [`Differentiation`](@ref)) by sampling from a contained [`Template`](@ref)
with a fixed `seed`.

# Specification

In JSON, a `RandomDifferentiation.Definition` is specified as a JSON object
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

Models.describe(definition::Definition) = Models.Label(" \
    'regulation/random-differentiation' definition \
    with seed '$(definition.seed)'\
")

@kwdef struct Node
    next::Union{Node, Float64}
    alternative::Union{Node, Float64}
end

assemble_differentiation(
    ratio::Float64;
    template,
    name,
    randomness::AbstractRNG
) = (;
    differentiation = V1.Gene(
        name = Symbol(name),
        base_rates = rand(randomness, template.differentiator_base_rates),
    ),
    ratio,
    differentiators = [Symbol(name)],
)

function assemble_differentiation(
    node::Node;
    template,
    name = "differentiator",
    trigger = nothing,  # We special-case the root node by providing trigger.
    randomness::AbstractRNG,
)
    next = assemble_differentiation(
        node.next,
        name = "$(name)0";
        template,
        randomness,
    )
    alternative = assemble_differentiation(
        node.alternative,
        name = "$(name)1";
        template,
        randomness,
    )
    ratio = next.ratio + alternative.ratio
    differentiator = @something(
        trigger,
        V1.Gene(
            name = Symbol(name),
            base_rates = rand(randomness, template.differentiator_base_rates),
        )
    )
    duration = rand(randomness, template.duration)
    timer = V1.Gene(
        name = Symbol(),
        base_rates = rand(randomness, template.timer_base_rates),
    )

    (;
        differentiation = Differentiation.Transient(
            ratio = next.ratio / ratio,
            next = next.differentiation,
            alternative = alternative.differentiation;
            differentiator,
            duration,
            timer,
            template.buffer,
            template.timer_deposit,
            template.buffer_deposit,
        ),
        ratio,
        differentiators = [
            trigger === nothing ? [differentiator.name] : Symbol[]
            next.differentiators
            alternative.differentiators
        ]
    )
end

function Base.rand(randomness::AbstractRNG, template::Template)
    # Determine target terminal state ratios:
    ratios = Float64.(rand(randomness, template.differentiation.ratios))
    isempty(ratios) && error("no terminal states defined")
    ratios ./= sum(ratios)

    # Assemble the states into a Huffman tree:
    queue = PriorityQueue()
    for ratio in ratios
        enqueue!(queue, Ref(nothing) => ratio)
    end
    while length(queue) > 1
        next, next_ratio = peek(queue)
        dequeue!(queue)
        alternative, alternative_ratio = peek(queue)
        dequeue!(queue)

        node = Node(
            next = something(next[], next_ratio),
            alternative = something(alternative[], alternative_ratio),
        )
        enqueue!(queue, Ref(node) => next_ratio + alternative_ratio)
    end
    root = dequeue!(queue)[]

    # Convert the tree to a Differentiation.Transient definition, collecting
    # the names of the created differeniators:
    (; differentiation, differentiators) = assemble_differentiation(
        root,
        template = template.differentiation;
        template.differentiation.trigger,
        randomness,
    )

    # Sample a peripheral (Kronecker-linked) network, and have its genes
    # regulated by differentiators chosen uniformly at random:
    peripheral =
        if template.peripheral !== nothing
            rand(randomness, something(template.peripheral))
        else
            V1.Definition()
        end
    for gene in peripheral.genes
        for _ in 1:rand(randomness, template.activation.count)
            push!(
                gene.activation.slots,
                V1.HillRegulator(
                    from = rand(randomness, differentiators),
                    at = rand(randomness, template.activation.at),
                    k = rand(randomness, template.activation.k),
                )
            )
        end
        for _ in 1:rand(randomness, template.repression.count)
            push!(
                gene.repression.slots,
                V1.HillRegulator(
                    from = rand(randomness, differentiators),
                    at = rand(randomness, template.repression.at),
                    k = rand(randomness, template.repression.k),
                )
            )
        end
    end

    deposit = Dict(
        template.differentiation.trigger =>
            template.differentiation.trigger_deposit
    )

    Differentiation.Definition(
        meta = Dict{Symbol, Any}(:ratios => ratios);
        differentiation,
        peripheral,
        deposit,
    )
end

"""
    build(specification::AbstractDict{Symbol})
    build(definition::Definition; method::Symbol = :default)

Construct a randomly differentiating `SciML.JumpModel` from a
[`Definition`](@ref).

When interpreting a JSON specification, this function (in its first form) is
called to construct a concrete regulation model on encountering a
`{"{regulation/random-differentiation}": {...}}` literal. It will first
destructure the parsed JSON into a `Definition` and then proceed from there.

The result is constructed by first making a concrete
`Differentiation.Definition` from [`definition.template`](@ref Template) to
obtain the corresponding `Model` and then wrapping that up in a
[`Models.Wrapped`](@ref) with `definition`. This will result in the following
stack of abstractions:
- [`SciML.JumpModel`](@ref Models.SciML.JumpModel), specified by a
- `Catalyst.ReactionSystem`, specified by a
- [`V1.Definition`](@ref), specified by a
- [`Differentiation.Definition`](@ref), specified by
- `definition::`[`RandomDifferentiation.Definition`](@ref), potentially
  specified by
- `specification`

Note that typically, using `Differentiation` models require preparation of the
system state; see [`Differentiation.bootstrap`](@ref), which should also be used
for `RandomDifferentiation` models.

# Specification

Randomly differentiating networks are specified in JSON as
`{"{regulation/random-differentiation}": <definition>}` where `<definition>`
specifies a [`Definition`](@ref) as described there.

For an example, see `examples/specification/random-differentiation.json`.
"""
function build end

build(specification::AbstractDict{Symbol}) = build(
    Definition(
        seed = specification[:seed],
        template = Specifications.cast(Template, specification),
    ),
    method = Symbol(get(specification, :method, "default")),
)

build(definition::Definition; method::Symbol = :default) = Models.Wrapped(
    model = Differentiation.build(
        # Deterministically fill in the template to create a concrete
        # Differentiation.Definition from it:
        rand(Xoshiro(definition.seed), definition.template);
        method,
    );
    definition,
)

Specifications.constructor(::Val{Symbol("regulation/random-differentiation")}) =
    build

end
