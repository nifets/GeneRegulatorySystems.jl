# Regulation

All currently implemented gene regulation models assume (discrete-space) pure jump process dynamics and (only) allow stochastically realizing concrete trajectories.
For this purpose, they are first instantiated as [`Models.SciML.JumpModel`](@ref)s containing `Catalyst.ReactionSystem`s that are assembled according to various kinds of model definitions.
This compilation process proceeds in multiple stages that recursively translate higher- to lower-level definitions, and each stage wraps its result in a [`Wrapped`](@ref Models.Wrapped) model along with its (intermediate) definition.

While most of this package assumes a dynamics interpretation in terms of jump processes and discrete counts, when the package is used as a library, the intermediate definitions are accessible (including the underlying `Catalyst.ReactionSystem` and `JumpProcesses.JumpSystem`).
It is therefore possible, for example, to use these objects for static analysis with other SciML tooling or to simulate trajectories from an ODE relaxation.
Given any `Wrapped` model object, its `definition` property holds just that, and accessing the `model` property peels off one layer.
For example, given a JSON model specification
```jldoctest onion; setup = :(using GeneRegulatorySystems), output = false
d = """
    {"{regulation/kronecker}": {
        "seed": "seed",
        "base_rates": {"\$": ["defaults", "gene", "base_rates"]},
        "activation": {
            "adjacency": {
                "initiator": [[0.7, 0.6], [0.4, 0.2]],
                "power": 3
            },
            "at": ["LogNormal", 2.0, 1.0]
        },
        "repression": {
            "adjacency": {
                "initiator": [
                    [0.7, 0.6],
                    [0.5, 0.2]
                ],
                "power": 3
            },
            "at": ["LogNormal", 2.0, 1.0]
        }
    }}
"""

nothing

# output

```
we have:

```jldoctest onion
julia> model = Models.parse(d);

julia> model isa Models.Wrapped
true

julia> typeof(model.definition)
GeneRegulatorySystems.Models.KroneckerNetworks.Definition

julia> typeof(model.model.definition)
GeneRegulatorySystems.Models.V1.Definition

julia> typeof(model.model.model.definition)
Catalyst.ReactionSystem{Catalyst.NetworkProperties{Int64, SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymbolicUtils.SymReal}}}

julia> typeof(model.model.model.model)
GeneRegulatorySystems.Models.SciML.JumpModel

```

Higher-order model blueprints such [Differentiation](@ref) are expressed in terms of a common first-iteration model of simplified gene regulation defined in the [`Models.V1`](@ref) module.
Some of these models (named "templates", currently [Kronecker-linked networks](@ref) and [Random differentiation](@ref)) are stochastic and only become concerete model definitions (i.e., their structure deterministic) by affixing a seed.

## V1 model

```@docs
Models.V1
Models.V1.build
Models.V1.Definition
Models.V1.Gene
Models.V1.EukaryoteBaseRates
Models.V1.ProkaryoteBaseRates
Models.V1.Activation
Models.V1.Repression
Models.V1.Proteolysis
Models.V1.knockout
Models.Reaction
Models.Reagents
```

## Kronecker-linked networks

```@docs
Models.KroneckerNetworks
Models.KroneckerNetworks.build
Models.KroneckerNetworks.Definition
Models.KroneckerNetworks.Template
Models.Sampling.BaseRatesTemplate
Models.Sampling.Nonnegative
Models.KroneckerNetworks.NetworkTemplate
Models.KroneckerNetworks.ActivationNetworkTemplate
Models.KroneckerNetworks.RepressionNetworkTemplate
Models.KroneckerNetworks.ProteolysisNetworkTemplate
```

## Differentiation

```@docs
Models.Differentiation
Models.Differentiation.build
Models.Differentiation.Definition
Models.Differentiation.Transient
Models.Differentiation.bootstrap
```

## Random differentiation

```@docs
Models.RandomDifferentiation
Models.RandomDifferentiation.build
Models.RandomDifferentiation.Definition
Models.RandomDifferentiation.Template
Models.RandomDifferentiation.DifferentiationTemplate
Models.RandomDifferentiation.InterRegulationTemplate
```
