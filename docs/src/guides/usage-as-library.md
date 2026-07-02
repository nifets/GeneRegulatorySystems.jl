# Usage as a library

Everything the CLI tools do can be replicated manually using the library interface, that is, `using GeneRegulatorySystems` from a Julia session or program.
This includes the construction of the actual (Catalyst.jl-based) regulation models that are simulated as part of an experiment `Schedule`.

The following describes how to obtain these models so that they can be studied in isolation using tools from the SciML ecosystem.
In principle, they could then also be used in larger constructions such as inference procedures.

## Setup

The package can simply be installed using Julia's package manager.
Since it is not published to the General registry yet, you need to reference the package repository directly.
For example, type
```
]add https://github.com/drostlab/GeneRegulatorySystems.jl
```
in the Julia REPL after activating your work environment, like a toy environment (`]activate --temp`).
You can then run the examples on this page.

### Definition of new model templates

The CLI tools do not support an explicit extension mechanism, so if you indend to [define new model templates](@ref "The `Model` interface") and actually make them available as part of the scheduling language, you need to modify the the Git repository (and ideally contribute any new models upstream eventually).

If you have [installed the CLI tools by cloning the Git repository](@ref "CLI setup"), it may make sense to `]`[`develop`](https://pkgdocs.julialang.org/v1/repl/#repl-develop) that path directly.
For example, in a Julia REPL after activating your work environment, type
```
]develop ~/src/grs
```
where you replace `~/src/grs` by the repository location.
Changes to that repository then affect both the CLI tools and the work environment.

### Pluto.jl

If you intend to work on a [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebook, because its built-in package manager only allows registered packages, you need to use a ["Pkg cell"](https://plutojl.org/en/docs/packages-advanced/#pattern-the-pkg-cell).
For example, you can define a cell
```julia
begin
    import Pkg
    Pkg.activate(mktempdir())

    Pkg.add(url = "https://github.com/drostlab/GeneRegulatorySystems.jl")
    using GeneRegulatorySystems

    # maybe add other packages here...
end
```
to add and load this package.

Like above, if you want to operate on a locally checked out repository (e.g. because you need access to new model templates that you have registered for the CLI tools), you can `]develop` that repository.
For example, you can define a cell
```julia
begin
    import Pkg
    Pkg.activate(mktempdir())

    Pkg.develop(path = joinpath(ENV["HOME"], "src", "grs"))
    using GeneRegulatorySystems
end
```
(where you may need to change the path definition if you installed the package somewhere else).
However, you may need to restart the notebook after making changes to the repository.

## Constructing models

```@example usage
using GeneRegulatorySystems
```

The lowest-level way to construct a regulation model `f!` is to call a model constructor with a corresponding definition.
For example,
```@example usage
base_rates = Models.V1.EukaryoteBaseRates(
    activation = 2.5,
    deactivation = 10.0,
    trigger = 6.6e-7,
    transcription = 0.001,
    processing = 0.02,
    translation = 2.5e-9,
    abortion = 0.01,
    premrna_decay = 0.001,
    mrna_decay = 0.001,
    protein_decay = 3.0e-10,
)

definition = Models.V1.Definition(genes = [
    Models.V1.Gene(name = :a; base_rates)
    Models.V1.Gene(
        name = :b,
        repression = Models.V1.Repression(slots = [
            Models.V1.HillRegulator(from = :a, at = 2.0)
        ]);
        base_rates,
    )
])

f! = Models.V1.build(definition)
nothing # hide
```
builds a simple 2-gene network where one gene represses the other.
This same model can also be constructed from a JSON specification, either in its `JSON.parse`d form
```@example usage
base_rates = Dict(
    :activation => 2.5,
    :deactivation => 10.0,
    :trigger => 6.6e-7,
    :transcription => 0.001,
    :processing => 0.02,
    :translation => 2.5e-9,
    :abortion => 0.01,
    :premrna_decay => 0.001,
    :mrna_decay => 0.001,
    :protein_decay => 3.0e-10,
)

specification = Dict(:genes => [
    Dict(:name => "a", :base_rates => base_rates)
    Dict(
        :name => "b",
        :base_rates => base_rates,
        :repression => [Dict(:from => "a", :at => 2.0)],
    )
])

f! = Models.V1.build(specification)
nothing # hide
```
or directly from raw JSON:
```@example usage
import JSON

base_rates_json = """
{
    "activation": 2.5,
    "deactivation": 10.0,
    "trigger": 6.6e-7,
    "transcription": 0.001,
    "processing": 0.02,
    "translation": 2.5e-9,
    "abortion": 0.01,
    "premrna_decay": 0.001,
    "mrna_decay": 0.001,
    "protein_decay": 3e-10
}
"""

specification_json = """
    {"genes": [
        {"name": "a", "base_rates": $base_rates_json},
        {
            "name": "b",
            "base_rates": $base_rates_json,
            "repression": [{"from": "a", "at": 2.0}]
        }
    ]}
"""

specification = JSON.parse(specification_json, dicttype = Dict{Symbol, Any})

f! = Models.V1.build(specification)
nothing # hide
```
Note that the JSON **must be parsed with `dicttype = Dict{Symbol, Any}`** and likewise any `Dict`/`Array`-specification as shown above must use `Dict{Symbol, Any}` for every `Dict`.
Finally, a `Model` can be constructed via [`Models.parse`](@ref) (or loaded from a model specification file using [`Models.load`](@ref))

```@example usage
f! = Models.parse("""
    {"{regulation/v1}": $specification_json}
""")
nothing # hide
```
which means that the JSON specification is interpreted as a literal of the [scheduling language](@ref "Scheduling").
All of the above examples produce the same regulation model.
The result `f!` is a `Wrapped` model that contains a `model::Model` tagged with an arbitrary `definition` object:

```@example usage
typeof(f!)
```
In this case, the model is twice-wrapped, where the outer definition
```@example usage
f!.definition
```
is the `V1.Definition` defining the model and the inner definition
```@example usage
f!.model.definition
```
is the `Catalyst.ReactionSystem` intermediate representation.
The inner model
```@example usage
typeof(f!.model.model)
```
is the actual [`JumpModel`](@ref Models.SciML.JumpModel) containing the prepared `system` (a `ModelingToolkit.JumpSystem`), `method` (a `JumpProcesses.AbstractAggregatorAlgorithm`) and `parameters`.

See [Instantiating models](@ref) and [Regulation](@ref) for information on how to construct other regulation models.

## Simulating models

Models can be simulated either using [The `Model` interface](@ref), or directly via the contained SciML objects.
The former is what the [Scheduling](@ref) machinery and CLI tools use, while [the latter is probably what you want](@ref "Using SciML directly") if you are using this package from Julia.
We will quickly start with the `Model` interface anyway for rough overview of what the [`experiment` tool](@ref) does internally.

### Using the `Model` interface

!!! note "Implementation note"
    The reason this interface exists in the first place is because the package and its first-iteration base regulation model were originally not written in terms of Catalyst.jl/SciML.
    [Scheduling](@ref) and, by extension, the [experiment runner](@ref "`experiment` tool"), can stitch together heterogeneous model types (besides SciML-based regulation models) along the time dimension.
    Another design goal was reproducibility and pervasive seed control in the presence of arbitrary intervention schedules and trajectory branching.

    Even though all regulation is currently implemented exclusively using Catalyst.jl and JumpProcesses.jl, and interventions could probably now be implemented in that way as well, the current scheduling subsystem is by now quite flexible and not all of its constructs map one-to-one to equivalents in SciML-land in a straightforward way.
    Further, SciML is quite the moving target and some isolation seemed helpful.

    However, in some situations this flexibility comes at the cost of performance, and duplicating complexity is always annoying, and further ModelingToolkit.jl seems to become more stable and complete, so I may still move the scheduling system closer to SciML eventually.
    In the meantime, as it is currently structured it still has lots of potential for optimization if performance becomes an issue in specific instances.

To simulate a trajectory, you need to first construct an initial state and then invoke the model `f!` as described in [Models](@ref):
```@example usage
# Initial state for V1 regulation requires some molecular species present;
# defaults are:
#     Dict(:polymerases => 5*10^5, :ribosomes => 2*10^6, :proteasomes => 10^6)
x = Models.FlatState(counts = Models.load_defaults()[:bootstrap])

# f! contains a JumpModel, which requires a JumpState; convert it:
x = Models.adapt!(x, f!)

x = f!(x, 100.0)  # step 100 time units

Models.FlatState(x)  # show as FlatState
```
This will however not record the trajectories but only keep the final counts.
The `JumpState` holds the SciML `x.integrator` and thus its solution object `x.integrator.sol`, but it will only contain values for the final timepoint:
```@example usage
x.integrator.sol
```

To retain intermediate jump events, the model can be invoked with the `record` keyword argument:
```@example usage
x = f!(x, 1.0, record = true)
x.integrator.sol
```
The contained integrator will then contain the trajectory segment for this invocation, (unfortunately) in a dense format.
The [`experiment` tool](@ref) accesses this segment via [`Models.each_event`](@ref), which calls back for each initial value and then for each jump event:
```@example usage
Models.each_event(x) do t, name, value
    @show (; t, name, value)
end
```
Between saving nothing and everything, other saving options from SciML are not supported here; they are all emulated via the [Scheduling system](@ref "Scheduling").
For example, to record the state only at a single time slice, we can advance to that timepoint without recording and then record without advancing:
```@example usage
x = f!(x, 99.0)
x = f!(x, 0.0, record = true)
x.integrator.sol
```
Note that each invocation of the `JumpModel` `f!` will initially clear the previously recorded trajectory segment, which therefore needs to be extracted immediately.

### Using SciML directly

As previously noted, the actual regulation `JumpModel` is defined by a `method` applied to a `system` along with a set of `parameters`.
These can be used directly to sample trajectories using the standard JumpProcesses.jl interface:
```@example usage
using Random
using JumpProcesses
using ModelingToolkit
using Plots

(; system, method, parameters) = f!.model.model

bootstrap = Models.load_defaults()[:bootstrap]
normalize_name = GeneRegulatorySystems.Models.SciML.normalize_name
u₀ = [
    s => Int(get(bootstrap, normalize_name(s), 0))
    for s in unknowns(system)
]
aggregator = method
rng = Xoshiro()
problem = JumpProblem(system, vcat(u₀, parameters), (0.0, 3e4); aggregator, rng)

solution = solve(problem, SSAStepper())

plot(solution, idxs = (:a₊proteins, :b₊proteins))
```

But since intermediate definitions are retained as described above, we also have access to the `ReactionSystem` and can thus use it with the [tools Catalyst.jl provides](https://docs.sciml.ai/Catalyst/stable/). For example:

```@example usage
using CairoMakie
using Catalyst
using GraphMakie
using NetworkLayout

reaction_system = f!.model.definition
plot_network(reaction_system)
```

```@example usage
conservationlaw_constants(flatten(reaction_system))
```

This allows us to easily construct an ODE approximation of the jump process dynamics

```@example usage
using OrdinaryDiffEqTsit5

approximation = ode_model(reaction_system)
```

and simulate from it:

```@example usage
u₀′ = map(x -> x.first => float(x.second), u₀)
problem = ODEProblem(complete(approximation), u₀′, (0.0, 3e4), parameters)
solution = solve(problem, Tsit5())
Plots.plot!(solution, idxs = (:a₊proteins, :b₊proteins))
```

## Reifying `Model`s from `Schedule`s

The regulation `JumpModel`s that a `Schedule` constructs during its execution can be plucked out and used in isolation.
A `Schedule` can be loaded from its JSON specification using the [`Models.load`](@ref) function:

```@example usage
source = "../../../examples/toy/repressilator.schedule.json"
schedule! = Models.load(source, seed = "(ignored)")
```

In principle, `Models.load` can load regulation model specifications directly as shown above, but since the specification is interpreted as an expression of the [scheduling language](@ref "Schedule specification"), it may contain references that need to be substituted in to make it a complete model specification.
The referenced values need to be provided as bindings to `Models.load` so that it can create a closed `Schedule`.
The schedule specification loaded in the preceding example
```@example usage
read(source, String) |> print
```
refers to the `seed`, which is why it needs to be provided as a keyword argument, even though we will not actually use it in the following.

The actual regulation model can now be [reified](@ref Scheduling.reify) from its path within the schedule:
```@example usage
repressilator! = Scheduling.reify(schedule!, "+.do")
reaction_system = repressilator!.model.model.definition
```
The specification path refers to a named binding `:do` defined within a sub-schedule (in this case, the one nested inside the top-level `Scope`/`{}`); reification will first create that nested `Schedule`, evaluating the definition in the process, and then look up the resulting binding by name.
For more on this, see [Paths and reification](@ref).

Compared to the directly constructed models [from above](@ref "Constructing models"), it is wrapped once more, now additionally recording its definition location within `schedule!`:
```@example usage
repressilator!.definition
```

And just [like above](@ref "Simulating models"), we can now use the contained `JumpModel`, for example to manually run simulations:
```@example usage
(; system, method, parameters) = repressilator!.model.model.model

start = merge(bootstrap, Dict(Symbol("1.proteins") => 100))
u₀ = [
    s => Int(get(start, normalize_name(s), 0))
    for s in unknowns(system)
]

aggregator = method
rng = Xoshiro("seed")
problem = JumpProblem(system, vcat(u₀, parameters), (0.0, 1e5); aggregator, rng)

solution = solve(problem, SSAStepper())

Plots.plot(
    solution,
    idxs = (Symbol("1₊proteins"), Symbol("2₊proteins"), Symbol("3₊proteins")),
    camera = (60, 60),
)

u₀′ = map(x -> x.first => float(x.second), u₀)
problem = ODEProblem(reaction_system, u₀′, (0.0, 1e5), parameters)
solution = solve(problem, Tsit5())
Plots.plot!(
    solution,
    idxs = (Symbol("1₊proteins"), Symbol("2₊proteins"), Symbol("3₊proteins"))
)
```