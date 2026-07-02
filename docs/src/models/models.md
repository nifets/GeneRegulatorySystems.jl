# Models

In this package, the common interface to advance simulations is to apply functions of the form `x = f!(x, Δt::Float64; parameters...)` to the current simulation state `x`, and we call the effects of a single model invocation a *simulation segment*.
The state `x` has no fixed representation and may even change type within `f!`, so while `f!` may modify `x` in place, the canonical result is given by its return value.
The time step `Δt ≥ 0.0` indicates the units of simulation time that `f!` is supposed (but not required) to advance.
We refer to functions that adhere to this contract as *models*, and all of those provided in this package are implemented as (immutable) functors of type `<: Model`.

Models are of two kinds: *dynamics models* that represent continuous-time behaviour and thus advance the simulation state along the time axis, and *instant models* that represent instantaneous adjustments and do not actually advance along the time axis.
The former includes all [gene regulation models](@ref "Regulation"), while the latter, defined as `<: Instant`, includes various forms of instant adjustments that may be used to construct [interventions](@ref "Plumbing") as well as [observation schemes](@ref "Extraction").
Besides those provided, it is straightforward to define new `Model`s and state representations and wire them up with the rest of the package; the following describes the interface that needs to be implemented to achieve this.

## State

Different models may require different representations of the simulation state to operate on. There is no common supertype of these state representations, but all of these `State` types must at least carry information about the current simulation time, the values of all modeled dynamic variables at that time, and the current state of a random number generator.
For the dynamics models implemented in this package, the state object also accumulates the full simulation trajectories of a simulation segment if their output is required.
This information is accessed using the following methods, which must be implemented for all state types:

```@docs
Models.t
Models.randomness
Models.each_event
```

After the trajectories have been read, they may be cleared to conserve memory:

```@docs
Models.empty_trajectory!
```

Most `Instant` models in this package use `FlatState` as their state representation, which simply contains the current species counts in a `Dict{Symbol, Int}` (i.e., not nested for subsystems) but does not retain full trajectories.
For `FlatState`s, `each_event` emits one event for each contained species at the current simulation time.

```@docs
Models.FlatState
```

For example:
```jldoctest; setup = :(using GeneRegulatorySystems)
julia> f! = Models.Plumbing.Wait();

julia> g! = Models.Plumbing.adder(Dict(:a => 10));

julia> x = Models.FlatState();

julia> x = f!(x, 100.0);

julia> x = g!(x, 0.0);

julia> Models.t(x)
100.0

julia> Models.each_event(x) do t, key, value
           println((; t, key, value))
       end
(t = 100.0, key = :a, value = 10)
```

## Instantiating models

Although model functors can be created directly like in the preceding example, in practice they will likely be constructed by loading a JSON definition:

```@docs
Models.parse
Models.load
```

## SciML-based jump process models

While the simulation interface in principle allows the definition of arbitrary dynamics, all of the [first-iteration gene regulation model templates](@ref "Regulation") currently included in this package are ultimately interpreted as sampling trajectories from (discrete-space) pure jump processes.
This choice is also reflected in the current implementation of instant adjustment and observation models as well as the command-line wrapper scripts and their shared [results format](@ref "Results format") (which records discrete state change events in long format).

Originally built using a custom Gillespie-type integrator, the included regulation models have subsequently been migrated onto the SciML ecosystem, specifically Catalyst.jl and JumpProcesses.jl.
This allows using those tools to assemble systems and to sample trajectories with fast off-the-shelf stochastic simulation algorithms such as `RSSACR`, but also enables further analysis of the intermediate model representations with the various SciML packages.

Common machinery to build and simulate such SciML-based models are defined in the `Models.SciML` module, which contains `JumpModel <: Model{JumpState}` and the corresponding `JumpState`.

```@docs
Models.SciML.JumpModel
Models.SciML.JumpState
```

For `JumpState`s, `each_event` emits events depending on whether the preceding model invocation was instructed to record the trajectory or not.
Without recording, `each_event` behaves similarly to the `FlatState` case, emitting one event for each contained species at the final simulation time.
With recording, `each_event` emits one event per species at the beginning of the trajectory, and subsequently one event for every recorded change.
The species names will be converted to `Symbol`s, replacing potential `₊`s (subsystem separators) by `.`s.

## Connecting model invocations

Consecutive invocations of models that operate on different state representations require a state conversion in between.
For this purpose, the `adapt!` function can be used: it transforms a state to the type required for a given `Model`, and it is predefined for all state-model combinations included in this package.

```@docs
Models.adapt!
```

This design (having a type-unstable state representation) generally assumes that the total number of simulation segments is not too large and that the majority of computation time is spent within the primitive model invocations.
A further caveat is that the construction of JumpProcesses.jl `JumpProblem`s and initializing its integrators is quite costly and therefore repeatedly `adapt!`ing state to `JumpModel`s, for example after executing periodic interventions, can be slow in our current implementation.

Defining a new `Model` or state type `State` typically also requires defining new `adapt!` methods.
However, if no specific method is found for a given state-model pair, `adapt!` defaults to retrying with a `FlatState` copy of the state.
This means that often, defining `FlatState(x::State)` and `adapt!(x::FlatState, f!::Model, _copy)` is sufficient ([see below](@ref "The `Model` interface")).

Model invocations both of instant and dynamics models can be organized in [`Schedule`](@ref "Scheduling")s, which are themselves `Model`s that execute multiple simulation segments when invoked.
The scheduling machinery internally uses `adapt!` to piece the model segments together.

Typically, `adapt!` may modify or alias its state parameter when preparing the result.
But since `Schedule` also supports trajectory branching, `adapt!` can be instructed (with its `copy` parameter) to make sure that the resulting state is independent of the input (besides always aliasing the used random number generator).

## The `Model` interface

```@docs
Models.Model
Models.Instant
Models.Wrapped
```
