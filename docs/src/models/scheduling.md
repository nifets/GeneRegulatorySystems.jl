# Scheduling

The scheduling system forms the core of this package as it ties the various included regulation and instant adjustment models together into single reproducible experiment definitions that, besides the gene regulatory systems to be simulated, may also include initial setup and seed control, simulated -omics extraction protocols, (optionally regularly repeating) interventions as well as simulation branching into independent samples sharing parts of their history.

The primary entity in this module is the `Schedule`, which is a type of `Model` that advances the simulation by organizing simulation segments that proceed by invoking other `Model`s.
Since the orchestrated `Model`s may themselves be `Schedule`s, it is possible to assemble complex sequences of simulation segments and conceptually integrate the accumulated state changes into single unified trajectories.

Each `Schedule` has a `Specification` that represents the exact instructions for how to advance.
The `Specification` also has an [alternative representation](@ref Specification) as a JSON document from which it can be conveniently constructed.

Details are given below, but broadly, the scheduling system supports
* construction of `Schedule`s and other `Model`s from their JSON representation,
* a simple templating mechanism allowing the definition of named bindings and their insertion into specifications before model construction,
* iteration both over explicit enumerations of specifications and over values to consecutively insert into a nested specification template,
* repetition of (finite-time) nested `Schedule`s,
* automatic conversion of the simulation state between the representations required by the consecutive segments' `Model`s,
* random seed control,
* hooks for output control and progress reporting, and
* trajectory branching.

To support this multitude of requirements, `Specification`s are assembled from reusable building blocks (which are themselves `Specification`s), effectively defining a domain-specific programming language.
Conceptually, the JSON representation is the source code of the scheduling language, its `Dict`/`Array` form after `JSON.parse`ing is an intermediate representation, its `Specification` form is an abstract syntax tree, and a `Schedule` containing it is an interpreter that is invoked whenever the simulation state should be advanced.
The [syntax of the scheduling language](@ref "Schedule specification") is therefore defined by the allowed combinations of `Specification`s, and the [associated semantics](@ref "Schedule semantics") defined in terms of which simulation segments are produced and executed.

The [experiment command line tool](@ref "`experiment` tool") wraps all of this functionality by loading a collection of JSON specifications, interpreting them and collecting the results in a [structured format](@ref "Results format"), along with an index of simulation segments and the specifications used to produce the data.
For reproducible simulation experiments, this is the recommended way to use the package.

Each simulation segment can be addressed by a path relative to the `Schedule` that produced it, and the index includes each segment's path.
Since the `Schedule` can be recreated from its `Specification` (and the associated bindings) and the sequence of segments it produces is deterministic, the scheduling system's iteration state can be [fully reified](@ref "Paths and reification") for any simulation segment, including simulation context such as the exact `Model` in effect.
This can be used to obtain that model for export or analysis without re-running the full experiment, and also to better understand or debug the scheduling subsystem.

## Schedule specification

A `Schedule`'s simulation plan is represented by its `specification`, which is effectively a tree of syntactic elements (all `<: Specification`) of the scheduling mini-language, with `Scope`, `List` and `Each` being inner nodes and `Template`, `Load` and `Slice` terminal nodes.

Besides its `specification`, each `Schedule` keeps track of a map of named value `bindings`.
`Template`s (and by extension, indirectly also other `Specification`s enclosing them) may include references to such named values that are to be inserted when stepping through the `Schedule`; these values may be defined by an enclosing `Scope` or `Each`, or they may be free and must then be injected into the `Schedule`'s bindings on construction.
`Schedule`s therefore close over their `specification`, which in that process will ensure that no references are left dangling.

The following are the potential elements of a `Specification`:

```@docs
Specifications.Template
Specifications.Slice
Specifications.Load
Specifications.Sequence
Specifications.List
Specifications.Each
Specifications.Scope
```

## JSON representation

While `Schedule`s and their `specification`s can be built by hand, users will typically construct them using the `Specification` function from their alternative representation as JSON documents:

```@docs
Specification
```

When parsing the alternative representation as a `Specification`, whenever a `Template` is constructed, the corresponding constructor, to be called when it is eventually `expand`ed, is looked up by calling `Specifications.constructor`. This allows the definition of (non-`String` and non-number) terminal values within the JSON specification language, including all actual `Model`s.

```@docs
Specifications.constructor(::Symbol)
```

To see which objects can be defined in the language using the `{"{...}": ...}` syntax, you can simply call `methods(Specifications.constructor)`. To register new types of objects and in this way support them in the language, you may define new methods of the form `constructor(::Val{:...})`.

## Schedule semantics

When a `Schedule` is invoked to advance the simulation state (i.e. by calling it as a functor, like any `Model`), its exact behavior is determined by the type of its (top-level) `specification`, which may involve constructing and recursively avancing on nested `Schedule`s until the recursion terminates on the non-`Schedule` (*primitive*) `Model`s to actually produce and execute simulation segments.
These terminal models are wrapped in `Primitive`s (which are also `<: Model`) that delegate simulation but add hooks for output handling and progress reporting and further automatically convert the simulation state to the representation required by the wrapped model.

```@docs
Scheduling.Primitive
```

As the interpreter descends on the `specification`, the constructed nested `Schedule`s or `Primitive`s obtain new or replaced bindings either from direct definition in the `specification` or from implicit built-in behavior (mostly related to output control), and they further keep track of their path in the recursion (see [Paths and reification](@ref)).

As a reminder, since each `Schedule` `f!` is a `Model`, it may advance the simulation state `x` by `Δt ≥ 0.0` units of simulation time when being called as a functor like `f!(x, Δt; ...)`.
This call dispatches on the specific type of the `Schedule`'s (top-level) specification, and the following describes the resulting behavior.

```@docs
Schedule
Scheduling.Merge
```

## Paths and reification

As a `Schedule` invocation descends on its `Specification`, it keeps track of its current path in that tree and includes it when constructing the terminal `Primitive` models that ultimately produce the simulation segments.
In this way, all segments generated by a `Schedule` and its top-level `Specification` are uniquely identified by their path, and path prefixes likewise address contiguous ranges of simulation segments associated with inner nodes of the specification.

Further, the current path is recorded for all definitions evaluated during schedule execution:

```@docs
Scheduling.Locator
```

Every path is a `String` that consists of segments, each describing a single step of descent:
* descending on a `Scope` with `branch` unset appends `"+"`,
* descending on a `Sequence` appends a `"/"` and the within-sequence index if a directly enclosing `Scope` has `branch` set, and a `"-"` and the within-sequence index otherwise, and
* evaluating a binding definition in a `Scope` appends a `"."` and the corresponding key.

While the terminal `Primitive` models may advance the simulation state stochastically, their construction and organization as part of `Schedule` execution is fully deterministic, and further independent between recursion branches at each inner `Specification` node.
This means that each object produced in the process of stepping through the `Schedule` can be reified exactly, given only the root `Schedule` (defined by `specification` and `bindings`) and the object's corresponding path.
This functionality is exposed through the `reify` function:

```@docs
Scheduling.reify
```

The [`experiment` tool](@ref) traces all simulation segments' `Primitive` paths and includes them in its results index so that they unambiguously identify their definition location and can also be reified if needed.

Reification can be useful to obtain a specific object or model that is defined within a JSON specification document, such as for export or further analysis. It can also be used for better understanding or debugging the scheduling mechanism. To assist with this, the [`reify` tool](@ref) provides a CLI wrapper script to the `reify` function that supports pretty printing and can be pointed either at an experiment results location or directly at a JSON specification file.
