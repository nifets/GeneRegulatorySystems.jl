module Scheduling

using ...GeneRegulatorySystems: GeneRegulatorySystems
using ..Models: Models, Model, Instant, FlatState
using ..Models.Plumbing: Wait
using ...Specifications:
    Specifications,
    Each,
    List,
    Load,
    Scope,
    Sequence,
    Slice,
    Specification,
    Template

using Logging: LogLevel, @logmsg
import Random

Progress = LogLevel(-2)

"""
    Locator

Contains a `path` to an object within a [`Schedule`](@ref).

As a `Schedule` is executed, `Locator`s will be bound (with names starting on
"^") alongside the explicitly defined bindings to record the path within the
`Schedule` where the definition was evaluated. These "source" bindings can
therefore be referenced in [`Template`](@ref)s, and are further used to wrap
evaluated [`Model`](@ref)s in [`Models.Wrapped`](@ref) to remember where they
were originally defined.
"""
struct Locator
    path::String
end

Specifications.paste(locator::Locator) =
    join("-$s" for s in split(locator.path, r"[+-/]") if length(s) > 0)

Models.describe(locator::Locator) =
    Models.Label("specification at '$(locator.path)'")

"""
Wraps a non-`Schedule` `Model` to be invoked in the process of executing a
`Schedule`, adding additional behavior around the forwarded invocation.

# Invocation

    (primitive!::Primitive)(
        x,
        Δt;
        path,
        trace = nothing,
        dryrun = nothing,
        consolidated_progress = nothing,
        context...,
    )

Delegate to another `Model` `primitive!.f!`, adding pre- and post-processing.

This produces a single simulation segment; it
- converts the simulation state to the representation required by the wrapped
  `Model`,
- potentially reseeds `x`'s randomness,
- reports progress via `@logmsg` or `consolidated_progress`, and
- if `trace` is given, calls it back with the new simulation state and appends
  various ancillary information in that call, including `into` to signal if and
  where results should be saved.

The reseeding happens if `seed` is set on construction, which is always the case
when the `Primitive` is executed as part of a `Schedule`. The seed is then
accumulated on descent by appending any iteration indices. On invocation, the
`Primitive`'s seed will further be suffixed by `"+"` and then the current model
time `Models.t(x)` (if it is nonzero) to ensure that every segment has a unique
seed even if the `Primitive` is nested within a repeating `Scope`, unless
otherwise specified by the user. There are two exceptions to this policy.

Firstly, as a convenience to enable iterative development of a schedule
specification, the first `Primitive`'s `seed` in any `Sequence` will be
automatically replaced by its parent's, potentially recursively (that is,
trailing `"-1"`s will be clipped); in this way, breaking up any `Primitive` into
a new `Sequence` will ensure that its first segment keeps its trajectory
unchanged up to the first new segment boundary.

Secondly, reseeding is suppressed in a directly enclosing repeating `Scope`,
except for its first invocation. This only applies to regularly sampled stepping
(i.e., `{"to" ⟨number⟩, "step": ⟨number⟩}`) and is meant to align the
trajectories with the contiguous case.

The wrapped models are expected to retain intermediate results for their last
invocation in the simulation state `x` if `into` is not `nothing` so they can be
saved during the `trace` callback, after which they will be dropped.

If `dryrun` is given, execution short-circuits by calling that back instead. If
`consolidated_progress` is given, progress logging is suppressed and the
argument is forwarded to `primitive!.f!` and `trace`, instructing them to only
report progress by calling back.
"""
@kwdef struct Primitive <: Model{Any}
    f!::Model
    skip::Float64 = 0.0
    path::String
    bindings::Dict{Symbol, Any}
    into::Union{String, Nothing} = sink(bindings)
    seed::Union{String, Nothing} = replace(bindings[:seed], r"(-1)+$" => "")
end

Models.unwrap(primitive!::Primitive) = Models.unwrap(primitive!.f!)
Models.describe(primitive!::Primitive) = Models.describe(primitive!.f!)

function (primitive!::Primitive)(
    x,
    Δt::Float64;
    path,
    trace = nothing,
    dryrun = nothing,
    consolidated_progress = nothing,
    done = 0.0,
    context...,
)
    verbose = consolidated_progress === nothing

    f! = Models.unwrap(primitive!)
    if primitive!.skip > 0.0
        Δt = min(Δt, primitive!.skip)
    end

    if dryrun !== nothing
        x′ = x isa FlatState ? x : FlatState(t = Models.t(x))
        if f! isa Models.Instant
            Δt = 0.0
        end
        dryrun(primitive!, x′, Δt; path, primitive!.into, context...)
        if isfinite(Δt)
            x′.t += Δt
        end
        return x′
    end

    verbose && @logmsg(
        Progress,
        :adapting,
        at = path,
        todo = "$(nameof(typeof(x))) to $(nameof(typeof(f!))) \
            ($(primitive!.path))",
    )
    from = Models.t(x)
    seed = primitive!.seed
    if seed !== nothing && done == 0.0
        Random.seed!(Models.randomness(x), from > 0.0 ? "$seed+$from" : seed)
    end
    x = Models.adapt!(x, f!)

    verbose && @logmsg Progress :advancing at = path
    record = primitive!.into !== nothing
    traced = (;
        path,
        from,
        primitive!,
        primitive!.into,
        consolidated_progress,
    )
    if trace === nothing
        x = f!(x, Δt; path, consolidated_progress, context...)
    elseif primitive!.skip > 0.0
        x = f!(x, Δt; path, consolidated_progress, context...)
        trace(x; traced..., into = nothing)
        record && trace(x; traced..., from = Models.t(x))
    else
        x = f!(x, Δt; path, consolidated_progress, context..., record)
        trace(x; traced...)
    end
    Models.empty_trajectory!(x)

    verbose && @logmsg Progress :done at = path
    x
end

struct Branched
    stem
    branches::Vector
end

Branched(x) = Branched(x, [])
Models.randomness(x::Branched) = Models.randomness(x.stem)
Models.t(x::Branched) = Models.t(x.stem)
each_event(callback::Function, x::Branched) = each_event(callback, x.stem)

Models.FlatState(x::Branched) = Models.FlatState(x.stem)
Models._adapt!(x::Branched, ::Model{Branched}, ::Val{false}) = x
Models._adapt!(x::Branched, f!::Model, copy::Val) =
    Models._adapt!(x.stem, f!, copy)

"""
    Schedule{S <: Specification} <: Model{Any}

A `Model` that advances the simulation by delegating to other `Model`s for a
sequence of simulation segments that is organized according to a `specification`
of type `S`.

This process may involve multiple levels of recursively constructing and
executing `Schedule`s, reflecting the potentially nested structure of
`specification`. When descending on the specification, the corresponding
`Schedule`s may accumulate a collection of named values (`bindings`) that may be
inserted in place of free references within the nested specifications before
their interpretation; in this way, `Schedule`s support a limited amount of
templating.

# Specification

See [Schedule specification](@ref).

# Invocation

    (f!::Schedule{Template})(x, Δt::Float64; path, context...)

Expand the template, convert it to a `Model` and forward the call to it.

Effectively, this means either a recursion to sub-schedules (if the template
evaluates to a `Specification`), or the execution of primitive simulation
segments, either recording everything (if the template evaluates to a `Model`)
or only at the last timepoint (if it evaluates to a number).

Specifically, depending on the expanded value,
- if it is already a `Model`, it first gets wrapped in a `Wrapped` (to tag it
  with `f!.path` for later reference), and then further in a `Primitive` (which
  adds output and progress reporting when invoked). The latter's `into` field
  is determined from `f!.bindings[:into]`; if it is `"{channels}"`, `into` is
  set to `f!.bindings[:channel]` instead. Otherwise,
- if it is a `Specification`, it gets baked into a new `Schedule` (closing over
  `f!.bindings`), or
- if it is a number, a model gets looked up at `f!.bindings[:do]` and placed
  into a `Primitive` with `skip` set to the expanded number (resulting in a
  segment without output advancing to that timepoint, and another instant
  segment with output). This is a shortcut provided to handle the common case of
  discretely sampling along the time axis by setting any `step::Specification`
  to the desired step size, but it requires defining the `Model` as `:do` in an
  enclosing `Scope`. If `:do` is not bound, fall back to
  [`Wait`](@ref Models.Plumbing.Wait)).

---
    (f!::Schedule{Slice})(x, Δt::Float64; context...)

Look up a `Model` in `f!.bindings[:do]` and forward the call to it.

Read this as "step in infinitesimal slices until the simulatation budget `Δt` is
exhausted". If `:do` is not bound, fall back to
[`Wait`](@ref Models.Plumbing.Wait)).

---
    (f!::Schedule{Load})(x, Δt::Float64; load, context...)

Load a `Specification` from a JSON file, turn it into a `Schedule` and forward
the call to it.

`f!.path` is passed to `load`, which needs to be a function that returns data
structures as they would be produced by calling `JSON.parse` with
`dicttype = Dict{Symbol, Any}`.

---
    (f!::Schedule{<:Sequence})(x, Δt::Float64; context...)

Iterate the `Model`s specified by `f!.specification` (a `List` or `Each`) and
invoke them in sequence.

The exact behavior depends on whether `f!.branch` was set (by the directly
enclosing `Scope`):
- If so, simulation will not advance the state `x`, but will advance copies
  instead, one for each item in the sequence. The copies will share the same
  `randomness` instance, so because they draw from that randomness in order,
  their trajectories will start to differ at the branch (copy) time point. All
  advanced copies will be returned together with the original `x` as a
  `Branched` state so that they can optionally be merged (see
  [`Merge`](@ref Models.Scheduling.Merge)), but typically the branched
  components will instead be dropped downstream. (Note that by this point, their
  trajectories likely already have been `trace`d in the respective `Primitive`
  invocations.)
- Otherwise, the items are invoked in turn on the same state `x`. After each
  step, the remaining simulation budget `Δt` will be decreased by the advanced
  time interval. This means that steps may be invoked with `Δt == 0.0`, which
  typically means that dynamic `Model`s will have no effect but `Instant` models
  will always be applied.

The nested `Model`s' `path`s will be suffixed by their iteration index,
separated either by `"/"` if `f!.branch` was set or by `"-"` otherwise.
Additionally, `f!.bindings[:channel]` and `f!.bindings[:seed]` will be suffixed
by `"-"` and the iteration index. This extension will be skipped for `:channel`
if the `Sequence` contains pure replicates (i.e., it is an `Each` that does not
bind its iteration variable.)

---
    (f!::Schedule{Scope})(x, Δt::Float64; context...)

Advance by constructing a nested `Model`, optionally evaluating and adding new
`bindings` to its context, and either invoking it exactly once or, if
`bindings[:to]` is set, repeatedly until that simulation budget is exhausted.

The new `bindings` are determined by merging (the prior) `f!.bindings` and new
entries obtained from `f!.specification.definitions`. If `f!.barrier` is set,
this will only include `:rootseed`, `:seed`, `:into`, `:channel` and `:defaults`
from `f!.bindings`. In either case, the definitions may contain references to
`f!.bindings`, and new definitions will shadow prior bindings of the same name.
(`f!.barrier` is currently only set when parsing a `Load`/`:<` literal from the
JSON representation.)

The so extended bindings are then used to construct a new `Model` (a `Schedule`
or `Primitive`) from the `Specification` in `f!.step` and invoke it. The new
`Model`'s `path` will be suffixed by `"+"` to signify descending on a `Scope`,
unless `f!.branch` is set (because then the information is redundant since
branching can only be specified in a `Scope` and the next path component is then
guaranteed to start with `"/"`).

If `f!.specification.definitions[:to]` is not set (i.e. defined directly in this
`Scope`), the call is just forwarded to that new `Model`. Otherwise, the
simulation time budget `Δt` is clipped to that value and the new `Model` is then
invoked repeatedly, each time deducting the actually advanced simulation time
from `Δt`, until it is exhausted. (The invocation will pass the full remaining
`Δt` each time, but the nested `Model` is allowed to advance less than that, for
example because it is a `Schedule` that itself defines `:to`.)
"""
@kwdef struct Schedule{S <: Specification} <: Model{Any}
    specification::S
    bindings::Dict{Symbol, Any} = Dict{Symbol, Any}(
        :rootseed => "",
        :seed => "",
        :into => "",
        :channel => "",
        :defaults => Dict{Symbol, Any}(),
    )
    branch::Bool = false
    path::String = ""
end

(f!::Schedule)(x = FlatState(); context...) = f!(x, Inf; context...)

Models.adapt!(x, ::Schedule, _copy::Val{false}) = x

function sink(bindings::Dict{Symbol, Any})
    into = get(bindings, :into, nothing)
    if into == "{channels}"
        into = bindings[:channel]
    end
    into
end

evaluate_bindings(f!::Schedule{Scope}) = merge(
    # from outer scope:
    if f!.specification.barrier
        # Implicitly retain seed, output control and defaults, even if this
        # scope has a barrier (because its step is a Load); these are always
        # available in the loaded specification.
        Dict{Symbol, Any}(
            keep => f!.bindings[keep]
            for keep in (:rootseed, :seed, :into, :channel, :defaults)
            if haskey(f!.bindings, keep)
        )
    else
        f!.bindings
    end,

    # evaluated definitions:
    Dict{Symbol, Any}(
        name => evaluate(specification, path = "$(f!.path).$name"; f!.bindings)
        for (name, specification) in f!.specification.definitions
    ),

    # evaluated definitions' paths
    Dict{Symbol, Any}(
        Symbol("^$name") => Locator(f!.path)
        for name in keys(f!.specification.definitions)
    )
)

track_model(model::Model; locator) = Models.Wrapped(definition = locator; model)
track_model(x; _...) = x

evaluate(template::Template; bindings, path) = track_model(
    Specifications.expand(template; bindings),
    locator = Locator(path),
)
evaluate(x; _...) = x

model(template::Template; bindings, branch, path) =
    model(evaluate(template; bindings, path); bindings, branch, path)

model(specification::Specification; bindings, branch, path) =
    Schedule(; specification, bindings, branch, path)

function model(f!::Model; bindings, branch, path)
    branch && error("cannot branch here: not a Sequence")
    Primitive(; f!, path, bindings)
end

function model(at::Real; bindings, branch, path)
    branch && error("cannot branch here: not a Sequence")
    f! = get(bindings, :do, Wait())
    path = "$(bindings[Symbol("^do")].path).do"
    Primitive(skip = Float64(at); f!, path, bindings)
end

function item(
    step;
    i,
    x = nothing,
    as = Symbol(""),
    isinline,
    bindings,
    path,
)
    redefinitions = Dict{Symbol, Any}(:seed => "$(bindings[:seed])-$i")
    if isinline || as != Symbol("")
        redefinitions[:channel] = "$(bindings[:channel])-$i"
    end
    if as != Symbol("")
        redefinitions[as] = evaluate(x, path = "$path$i"; bindings)
    end

    bindings = merge(bindings, redefinitions)
    model(step; bindings, branch = false, path = "$path$i")
end

models(list::List; bindings, path) = (
    item(specification; i, bindings, path, isinline = true)
    for (i, specification) in enumerate(list.items)
)

models(each::Each; bindings, path) = (
    item(each.step; i, x, each.as, bindings, path, isinline = false)
    for (i, x) in enumerate(evaluate(each.items; bindings, path))
)

load_schedule(f!::Schedule{Load}; load) = Schedule(
    specification = Specification(
        load(f!.specification.path),
        bound = Set(keys(f!.bindings)),
    );
    f!.bindings,
    f!.branch,
    f!.path,
)

consolidate_progress(::Model) = false
consolidate_progress(f!::Primitive) = !(Models.unwrap(f!.f!) isa Instant)
consolidate_progress(f!::Schedule{Scope}) =
    f!.specification.step isa Slice &&
    !haskey(f!.specification.definitions, :do)

function (f!::Schedule{Slice})(x, Δt::Float64; context...)
    path =
        if haskey(f!.bindings, :do)
            "$(f!.bindings[Symbol("^do")].path).do"
        else
            f!.path
        end

    # When this Slice is nested in a Scope (which is the normal case), context
    # contains path, which overrides the following empty path:
    model(
        get(f!.bindings, :do, Wait());
        f!.bindings,
        f!.branch,
        path,
    )(x, Δt, path = ""; context...)
end

(f!::Schedule{Template})(x, Δt::Float64; path, context...) =
    model(
        evaluate(f!.specification; f!.bindings, path);
        f!.bindings,
        f!.branch,
        path,
    )(x, Δt; path, context...)

function (f!::Schedule{Scope})(
    x,
    Δt::Float64;
    trace = nothing,
    consolidated_progress = nothing,
    context...,
)
    f!.branch && error("cannot branch at '$(f!.path)': not a Sequence")
    verbose = consolidated_progress === nothing

    verbose && @logmsg Progress :preparing at = f!.path
    x = Models.adapt!(x, f!)  # potentially unwrap Branched
    bindings = evaluate_bindings(f!)
    path = "$(f!.path)$(f!.specification.branch ? '/' : '+')"
    step! = model(
        f!.specification.step;
        bindings,
        f!.specification.branch,
        path,
    )

    if haskey(bindings, :to) && bindings[Symbol("^to")].path == f!.path
        Δt = min(Δt, bindings[:to])
        verbose && @logmsg Progress :preparing at = f!.path todo = Δt
        if verbose && consolidate_progress(step!)
            consolidated_progress =
                (message; todo = nothing, done = 0) ->
                    @logmsg Progress message at = f!.path todo done
        end
        done = 0.0
        while 0.0 < Δt
            current = Models.t(x)
            x = step!(
                x,
                Δt;
                trace,
                consolidated_progress,
                context...,
                path,
                done,
            )
            advance = Models.t(x) - current
            0.0 < advance || error("cannot progress")
            Δt -= advance
            done += advance
            if consolidated_progress === nothing
                @logmsg Progress :repeating at = f!.path done
            end
        end
    else
        verbose && @logmsg Progress :descending at = f!.path
        x = step!(x, Δt; trace, context..., path)
    end

    if haskey(f!.specification.definitions, :flush) && trace !== nothing
        trace(into = sink(bindings), flush = bindings[:flush]; f!.path)
    end
    verbose && @logmsg Progress :done at = f!.path

    x
end

function (f!::Schedule{<:Sequence})(x, Δt::Float64; context...)
    path = f!.branch ? f!.path : "$(f!.path)-"
    @logmsg Progress :preparing at = path

    steps = models(f!.specification; f!.bindings, path)

    @logmsg Progress :iterating at = path todo = length(steps)
    if f!.branch
        x = Branched(x)
        for (i, step!) in enumerate(steps)
            x′ = Models.adapt!(x.stem, step!, copy = true)
            x′ = step!(x′, Inf; context..., path = "$path$i")
            # In case the branch itself ended with a Branched state: it had its
            # chance for Merge... the following will implicitly unwrap it.
            push!(x.branches, FlatState(x′))
            @logmsg Progress :iterating at = path done = i
        end
    else
        for (i, step!) in enumerate(steps)
            current = Models.t(x)
            x = step!(x, Δt; context..., path = "$path$i")
            Δt -= Models.t(x) - current
            @logmsg Progress :iterating at = path done = i
        end
    end

    @logmsg Progress :done at = path

    x
end

(f!::Schedule{Load})(x, Δt::Float64; load, context...) =
    load_schedule(f!; load)(x, Δt; load, context..., f!.path)

"""
    reify(x, path; load = nothing)

Recreate an object by repeatedly descending on the definition object `x`, as
selected by `path`, expanding the required definitions along the way.

When called directly, `x` will typically be a `Schedule`, but it doesn't have to
be: As a convenience, `reify` can index into `AbstractVector`s, `AbstractDict`s
and other objects (by accessing their indices, keys or properties).

Reification will follow the same rules of descent though the definition object
as the corresponding direct invocation, but instead of walking the full tree
will only descend on one branch per inner node, as selected by `path`,
implicitly reifying further definition objects along the way.

If any of the intermediate definition objects are of type `Schedule{Load}`, the
`load` keyword must be given, analogously to invoking the [`Schedule`](@ref), so
that `reify` knows how to execute the `Load`.
"""
function reify end

reify(x; context...) = reify(x, ""; context...)

reify(x, path::AbstractString; _...) =
    isempty(path) ? x : Specifications.pluck(x, split(path, '.'))

reify(f!::Schedule{Slice}, path::AbstractString; context...) = reify(
    model(
        get(f!.bindings, :do, Wait());
        f!.bindings,
        f!.branch,
        path,
    ),
    path;
    context...,
)

reify(f!::Schedule{Template}, path::AbstractString; context...) = reify(
    model(
        evaluate(f!.specification; f!.bindings, f!.path);
        f!.bindings,
        f!.branch,
        f!.path,
    ),
    path;
    context...,
)

function reify(f!::Schedule{Scope}, path::AbstractString; context...)
    isempty(path) && return f!
    token = path[1]
    tail = path[2:end]
    bindings = evaluate_bindings(f!)
    if token == '.'
        reify(bindings, tail)
    elseif token == '+' || token == '/'
        step! = model(
            f!.specification.step;
            bindings,
            f!.specification.branch,
            path = "$(f!.path)$token",
        )
        reify(step!, tail; context...)
    else
        error("cannot descend to '$path' in Scope at '$(f!.path)'")
    end
end

function reify(f!::Schedule{<:Sequence}, path::AbstractString; context...)
    isempty(path) && return f!
    m = match(r"^(-?)(\d+)(.*)", path)
    if m !== nothing
        prefix, head, tail = m
        reify(f!, prefix, parse(Int, head), tail; context...)
    else
        reify(f!.bindings, path; context...)
    end
end

function reify(
    f!::Schedule{List},
    prefix::AbstractString,
    i::Int,
    tail::AbstractString;
    context...,
)
    path = "$(f!.path)$prefix"
    step = f!.specification.items[i]
    step! = item(step; i, isinline = true, f!.bindings, path)
    reify(step!, tail; context...)
end

function reify(
    f!::Schedule{Each},
    prefix::AbstractString,
    i::Int,
    tail::AbstractString;
    context...,
)
    each = f!.specification
    step! = item(
        each.step;
        i,
        x = evaluate(each.items; f!.bindings, f!.path)[i],
        each.as,
        isinline = false,
        f!.bindings,
        path = "$(f!.path)$prefix",
    )
    reify(step!, tail; context...)
end

reify(f!::Schedule{Load}, path::AbstractString; load, context...) =
    reify(load_schedule(f!; load), path; load, context...)

function reify(primitive!::Primitive, path::AbstractString; _...)
    isempty(path) && return primitive!
    token = path[1]
    tail = path[2:end]
    if token == '+'
        reify(primitive!.f!, tail)
    elseif token == '.'
        reify(primitive!.bindings, tail)
    else
        error("cannot descend to '$path' in Primitive at '$(primitive!.path)'")
    end
end

Specifications.locate(f!::Schedule, path::AbstractString) = Specifications.locate(f!.specification, path)

Specifications.set(f!::Schedule, path::AbstractString, new::Specification) = Schedule(;
    specification = Specifications.set(f!.specification, path, new),
    f!.bindings, f!.branch, f!.path
)

prebuild(f!::Schedule, path::AbstractString) =
    Specifications.set(f!, path, Template(reify(f!, path)))

"""
    Merge <: Instant{Branched}

Instantly collapse a `Branched` state to a `FlatState` by replacing the `stem`
with aggregated counts from the `branches`.

# Specification

Currently, only `+` is supported as aggregation. Specified in JSON as
`{"{merge}": "+"}`.

# Invocation

    (f!::Merge)(x::Branched, _Δt::Float64; _...)

Use the function `f!.merge` to aggregate all of the `x.branches` and return a
new `FlatState` with the aggregated counts, but retaining `x.stem.t` and
`x.stem.randomness`.
"""
struct Merge <: Instant{Branched}
    merge::Function
end

Merge(operation::AbstractString) = operation |> Symbol |> Val |> Merge
Merge(::Val{:+}) = Merge(+)

Specifications.constructor(::Val{:merge}) = Merge

function (f!::Merge)(x::Branched, _Δt::Float64; _...)
    accumulator = FlatState(
        t = Models.t(x.stem),
        randomness = Models.randomness(x.stem),
    )
    for b in x.branches
        mergewith!(f!.merge, accumulator.counts, b.counts)
    end
    accumulator
end

end
