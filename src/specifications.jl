module Specifications

using Base: Fix2

range_(x) = range(; x...)
logrange_(x) = logrange(x...)

"""
    constructor(name::Symbol)

Select by `name` and return a function that accepts a substituted `Template`
value and constructs an object of the selected kind from it.

This function's methods define which object literals can be used in the JSON
specification language and how they should be interpreted after template
substitution. Each method must accept a single argument of a type as it would
be produced by calling `JSON.parse(..., dicttype = Dict{Symbol, Any})`.
"""
constructor(name::Symbol) = constructor(Val(name))
constructor(name::AbstractString) = constructor(Symbol(name))
constructor(::Val{Symbol("")}) = Dict{String, Any}
constructor(::Val{:range}) = range_
constructor(::Val{:logrange}) = logrange_

"""
Abstract supertype of all syntactic elements of the scheduling language.

# Construction

    Specification(x; bound::Set{Symbol} = Set{Symbol}(), as::Symbol = :step)

Construct a `Specification` from nested `Dict`/`Vector` objects such as they are
obtained by loading JSON via `JSON.parse(..., dicttype = Dict{Symbol, Any})`.

`Specification` recursively interprets a JSON document (or part thereof) and
returns a `Specification` subtype depending on the type and shape of the JSON.
It will be interpreted as either a `:step`, a `:value` or `:items`; the
top-level document will be interpreted as a `:step`. Specifically:
- When expecting a `:step`:
  - If `x` is a `Vector` (JSON Array), it will be parsed as a [`List`](@ref) of
    `:step` `Specification`s.
  - If `x` is a `Dict` (JSON Object), the earliest matching rule of the
    following applies:
    1. (*reference literal*) If `x` contains a `:\$` key (JSON name `"\$"`), it
        is parsed as a [`Template`](@ref) expanding to the binding referenced by
        the corresponding value.
    2. (*load literal*) If `x` contains a `:<` key (JSON name `"<"`), it is
        parsed as a [`Load`](@ref) of the file referenced by the corresponding
        value. This `Load` will be wrapped in a [`Scope`](@ref) that has
        `barrier` set and collects all the other mappings in `x` as
        `definitions`, interpreting the mapped values as `:value`
        `Specification`s.
    3. (*template literal*) If `x` contains a single mapping, and that mapping's
        key is enclosed in braces (JSON names `"{...}"`), it is parsed as a
        [`Template`](@ref) that is expanded by transforming the substituted
        value using a function returned by passing the key (without the braces)
        to [`Specifications.constructor`](\
          @ref GeneRegulatorySystems.Specifications.constructor\
        ).
    4. (*each*) If `x` contains an `:each` key (JSON name `"each"`), it is
        parsed as an [`Each`](@ref). The iterable `items` are defined by
        `x[:each]` interpreted as an `:items` `Specification`; the `step` is
        defined by `x[:step]` interpreted as a `:step` `Specification`, and the
        index variable name is optionally defined by `x[:as]`. If there are any
        other mappings in `x`, they will be collected and the `Each` wrapped in
        a [`Scope`](@ref) using these definitions. In other words, the
        corresponding definitions are available in the `items` and `step`
        definitions, but cannot refer to the index variable.
    5. (*scope*) If `x` is not empty, it is parsed as a [`Scope`](@ref), with
        the `step` defined by `x[:step]` interpreted as a `:step`
        `Specification` (defaulting to `Slice()`) and `branch` optionally set by
        `x[:branch]`. All the other mappings in `x` are collected as
        `definitions`, interpreting the mapped values as `:value`
        `Specification`s.
    6. (*slice*) Otherwise `x` is empty and is parsed as
        [`Slice()`](@ref Slice).
  - Otherwise, `x` is parsed as a [`Template`](@ref) expanding to `x`.
- When expecting `:items`:
  - If `x` is a `Dict` (JSON Object), the earliest matching rule of the
    following applies:
    1. (*reference literal*) as above in the `:step` case
    2. (*template literal*) as above in the `:step` case
    3. (*range literal*) Otherwise `x` is parsed as a [`Template`](@ref) that is
        expanded by calling Julia's `range` function, splatting the substituted
        value as keyword arguments.
  - Otherwise, `x` is parsed as a [`Template`](@ref) expanding to `x`.
- When expecting a `:value`:
  - If `x` is a `Dict` (JSON Object), the earliest matching rule of the
    following applies:
    1. (*reference literal*) as above in the `:step` case
    2. (*load literal*) as above in the `:step` case
    3. (*template literal*) as above in the `:step` case
    4. Otherwise, `x` is parsed as a [`Template`](@ref) expanding to `x`.
  - Otherwise, `x` is parsed as a [`Template`](@ref) expanding to `x`.
"""
abstract type Specification end

"""
    Slice <: Specification

Empty singleton that represents an infinitesimal-time step in the simulation.

It acts as a sentinel element in the scheduling language and roughly has the
role of `Nothing`.
"""
struct Slice <: Specification end

function references(s::AbstractString, bound::Set{Symbol})
    found = Set(map(Symbol ∘ first, eachmatch(r"\$\{(\^?[^\$\{\}\\]+)\}", s)))
    found ⊆ bound || error("undefined references: ", setdiff(found, bound))
    found
end

function references(x::AbstractDict{Symbol}, bound::Set{Symbol})
    found = mapreduce(
        Fix2(references, bound),
        union,
        values(filter(!=(:$) ∘ first, x)),
        init = Set{Symbol}(),
    )

    if haskey(x, :$)
        direct = Symbol(x[:$] isa AbstractVector ? first(x[:$]) : x[:$])
        found = union(found, Set([direct]))
    end

    found ⊆ bound || error("undefined references: ", setdiff(found, bound))
    found
end

references(xs::AbstractVector, bound::Set{Symbol}) =
    mapreduce(Fix2(references, bound), union, xs, init = Set{Symbol}())

references(_x, _bound::Set{Symbol}) = Set{Symbol}()

@doc raw"""
    Template <: Specification

Contains instructions for instantiating (*"expanding"*) a value from a
definition.

When stepping through a `Schedule`, `expand`ing `Template`s produces all
non-`Specification` values that influence the schedule's behavior, which
includes all primitive `Model`s.

The definition may contain references to named bindings. When `expand`ed, these
references are first replaced by their values, and then the function held in the
`Template`'s `constructor` field is called on the result. References come in two
forms:
- Any `Dict` `x` that contains a `:$` key will be replaced by the object
  addressed by `x[:$]`. If that value is a `String`, it refers to the binding
  of that same name. If it is a `Vector`, it refers to a nested object addressed
  by its items interpreted as path components; each item will in turn descend
  by accessing the respective key, index or property. If the found value `x′` is
  also a `Dict` and there are other mappings in `x`, they are merged into (a
  shallow copy of) `x′`, overriding previously existing mappings.
- Any `String` `s` that contains substrings of the form `"${binding}"` will be
  replaced by a `String` with that reference substituted by the respective
  binding's `repr` for `String`s and `Number`s, or the literal `"__omitted__"`
  otherwise.

Substitution will in geneneral not return independent objects but rather alias
intermediate `Dict`s, `Vector`s and other objects into the substituted objects
if they contain no substitutions of of their own. In other words, the produced
data structures are treated as *persistent* (and therefore *immutable*) during
expansion.
"""
@kwdef struct Template <: Specification
    value
    constructor
    free::Set{Symbol}

    Template(value, constructor, free::Set{Symbol}) =
        # TODO: Perhaps we should add a (single element by default, but
        # configurable) cache to Template? This would alleviate the need for
        # this distinction:
        if isempty(free)
            # Eagerly construct the payload since value is already fully
            # closed; this will prevent unnecessary reconstructions down the
            # line for this common case.
            new(constructor(value), identity, free)
        else
            new(value, constructor, free)
        end
end

Template(value; constructor = identity, bound = Set{Symbol}()) =
    Template(free = references(value, bound); value, constructor)

reference(name::Symbol; bound = Set([name])) =
    Template(Dict(:$ => name); bound)

"""
    Scope <: Specification

Contains named value `definitions` (mostly `Template`s) that apply to a nested
context specified by `step` (which is a `Specification`).

In that sense it is equivalent to a lexical scope in any programming language,
but it can additionally be thought of conceptually as applying to a range on the
time axis during simulation, either filling its parent range or, if
`definitions[:to]` is set, limited to that duration.

`Template` values in `definitions` may contain references to bindings from a
surrounding scope (`Scope` or `Each`, or by inclusion on `Schedule`
construction).

The `definitions` shadow bindings of the same name from a surrounding scope. If
the `barrier` flag is set, any bindings from a surrounding scope will not be
available in the nested `step`; only new bindings from `definitions` will be
included (with some exceptions, see
[`Schedule{Scope}`](@ref GeneRegulatorySystems.Models.Scheduling.Schedule)).

If the `branch` flag is set, the `step` must be a `Sequence` specification and
should then be interpreted as specifying independent simulation branches instead
of the default behavior of acting on the same simulation state one after the
other.
"""
@kwdef struct Scope <: Specification
    definitions::Dict{Symbol, Specification}
    step::Specification = Slice()
    barrier::Bool = false
    branch::Bool = false
    free::Set{Symbol} = union(
        mapreduce(free, union, values(definitions), init = Set{Symbol}()),
        barrier ? Set{Symbol}() : setdiff(free(step), keys(definitions)),
    )
end

"""
    Sequence <: Specification

Abstract supertype of specifications that can be iterated.

The meaning of the specified `Sequence` (to be interpreted when executing a
schedule) depends on whether a directly enclosing `Scope` has its `branch` flag
set.
"""
abstract type Sequence <: Specification end

"""
    List <: Sequence

Represents a static list of specifications.
"""
@kwdef struct List <: Sequence
    items::Vector{Specification}
    free::Set{Symbol} = mapreduce(free, union, items, init = Set{Symbol}())
end

List(items) = List(; items)

"""
    Each <: Sequence

Represents a sequence of specifications defined implicitly by setting a named
binding, in turn, to each item from an ordered collection of `item`s and
evaluating a nested `Specification` `step` in the resulting context.

The iteration variable, defined by `as`, does not necessarily need to be named
or used; if it is not, the `Each` effectively represents a repetition of the
same (nested) specification.
"""
@kwdef struct Each <: Sequence
    items::Template
    as::Symbol = Symbol("")
    step::Specification
    free::Set{Symbol} = union(
        free(items),
        setdiff(free(step), Set(as == Symbol("") ? Symbol[] : [as])),
    )
end

"""
    Load <: Specification

Represents an instruction to load, parse and insert a `Specification` from a
file.

Its `path` is relative and given context when invoking the containing
[`Schedule{Load}`](@ref GeneRegulatorySystems.Models.Scheduling.Schedule) via
the `load` function argument.
"""
struct Load <: Specification
    path::AbstractString
end

free(::Slice) = Set{Symbol}()
free(::Load) = Set{Symbol}()
free(s::Specification) = s.free

"""
    locate(x, path) -> Specification

    Return the specification in `x` at `path`.
"""
function locate end

locate(s::Specification, path::AbstractString) =
    isempty(path) ? s :
    error("cannot descend to '$path' in $(nameof(typeof(s)))")

function locate(s::Scope, path::AbstractString)
    isempty(path) && return s
    token = path[1]
    tail = path[2:end]
    if token == '.'
        m = match(r"^([^+/.]+)(.*)", tail)
        locate(s.definitions[(Symbol(m[1]))], m[2])
    elseif token == '+' || token == '/'
        locate(s.step, tail)
    else
        error("cannot descend to '$path' in Scope")
    end
end

function locate(list::List, path::AbstractString)
    isempty(path) && return list
    m = match(r"^(-?)(\d+)(.*)", path)
    m === nothing && error("cannot descend to '$path' in List")
    _, head, tail = m
    locate(list.items[parse(Int, head)], tail)
end

function locate(each::Each, path::AbstractString)
    isempty(path) && return each
    m = match(r"^(-?)(\d+)(.*)", path)
    m === nothing && error("cannot descend to '$path' in Each")
    _, _, tail = m
    locate(each.step, tail)
end

"""
    set(x, path::AbstractString, new::Specification)

    Substitute the specification in `x` at `path` with `new`.
"""
function set end

set(::Specification, path::AbstractString, new::Specification) =
    isempty(path) ? new : error("cannot descend to '$path' in $(nameof(typeof(s)))")

function set(s::Scope, path::AbstractString, new::Specification)
    isempty(path) && return new
    token = path[1]
    tail = path[2:end]
    if token == '.'
        m = match(r"^([^+/.]+)(.*)", tail)
        name = Symbol(m[1])
        new_definitions = copy(s.definitions)
        new_definitions[name] = set(s.definitions[name], m[2], new)
        Scope(; definitions = new_definitions, s.step, s.barrier, s.branch)
    elseif token == '+' || token == '/'
        Scope(; s.definitions, step = set(s.step, tail, new), s.barrier, s.branch)
    else
        error("cannot descend to '$path' in Scope")
    end
end

function set(s::List, path::AbstractString, new::Specification)
    isempty(path) && return new
    m = match(r"^(-?)(\d+)(.*)", path)
    m === nothing && error("cannot descend to '$path' in List")
    _, head, tail = m
    i = parse(Int, head)
    items = copy(s.items)
    items[i] = set(s.items[i], tail, new)
    List(; items)
end

function set(s::Each, path::AbstractString, new::Specification)
    isempty(path) && return new
    m = match(r"^(-?)(\d+)(.*)", path)
    m === nothing && error("cannot descend to '$path' in Each")
    _, _, tail = m
    Each(; s.items, s.as, step = set(s.step, tail, new))
end

paste(x::AbstractString) = x
paste(x::Number) = repr(x)
paste(::Any) = "__omitted__"

pluck(x, path::AbstractVector) = foldl(pluck, path, init = x)
pluck(x::AbstractDict{Symbol}, key::AbstractString) = pluck(x, Symbol(key))
pluck(x::AbstractDict{Symbol}, key::Symbol) = key == Symbol("") ? x : x[key]
pluck(xs::AbstractVector, key::AbstractString) = pluck(xs, parse(Int, key))
pluck(xs::AbstractVector, i::Int) = xs[i]
pluck(x, key::AbstractString) =
    isempty(key) ? x : getproperty(x, Symbol(key))

function substituted(target::AbstractString; bindings::AbstractDict{Symbol})
    new = replace(
        target,
        ("\${$key}" => paste(value) for (key, value) in bindings)...
    )

    new == target ? nothing : Some(new)
end

substituted(_target; _...) = nothing

function substituted(target::AbstractVector; bindings::AbstractDict{Symbol})
    changed = Dict(
        i => something(value)
        for (i, value) in zip(
            LinearIndices(target),
            substituted.(target; bindings)
        )
        if !isnothing(value)
    )

    isempty(changed) ? nothing : Some([
        get(changed, i, old) for (i, old) in pairs(IndexLinear(), target)
    ])
end

function substituted(
    target::AbstractDict{Symbol};
    bindings::AbstractDict{Symbol},
)
    if haskey(target, :$)
        specializations = filter(!=(:$) ∘ first, target)
        prototype = pluck(bindings, target[:$])  # fail on undefined reference
        if isempty(specializations)
            Some(prototype)
        else
            Some(merge(prototype, substitute(specializations; bindings)))
            #    ^ fail unless template is a Dict
        end
    else
        changed = Dict(
            key => something(value)
            for (key, value) in zip(
                keys(target),
                substituted.(values(target); bindings)
            )
            if !isnothing(value)
        )

        isempty(changed) ? nothing : Some(Dict{Symbol, Any}(
            key => get(changed, key, old) for (key, old) in target
        ))
    end
end

substitute(target; bindings::AbstractDict{Symbol}) =
    something(substituted(target; bindings), Some(target))

function expand(
    template::Template;
    bindings::AbstractDict{Symbol}
)
    bindings = Dict{Symbol, Any}(k => bindings[k] for k in template.free)
    template.constructor(substitute(template.value; bindings))
end

maybe_reference(x::AbstractDict{Symbol}; bound) =
    haskey(x, :$) ? Template(x; bound) : nothing

maybe_load(x::AbstractDict{Symbol}; bound) =
    if haskey(x, :<)
        Scope(
            definitions = Dict{Symbol, Specification}(
                key => Specification(value, as = :value; bound)
                for (key, value) in x
                if key != :<
            ),
            step = Load(x[:<]),
            barrier = true,
        )
    else
        nothing  # x is not a Load literal
    end

function maybe_tagged_template(x::AbstractDict{Symbol}; bound)
    length(x) == 1 || return
    key, value = only(x)
    m = match(r"\{(?<kind>.*)\}", String(key))
    m !== nothing || return

    Template(value, constructor = constructor(m[:kind]); bound)
end

function maybe_scope(x::AbstractDict{Symbol}; bound)
    isempty(x) && return

    definitions = Dict{Symbol, Specification}(
        key => Specification(value, as = :value; bound)
        for (key, value) in x
        if key != :step && key != :branch
    )
    bound = union(
        bound,
        keys(definitions),
        Set(Symbol("^$k") for k in keys(definitions)),
    )
    step = haskey(x, :step) ? Specification(x[:step]; bound) : Slice()
    branch = get(x, :branch, false)

    Scope(; definitions, step, branch)
end

function maybe_each(x::AbstractDict{Symbol}; bound)
    haskey(x, :each) || return

    definitions = Dict{Symbol, Specification}(
        key => Specification(value, as = :value; bound)
        for (key, value) in x
        if key ∉ [:each, :as, :step, :branch]
    )
    bound = union(
        bound,
        keys(definitions),
        Set(Symbol("^$k") for k in keys(definitions)),
    )
    items = Specification(x[:each], as = :items; bound)

    if haskey(x, :as)
        as = Symbol(x[:as])
        bound = union(bound, [as, Symbol("^$as")])
    else
        as = Symbol("")
    end

    step = haskey(x, :step) ? Specification(x[:step]; bound) : Slice()
    each = Each(; items, as, step)

    if isempty(definitions) && !haskey(x, :branch)
        each  # elide unnecessary Scope
    else
        Scope(step = each, branch = get(x, :branch, false); definitions)
    end
end

Specification(x; bound::Set{Symbol} = Set{Symbol}(), as::Symbol = :step) =
    Specification(x, Val(as); bound)

Specification(value, ::Val; bound::Set{Symbol}) = Template(value; bound)

Specification(xs::AbstractVector, ::Val{:step}; bound::Set{Symbol}) =
    List(items = Specification.(xs; bound))

Specification(
    x::AbstractDict{Symbol},
    ::Val{:value};
    bound::Set{Symbol}
) = @something(
    maybe_reference(x; bound),
    maybe_load(x; bound),
    maybe_tagged_template(x; bound),
    Template(x; bound),
)

Specification(
    x::AbstractDict{Symbol},
    ::Val{:items};
    bound::Set{Symbol}
) = @something(
    maybe_reference(x; bound),
    maybe_tagged_template(x; bound),
    Template(x, constructor = range_; bound),
)

Specification(
    x::AbstractDict{Symbol},
    ::Val{:step};
    bound::Set{Symbol}
) = @something(
    maybe_reference(x; bound),
    maybe_load(x; bound),
    maybe_tagged_template(x; bound),
    maybe_each(x; bound),
    maybe_scope(x; bound),
    Slice(),
)

cast(::Type{T}, x::T; _...) where {T} = x
cast(T::Type, x::AbstractDict{Symbol}, ::Val{K}; context) where {K} =
    cast(fieldtype(T, K), x[K]; context)
cast(::Type{Symbol}, x; _...) = Symbol(x)
cast(T::Type{<:Real}, x::Real; _...) = convert(T, x)
cast(T::Type{<:Real}, x::AbstractString; _...) = parse(T, x)
cast(::Type{Vector{T}}, xs::Vector{T}; _...) where {T} = xs
cast(::Type{Vector{T}}, xs::AbstractVector; context = nothing) where {T} =
    cast.(T, xs; context)
cast(::Type{Dict{Symbol, T}}, x::AbstractDict{Symbol}; context) where {T} =
    Dict(key => cast(T, value; context) for (key, value) in x)
cast(T::Type, x::AbstractDict{Symbol}; context = x) = T(; (
    key => cast(T, x, Val(key); context)
    for key in keys(x)
    if hasfield(T, key)
)...)
cast(::Type{Union{Some{T}, Nothing}}, ::Nothing; _...) where {T} = nothing
cast(::Type{Union{Some{T}, Nothing}}, x; context) where {T} =
    Some(cast(T, x; context))
cast(  # disambiguation
    ::Type{Union{Some{T}, Nothing}},
    x::AbstractDict{Symbol};
    context,
) where {T} = Some(cast(T, x; context))

representation(x::Dict) = x
representation(xs::AbstractVector) = representation.(xs)
representation(x::Integer) = x
representation(x::AbstractFloat) = x
representation(x::AbstractString) = x
representation(x::Symbol) = string(x)
representation(x; simple = false, rest...) =
    representation(x, Val(simple); rest...)
function representation(x, ::Val{true}; omit_defaults = Pair{Symbol, Any}[])
    result = Dict{Symbol, Any}(
        key => representation(getproperty(x, key))
        for key in propertynames(x)
    )
    for (key, default) in omit_defaults
        if result[key] == default
            delete!(result, key)
        end
    end
    result
end

end
