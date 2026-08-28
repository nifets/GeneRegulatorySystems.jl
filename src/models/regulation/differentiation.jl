"""
Contains components to build `Models.V1`-based `SciML.JumpModel`s for gene
regulation that exhibit trajectories following a predefined differentiation
tree.

These models can be constructed by [`build`](@ref) from a [`Definition`](@ref)
or from the corresponding JSON specification.

They include a *core network* controlling differentiation according to a binary
tree with specified developmental timing and split ratios. For each node in the
tree, the resulting model will contain a signalling gene which indicates (by
high expression) that its branch was taken. In addition, an arbitrary `V1`
*peripheral network* may be specified and regulated from the core network.

The differentiation machinery requires setting up initial state (in addition to
the polymerases, ribosomes and proteasomes required anyway). The required
initial deposits can be triggered by first invoking an `Instant` bootstrap
model -- constructed by [`bootstrap`](@ref) -- before invoking the actual
regulation model.
"""
module Differentiation

using ..Models: Models, SciML, V1, Plumbing
import ..Specifications: constructor, cast, representation

# To reproduce the parameters in the following two functions, see:
#
#     contrib/differentiation_calibration/README.txt

proportion_rate(target::Float64) =
    0.0001003483 * (1.040447 / (1.011666 - target) - 1.0)^-0.024953787457163392

timing_factor(duration::Float64) =
    0.00003167 * max(600.0, duration)^-1.031

"""
    Transient

Defines a module of inter-regulated genes that engineer a time-controlled
bifurcation in the simulated trajectory.

Such a module consists of a `differentiator` gene that has higher expression
whenever the corresponding transient state has been reached or surpassed, and
a `timer` gene that controls when the trajectory should bifurcate into one of
the defined downstream states. The `Transient` may also contain a dimerization
`buffer` that can make timing more robust.

The (two) downstream states `next` and `alternative` can either be terminal
(differentiator) genes or `Transient`s of their own, and the target probability
of moving into the `next` state is given by `ratio`.

In the current implementation, initial deposits of molecules must be made for
the `timer` (and potentially the `buffer`) before the differentiation model is
invoked for the assembled machinery to work properly. The timings and split
proportions are calibrated for the default deposit amounts as defined in
`defaults.specification.json`, so for now these values should always be used.

# Specification

In JSON, the [`Definition`](@ref) containing a `Transient` (perhaps indirectly)
is usually specified as part of a [`Schedule`](@ref Models.Scheduling.Schedule)
and therefore has `"defaults"` available. While many `Transient` properties can
be set by the specification, most should typically be left at their defaults.
In that simple case, the `Transient` is specified as a JSON object
```
{
    "\$": ["defaults", "differentiation"],
    "duration": <duration>,
    "ratio": <ratio>,
    "next": <next>,
    "alternative": <alternative>
}
```
to mostly use the differentiation defaults and only override node-specific
attributes. Here, `<duration>` must be a JSON number specifying the delay (in
simulation time units) after entering this `Transient` when the trajectory
should start bifurcating to the downstream states `<next>` and `<alternative>`,
which specify either nested `Transient`s for transient states or
[`V1.Gene`](@ref)s for terminal states. Optionally specifying `<ratio>` as a
JSON number (in the unit range) changes the target probability to take the
`next` branch from its default value of `0.5`.

A special requirement applies to the initial (top-level) `Transient` in a
`Differentiation`: if its `differentiator` represents a changing quantity in
the system such a gene (with protein production and decay), the timing for that
transient state will be inaccurate. The easiest way to prevent this is to set
e.g. `"differentiator": "differentiator"` on that intial `Transient` to have it
refer to a molecule species `"differentiator"`, and to initially add a
sufficient amount of it when the simulation state is set up anyway, for example
like
```
{"{add}": {
    "\$": ["defaults", "bootstrap"],
    "differentiator": 500
}}
```
in the first step of a `Schedule`.

In cases where the defaults are not available, or where more control is needed,
the set of user-facing parameters that make up a `Transient` definition can be
specified in more detail as a JSON object
```
{
    "differentiator": <differentiator>,
    "duration": <duration>,
    "timer": <timer>,
    "buffer": [<forward>, <reverse>],
    "ratio": <ratio>,
    "next": <next>,
    "alternative": <alternative>,
    "timer_deposit": <timer_deposit>,
    "buffer_deposit": <buffer_deposit>,
}
```
where:
- `<differentiator>` specifies the [`V1.Gene`](@ref) whose increased expression
  should indicate that this `Transient` state has been reached. This may either
  specify the gene directly or alternatively a JSON string that names an already
  existing gene (e.g. in the peripheral network), which will in this case
  participate in the new reactions controlling the differentiation. As a special
  case, only for the initial (top-level) `Transient` in a `Differentiation`,
  `<differentiator>` as a JSON string may refer to an arbitrary (non-gene)
  molecular species to use as the upstream switch. The whole `Differentiation`
  is then enabled by the presence of this factor.
- `<duration>` must be a JSON number specifying the target time to remain in
  the transient state as described above.
- `<timer>` specifies a [`V1.Gene`](@ref) to represent the remaining time in
  this transient state.
- `<forward>` and `<reverse>` must be JSON numbers specifying the reaction rates
  of the timer dimerization buffer. If the containing JSON array is not set or
  empty, `timer` will have no such buffer.
- `<ratio>` must be a JSON number (defaulting to `0.5`) defining the trajectory
  split proportion as described above.
- `<next>` and `<alternative>` specify the follow-on stages that the dynamics
  should bifurcate into. They may either be specified as another `Transient`,
  as a [`V1.Gene`](@ref) to signify a terminal differentiator gene as described
  above, or as a JSON string referring to such a gene (e.g. in the peripheral
  network). The downstream differentiators will participate in the synthesized
  reactions controlling the differentiation. Unnamed differentiators will be
  automatically named depending on their parent differentiator, appending `"0"`
  for `<next>` and `"1"` for `<alternative>`.
- `<timer_deposit>` must be a JSON object mapping any of `"elongations"`,
  `"premrnas"`, `"mrnas"` and `"proteins"` to JSON (integer) numbers that the
  [`bootstrap`](@ref) model should initialize this `Transient`'s timer gene
  species to.
- `<buffer_deposit>` must be a JSON (integer) number that the
  [`bootstrap`](@ref) model should initialize the `buffer` to.
"""
@kwdef struct Transient
    differentiator::Union{V1.Gene, Symbol}
    duration::Float64
    timer::V1.Gene
    buffer::Vector{Float64} = Float64[]
    ratio::Float64 = 0.5
    next::Union{Transient, V1.Gene, Symbol}
    alternative::Union{Transient, V1.Gene, Symbol}

    # These values will be picked up when using the (optional) bootstrap
    # mechanism via the "{bootstrap/differentiation}" instant Model:
    timer_deposit::Dict{Symbol, Int} = Dict{Symbol, Int}()
    buffer_deposit::Int = 0

    # The following defaults can be overridden (e.g. in a specification), but
    # that may break the timing and ratio of differentiation, the calculations
    # of which currently assume these fixed values:
    timer_trigger_at::Float64 = 0.5
    timer_brake_at::Float64 = 0.5
    timer_repression::Float64 = 2.0
    timer_proteolysis::Float64 = 0.000002
    differentiator_self_activation::Float64 = 2.0
    differentiator_proteolysis::Float64 = 0.0001
end

"""
    Definition

Defines how to construct a [`SciML.JumpModel`](@ref) (via [`V1`](@ref)) from a
specific differentiation tree and peripheral network.

A `Definition` contains instructions for [`build`](@ref) to assemble a
(lower-level) `V1.Definition` that exhibits differentiating behavior; it
consists of definitions for the `differentiation` tree, the `peripheral` network
to be regulated downstream, and the `deposit` of chemical species required
initially to boostrap the system.

The `differentiation` is a nested structure of [`Transient`](@ref)s that each
represent an inner node of the tree and correspond to transient cell states
during differentiation.

# Specification

In JSON, a `Differentiation.Definition` is specificed as a JSON object
specifying a [`V1.Definition`](@ref) (as described there) that is taken as the
downstream peripheral network and further contains an additional mapping
`"differentiation": <differentiation>` that defines the differentiation tree.
Here, `<differentiation>` must specify a [`Transient`](@ref). (The `deposit`
definitions are collected automatically from that differentiation.)
```
"""
@kwdef struct Definition
    differentiation::Transient
    peripheral::V1.Definition
    deposit::Dict{Symbol, Int} = Dict{Symbol, Int}()
    meta::Dict{Symbol, Any} = Dict{Symbol, Any}()
end

representation(x::Transient) = representation(
    x,
    simple = true,
    omit_defaults = [
        :buffer => []
        :timer_deposit => Dict()
        :buffer_deposit => 0
        :timer_trigger_at => 0.5
        :timer_brake_at => 0.5
        :timer_repression => 2.0
        :differentiator_self_activation => 2.0
        :differentiator_proteolysis => 0.0001
    ]
)
representation(x::Definition) = Dict{Symbol, Any}(
    Symbol("{regulation/differentiation}") => merge(
        only(values(representation(x.peripheral))),
        Dict{Symbol, Any}(:differentiation => representation(x.differentiation))
    )
)

Models.describe(::Definition) =
    Models.Label("'regulation/differentiation' definition")

function cast(
    ::Type{Transient},
    x::AbstractDict{Symbol},
    ::Val{:differentiator};
    context,
)
    d = x[:differentiator]
    if d isa AbstractString
        Symbol(d)
    else
        cast(V1.Gene, merge(Dict(:name => ""), d); context)
    end
end

cast(::Type{Transient}, x::AbstractDict{Symbol}, ::Val{:timer}; context) =
    cast(V1.Gene, merge(Dict(:name => ""), x[:timer]); context)

cast(::Type{Union{Symbol, Transient, V1.Gene}}, x::AbstractString; _...) =
    Symbol(x)
cast(
    ::Type{Union{Symbol, Transient, V1.Gene}},
    x::AbstractDict{Symbol};
    context,
) =
    if haskey(x, :next)
        cast(Transient, x; context)
    else
        cast(V1.Gene, merge(Dict(:name => ""), x); context)
    end

cast(::Type{Definition}, x::AbstractDict{Symbol}; context = x) = Definition(
    differentiation = cast(Transient, x[:differentiation]; context),
    peripheral = cast(V1.Definition, x; context),
)

function make_timer!(
    gene::V1.Gene;
    default_name,
    duration,
    buffer_rates,
    genes,
    reactions,
)
    name = gene.name == Symbol() ? Symbol(default_name) : gene.name
    timing = timing_factor(duration)
    timer = V1.Gene(;
        name,
        base_rates = typeof(gene.base_rates)(;
            (
                name => getproperty(gene.base_rates, name)
                for name in propertynames(gene.base_rates)
            )...,
            transcription = timing * gene.base_rates.transcription,
            protein_decay = timing * gene.base_rates.protein_decay,
        ),
        activation = deepcopy(gene.activation),
        repression = deepcopy(gene.repression),
        proteolysis = deepcopy(gene.proteolysis),
    )

    # If specified, add a dimerization buffer reaction for the timer gene's
    # proteins.
    if !isempty(buffer_rates)
        if length(buffer_rates) == 1
            k⁺ = k⁻ = only(buffer_rates)
        elseif length(buffer_rates) == 2
            k⁺, k⁻ = buffer_rates
        else
            error("timer dimerization buffer reaction has too many rates")
        end
        from = Models.Reagents(Dict(timer.name => 2))
        to = Models.Reagents(Dict(Symbol("$(timer.name)_buffer") => 1))
        push!(reactions, Models.Reaction(;
            name = Symbol("$(timer.name)_buffer"), from, to, k⁺, k⁻,
        ))
    end

    genes[timer.name] = timer

    timer
end

function obtain_differentiator!(gene::V1.Gene; default_name, genes)
    # register gene defined inline
    result =
        if gene.name == Symbol()
            V1.Gene(gene; name = Symbol(default_name))
        else
            gene
        end
    genes[result.name] = deepcopy(result)
end

function obtain_differentiator!(reference::Symbol; default_name, genes)
    # replace referenced gene with a copy
    result = pop!(genes, reference)
    genes[result.name] = deepcopy(result)
end

obtain_differentiator!(transient::Transient; default_name, genes) =
    obtain_differentiator!(transient.differentiator; default_name, genes)

descend!(::Any; _...) = nothing
function descend!(
    transient::Transient;
    trigger,
    brake = nothing,
    genes,
    reactions,
    deposit,
)
    # Obtain or create genes for timing and downstream differentiators:
    timer = make_timer!(
        transient.timer;
        default_name = "$(trigger)_timer",
        transient.duration,
        buffer_rates = transient.buffer,
        genes,
        reactions,
    )
    next = obtain_differentiator!(
        transient.next,
        default_name = "$(trigger)0";
        genes,
    )
    alternative = obtain_differentiator!(
        transient.alternative,
        default_name = "$(trigger)1";
        genes,
    )

    # Add regulation triggering timer decay:
    push!(
        timer.repression.slots,
        V1.HillRegulator(from = trigger, at = transient.timer_trigger_at),
    )

    # The alternative upstream differentiator prevents timer decay:
    brake === nothing || push!(
        timer.activation.slots,
        V1.HillRegulator(from = brake, at = transient.timer_brake_at),
    )

    # Control the switch proportion by setting mutually repressive proteolysis
    # between the differentiators:
    let
        0.0 ≤ transient.ratio ≤ 1.0 ||
            error("differentiation ratio is not in the range [0, 1]")
        k1, k2 = proportion_rate.((transient.ratio, 1.0 - transient.ratio))
        push!(
            next.proteolysis.slots,
            V1.DirectRegulator(from = alternative.name, k = k1),
        )
        push!(
            alternative.proteolysis.slots,
            V1.DirectRegulator(from = next.name, k = k2),
        )
    end

    # Add self-activation for the differentiators:
    let at = transient.differentiator_self_activation
        push!(
            next.activation.slots,
            V1.HillRegulator(from = next.name; at),
        )
        push!(
            alternative.activation.slots,
            V1.HillRegulator(from = alternative.name; at),
        )
    end

    # Proteolytic repression from an undecayed timer prevents differentiation:
    let k = transient.differentiator_proteolysis, from = timer.name
        push!(next.proteolysis.slots, V1.DirectRegulator(; from, k))
        push!(alternative.proteolysis.slots, V1.DirectRegulator(; from, k))
    end

    # Repression from any downstream differentiator keeps the timer depleted
    # in the differentiated state:
    let at = transient.timer_repression, k = transient.timer_proteolysis
        push!(
            timer.repression.slots,
            V1.HillRegulator(from = next.name; at),
        )
        push!(
            timer.repression.slots,
            V1.HillRegulator(from = alternative.name; at),
        )
        push!(
            timer.proteolysis.slots,
            V1.DirectRegulator(from = next.name; k),
        )
        push!(
            timer.proteolysis.slots,
            V1.DirectRegulator(from = alternative.name; k),
        )
    end

    # Extend the initial deposit to let the (optional) separate bootstrap stage
    # ("{boostrap/differentiation}" instant model) know what to initialize:
    for (kind, value) in transient.timer_deposit
        deposit[Symbol("$(timer.name).$kind")] = value
    end
    if !isempty(transient.buffer)
        deposit[Symbol("$(timer.name)_buffer")] = transient.buffer_deposit
    end

    # Recurse for any nested differentiations:
    descend!(
        transient.next,
        trigger = next.name,
        brake = alternative.name;
        genes,
        reactions,
        deposit,
    )
    descend!(
        transient.alternative,
        trigger = alternative.name,
        brake = next.name;
        genes,
        reactions,
        deposit,
    )
end

"""
    build(specification::AbstractDict{Symbol})
    build(definition::Definition; method::Symbol = :default)

Construct a differentiating `SciML.JumpModel` from a [`Definition`](@ref).

When interpreting a JSON specification, this function (in its first form) is
called to construct a concrete regulation model on encountering a
`{"{regulation/differentiation}": {...}}` literal. It will first destructure the
parsed JSON into a `Definition` and then proceed from there.

The result is constructed by first interpreting the peripheral network
specification as a `V1` model and extending that by the differentiation
specification, and then wrapping that up in a [`Models.Wrapped`](@ref) with the
model definition. This will result in the following stack of abstractions:
- [`SciML.JumpModel`](@ref Models.SciML.JumpModel), specified by a
- `Catalyst.ReactionSystem`, specified by a
- [`V1.Definition`](@ref), specified by
- `definition::`[`Differentiation.Definition`](@ref), potentially specified by
- `specification`

Note that typically, using `Differentiation` models requires preparation of the
system state; see [`bootstrap`](@ref).

# Specification

Differentiated models are specified in JSON as
`{"{regulation/differentiation}": <definition>}` where `<definition>` specifies
a [`Definition`](@ref) as described there.

For an example, see `examples/specification/differentiation.schedule.json`.
"""
function build end

build(specification::AbstractDict{Symbol}) = build(
    cast(Definition, specification),
    method = Symbol(get(specification, :method, "default"))
)

function build(definition::Definition; method::Symbol)
    # Shallow-copy genes, reactions and deposit:
    genes = Dict(gene.name => gene for gene in definition.peripheral.genes)
    reactions = copy(definition.peripheral.reactions)
    deposit = copy(definition.deposit)

    # Extend them according to the differentiation definition:
    root = definition.differentiation.differentiator
    trigger =
        if root isa V1.Gene
            obtain_differentiator!(
                root,
                default_name = "differentiator";
                genes,
            ).name
        else
            # Exclusively for the root differentiator, we allow it to to alias
            # a non-gene species. (All consequent differentiators must be genes
            # because they will be transcriptionally regulated.)
            root
        end
    descend!(definition.differentiation; trigger, genes, reactions, deposit)

    # Compile down to a V1 model:
    model = V1.build(
        V1.Definition(
            genes = collect(values(genes));
            definition.peripheral.polymerases,
            definition.peripheral.ribosomes,
            definition.peripheral.proteasomes,
            reactions,
        );
        method,
    )

    # Replace the original Differentiation.Definition by an extended variant
    # defining which species to deposit for the newly created timers when we
    # later bootstrap their states (if requested):
    definition = Definition(;
        definition.differentiation,
        definition.peripheral,
        deposit,
        definition.meta,
    )

    Models.Wrapped(; definition, model)
end

constructor(::Val{Symbol("regulation/differentiation")}) = build

"""
    bootstrap(model)

Construct an `Instant` model that deposits species required by the
differentiation `model` before regulation can start.

# Specification

A bootstrap model is specified in JSON by referring to the differentiation
model that should be prepared for regulation. This typically means that they are
both specified as part of a `Schedule` where the differentiation is defined in
an outer scope (e.g. bound to `"do"`) and the bootstrap model is then specified
as `{"{bootstrap/differentiation}": {"\$": "do"}}`.

For an example, see `examples/specification/differentiation.schedule.json`.
"""
function bootstrap end
bootstrap(model::Models.Wrapped) = bootstrap(model.definition, model.model)
bootstrap(d::Definition, ::Models.Model) = Plumbing.setter(d.deposit)
bootstrap(::Any, model::Models.Model) = bootstrap(model)

constructor(::Val{Symbol("bootstrap/differentiation")}) = bootstrap

end
