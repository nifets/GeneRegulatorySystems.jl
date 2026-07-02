module InspectTool

include("$(@__DIR__)/../../common.jl")
include("$(@__DIR__)/visualization.jl")

using .Common: Dimension, artifact, warn_incompatible_versions
using .Visualization: Catenation

import Arrow
using Chain
using DataFrames
using GeneRegulatorySystems
using GLMakie
using PrecompileTools

const LIMITS = (;
    catenations = 500,
    groups = 32,
)

@kwdef struct AdjacentPrefixes
    parent::String
    next::String
    previous::String
    firstborn::String
end

@kwdef struct PreparedData
    index::DataFrame
    events::Union{Dict{Symbol, Dict{Int, Catenation}}, Nothing}
    model::Union{Models.Description, Nothing}
    groups::Union{Vector{String}, Nothing}
    group_colors::Visualization.GroupColors
    adjacents::AdjacentPrefixes
end

@kwdef struct Selection
    channel::String = ""
    items_prefix::String = ""
    label_pattern::String = ""
end

const BRANCH_PATTERN = r"^(?:.*/\d+)?"

branch(path) = match(BRANCH_PATTERN, path).match

function arrange(index)
    result = Int[]
    track = 0
    tracks = Dict{String, Int}()
    watermarks = [0.0]
    for segment in eachrow(index)
        if segment.branch != branch
            track = get!(tracks, segment.branch) do
                if track > 0
                    watermarks[track] = Inf
                end
                @something(
                    findfirst(≤(segment.from), watermarks),
                    lastindex(push!(watermarks, NaN)),
                )
            end
        end
        watermarks[track] = segment.to
        push!(result, track)
    end
    result
end

function load_index(location)
    result = @chain begin
        artifact(:index; prefix = location)
        Arrow.Table
        DataFrame
    end

    result.branch = branch.(result.path)
    result.track = arrange(result)

    result
end

function filter(index; selection)
    template = Dict(
        :channel => channel -> :into => ByRow(contains(channel)),
        :items_prefix => prefix -> (
            :path => ByRow(contains(Regex("^\\Q$prefix\\E(?=[/+-]|\$)")))
        ),
        :label_pattern =>
            label_pattern -> :label => ByRow(contains(Regex(label_pattern))),
    )

    criteria = [
        template[name](getproperty(selection, name))
        for name in propertynames(selection)
        if !isempty(getproperty(selection, name))
    ]

    subset(index, :count => ByRow(>(0)), criteria...)
end

function cut(filtered)
    result = Catenation[]

    front = 0
    back = 0
    branch = nothing
    for segment in eachrow(filtered)
        segment.count > 0 || continue  # Ignore segments without recorded data.

        # Terminate and emit an ongoing run when encountering a non-instant
        # segment or on switching branches:
        if front > 0 && (segment.branch != branch || segment.from < segment.to)
            push!(result, Catenation(; front, back))
            front = back = 0
        end
        branch = segment.branch

        current = rownumber(segment)
        if segment.from < segment.to
            # Emit a single-segment run:
            push!(result, Catenation(front = current, back = current))
        else
            # Initiate or extend an ongoing run with the current segment:
            back = current
            if front == 0
                front = current
            end
        end
    end

    # Emit the potentially ongoing run:
    front > 0 && push!(result, Catenation(; front, back))

    result
end

function place!(catenation::Catenation, event)
    dimension = Dimension(event.name)
    series = get!(catenation.series, dimension) do
        Visualization.seriestype(dimension)()
    end

    push!(series.ts, event.t)
    push!(series.ys, event.value)

    nothing
end

augment!(catenation::Catenation) =
    for (dimension, active) in catenation.series
        dimension.kind == :active || continue
        inactive_dimension = Dimension(:inactive, dimension.group)
        catenation.series[Dimension(:activity, dimension.group)] =
            if haskey(catenation.series, inactive_dimension)
                Visualization.FractionSeries(
                    active,
                    catenation.series[inactive_dimension],
                )
            else
                Visualization.FractionSeries(
                    ys = Float64.(active.ys);
                    active.ts,
                )
            end
    end

function load_events(filtered; location)
    # Since Makie does not handle large numbers of plot objects well, we
    # optimize the common case where the simulation state was sampled only at
    # equidistant time steps by concatenating these consecutive slices into
    # longer time series that can later be plotted (in a different style) as
    # single Makie plot objects.

    # We first bracket the segments into these potentially composite ranges,
    # each with an initially empty collection of time series attached.
    catenations = cut(filtered)
    length(catenations) ≤ LIMITS.catenations || return nothing
    catenations_index = Dict(
        filtered[k, :i] => j
        for j in LinearIndices(catenations)
        for k in (catenations[j].front : catenations[j].back)
        if filtered[k, :count] > 0
    )

    # Next we load all linked event streams and sort the events into their
    # respective timeseries.
    for channel in unique(subset(filtered, :count => ByRow(>(0))).into)
        events = @chain location begin
            dirname
            joinpath(channel)
            Arrow.Table
            DataFrame
        end
        for event in eachrow(events)
            haskey(catenations_index, event.i) || continue
            catenation = catenations[catenations_index[event.i]]
            place!(catenation, event)
        end
    end

    # Finally we regroup the timeseries (i.e. break up the catenations) by
    # kind, each catenation indexed by its final segment so that we can later
    # look them up from backlinks to connect them.
    result = Dict{Symbol, Dict{Int, Catenation}}()
    for catenation in catenations
        augment!(catenation)

        by_kind = Dict{Symbol, Catenation}()
        for (dimension, series) in catenation.series
            catenation′ = get!(by_kind, dimension.kind) do
                Catenation(; catenation.front, catenation.back)
            end
            catenation′.series[dimension] = series
        end

        for (kind, catenation′) in by_kind
            catenations′ = get!(Dict{Int, Catenation}, result, kind)
            catenations′[catenation′.back] = catenation′
        end
    end

    result
end

function backlinks(index)
    result = Int[]
    tips = [0]
    for segment in eachrow(index)
        local tip, branch
        while true
            tip = last(tips)
            branch = iszero(tip) ? "" : index[tip, :branch]
            startswith(segment.branch, branch) && break
            pop!(tips)
        end
        push!(result, tip)
        segment.branch == branch || push!(tips, tip)
        tips[end] = rownumber(segment)
    end
    result
end

function prepare(index; selection = Selection(), location)
    parent, current = let
        m = match(r"(.*)(\+|[/-]\d+)$", selection.items_prefix)
        m === nothing ? ["", ""] : m
    end

    index = filter(
        index,
        selection = Selection(
            items_prefix = parent;
            selection.channel,
            selection.label_pattern,
        ),
    )

    segment_pattern = r"^(\+|[/-]?\d+).*"
    siblings = unique(
        match(segment_pattern, chopprefix(path, parent))[1]
        for path in index.path
    )
    previous, next = let
        n = length(siblings)
        i = findfirst(==(current), siblings)
        if n < 2 || i === nothing
            "", ""
        else
            (
                parent * siblings[mod1(i - 1, n)],
                parent * siblings[mod1(i + 1, n)],
            )
        end
    end

    index = filter(index, selection = Selection(; selection.items_prefix))
    index.previous = backlinks(index)

    firstborn =
        if isempty(index)
            ""
        else
            prefix = selection.items_prefix
            first_tail = chopprefix(first(index.path), prefix)
            m = match(segment_pattern, first_tail)
            m === nothing ? "" : prefix * m[1]
        end

    adjacents = AdjacentPrefixes(; previous, next, firstborn, parent)

    events = load_events(index; location)
    if events === nothing
        groups = nothing
    else
        groups = unique(
            dimension.group
            for catenations in values(events)
            for catenation in values(catenations)
            for dimension in keys(catenation.series)
            if !isempty(dimension.group)
        )
        if length(groups) > 32
            groups = nothing
        else
            sort!(groups)
        end
    end

    group_colors = Visualization.GroupColors(groups)

    model_locators = unique(subset(index, [:from, :to] => ByRow(<)).model)
    if isempty(model_locators)
        model_locators = unique(index.model)
    end
    model =
        if length(model_locators) == 1
            Models.describe(
                # The seed will never be accessed because it will be overridden
                # in the loaded experiment specification.
                Common.reify(only(model_locators), seed = nothing; location)
            )
        else
            nothing
        end

    PreparedData(; index, events, groups, group_colors, model, adjacents)
end

function attach_display!(figure, ::Val{:selector}; data, selection, _...)
    items_prefix =
        isempty(selection[].items_prefix) ? " " : selection[].items_prefix
    label_pattern =
        isempty(selection[].label_pattern) ? " " : selection[].label_pattern

    navigation = Dict{Symbol, Button}()
    for (; name, label) in (
        (name = :parent, label = "↰"),
        (name = :previous, label = "▲"),
        (name = :next, label = "▼"),
        (name = :firstborn, label = "↳"),
    )
        value = getproperty(data.adjacents, name)
        is_active = !isempty(value) || name == :parent && items_prefix != " "
        navigation[name] = Button(
            figure,
            buttoncolor =
                is_active ? RGBf(0.92, 0.92, 0.92) : RGBf(0.97, 0.97, 0.97),
            labelcolor =
                is_active ? RGBf(0.0, 0.0, 0.0) : RGBf(0.7, 0.7, 0.7),
            tellheight = false;
            label,
        )
        is_active && on(navigation[name].clicks) do n
            selection[] = Selection(;
                selection[].channel,
                items_prefix = value,
                selection[].label_pattern,
            )
        end
    end

    items_prefix_textbox = Textbox(
        figure,
        stored_string = items_prefix,
        placeholder = "⟨prefix⟩",
    )
    on(items_prefix_textbox.stored_string) do items_prefix
        selection[] = Selection(;
            selection[].channel,
            items_prefix = strip(items_prefix),
            selection[].label_pattern,
        )
    end

    label_pattern_textbox = Textbox(
        figure,
        stored_string = label_pattern,
        placeholder = "⟨pattern⟩",
    )
    on(label_pattern_textbox.stored_string) do label_pattern
        selection[] = Selection(;
            selection[].channel,
            selection[].items_prefix,
            label_pattern = strip(label_pattern),
        )
    end

    navigation[:previous].height[] = Relative(1)
    navigation[:next].height[] = Relative(1)
    stepper = GridLayout(default_rowgap = 0, tellheight = false)
    stepper[1:2, 1] = [navigation[:previous], navigation[:next]]

    widgets = [
        Label(figure, "prefix:")
        navigation[:parent]
        items_prefix_textbox
        stepper
        navigation[:firstborn]
        Label(figure, "label:")
        label_pattern_textbox
    ]

    selector = GridLayout(default_colgap = 12, tellwidth = false)
    selector[1, 1:length(widgets)] = widgets
    colgap!(selector, 3, 8)
    colgap!(selector, 5, 32)

    selector
end

function attach_display!(figure, ::Val{:trajectory}; data, kinds, _...)
    if data.events === nothing
        Label(
            figure,
            "(>$(LIMITS.catenations) catenations)",
            tellheight = false,
        )
    elseif isempty(data.events)
        Label(figure, "(no data)", tellheight = false)
    elseif isempty(kinds)
        Label(figure, "(no kinds selected for display)", tellheight = false)
    else
        Visualization.attach_trajectory!(
            figure;
            data.index,
            data.events,
            kinds,
            data.group_colors,
        )
    end
end

function attach_display!(figure, ::Val{:model}; data, _...)
    if data.model === nothing
        Label(
            figure,
            "(model not unique among non-instant segments)",
            tellheight = false,
        )
    else
        Visualization.attach_model!(figure, data.model; data.group_colors)
    end
end

function attach_display!(figure, ::Val{:legend}; data, _...)
    if data.groups === nothing
        Label(figure, "(>32 groups)", tellwidth = false)
    else
        groups = sort(data.groups)
        Legend(
            figure,
            map(groups) do group
                MarkerElement(
                    marker = :circle,
                    markersize = 32,
                    color = data.group_colors[group],
                    strokewidth = 2,
                    strokecolor = :black,
                )
            end,
            groups,
            orientation = :horizontal,
            tellwidth = false,
            tellheight = true,
        )
    end
end

function attach_display!(figure, ::Val{:info}; data, _...)
    events_count = sum(data.index.count)
    message =
        if data.events === nothing
            "$events_count events \
                (>$(LIMITS.catenations) catenations, not loaded)"
        else
            groups_count =
                if data.groups === nothing
                    ">$(LIMITS.groups)"
                else
                    "$(length(data.groups))"
                end
            "$events_count events of $(length(data.events)) kinds \
                in $groups_count groups × $(nrow(data.index)) segments"
        end
    Label(figure, message, tellwidth = false)
end

function build_figure(;
    data::PreparedData,
    displays::AbstractSet{Symbol},
    kinds::AbstractVector{Symbol},
    selection::Observable{Selection},
    size::Tuple{Float64, Float64},
)
    figure = Figure(; size)
    subplots = Dict(
        name => attach_display!(
            figure,
            Val(name);
            data,
            kinds,
            data.group_colors,
            selection,
        )
        for name in displays
    )

    main_plots = collect(skipmissing((
        get(subplots, :trajectory, missing),
        get(subplots, :model, missing),
    )))
    main_grid = GridLayout(tellheight = false)
    main_grid[1, 1:length(main_plots)] = main_plots

    root_plots = collect(skipmissing((
        get(subplots, :selector, missing),
        main_grid,
        get(subplots, :legend, missing),
        get(subplots, :info, missing),
    )))
    figure[1:length(root_plots), 1] = root_plots

    figure
end

function main(;
    location,
    channel,
    items_prefix,
    label_pattern,
    displays,
    kinds,
    size,
    ppi,
    out,
    wait_for_close,
)
    warn_incompatible_versions(location)

    selection = Observable(
        Selection(
            channel = something(channel, ""),
            items_prefix = something(items_prefix, ""),
            label_pattern = something(label_pattern, ""),
        )
    )
    displays = Set(Symbol.(split(displays, ',', keepempty = false)))
    kinds = Symbol.(split(kinds, ',', keepempty = false))
    size = Tuple(parse.(Float64, split(size, 'x')))

    index = load_index(location)

    GLMakie.activate!()
    screen = GLMakie.Screen()

    if out !== nothing
        wait_for_close = false
    end

    on(selection) do selected
        empty!(screen)
        data = prepare(index, selection = selected; location)
        figure = build_figure(; data, displays, kinds, selection, size)
        if out === nothing
            display(screen, figure)
        else
            save(out, figure, px_per_unit = ppi / 96)
        end
    end

    notify(selection)
    if wait_for_close
        wait(screen)
    else
        close(screen)
    end
end

@setup_workload begin
    get(ENV, "JULIA_PKG_PRECOMPILE_AUTO", "1") == "0" &&
        error("Precompilation triggered implicitly; this should not happen.")

    mktempdir() do temporary
        location = "$temporary/"
        Arrow.write(
            artifact(:index, prefix = location),
            DataFrame(
                i = 1:3,
                path = ["-1", "-1", "-2"],
                from = [0.0, 1.0, 1.0],
                to = [0.0, 1.0, 2.0],
                model = ["-1", "-1", "-2"],
                label = ["", "", "SKG"],
                count = [2, 2, 1],
                into = fill(
                    basename(artifact(:events, "-1", prefix = location)),
                    3,
                ),
            )
        )

        Arrow.write(
            artifact(:events, "-1", prefix = location),
            DataFrame(
                i = [1, 1, 2, 2, 3],
                t = [0.0, 0.0, 1.0, 1.0, 1.5],
                name = [
                    Symbol("1.mrnas"),
                    Symbol("1.proteins"),
                    Symbol("1.mrnas"),
                    Symbol("1.proteins"),
                    Symbol("2.active"),
                ],
                value = [1, 0, 2, 1, 1],
            )
        )

        cp(
            "$(@__DIR__)/precompile.schedule.json",
            artifact(:specification, prefix = location),
        )

        @compile_workload begin
            main(
                channel = nothing,
                items_prefix = "-1",
                label_pattern = nothing,
                displays = "selector,trajectory,model,legend,info",
                kinds = "log10(mrnas),proteins",
                size = "1280x720",
                ppi = 96.0,
                out = nothing,
                wait_for_close = false;
                location,
            )
            main(
                channel = nothing,
                items_prefix = nothing,
                label_pattern = "SKG",
                displays = "selector,trajectory,model,legend,info",
                kinds = "activity",
                size = "1280x720",
                ppi = 300.0,
                out = nothing,
                wait_for_close = false;
                location,
            )
        end
    end
end

end
