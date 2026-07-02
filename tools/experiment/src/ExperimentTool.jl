module ExperimentTool

include("$(@__DIR__)/../../common.jl")
using .Common: repository_version, artifact

import Arrow
using Chain
using DataFrames
using GeneRegulatorySystems
import JSON
using LoggingExtras
using PrecompileTools
import ProgressLogging
using TerminalLoggers: TerminalLogger
using UUIDs

import Dates

abstract type ProgressLogger <: AbstractLogger end

Logging.min_enabled_level(::ProgressLogger) = Scheduling.Progress

Logging.shouldlog(::ProgressLogger, level, module_, _group, _id) =
    level == Scheduling.Progress && (
        module_ == Scheduling ||
        module_ == Models.SciML ||
        module_ == ExperimentTool
    )

struct SimpleProgressLogger <: ProgressLogger end

function Logging.handle_message(
    ::SimpleProgressLogger,
    _level,
    message::Symbol,
    _rest...;
    at,
    todo = nothing,
    done = nothing,
)
    if message == :done
        @info "`$at` done"
    elseif message == :saved
        @info "Saved $done events into `$at`"
    end
end

@kwdef struct BarProgressLogger <: ProgressLogger
    ids::Dict{String, UUID} = Dict{String, UUID}()
    todo::Dict{String, Real} = Dict{String, Real}()
end

function Logging.handle_message(
    logger::BarProgressLogger,
    _level,
    message::Symbol,
    _rest...;
    at,
    todo = nothing,
    done = 0,
)
    if message == :saved
        @info "Saved $done events into `$at`"
        return
    end

    id = get!(uuid4, logger.ids, at)

    if message == :done || message == :descending
        @info ProgressLogging.Progress(id, done = true, name = "$at $message")
        delete!(logger.ids, at)
        delete!(logger.todo, at)
        return
    end

    if todo isa Real && isfinite(todo)
        logger.todo[at] = todo
    end

    if todo isa Real && !isfinite(todo) && done <= 0
        # To unclutter the bar stack we will pretend we are done with this path
        # when we freshly step into a Scope. If we make a second step in this
        # scope, we will just spawn a new bar with that log message.
        @info ProgressLogging.Progress(id, done = true, name = "$at ...")
        return
    end

    todo = @something(todo, get(logger.todo, at, nothing), Some(nothing))
    current =
        if message == :iterating
            done + 1
        elseif done isa AbstractFloat
            round(done, digits = 1)
        else
            done
        end
    fraction = nothing
    if todo === nothing
        details = ""
    elseif todo isa Real
        if isfinite(todo)
            fraction = done / todo
            details = "($current/$todo)"
        else
            details = "($current)"
        end
    else
        details = "($todo)"
    end
    action = lpad(message, 10)
    description = join(filter(!isempty, [action, at, details]), ' ')
    @info ProgressLogging.Progress(id, fraction, name = description)
end

function map_artifacts(paths)
    preliminary_map = Dict(p => basename(p) for p in paths)
    if allunique(vcat(artifact(:specification), values(preliminary_map)...))
        return preliminary_map
    else
        return Dict(p => "_$(i)_$(basename(p))" for (i, p) in enumerate(paths))
    end
end

function assert_compatible_versions(specification)
    if (specification[:_julia_version] != "v$VERSION")
        @error(
            "Experiment was prepared with a different Julia version.",
            specification.bindings[:_julia_version],
            "v$VERSION",
        )
        @info "This is disallowed to ensure reproducibility."
        throw(:help)
    end

    if (specification[:_version] != repository_version())
        @error(
            "Experiment was prepared with a different " *
            "GeneRegulatorySystems.jl version.",
            specification[:_version],
            repository_version(),
        )
        @info "This is disallowed to ensure reproducibility."
        throw(:help)
    end
end

function prepare!(; location, specifications, seed)
    specification_path = artifact(:specification; prefix = location)
    if ispath(specification_path)
        @error(
            "Cannot prepare experiment, specification already exists.",
            specification_path,
        )
        throw(:help)
    elseif isempty(specifications)
        @error "No specifications given to be prepared."
        throw(:help)
    else
        paths = realpath.(specifications)
        artifacts = map_artifacts(paths)
        wrapped = map(paths) do p
            name = replace(artifacts[p], r"(\.schedule)?(\.json)$" => "")
            filename = basename("$(location)$(artifacts[p])")
            Dict(
                :< => filename,
                :seed => "\${rootseed}-$name",
                :into => "-$name",
                :flush => true,
            )
        end

        mkpath(dirname(specification_path))
        for (source, target) in artifacts
            cp(source, "$(location)$target"; follow_symlinks = true)
        end
        open(specification_path, "w") do file
            JSON.print(
                file,
                merge(
                    Dict(
                        :rootseed => seed,
                        :_version => repository_version(),
                        :_julia_version => "v$VERSION",
                    ),
                    if length(wrapped) == 1
                        s = only(wrapped)
                        Dict(:step => Dict(:< => s[:<], :seed => s[:seed]))
                    else
                        Dict(:step => wrapped, :branch => true)
                    end,
                ),
                4,
            )
        end
    end
end

@kwdef struct Channel
    is::Vector{Int64} = Int64[]
    ts::Vector{Float64} = Float64[]
    names::Vector{Symbol} = Symbol[]
    values::Vector{Int64} = Int64[]
end

@kwdef mutable struct Sink
    location::String
    i::Int = 0
    index = []
    threshold::Int = 1000000
    channels::Dict{String, Channel} = Dict{String, Channel}()
end

function flush!(sink::Sink; matching = nothing, finalize = false)
    sink.i > 0 || return

    for into in keys(sink.channels)
        if matching === nothing || startswith(into, matching)
            flush!(sink, into)
        end
    end

    finalize || return
    index = Tables.columntable(sink.index)
    Arrow.write(
        artifact(:index; prefix = sink.location),
        (;
            index.i,
            index.path,
            index.from,
            index.to,
            model = Arrow.DictEncode(index.model),
            label = Arrow.DictEncode(index.label),
            index.count,
            into = Arrow.DictEncode(index.into),
        )
    )
end

function flush!(sink::Sink, into::AbstractString)
    channel = pop!(sink.channels, into)
    filename = artifact(:events, into, prefix = sink.location)
    events = (;
        i = channel.is,
        t = channel.ts,
        name = channel.names,
        value = channel.values,
    )
    count = length(events.t)
    if isfile(filename)
        Arrow.append(filename, events)
    else
        Arrow.write(filename, events, file = false)
    end
    @logmsg(Scheduling.Progress, :saved, at = filename, done = count)
end

function (sink::Sink)(; path, into = nothing, flush = nothing, _...)
    @logmsg Scheduling.Progress :flushing at = path
    matching = flush === true ? into : flush
    matching !== nothing && flush!(sink; matching)
end

function (sink::Sink)(
    state;
    path,
    from,
    primitive!,
    into = nothing,
    consolidated_progress = nothing,
    _...
)
    sink.i += 1

    to = Models.t(state)
    model = primitive!.path
    label = get(primitive!.bindings, :label, "")

    if into === nothing
        push!(
            sink.index,
            (; sink.i, path, from, to, model, label, count = 0, into = "")
        )
        return
    end

    if consolidated_progress === nothing
        @logmsg Scheduling.Progress :collecting at = path todo = "into $into"
    else
        consolidated_progress(:collecting, done = to)
    end
    channel = get!(Channel, sink.channels, into)
    count = 0
    Models.each_event(state) do t::Float64, name::Symbol, value::Int64
        if length(channel.values) ≥ sink.threshold
            flush!(sink, into)
            channel = sink.channels[into] = Channel()
        end
        push!(channel.is, sink.i)
        push!(channel.ts, t)
        push!(channel.names, name)
        push!(channel.values, value)
        count += 1
    end

    filename = basename(artifact(:events, into, prefix = sink.location))
    push!(
        sink.index,
        (; sink.i, path, from, to, model, label, count, into = filename),
    )
end

function with_progress(run!, progress::Symbol)
    if progress == :bars && stdout isa Base.TTY
        with_logger(run!, TeeLogger(current_logger(), BarProgressLogger()))
    elseif progress == :simple
        with_logger(run!, TeeLogger(current_logger(), SimpleProgressLogger()))
    else
        run!()
    end
end

function simulate!(; location, progress, dry)
    load(path) = JSON.parsefile(
        "$(dirname(location))/$path";
        dicttype = Dict{Symbol, Any}
    )

    function dryrun(primitive!, x, Δt; path, into, _...)
        modelname = nameof(typeof(Models.unwrap(primitive!.f!)))
        maybe_source = primitive!.path == path ? "" : " ($(primitive!.path))"
        isinterval = isfinite(Δt) && 0.0 < Δt
        interval = isinterval ? "from $(x.t) to $(x.t + Δt)" : "at $(x.t)"
        maybe_into =
            if into === nothing
                ""
            elseif isinterval && 0.0 < primitive!.skip
                ", record final into '$into'"
            else
                ", record into '$into'"
            end

        @info "$path: $modelname$maybe_source $interval$maybe_into"
    end

    specification = load(basename(artifact(:specification; prefix = location)))
    assert_compatible_versions(specification)
    schedule! = Model(
        specification,
        bindings = Dict(
            :into => "",
            :channel => "",
            :seed => "",
            :defaults => Models.load_defaults(),
        ),
    )

    if dry
        progress = :none
    end

    with_progress(progress) do
        state = Models.FlatState()
        sink = Sink(; location)
        schedule!(state; load, trace = sink, dryrun = dry ? dryrun : nothing)
        flush!(sink, finalize = true)
    end
end

function main(;
    location,
    prepare,
    simulate,
    progress,
    dry,
    seed,
    specifications,
)
    global_logger(TerminalLogger())

    timestamp = Dates.now()
    location = replace(location, "{TIMESTAMP}" => timestamp)

    if !any([prepare, simulate])
        prepare = true
        simulate = true
    end

    prepare && prepare!(; location, specifications, seed)
    simulate && simulate!(; location, progress, dry)

    nothing
end

@setup_workload begin
    get(ENV, "JULIA_PKG_PRECOMPILE_AUTO", "1") == "0" &&
        error("Precompilation triggered implicitly; this should not happen.")

    specification = "$(@__DIR__)/precompile.schedule.json"
    mktempdir() do location
        @compile_workload begin
            main(
                location = "$location/$(Dates.now())/",
                prepare = true,
                simulate = true,
                progress = :none,
                dry = false,
                seed = "seed",
                specifications = [specification],
            )
        end
    end
end

end
