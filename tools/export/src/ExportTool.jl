module ExportTool

include("$(@__DIR__)/../../common.jl")
using .Common: artifact, warn_incompatible_versions

import Arrow
using Chain
import CSV
using DataFrames
import JSON
import Muon
using PrecompileTools

ensure_binary(sink; format) =
    if sink isa Base.TTY
        @error "Refusing to write `$format` binary data to terminal."
        @info "Either provide `sink`, redirect output, or use another format."
        throw(:help)
    end

ensure_directory(sink; format) =
    if !(sink isa AbstractString)
        @error "Cannot write `$format` data: `$sink` is not a directory."
        throw(:help)
    elseif !isdirpath(sink)
        @error "Cannot write `$format` data: `$sink` does not name a directory."
        @info "Perhaps you meant `$sink/`?"
        throw(:help)
    end

branch(path; pattern = r"^(?:.*/\d+)?") = match(pattern, path).match

function matchfirst(name; kinds)
    for kind in kinds
        m = match(kind, string(name))
        m !== nothing && return m
    end
    return nothing
end

renamed(name, match) =
    isempty(match.captures) ? String(name) : first(match.captures)

function getlayer(match)
    d = Dict(match)
    isempty(d) && return ""  # no match groups: "" is the base layer
    haskey(d, 1) && return ""  # unnamed match group: also base layer
    return first(only(d))  # single match group's name determines layer
end

function load_wide(index; source, kinds)
    is = unique(index.i)
    channels = unique(index.into)
    @chain begin
        mapreduce(vcat, channels) do channel
            @chain source begin
                dirname
                joinpath(channel)
                Arrow.Table
                DataFrame
                subset(:i => ByRow(in(is)))
                transform(:name => (ns -> matchfirst.(ns; kinds)) => :match)
                subset(:match => ByRow(!isnothing))
                transform(
                    [:name, :match] => ByRow(renamed) => :dimension,
                    :match => ByRow(getlayer) => :layer,
                )
            end
        end
        innerjoin(select(index, :i, :label), on = :i)
        unstack([:layer, :dimension], :label, :value, combine = last, fill = 0)
        sort([:layer, :dimension])
    end
end

fuse(dimension, layer) = isempty(layer) ? dimension : "$dimension.$layer"

function main(;
    source,
    sink,
    format,
    path,
    model,
    labeled,
    relabel_with_segment_id,
    relabel_with_timestamp,
    events,
    kinds,
    dry
)
    warn_incompatible_versions(source)

    path = [Regex("^\\Q$prefix\\E(?=[/+-]|\$)") for prefix in path]
    model = [Regex("^\\Q$prefix\\E(?=[/+-]|\$)") for prefix in model]
    isempty(labeled) && push!(labeled, r".")  # Only labeled segments by default

    criteria = Dict{Symbol, Base.Callable}(
        column => ByRow(x -> any(contains(x, r) for r in rs))
        for (column, rs) in [
            (:path, path)
            (:model, model)
            (:label, labeled)
            (:into, events)
        ]
        if !isempty(rs)
    )

    isempty(kinds) && push!(kinds, r"")  # Include all by default

    index = @chain begin
        artifact(:index; prefix = source)
        Arrow.Table
        DataFrame
        subset(:count => ByRow(>(0)))
        subset(criteria...)
    end

    if relabel_with_segment_id
        index.label .= index.label .* "_" .* repr.(index.i)
    end

    by_label = groupby(index, :label)
    for group in by_label
        branches = unique(branch.(group.path))
        if length(branches) > 1
            @error(
                "Multiple branches share the same label.",
                first(group.label),
                branches,
            )
            @info "This is disallowed because it easy to get wrong as only the \
                final segment is retained for each label to create any output \
                observation. You can use the `--relabel-with-segment-id` \
                option to make the labels unique, but be aware that this might \
                create a large number of output observations."
            throw(:help)
        end
    end
    index = combine(by_label, last)

    if relabel_with_timestamp
        index.label .= index.label .* "_" .* repr.(index.to)
    end

    if dry
        @show select(
            index,
            :i,
            :path,
            :to => :at,
            :model,
            :label,
            :count,
            :into => :events,
        )
        return nothing
    end

    result = load_wide(index; source, kinds)

    sink = @something(sink, stdout)
    if format === nothing
        if sink isa AbstractString
            chopped, extension = splitext(sink)
            if isempty(extension) && isdirpath(chopped)
                chopped, extension = splitext(chop(chopped))
            end
            if startswith(extension, '.')
                format = Symbol(extension[2:end])
            end
        else
            format = :csv
        end
    end

    if format == :h5ad  # (layered format)
        layers = groupby(result, :layer)
        for layer in layers
            sort!(layer, :dimension)
        end
        if !allequal(l -> l.dimension, layers)
            @error("Extracted AnnData layers are not homogeneous.")
            @info "This means that the layers describe distinct sets of genes."
            @info "Thus they cannot be stored in a single AnnData object."
            throw(:help)
        end

        baselayer = get(layers, ("",), first(layers))
        h5ad = Muon.AnnData(
            X = Array(select(baselayer, Not(:dimension, :layer)))'
        )
        h5ad.obs_names .= names(baselayer, Not(:dimension, :layer))
        h5ad.var_names .= String.(baselayer.dimension)

        for (groupkey, layer) in pairs(layers)
            layername = first(groupkey)
            layername == "" && continue
            h5ad.layers[layername] =
                Array(select(layer, Not(:dimension, :layer)))'
        end

        if !haskey(layers, ("",))
            # Base layer had no data so a provisional arbitrary layer was used.
            # Now replace its values.
            h5ad.X = sum(values(h5ad.layers))
        end

        if sink isa AbstractString
            Muon.writeh5ad(sink, h5ad)
        else
            ensure_binary(sink; format)
            mktemp() do temporary, _io
                Muon.writeh5ad(temporary, h5ad)
                write(sink, read(temporary))
            end
        end
    else  # (flat formats)
        select!(
            result,
            [:dimension, :layer] => ByRow(fuse) => :name,
            Not(:dimension, :layer),
        )
        segments_as_columns = result
        segments_as_rows = permutedims(result, :name)

        if format == :arrow
            ensure_binary(sink; format)
            Arrow.write(sink, segments_as_rows)
        elseif format == :csv
            CSV.write(sink, segments_as_rows)
        elseif format == :tsv
            CSV.write(sink, segments_as_rows, delim = '\t')
        elseif format == :bonsai
            ensure_directory(sink; format)
            mkpath(sink)
            CSV.write(
                "$(sink)features.txt",
                select(segments_as_columns, Not(:name)),
                delim = '\t',
                writeheader = false,
            )
            CSV.write(
                "$(sink)geneID.txt",
                select(segments_as_columns, :name),
                writeheader = false,
            )
            CSV.write(
                "$(sink)cellID.txt",
                select(segments_as_rows, :name),
                writeheader = false,
            )
        else
            @error "Cannot export to unknown format `$format`."
            throw(:help)
        end
    end

    nothing
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
                label = ["A", "B", "C"],
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

        write(artifact(:specification, prefix = location), "[]")
        defaults = (;
            source = location,
            format = nothing,
            path = String[],
            model = String[],
            labeled = Regex[],
            relabel_with_segment_id = false,
            relabel_with_timestamp = false,
            events = Regex[],
            kinds = Regex[],
            dry = false,
        )

        @compile_workload begin
            main(; defaults..., sink = "$location/export.csv")
            main(; defaults..., sink = "$location/export.arrow")
            main(; defaults..., sink = "$location/export.h5ad")
            main(;
                defaults...,
                sink = "$location/export.anndata",
                format = :h5ad,
                kinds = [r"(\d+)\.proteins", r"(?<transcripts>\d+)\.mrnas"],
            )
        end
    end
end

end
