module Common

using GeneRegulatorySystems

import JSON

import LibGit2

repository_version() = LibGit2.format(
    LibGit2.GitDescribeResult(
        LibGit2.GitRepo(@__DIR__() |> dirname)
    );
    options = LibGit2.DescribeFormatOptions(;
        dirty_suffix = Base.unsafe_convert(
            Cstring,
            Base.cconvert(Cstring, "-dirty")
        )
    )
)

artifact(kind::Symbol, name = nothing; prefix = "") =
    "$prefix$(artifact(Val(kind), name))"
artifact(::Val{:specification}, ::Nothing) = "experiment.schedule.json"
artifact(::Val{:index}, ::Nothing) = "index.arrow"
artifact(::Val{:events}, into) = "events$into.stream.arrow"

function warn_incompatible_versions(location)
    specification = JSON.parsefile(artifact(:specification; prefix = location))
    specification isa Dict || return
    haskey(specification, "_version") || return

    warned = false
    if (specification["_julia_version"] != "v$VERSION")
        @warn(
            "Experiment was prepared with a different Julia version.",
            specification["_julia_version"],
            "v$VERSION",
        )
        warned = true
    end
    if (specification["_version"] != repository_version())
        @warn(
            "Experiment was prepared with a different " *
            "GeneRegulatorySystems.jl version.",
            specification["_version"],
            repository_version(),
        )
        warned = true
    end
    warned && @info "Reified objects may be incorrect."

    nothing
end

reify(path; location, seed) =
    reify(artifact(:specification, prefix = location), path; seed)
reify(file, path; seed) = Scheduling.reify(
    Models.load(file; seed),
    path;
    load = filename -> JSON.parsefile(
        joinpath(dirname(file), filename),
        dicttype = Dict{Symbol, Any},
    )
)

end
