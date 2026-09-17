module GRSApp

using PrecompileTools

using CytoscapeJS
using GeneRegulatorySystems
using GeneRegulatorySystems: Models, SPECIFICATION_EXAMPLES
using Pluto

const Vis = GeneRegulatorySystems.Visualisation

include("trajectories/trajectories.jl")
include("pagestate.jl")
include("JSONEditor.jl")

export Trajectories, PageState, JSONEditor

function serve(; port = nothing, docs_port = 8001)
    app = pkgdir(GRSApp)
    root = normpath(joinpath(app, ".."))
    docs = joinpath(root, "docs")
    build = joinpath(docs, "build")

    code = isdir(build) ?
        "using LiveServer; serve(dir=$(repr(build)), port=$(docs_port))" :
        "using LiveServer; cd($(repr(root))) do; servedocs(port=$(docs_port)); end"

    docs_process = run(
        Cmd(`$(Base.julia_cmd()) --project=$(docs) -e $(code)`; detach = true);
        wait = false,
    )

    try
        notebook = joinpath(app, "gene_regulatory_systems.jl")
        isnothing(port) ?
            Pluto.run(; notebook) :
            Pluto.run(; notebook, port)
    finally
        process_running(docs_process) && kill(docs_process)
    end
end

function (@main)(args)
    serve()
    return 0
end

@compile_workload begin
    for filename in ("minimal.schedule.json", "differentiation.schedule.json")
        schedule! = Models.load("$SPECIFICATION_EXAMPLES/$filename")

        network = Vis.Network(schedule!)
        groups = string.(network.groups)
        colors = Vis.group_colors(network.groups).colors
        CytoscapeJS.Cytoscape(network; group_colors = colors)

        sink = Trajectories.Sink()
        schedule!(; trace = sink)

        trace = Trajectories.catenate(sink)
        levels = Trajectories.lod(trace)
        for path in ("", first(Vis.paths(network)))
            Trajectories.select(levels, path)
        end

        snapshot = Trajectories.snapshots(
            Trajectories.select(trace, ""), groups, "proteins")
        for components in (2, 3), coloring in (:genes, :time)
            Trajectories.project(
                snapshot; components, group_colors = colors, coloring)
        end
    end
end

end
