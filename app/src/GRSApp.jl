module GRSApp

using PrecompileTools

using CytoscapeJS
using GeneRegulatorySystems
using GeneRegulatorySystems: Models, SPECIFICATION_EXAMPLES
using Pluto
import Bonito

const Vis = GeneRegulatorySystems.Visualisation

include("trajectories/trajectories.jl")
include("pagestate.jl")
include("JSONEditor.jl")
include("DarkMode.jl")

export Trajectories, PageState, JSONEditor, DarkMode

function serve(; port = nothing, docs_port = 8001)
    app = pkgdir(GRSApp)
    root = normpath(joinpath(app, ".."))
    build = joinpath(root, "docs", "build")

    isdir(build) || run(
        Cmd(`$(Base.julia_cmd()) --project=docs docs/make.jl`; dir = root);
        wait = false,
    )

    docs = try
        server = Bonito.Server("127.0.0.1", docs_port)
        Bonito.route!(server, r".*" => Bonito.FolderServer(build))
        server
    catch exception
        @warn "Could not serve the docs on port $docs_port." exception
        nothing
    end
    isnothing(docs) || (ENV["GRS_DOCS_URL"] = "http://localhost:$(docs.port)")

    try
        notebook = joinpath(app, "gene_regulatory_systems.jl")
        isnothing(port) ?
            Pluto.run(; notebook) :
            Pluto.run(; notebook, port)
    finally
        isnothing(docs) || close(docs)
    end
end

function (@main)(args)
    serve()
    return 0
end

@compile_workload begin
    for (method, genes) in ((:pca, 10), (:umap, 60)), components in (2, 3)
        Trajectories.project(Val(method), rand(genes, 200), components)
    end

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
        projections = [
            Trajectories.project(snapshot; components, group_colors = colors, coloring)
            for components in (2, 3) for coloring in (:genes, :time, :schedule)
        ]
        samples = Trajectories.counts(Trajectories.select(trace, ""), groups, "proteins")
        instants = Trajectories.instants(sink; path = "", channel = nothing)

        for theme in (Trajectories.WGLMakie.Theme(), Trajectories.WGLMakie.theme_dark())
            selection = Trajectories.WGLMakie.Observable("")
            apps = [
                [
                    Trajectories.embed(theme) do
                        Trajectories.render(
                            Trajectories.select(levels, "");
                            tracks = ["activity", "mrnas", "proteins"],
                            selected_genes = groups,
                            group_colors = colors,
                            aggregate_mode,
                            instants,
                        )
                    end
                    for aggregate_mode in (:raw, :aggregate)
                ];
                [
                    Trajectories.embed(theme, PageState.link(PageState.Shared("path"), selection)) do
                        Trajectories.render(projection; selection)
                    end
                    for projection in projections
                ];
                Trajectories.embed(theme) do
                    Trajectories.render(samples; group_colors = colors)
                end
            ]
            Bonito.force_connection(Bonito.PlutoConnection) do
                Bonito.force_asset_server(Bonito.NoServer) do
                    foreach(app -> show(IOBuffer(), MIME"text/html"(), app), apps)
                end
            end
        end
    end
    Bonito.cleanup_globals()
    Trajectories.WGLMakie.Makie.cleanup_globals()
end

end
