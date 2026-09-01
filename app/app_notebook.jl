### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 40d24905-5708-4580-9ae5-71e1c1de01a9
begin
    import Pkg
    Pkg.activate(".")
    using Base64
    using HypertextLiteral
end

# ╔═╡ 34a8deb2-0db9-456c-ba38-eabc3254ab50
md"""
## Logic
"""

# ╔═╡ 76de62b9-be88-4b96-b801-b40872dfa0be
md"""
### Schedule
"""

# ╔═╡ 16d523e4-ebec-4d6e-b5e1-7aaf949ec5fc
md"""
### Network
"""

# ╔═╡ c7cd57d6-464d-4dbc-b955-699802e60be7
md"""
### Run simulation
"""

# ╔═╡ 35e4e15a-e1f7-42c1-8ded-5714234811c6
md"""
## Styling
"""

# ╔═╡ 4a29ccba-5cdc-47e1-953c-5eee6d1edbe7
@bind dark_mode @htl("""
<span>
<script>
    const root = currentScript.parentElement
    const media = window.matchMedia("(prefers-color-scheme: dark)")

    const update = () => {
        root.value = media.matches
        root.dispatchEvent(new CustomEvent("input"))
    }

    media.addEventListener("change", update)
    update()

    invalidation.then(() =>
        media.removeEventListener("change", update)
    )
</script>
</span>
""")

# ╔═╡ 375349c2-33bb-45fb-b412-672c678b93b3
@htl("""
 <style id="grs-dashboard-styles">
    .grs-dashboard {
        box-sizing: border-box;
        width: 100%;
        overflow: hidden;
        font-family: Montserrat, sans-serif;
    }
    .dashboard-option {
        display: flex;
        align-items: center;
        gap: 0.5rem;
        font-family: Montserrat, sans-serif;
        font-size: 0.8rem;
        font-weight: 300;
    }
    .dashboard-option select,
    .dashboard-option option {
        font-family: Montserrat, sans-serif;
    }
    .dashboard-option select {
        font-size: 0.7rem;
    }
    .schedule-panel {
        display: flex;
        flex-direction: column;
        min-width: 0;
        min-height: 0;
        gap: 0.5rem;
    }
    .schedule-error {
        max-height: 4rem;
        overflow: auto;
        padding: 0.6rem 0.75rem;
        color: #ffb4ab;
        background: rgba(255, 80, 80, 0.12);
        border: 1px solid rgba(255, 100, 100, 0.35);
        border-radius: 6px;
        font-size: 0.7rem;
        white-space: pre-wrap;
    }
    .schedule-error:empty {
        display: none;
    }

    body:has(.grs-dashboard):not(:has(.notebook-toggle input:checked))
        pluto-cell:not(:has(.grs-dashboard)) {
        display: none !important;
    }
    body:not(:has(.notebook-toggle input:checked))
        pluto-cell:has(.grs-dashboard) > pluto-shoulder,
    body:not(:has(.notebook-toggle input:checked))
        pluto-cell:has(.grs-dashboard) > button.add_cell {
        display: none !important;
    }
    body:not(:has(.notebook-toggle input:checked))
        pluto-cell:has(.grs-dashboard) > pluto-trafficlight {
        display: none !important;
    }
    body:not(:has(.notebook-toggle input:checked))
        pluto-editor > main > preamble {
        display: none !important;
    }
    body:not(:has(.notebook-toggle input:checked)) #helpbox-wrapper {
        display: none !important;
    }
    /* hide the dashboard cell's source code outside dev mode. */
    body:not(:has(.notebook-toggle input:checked))
        pluto-cell:has(.grs-dashboard) > pluto-input,
    body:not(:has(.notebook-toggle input:checked))
        pluto-cell:has(.grs-dashboard) > pluto-runarea {
        display: none !important;
    }
    body:not(:has(.notebook-toggle input:checked)) footer {
        display: none !important;
    }
</style>
""")

# ╔═╡ 0187c445-f7a0-4eae-b55d-2bd955df2c2c
begin
    logo_data = "data:image/png;base64," *
        base64encode(read(joinpath(@__DIR__, "assets", "logo.png")))

    splash_screen = @htl("""
<div class="grs-splash">
    <div class="grs-splash-content">
        <div class="grs-splash-title">Gene Regulatory Systems</div>
        <div class="grs-splash-logo">
            <img src=$(logo_data)>
        </div>
        <div class="grs-splash-loading">loading</div>
    </div>
</div>

<style>
    body:not(.grs-ready) pluto-cell {
        visibility: hidden !important;
    }

    .grs-splash {
        visibility: visible;
        position: fixed;
        inset: 0;
        z-index: 99999;
        display: grid;
        place-items: center;
        background: var(--main-bg-color);
    }

    .grs-splash-content {
        display: grid;
        justify-items: center;
        gap: 1.5rem;
        padding: 2.5rem 3rem;
        border-radius: 1rem;
        background: var(--main-bg-color);
    }

    .grs-splash-title {
        font-size: 1.5rem;
        font-family: Montserrat, sans-serif;
    }

    .grs-splash-logo {
        filter: drop-shadow(-4.2px 6px 6px rgba(150, 150, 150, 0.3));
    }

    .grs-splash-logo img {
        display: block;
        width: 120px;
        animation: grs-spin 2s linear infinite;
    }

    .grs-splash-loading {
        font-size: 0.8rem;
        font-family: Montserrat, sans-serif;
    }

    body.grs-ready .grs-splash {
        display: none;
    }

    @keyframes grs-spin {
        0%   { transform: rotate(0); }
        80%  { transform: rotate(2turn); }
        100% { transform: rotate(2turn); }
    }
    .grs-splash-loading::after {
        content: "";
        display: inline-block;
        width: 1.5em;
        text-align: left;
        animation: grs-dots 2s steps(1, end) infinite;
    }
    
    @keyframes grs-dots {
        0%   { content: "."; }
        27%  { content: ".."; }
        54%  { content: "..."; }
        80%  { content: ""; }
        100% { content: ""; }
    }
</style>
""")
end

# ╔═╡ 12e75652-a3fd-11f1-9335-a986ee44237f
begin
    splash_screen

    using Revise
    using GeneRegulatorySystems
    using CytoscapeJS
    using Bonito
    using PlutoUI
    const Layout = PlutoUI.ExperimentalLayout
    import JSON
    import Colors
    const GRS = GeneRegulatorySystems
    const Vis = GRS.Visualisation
    using WGLMakie
    WGLMakie.activate!()
    Revise.includet(@__MODULE__, "src/JSONEditor.jl")
    Revise.includet(@__MODULE__, "src/ValueMultiSelect.jl")
end;

# ╔═╡ a2866674-5e4e-4738-b18e-90d122d2d161
schedule_options = let
	dir = GRS.SPECIFICATION_EXAMPLES
	[
		file => replace(basename(file), ".schedule.json" => "")
		for file in sort(readdir(dir; join=true))
		if endswith(file, ".schedule.json")
	]
end;

# ╔═╡ 7abbbe81-3d0b-4d94-a69f-32c66514687f
schedule_header = @htl("""
<div style="display: flex; align-items: center; gap: 1rem; margin-bottom: 0.75rem;">
<label class="dashboard-option">
    <span>schedule:</span>
    $(@bind selected_schedule PlutoUI.Select(schedule_options))
</label>
<label class="dashboard-option">
    <span>simulate: </span>
    $(@bind run_simulation PlutoUI.Switch(default=false))
</label>
</div>
""")

# ╔═╡ 8392f056-369c-45fb-9c00-719172e82616
schedule_editor = @bind schedule_json JSONEditor(read(selected_schedule, String); height="800px")

# ╔═╡ 69993fb8-daca-4b61-a711-907009dbee24
begin
    gene_colours(_) = Dict{String,String}()

    gene_colours(items::AbstractVector) =
        mapreduce(gene_colours, merge, items; init=Dict{String,String}())

    function gene_colours(node::AbstractDict)
        colours = mapreduce(
            gene_colours,
            merge,
            values(node);
            init=Dict{String,String}(),
        )

        v1 = get(node, Symbol("{regulation/v1}"), nothing)
        v1 isa AbstractDict || return colours

        genes = get(v1, :genes, [])
        digits = ndigits(length(genes))

        merge!(colours, Dict(
            string(get(gene, :name, lpad(index, digits, '0'))) => let
                col = string(gene[:color])
                parse(Colors.Colorant, col)
                col
            end
            for (index, gene) in enumerate(genes)
            if gene isa AbstractDict &&
               get(gene, :color, nothing) isa AbstractString
        ))
    end
end

# ╔═╡ ea610eb8-caa3-4639-8d25-58904e600e52
schedule!, gene_colors, network, schedule_error = try
	spec = JSON.parse(schedule_json; dicttype=Dict{Symbol, Any})
	colors = gene_colours(spec)
	schedule = GRS.Models.parse(schedule_json)
	schedule, colors, Vis.Network(schedule), nothing
catch exception
	nothing, nothing, Vis.Network(), sprint(showerror, exception)
end;

# ╔═╡ cdbda6eb-b5b9-4e9e-bac0-f334c19ae28e
schedule_feedback = @htl("""
<div class="schedule-error">$(something(schedule_error, ""))</div>
""")

# ╔═╡ 50fb4e60-7911-4f68-bd4a-166be8d4d44b
schedule_panel = Layout.DOMElement(
	tag="section",
	attributes=Dict("class" => "schedule-panel"),
		children=[
			schedule_editor,
			schedule_feedback
		]
)

# ╔═╡ 180bf479-2269-4d94-a9f6-b3c060c17ab9
species_options = sort!(unique(
    Symbol(last(split(string(node.name), '.')))
    for node in network.nodes
    if node.kind === :species && !isnothing(node.parent)
))

# ╔═╡ 8bff584c-2987-4584-aaed-b1ce4104d891
group_colors = merge(
    Vis.group_colors(network.groups).colors,
    something(gene_colors, Dict())
)

# ╔═╡ 4d8d901f-d3c2-4b0a-b2c7-4ed8b2a5280d
cytoscape_graph = CytoscapeJS.Cytoscape(network; group_colors, height="600px");

# ╔═╡ 336d7ff6-065e-43b6-8733-eac22ee17454
selected_genes = cytoscape_graph.selection

# ╔═╡ aa943167-f20e-4349-a14b-d512c8005ab0
network_view = Bonito.App(cytoscape_graph);

# ╔═╡ 5f81c7aa-273c-49e2-ae61-f473f810d6c3
selected_species = Bonito.Observable(string.(filter(in((:active, :mrnas, :proteins)), species_options)))

# ╔═╡ d66ca2cc-2d24-4f3e-86dc-b24ff4cdabbe
selection_header = @htl("""
<div style="display: flex; align-items: center; gap: 1rem;">
<label class="dashboard-option">
    <span>model: </span>
    $(@bind selected_model PlutoUI.Select(Vis.paths(network), default=nothing))
</label>
<label class="dashboard-option">
    <span>genes: </span>
    $(value_multiselect(cytoscape_graph.selection, string.(network.groups), size=2))
</label>
<label class="dashboard-option">
    <span>species: </span>
    $(value_multiselect(selected_species, string.(species_options); size=2))
</label>
</div>
""");

# ╔═╡ 27655a09-73a0-4370-933b-2e4fe2e2ecf3
CytoscapeJS.set_filter!(cytoscape_graph, :presentIn, selected_model);

# ╔═╡ 8ab6e516-5c63-4879-8807-2f40b02b782d
simulation = let
	events = NamedTuple[]
	segments = NamedTuple[]
	error = nothing

	trace = function (state; path, from, into=nothing, _...)
		i = length(segments) + 1
		push!(segments, (; i, path, from, to=GRS.Models.t(state), into))
		isnothing(into) || GRS.Models.each_event(state) do t, name, value
			push!(events, (; i, t, name, value))
		end
	end
	
	if run_simulation && !isnothing(schedule!) 
		try
			schedule!(; trace)
		catch exception
			error=sprint(showerror, exception)
		end
	end
	(; events, segments, error)
end

# ╔═╡ d8d20437-7e9d-4d0c-a179-57bb8c50c773
trajectory_view = with_theme(dark_mode ? theme_dark() : Theme()) do
    species = selected_species[]
    genes = string.(network.groups)

    fig = Figure(size=(900, 180 * max(1, length(species))))
    axes = Axis[]

    for (row, kind) in enumerate(species)
        axis = Axis(
            fig[row, 1];
            ylabel=kind,
            xlabel=row == length(species) ? "time" : "",
            backgroundcolor=:transparent,
            xzoomlock=false,
            yzoomlock=true,
            xpanlock=false,
            ypanlock=true,
        )
        push!(axes, axis)

        for gene in genes
            data = filter(
                event -> event.name == Symbol(gene, ".", kind),
                simulation.events,
            )
            isempty(data) && continue

            order = sortperm(getproperty.(data, :t))
            stairs!(
                axis,
                getproperty.(data, :t)[order],
                getproperty.(data, :value)[order];
                color=group_colors[gene],
                step=:post,
                visible=lift(selected_genes) do selected
                    gene in selected
                end,
            )
        end

        autolimits!(axis)
    end

    length(axes) > 1 && linkxaxes!(axes...)

    if !isempty(simulation.events)
        tmin, tmax = extrema(getproperty.(simulation.events, :t))
        xlims!(axes[1], tmin, tmax)
    end

    fig
end;

# ╔═╡ 1d72c85a-a4c4-4b3d-b69b-b7354cf8df0a
trajectory_panel = @htl("""
<div class="trajectory-panel">
    $(trajectory_view)

    <style>
        .trajectory-panel {
            width: 100%;
            min-width: 0;
            box-sizing: border-box;
            overflow: hidden;

            border: 1px solid #d4d4d8;
            border-radius: 6px;
        }

        .trajectory-panel > * {
            width: 100% !important;
            max-width: 100% !important;
        }

        .trajectory-panel canvas,
        .trajectory-panel svg {
            display: block;
            width: 100% !important;
            height: auto !important;
        }

        @media (prefers-color-scheme: dark) {
            .trajectory-panel {
                border-color: #6b7280;
            }
        }
    </style>
</div>
""")

# ╔═╡ 8eb4f0d7-ef76-4ca7-a117-a48685c12667
show_notebook_control =
    @bind show_notebook PlutoUI.Switch(default=false)

# ╔═╡ a4bcea34-2ee9-48ee-bf1a-6d1fd6256c91
begin
    favicon = "data:image/png;base64," *
    base64encode(read(joinpath(@__DIR__, "assets", "logo.png")))
    docs_url = "http://localhost:8001"
    app_header = @htl("""
<header class="grs-header">
    <div class="logo-wrapper">
        $(PlutoUI.LocalResource(joinpath(@__DIR__, "assets", "logo.png")))
    </div>
    <div>
        <div class="grs-title">Gene Regulatory Systems</div>
    </div>
    <a class="docs-link" href="$(docs_url)" target="_blank">
        docs <span class="external-link">↗</span>
    </a>
    <label class="notebook-toggle" title="show notebook internals">
        <span>dev</span>
        $(show_notebook_control)
    </label>
</header>

<style>
    .grs-header {
        box-sizing: border-box;
        display: flex;
        align-items: center;
        gap: 0.85rem;
        width: 100%;
        height: 60px;
        padding: 0 1rem;
        margin-bottom: 0.75rem;
        border-bottom: 1px solid color-mix(
            in srgb,
            currentColor 15%,
            transparent
        );
        font-family: Montserrat, sans-serif;
    }
    .logo-wrapper {
        filter: drop-shadow(-1.4px 2px 2px rgba(150, 150, 150, 0.3));
    }
    .logo-wrapper img {
        width: 40px;
        height: 40px;
        object-fit: contain;
    }
    @keyframes grs-logo-spin {
        to {
            transform: rotate(-360deg);
        }
    }
    pluto-editor:has(pluto-cell.running, pluto-cell.queued)
        .grs-header img {
        animation: grs-logo-spin 0.3s linear infinite;
    }

    .grs-title {
        font-size: 1.15rem;
        font-weight: 400;
        line-height: 1.2;
    }

    .docs-link {
        margin-left: auto;
        color: inherit;
        transition: opacity 120ms ease;
    }

    body:not(:has(.notebook-toggle input:checked)) header#pluto-nav {
        display: none !important;
    }
    .notebook-toggle {
        display: flex;
        gap: 0.35rem;
        opacity: 0.3;
        font-family: Montserrat, sans-serif;
        font-size: 0.7rem;
        cursor: pointer;
        transition: opacity 120ms ease;
    }

    .notebook-toggle:hover,
    .notebook-toggle:has(input:checked) {
        opacity: 0.8;
    }

</style>
<script>
    let icon = document.head.querySelector('link[rel~="icon"]')

    if (!icon) {
        icon = document.createElement("link")
        icon.rel = "icon"
        document.head.appendChild(icon)
    }

    icon.type = "image/png"
    icon.href = $(favicon)
    document.title = "Gene Regulatory Systems"
    document.body.classList.add("grs-ready")
</script>
""");
end;

# ╔═╡ c7d9d65d-10d4-4e88-a034-cd258438f54c
let
    column(children...) = Layout.DOMElement(
        tag="div",
        attributes=Dict(
            "style" => "display: grid; gap: 0.75rem; min-width: 0;",
        ),
        children=collect(children),
    )

    columns = Layout.DOMElement(
        tag="div",
        attributes=Dict(
            "style" => """
                display: grid;
                grid-template-columns: minmax(300px, 1fr) minmax(0, 2fr);
                gap: 1rem;
                align-items: start;
                width: 100%;
                min-width: 0;
            """,
        ),
        children=[
            column(schedule_panel),
            column(network_view, selection_header, trajectory_panel),
        ],
    )

    dashboard = Layout.DOMElement(
        tag="div",
        attributes=Dict("class" => "grs-dashboard"),
        children=[app_header, schedule_header, columns],
    )

    PlutoUI.WideCell(dashboard; max_width=1600)
end

# ╔═╡ Cell order:
# ╠═c7d9d65d-10d4-4e88-a034-cd258438f54c
# ╟─34a8deb2-0db9-456c-ba38-eabc3254ab50
# ╠═40d24905-5708-4580-9ae5-71e1c1de01a9
# ╠═12e75652-a3fd-11f1-9335-a986ee44237f
# ╟─76de62b9-be88-4b96-b801-b40872dfa0be
# ╠═a2866674-5e4e-4738-b18e-90d122d2d161
# ╠═7abbbe81-3d0b-4d94-a69f-32c66514687f
# ╠═8392f056-369c-45fb-9c00-719172e82616
# ╠═cdbda6eb-b5b9-4e9e-bac0-f334c19ae28e
# ╠═50fb4e60-7911-4f68-bd4a-166be8d4d44b
# ╟─69993fb8-daca-4b61-a711-907009dbee24
# ╠═ea610eb8-caa3-4639-8d25-58904e600e52
# ╟─16d523e4-ebec-4d6e-b5e1-7aaf949ec5fc
# ╠═8bff584c-2987-4584-aaed-b1ce4104d891
# ╠═4d8d901f-d3c2-4b0a-b2c7-4ed8b2a5280d
# ╠═27655a09-73a0-4370-933b-2e4fe2e2ecf3
# ╠═aa943167-f20e-4349-a14b-d512c8005ab0
# ╟─c7cd57d6-464d-4dbc-b955-699802e60be7
# ╠═336d7ff6-065e-43b6-8733-eac22ee17454
# ╠═180bf479-2269-4d94-a9f6-b3c060c17ab9
# ╠═5f81c7aa-273c-49e2-ae61-f473f810d6c3
# ╠═d66ca2cc-2d24-4f3e-86dc-b24ff4cdabbe
# ╠═8ab6e516-5c63-4879-8807-2f40b02b782d
# ╠═d8d20437-7e9d-4d0c-a179-57bb8c50c773
# ╠═1d72c85a-a4c4-4b3d-b69b-b7354cf8df0a
# ╟─35e4e15a-e1f7-42c1-8ded-5714234811c6
# ╠═4a29ccba-5cdc-47e1-953c-5eee6d1edbe7
# ╠═8eb4f0d7-ef76-4ca7-a117-a48685c12667
# ╠═a4bcea34-2ee9-48ee-bf1a-6d1fd6256c91
# ╠═375349c2-33bb-45fb-b412-672c678b93b3
# ╠═0187c445-f7a0-4eae-b55d-2bd955df2c2c
