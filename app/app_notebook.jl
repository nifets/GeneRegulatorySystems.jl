### A Pluto.jl notebook ###
# v1.0.3

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

# ╔═╡ 12e75652-a3fd-11f1-9335-a986ee44237f
begin
    import Pkg
    Pkg.activate(".")
    using Revise

    using GeneRegulatorySystems
    using CytoscapeJS
    using Bonito
    using PlutoUI
    const Layout = PlutoUI.ExperimentalLayout
    using HypertextLiteral
    using Base64
    import JSON
    import Colors
    const GRS = GeneRegulatorySystems
    const Vis = GRS.Visualisation
    Revise.includet(@__MODULE__, "src/JSONEditor.jl")
end;

# ╔═╡ 34a8deb2-0db9-456c-ba38-eabc3254ab50
md"""
## Logic
"""

# ╔═╡ 76de62b9-be88-4b96-b801-b40872dfa0be
md"""
### Schedule
"""

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
schedule_selector = @htl("""
<label class="dashboard-option">
    <span>schedule:</span>
    $(@bind selected_schedule PlutoUI.Select(schedule_options))
</label>
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

# ╔═╡ 58addff4-9608-48eb-9742-64aa06a2c6a5
schedule!, gene_colors, schedule_error = try
	spec = JSON.parse(schedule_json; dicttype=Dict{Symbol, Any})
	colors = gene_colours(spec)
	GRS.Models.parse(schedule_json), colors, nothing
catch exception
	nothing, nothing, sprint(showerror, exception)
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

# ╔═╡ 16d523e4-ebec-4d6e-b5e1-7aaf949ec5fc
md"""
### Network
"""

# ╔═╡ bd8d124f-6362-43ff-a4cf-44e6d41b22af
network = isnothing(schedule!) ? Vis.Network() : Vis.Network(schedule!);

# ╔═╡ 8bff584c-2987-4584-aaed-b1ce4104d891
group_colors = merge(
    Vis.group_colors(network.groups).colors,
    something(gene_colors, Dict())
)

# ╔═╡ d66ca2cc-2d24-4f3e-86dc-b24ff4cdabbe
model_selector = @htl("""
<label class="dashboard-option">
    <span>model: </span>
    $(@bind selected_model PlutoUI.Select(Vis.paths(network), default=nothing))
</label>
""")

# ╔═╡ 4d8d901f-d3c2-4b0a-b2c7-4ed8b2a5280d
cytoscape_graph = CytoscapeJS.Cytoscape(network; group_colors, height="600px");

# ╔═╡ 27655a09-73a0-4370-933b-2e4fe2e2ecf3
CytoscapeJS.set_filter!(cytoscape_graph, :presentIn, selected_model);

# ╔═╡ aa943167-f20e-4349-a14b-d512c8005ab0
network_view = Bonito.App(cytoscape_graph);

# ╔═╡ 8eb4f0d7-ef76-4ca7-a117-a48685c12667
show_notebook_control =
    @bind show_notebook PlutoUI.Switch(default=false)

# ╔═╡ 35e4e15a-e1f7-42c1-8ded-5714234811c6
md"""
## Styling
"""

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
</script>
""");
end;

# ╔═╡ c7d9d65d-10d4-4e88-a034-cd258438f54c
let 
	dashboard = Layout.DOMElement(
		tag="div",
		attributes=Dict("class" => "grs-dashboard"),
		children = [
			app_header,
			Layout.grid(
				[
					schedule_selector model_selector
					schedule_panel network_view
				];
	            column_gap="1rem",
	            row_gap="0.75rem",
	            style=Dict(
	                "grid-template-columns" =>
	                    "minmax(300px, 1fr) minmax(0, 2fr)",
	                "align-items" => "start",
	            ),
			),
		],
	)
	PlutoUI.WideCell(dashboard; max_width=1600)
end

# ╔═╡ 375349c2-33bb-45fb-b412-672c678b93b3
@htl("""
 <style id="grs-dashboard-styles">
    .grs-dashboard {
        box-sizing: border-box;
        width: 100%;
        height: 100vh;
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

# ╔═╡ Cell order:
# ╠═c7d9d65d-10d4-4e88-a034-cd258438f54c
# ╟─34a8deb2-0db9-456c-ba38-eabc3254ab50
# ╠═12e75652-a3fd-11f1-9335-a986ee44237f
# ╟─76de62b9-be88-4b96-b801-b40872dfa0be
# ╠═a2866674-5e4e-4738-b18e-90d122d2d161
# ╠═7abbbe81-3d0b-4d94-a69f-32c66514687f
# ╠═8392f056-369c-45fb-9c00-719172e82616
# ╠═58addff4-9608-48eb-9742-64aa06a2c6a5
# ╠═cdbda6eb-b5b9-4e9e-bac0-f334c19ae28e
# ╠═50fb4e60-7911-4f68-bd4a-166be8d4d44b
# ╠═69993fb8-daca-4b61-a711-907009dbee24
# ╟─16d523e4-ebec-4d6e-b5e1-7aaf949ec5fc
# ╠═bd8d124f-6362-43ff-a4cf-44e6d41b22af
# ╠═8bff584c-2987-4584-aaed-b1ce4104d891
# ╠═d66ca2cc-2d24-4f3e-86dc-b24ff4cdabbe
# ╠═4d8d901f-d3c2-4b0a-b2c7-4ed8b2a5280d
# ╠═27655a09-73a0-4370-933b-2e4fe2e2ecf3
# ╠═aa943167-f20e-4349-a14b-d512c8005ab0
# ╠═8eb4f0d7-ef76-4ca7-a117-a48685c12667
# ╟─35e4e15a-e1f7-42c1-8ded-5714234811c6
# ╠═a4bcea34-2ee9-48ee-bf1a-6d1fd6256c91
# ╠═375349c2-33bb-45fb-b412-672c678b93b3
