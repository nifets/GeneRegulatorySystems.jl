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

# ╔═╡ 12e75652-a3fd-11f1-9335-a986ee44237f
begin
    import Pkg
    Pkg.activate(".")
    using Revise

    using GeneRegulatorySystems
    using CytoscapeJS
    using Bonito
    using PlutoUI
    using HypertextLiteral
    using Base64
    const GRS = GeneRegulatorySystems
    const Vis = GRS.Visualisation
end;

# ╔═╡ ade08263-48f7-440e-8ec2-ac6321ad6354
function json_editor(contents; height="400px")
    @htl("""
    <div class="json-editor" style="--editor-height: $height">
        <textarea hidden>$contents</textarea>
        <div class="editor"></div>

        <style>
            .json-editor .cm-editor {
                height: var(--editor-height);
                border: 1px solid #d4d4d8;
                border-radius: 6px;
                font-size: 12px;
            }

            .json-editor .cm-scroller {
                overflow: auto;
            }

            @media (prefers-color-scheme: dark) {
                .json-editor .cm-editor {
                    border-color: #6b7280;
                }

                .json-editor .cm-activeLine,
                .json-editor .cm-activeLineGutter {
                    background-color: rgba(255, 255, 255, 0.06) !important;
                }
            }

            @media (prefers-color-scheme: light) {
                .json-editor .cm-activeLine,
                .json-editor .cm-activeLineGutter {
                    background-color: rgba(0, 0, 0, 0.06) !important;
                }
            }
        </style>

        <script>
            const root = currentScript.parentElement
            const parent = root.querySelector(".editor")
            const initialValue = root.querySelector("textarea").value

            const { basicSetup, EditorView } =
                await import("https://esm.sh/codemirror@6.0.2")
            const { json } =
                await import("https://esm.sh/@codemirror/lang-json@6.0.2")
            const { oneDark } =
                await import("https://esm.sh/@codemirror/theme-one-dark@6.1.3")

            const colorScheme =
                window.matchMedia("(prefers-color-scheme: dark)")

            let timeout
            let view

            const createEditor = contents => new EditorView({
                doc: contents,
                parent,
                extensions: [
                    basicSetup,
                    json(),
                    colorScheme.matches ? oneDark : [],
                    EditorView.updateListener.of(update => {
                        if (!update.docChanged) return

                        clearTimeout(timeout)
                        timeout = setTimeout(() => {
                            root.value = update.state.doc.toString()
                            root.dispatchEvent(new CustomEvent("input"))
                        }, 250)
                    }),
                ],
            })

            const updateTheme = () => {
                const contents = view.state.doc.toString()
                view.destroy()
                view = createEditor(contents)
            }

            root.value = initialValue
            view = createEditor(initialValue)

            colorScheme.addEventListener("change", updateTheme)

            invalidation.then(() => {
                clearTimeout(timeout)
                colorScheme.removeEventListener("change", updateTheme)
                view.destroy()
            })
        </script>
    </div>
    """)
end

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
schedule_editor = @bind schedule_json json_editor(read(selected_schedule, String); height="800px")

# ╔═╡ 5f174839-b6f2-4e9b-b2a4-ffbb525a3f11
schedule! = GRS.Models.parse(schedule_json);

# ╔═╡ bd8d124f-6362-43ff-a4cf-44e6d41b22af
network = Vis.Network(schedule!);

# ╔═╡ d66ca2cc-2d24-4f3e-86dc-b24ff4cdabbe
model_selector = @htl("""
<label class="dashboard-option">
    <span>model: </span>
    $(@bind selected_model PlutoUI.Select(Vis.paths(network)))
</label>
""")

# ╔═╡ 4d8d901f-d3c2-4b0a-b2c7-4ed8b2a5280d
cytoscape_graph = CytoscapeJS.Cytoscape(network; height="600px");

# ╔═╡ 27655a09-73a0-4370-933b-2e4fe2e2ecf3
CytoscapeJS.set_filter!(cytoscape_graph, :presentIn, selected_model);

# ╔═╡ aa943167-f20e-4349-a14b-d512c8005ab0
network_view = Bonito.App(cytoscape_graph)

# ╔═╡ 8eb4f0d7-ef76-4ca7-a117-a48685c12667
show_notebook_control =
    @bind show_notebook PlutoUI.Switch(default=false)

# ╔═╡ a4bcea34-2ee9-48ee-bf1a-6d1fd6256c91
begin
    favicon = "data:image/png;base64," *
    base64encode(read(joinpath(@__DIR__, "assets", "logo.png")))
    app_header = @htl("""
<header class="grs-header">
    $(PlutoUI.LocalResource(joinpath(@__DIR__, "assets", "logo.png")))

    <div>
        <div class="grs-title">Gene Regulatory Systems</div>
    </div>
    <label class="notebook-toggle" title="Show notebook internals">
        $(show_notebook_control)
        <span>Dev</span>
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

    .grs-header img {
        width: 40px;
        height: 40px;
        object-fit: contain;
    }

    .grs-title {
        font-size: 1.15rem;
        font-weight: 400;
        line-height: 1.2;
    }

    body:not(:has(.notebook-toggle input:checked)) header#pluto-nav {
        display: none !important;
    }
    .notebook-toggle {
        display: flex;
        align-items: center;
        gap: 0.35rem;
        margin-left: auto;
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

# ╔═╡ 8a1cc60e-8376-495f-af9d-8bada41627a0
PlutoUI.WideCell(
    @htl("""
    <div class="grs-dashboard">
        $(app_header)

        $(PlutoUI.ExperimentalLayout.grid(
            [
                schedule_selector model_selector
                schedule_editor   network_view
            ];
            column_gap="1rem",
            row_gap="0.75rem",
            style=Dict(
                "grid-template-columns" =>
                    "minmax(300px, 1fr) minmax(0, 2fr)",
                "align-items" => "start",
            ),
        ))
    </div>

    <style>
        .grs-dashboard {
            width: 100%;
            font-family: Montserrat, sans-serif;
        }
        body:not(:has(.notebook-toggle input:checked))
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
        
        body:not(:has(.notebook-toggle input:checked)) #helpbox-wrapper {
            display: none !important;
        }
        /* Hide the dashboard cell's source code outside dev mode. */
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
    """);
    max_width=1600,
)

# ╔═╡ Cell order:
# ╠═8a1cc60e-8376-495f-af9d-8bada41627a0
# ╠═12e75652-a3fd-11f1-9335-a986ee44237f
# ╟─ade08263-48f7-440e-8ec2-ac6321ad6354
# ╠═a2866674-5e4e-4738-b18e-90d122d2d161
# ╠═7abbbe81-3d0b-4d94-a69f-32c66514687f
# ╠═8392f056-369c-45fb-9c00-719172e82616
# ╠═5f174839-b6f2-4e9b-b2a4-ffbb525a3f11
# ╠═bd8d124f-6362-43ff-a4cf-44e6d41b22af
# ╠═d66ca2cc-2d24-4f3e-86dc-b24ff4cdabbe
# ╠═4d8d901f-d3c2-4b0a-b2c7-4ed8b2a5280d
# ╠═27655a09-73a0-4370-933b-2e4fe2e2ecf3
# ╠═aa943167-f20e-4349-a14b-d512c8005ab0
# ╟─8eb4f0d7-ef76-4ca7-a117-a48685c12667
# ╠═a4bcea34-2ee9-48ee-bf1a-6d1fd6256c91
