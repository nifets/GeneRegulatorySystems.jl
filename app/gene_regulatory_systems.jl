### A Pluto.jl notebook ###
# v1.0.1

#> [frontmatter]
#> title = "Gene Regulatory Systems"
#> description = "Interactive dashboard for simulating and exploring gene regulatory network models."

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
    Pkg.activate(@__DIR__)
    using Base64
    using HypertextLiteral
end

# ╔═╡ 29eb1509-3cea-4c50-9d8e-24bcd8a54bbb
md"""
## Layout
"""

# ╔═╡ 34a8deb2-0db9-456c-ba38-eabc3254ab50
md"""
## Logic
"""

# ╔═╡ 76de62b9-be88-4b96-b801-b40872dfa0be
md"""
### Schedule
"""

# ╔═╡ c81d5f60-2a47-4e93-b6d1-5f39c284ae71
schedule_dir = let
    dir = joinpath(@__DIR__, "schedules")
    isdir(dir) || mkpath(dir)
    dir
end;

# ╔═╡ 16d523e4-ebec-4d6e-b5e1-7aaf949ec5fc
md"""
### Network
"""

# ╔═╡ 8d1f0c53-6a24-4e91-b7d5-2c9a04f61b38
md"""
### Gene selection
"""

# ╔═╡ c7cd57d6-464d-4dbc-b955-699802e60be7
md"""
### Run simulation
"""

# ╔═╡ a1d5e8f2-4b73-4c96-8e05-71fc32a9b6d8
function fraction_reporter(emit, total; interval = 0.1)
    advanced = Dict{String, Float64}()
    guard = ReentrantLock()
    last = Ref(0.0)
    (stage; done = nothing, path = nothing, _...) -> begin
        done === nothing && return
        fraction = lock(guard) do
            advanced[something(path, "")] = done
            time() - last[] < interval && return nothing
            last[] = time()
            min(sum(values(advanced)) / max(total, eps()), 1.0)
        end
        fraction === nothing || emit(fraction)
    end
end

# ╔═╡ 525c8855-9758-4d61-a8da-a616396be848
md"""
### Trajectory View
"""

# ╔═╡ 76d9a00b-10fa-4248-b460-6c5abb2424dd
md"""
### Phase space
"""

# ╔═╡ 35e4e15a-e1f7-42c1-8ded-5714234811c6
md"""
## Styling
"""

# ╔═╡ 4a29ccba-5cdc-47e1-953c-5eee6d1edbe7
dark_mode_probe = @bind dark_mode @htl("""
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

# ╔═╡ 6a3d1f28-04bc-4e7a-9c15-8b2f70d6ae41
hide_notebook_styles = @htl("""
<style id="grs-hide-notebook">
    body:not(:has(.notebook-toggle input:checked)) {
        & header#pluto-nav,
        & #helpbox-wrapper,
        & footer,
        & pluto-editor > main > preamble {
            display: none !important;
        }

        & pluto-cell:has(.dashboard) {
            & > pluto-shoulder,
            & > button.add_cell,
            & > pluto-trafficlight,
            & > pluto-input,
            & > pluto-runarea {
                display: none !important;
            }
        }
    }

    body:has(.dashboard):not(:has(.notebook-toggle input:checked)) {
        & pluto-cell:not(:has(.dashboard)):not(:has(.grs-splash)):not(:has(pluto-logs-container)) {
            display: none !important;
        }

        & pluto-cell:has(pluto-logs-container) > *:not(pluto-logs-container) {
            display: none;
        }
    }
</style>
""")

# ╔═╡ 375349c2-33bb-45fb-b412-672c678b93b3
dashboard_styles = @htl("""
<style id="grs-dashboard-styles">
    body:not(:has(.notebook-toggle input:checked)) {
        & pluto-notebook {
            display: grid;
            grid-template-columns:
                minmax(160px, 0.75fr) minmax(0, 0.25fr) minmax(0, 2fr);
            grid-template-areas:
                "header          header   header"
                "schedule-header progress network-header"
                "left            left     network"
                "left            left     trajectory";
            grid-template-rows: auto auto auto 1fr;
            align-items: start;
            gap: 1rem;
            max-width: none;
            width: 100%;
            margin: 0;
            padding: 0 1.5rem 2rem;
            box-sizing: border-box;
        }

        & pluto-editor > main {
            max-width: none;
            width: 100%;
            margin: 0;
            padding: 0;
            align-self: stretch;
        }

        & pluto-cell:has(.dashboard) {
            margin: 0;
            min-height: 0;
            min-width: 0;
        }
    }

    @media (max-width: 900px) {
        body:not(:has(.notebook-toggle input:checked)) pluto-notebook {
            grid-template-columns: minmax(0, 1fr);
            grid-template-areas:
                "header" "schedule-header" "progress"
                "network-header" "network" "trajectory" "left";
        }
    }

    pluto-cell:has(.area-header) { grid-area: header; }
    pluto-cell:has(.area-left) { grid-area: left; }
    pluto-cell:has(.area-schedule-header) { grid-area: schedule-header; }
    pluto-cell:has(.area-network-header) { grid-area: network-header; }
    pluto-cell:has(.area-network) { grid-area: network; }
    pluto-cell:has(.area-trajectory) { grid-area: trajectory; }

    pluto-cell:has(pluto-logs-container) {
        grid-area: progress;
        padding: 0;
        margin-top: 0;
        min-height: 0;

        & > pluto-output {
            padding: 0 10px;
        }

        & pluto-logs-container:not(:empty) {
            background: none;
            padding: 0;
            margin: 0 1.3rem 0 0;
        }

        & pluto-log-dot-positioner {
            background: none;
            margin: 0;
        }

        & pluto-log-dot-positioner:not(.Progress) {
            display: none;
        }

        & pluto-log-dot {
            padding: 0;
        }

        & pluto-log-icon {
            display: none;
        }

        & pluto-progress-bar-container {
            flex: 1;
            outline: 1px solid color-mix(in srgb, currentColor 15%, transparent);
            outline-offset: -1px;
        }
    }

    .dashboard {
        display: grid;
        gap: 0.75rem;
        min-width: 0;
        font-family: Montserrat, sans-serif;
    }

    .dashboard-option {
        display: flex;
        align-items: center;
        gap: 0.5rem;
        font-family: Montserrat, sans-serif;
        font-size: 0.8rem;
        font-weight: 300;

        & select,
        & input,
        & option {
            font-family: Montserrat, sans-serif;
        }

        & select,
        & input {
            font-size: 0.7rem;
        }

    }

    .area-schedule-header .dashboard-header {
        flex-wrap: nowrap;
    }

    .area-schedule-header .dashboard-option select {
        min-width: 0;
        max-width: 8rem;
        text-overflow: ellipsis;
    }


    .dashboard-header {
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: 1rem;
    }

    .dashboard-header.stacked {
        align-items: flex-start;
    }

    .dashboard-option.stacked {
        flex-direction: column;
        align-items: flex-start;
    }

    .dashboard-controls {
        display: flex;
        flex-direction: column;
        gap: 0.35rem;
    }

    .network-provenance {
        margin-left: auto;
        text-align: right;
        font-size: 0.7rem;
        opacity: 0.6;
        line-height: 1.3;
    }

    .panel {
        box-sizing: border-box;
        width: 100%;
        min-width: 0;
        overflow: hidden;
        clip-path: inset(0 round 6px);
        border: 1px solid color-mix(in srgb, currentColor 15%, transparent);
        border-radius: 6px;

        & > * {
            width: 100% !important;
            height: 100% !important;
            max-width: 100% !important;
        }
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
</style>
""")


# ╔═╡ b3f27a51-6c84-4d19-9f52-8e1c47a0d6b2
reload_control = @bind reload_count @htl("""
<span class="grs-reload">
<input type="button" value="↻" title="reload">
<style>
    .grs-reload input {
        background: none;
        border: none;
        padding: 0;
        color: inherit;
        font-size: 1rem;
        line-height: 1;
        cursor: pointer;
        opacity: 0.6;
        transition: opacity 0.15s;
    }
    .grs-reload input:hover {
        opacity: 1;
    }
</style>
<script>
    const root = currentScript.parentElement
    const button = root.querySelector("input")
    let count = 0
    root.value = count

    button.addEventListener("click", () => {
        root.value = ++count
        root.dispatchEvent(new CustomEvent("input"))
    })
</script>
</span>
""");

# ╔═╡ 0187c445-f7a0-4eae-b55d-2bd955df2c2c
begin
    logo_data = "data:image/png;base64," *
        base64encode(read(joinpath(@__DIR__, "assets", "logo.png")))

    splash_screen = @htl("""
<div class="grs-splash">
    <button class="grs-splash-dismiss" type="button" title="dismiss"
            onclick="document.body.classList.add('grs-ready')">&times;</button>
    <div class="grs-splash-content">
        <div class="grs-splash-title">Gene Regulatory Systems</div>
        <div class="grs-splash-logo">
            <img src=$(logo_data)>
        </div>
        <div class="grs-splash-loading">loading</div>
    </div>
</div>

<script>
    let icon = document.head.querySelector('link[rel~="icon"]')

    if (!icon) {
        icon = document.createElement("link")
        icon.rel = "icon"
        document.head.appendChild(icon)
    }

    icon.type = "image/png"
    icon.href = $(logo_data)
    document.title = "Gene Regulatory Systems"

    let started = false

    const check = () => {
        if (document.querySelector("pluto-cell.running, pluto-cell.queued")) {
            started = true
        } else if (started) {
            document.body.classList.add("grs-ready")
            observer.disconnect()
        }
    }

    const observer = new MutationObserver(check)
    observer.observe(document.body, {
        subtree: true,
        childList: true,
        attributes: true,
        attributeFilter: ["class"],
    })
    check()

    invalidation.then(() => observer.disconnect())
</script>

<style>
    body:not(.grs-ready) pluto-cell {
        visibility: hidden !important;
    }

    body:not(.grs-ready) loading-bar {
        z-index: 100000;
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

    .grs-splash-dismiss {
        position: absolute;
        top: 1rem;
        right: 1rem;
        border: none;
        background: transparent;
        color: inherit;
        font-size: 1.25rem;
        opacity: 0.4;
        cursor: pointer;
    }

    .grs-splash-dismiss:hover {
        opacity: 1;
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

    .grs-splash-logo:hover img {
        animation-duration: 1s;
    }

    .grs-splash-loading {
        font-size: 0.8rem;
        font-family: Montserrat, sans-serif;
    }

    body.grs-ready .grs-splash,
    body.grs-ready pluto-cell:has(.grs-splash) {
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
    reload_count

    using Revise

    using GeneRegulatorySystems
    const GRS = GeneRegulatorySystems
    const Vis = GRS.Visualisation

    using CytoscapeJS
    using PlutoUI
    const Layout = PlutoUI.ExperimentalLayout
    import JSON
    import Colors
    import ProgressLogging

    using Bonito
    using WGLMakie
    WGLMakie.activate!()
    using GRSApp: Trajectories, PageState, JSONEditor
end;

# ╔═╡ 032084a8-991a-4eeb-a5e9-78110b2b7cbe
dashboard_panel(child; height=nothing) = isnothing(child) ? nothing : Layout.DOMElement(
        tag="div",
        attributes=isnothing(height) ?
            Dict("class" => "panel") :
            Dict("class" => "panel", "style" => "height: $(height)px;"),
        children=[child],
    )

# ╔═╡ 116c3c09-d5b1-4c86-8c08-c470df26099b
dashboard_area(name, children...) = Layout.DOMElement(
        tag="div",
        attributes=Dict("key" => name, "class" => "dashboard area-$(name)"),
        children=collect(Iterators.filter(!isnothing, children)),
    )

# ╔═╡ 8eb4f0d7-ef76-4ca7-a117-a48685c12667
show_notebook_control =
    @bind show_notebook PlutoUI.Switch(default=false)

# ╔═╡ a2866674-5e4e-4738-b18e-90d122d2d161
schedule_options = let
    entries(dir) = [
        file => replace(basename(file), ".schedule.json" => "")
        for file in sort(readdir(dir; join=true))
        if endswith(file, ".schedule.json")
    ]
    (;
        examples=entries(GRS.SPECIFICATION_EXAMPLES),
        user=entries(schedule_dir),
    )
end;

# ╔═╡ d4a916c8-7b52-49f6-a8e3-1c06d5b74f29
schedule_control = @htl("""
<span>
$(@bind selected_schedule PlutoUI.Select([
    schedule_options.examples; schedule_options.user
]))
<script>
    const select = currentScript.parentElement.querySelector("select")
    const split = $(length(schedule_options.examples))
    const options = Array.from(select.options)
    const chosen = select.value

    const group = (label, members) => {
        if (!members.length) return
        const holder = document.createElement("optgroup")
        holder.label = label
        members.forEach(option => holder.appendChild(option))
        select.appendChild(holder)
    }

    group("examples", options.slice(0, split))
    group("user", options.slice(split))
    select.value = chosen
</script>
</span>
""");


# ╔═╡ 7abbbe81-3d0b-4d94-a69f-32c66514687f
schedule_header = @htl("""
<div class="dashboard-header">
<label class="dashboard-option">
    <span>schedule:</span>
    $(schedule_control)
</label>
<label class="dashboard-option">
    <span>simulate:</span>
    $(@bind run_simulation PlutoUI.Switch(default=false))
</label>
</div>
""");

# ╔═╡ b4f7c2d9-1e58-4a63-9c07-8d5e21f4a3b6
dashboard_area("schedule-header", schedule_header)

# ╔═╡ 8392f056-369c-45fb-9c00-719172e82616
schedule_editor = @bind schedule_json JSONEditor(read(selected_schedule, String); height="800px");

# ╔═╡ e57c02d3-9f18-4a6b-b30c-6d2e81f45c93
autosave = if startswith(selected_schedule, schedule_dir) &&
        !ismissing(schedule_json) &&
        read(selected_schedule, String) != schedule_json
    write(selected_schedule, schedule_json)
end;

# ╔═╡ 69993fb8-daca-4b61-a711-907009dbee24
begin
    declared_colors(_) = Dict{String,String}()

    declared_colors(items::AbstractVector) =
        mapreduce(declared_colors, merge, items; init=Dict{String,String}())

    function declared_colors(node::AbstractDict)
        colours = mapreduce(
            declared_colors,
            merge,
            values(node);
            init=Dict{String,String}(),
        )

        genes = get(node, :genes, nothing)
        genes isa AbstractVector || return colours

        digits = ndigits(length(genes))

        merge!(colours, Dict(
            string(get(gene, :name, lpad(index, digits, '0'))) => let
                col = string(get(gene, :color, get(gene, :colour, "")))
                parse(Colors.Colorant, col)
                col
            end
            for (index, gene) in enumerate(genes)
            if gene isa AbstractDict &&
               get(gene, :color, get(gene, :colour, nothing)) isa AbstractString
        ))

        colours
    end
end


# ╔═╡ ea610eb8-caa3-4639-8d25-58904e600e52
schedule!, gene_colors, network, schedule_error = try
	spec = JSON.parse(schedule_json; dicttype=Dict{Symbol, Any})
	colors = declared_colors(spec)
	schedule = GRS.Models.parse(schedule_json)
	schedule, colors, Vis.Network(schedule), nothing
catch exception
	nothing, nothing, Vis.Network(), sprint(showerror, exception)
end;

# ╔═╡ b7e4c1a0-5d92-4f3e-8a61-2c9f0d7b4e83
provenance_lines = isnothing(schedule!) ?
    Dict{String, Vector{String}}() : Vis.provenance(schedule!)

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
);

# ╔═╡ e2b7f4a3-8c15-4d67-9e02-5fa31c7b8d46
total_duration = isnothing(schedule!) ? 0.0 : let total = Ref(0.0)
    schedule!(
        GRS.Models.FlatState();
        dryrun = (primitive!, x, dt; _...) ->
            (isfinite(dt) && dt > 0 && (total[] += dt); nothing),
        parallel = false,
    )
    total[]
end

# ╔═╡ 180bf479-2269-4d94-a9f6-b3c060c17ab9
track_options = intersect([:activity, :elongations, :premrnas, :mrnas, :proteins, :active, :inactive], unique(vcat(
    [:activity],
    [
        Symbol(last(split(string(node.name), '.')))
        for node in network.nodes
        if node.kind === :species && !isnothing(node.parent)
    ],
)))

# ╔═╡ 5f81c7aa-273c-49e2-ae61-f473f810d6c3
default_tracks = string.(filter(in((:activity, :mrnas, :proteins)), track_options))

# ╔═╡ 8bff584c-2987-4584-aaed-b1ce4104d891
group_colors = merge(
    Vis.group_colors(network.groups).colors,
    something(gene_colors, Dict())
);

# ╔═╡ 5e2c8a17-4b39-4f6d-9c81-7d0a3b5e62f4
model_options = Vis.paths(network)

# ╔═╡ 6c1e4b90-2a77-4d3e-9f58-31b0c7ea52d4
shared_species_control = @bind show_shared_species PlutoUI.Switch(default=false);

# ╔═╡ d9a6e3c2-7fb4-4156-ac83-4eb21f9d6a05
model_control = @bind selected_model PlutoUI.Select(
    model_options;
    default=length(model_options) == 1 ? only(model_options) : nothing,
);

# ╔═╡ c8f5d2b1-6ea3-4045-9b72-3da10e8c5f94
network_description = let
    lines = get(provenance_lines, something(selected_model, ""), String[])
    @htl("""
    <div class="network-provenance">
        $((@htl("<div>$(line)</div>") for line in lines))
    </div>
    """)
end;

# ╔═╡ 0d5a9f31-8e46-4c02-b7d1-9a2f6c48e713
network_header = @htl("""
<div class="dashboard-header stacked">
    <div class="dashboard-controls">
        <label class="dashboard-option">
            <span>model:</span>
            $(model_control)
        </label>
        <label class="dashboard-option">
            <span>shared species:</span>
            $(shared_species_control)
        </label>
    </div>
    $(network_description)
</div>
""");

# ╔═╡ 4d802163-af7b-4255-a6a9-3b17d52a842f
dashboard_area(
    "network-header",
    network_header
)

# ╔═╡ 4d8d901f-d3c2-4b0a-b2c7-4ed8b2a5280d
cytoscape_graph = let
    graph = CytoscapeJS.Cytoscape(
        network;
        group_colors,
        height="600px",
        include_shared=show_shared_species,
    )
    graph
end;

# ╔═╡ 27655a09-73a0-4370-933b-2e4fe2e2ecf3
CytoscapeJS.set_filter!(cytoscape_graph, :presentIn, selected_model);

# ╔═╡ 5f0c8a71-93de-4b02-a6c4-1e78d3520fb9
gene_selection = PageState.Shared("genes")

# ╔═╡ aa943167-f20e-4349-a14b-d512c8005ab0
network_view = PageState.sync(gene_selection, Bonito.App(cytoscape_graph));

# ╔═╡ 769e4fd2-71c3-418b-8838-0c8658dad724
dashboard_area(
    "network",
    dashboard_panel(network_view; height=600),
)

# ╔═╡ d1907f4c-3b26-4e85-9a70-52c8f31be6d4
gene_selection_bridge = @bind selected_gene_names PageState.bridge(gene_selection)

# ╔═╡ 6b2e94f1-8d05-4c73-b1a8-97f2e6c04a3d
selected_genes = selected_gene_names isa AbstractVector ?
    String.(selected_gene_names) : String[]

# ╔═╡ c5a70e39-4b82-4d16-9f38-e07a2c61b845
visible_genes = isempty(selected_genes) ? string.(network.groups) : selected_genes

# ╔═╡ 2b8e4ad5-6f93-4c71-9d02-af57b318ce03
trajectory_genes_control = PageState.picker(gene_selection, string.(network.groups))

# ╔═╡ a8d4dcb7-33b0-44ef-8736-8cce40bd6eb7
simulation = let
    sink = Trajectories.Sink()
    error = nothing
    calls = Ref(0)
    fractions = Float64[]
    started = time_ns();
    if run_simulation && !isnothing(schedule!)
        try
            ProgressLogging.@withprogress begin
                schedule!(;
                    trace=sink,
                    consolidated_progress = fraction_reporter(total_duration) do fraction
                        calls[] += 1
                        push!(fractions, fraction)
                        ProgressLogging.@logprogress fraction
                    end,
                )
            end
        catch exception
            error = sprint(showerror, exception)
        end
    end

    (; sink, error, calls = calls[], fractions, elapsed = (time_ns() - started) / 1e9)
end

# ╔═╡ fc7d5237-870d-4549-9b95-d6eb7d508203
path_options = let
    labels = Vis.path_labels(simulation.sink.records)
    [
        "" => "all paths"
        (
            path => Vis.describe(path, get(labels, path, ""))
            for path in Vis.paths(simulation.sink.records)
        )...
    ]
end

# ╔═╡ 1f6b30ca-8d47-4e52-9a03-6c8e1d47b259
distribution_path_control = @bind secondary_path PlutoUI.Select(
    [nothing => "none"; path_options];
    default=nothing,
);

# ╔═╡ 1a7f39c4-5d82-4e60-b3a1-8c46f207de92
trajectory_path_control = @bind selected_path PlutoUI.Select(path_options);

# ╔═╡ 7afe83d0-f91c-4d6d-80bb-82d634af59ed
trajectory_display_control = @bind aggregate_mode PlutoUI.Select(
	[
        :raw => "individual branches",
        :aggregate => "branch aggregate",
    ];
    default=:aggregate,
);

# ╔═╡ 3c9f5be6-70a4-4d82-8e13-b0681c429df1
trajectory_tracks_control = @bind selected_tracks PlutoUI.MultiSelect(
    string.(track_options);
    default=default_tracks,
    size=4,
);

# ╔═╡ 5b1e9c47-3a82-4d0f-9e61-7c2f8a4d6b30
trace = Trajectories.catenate(simulation.sink)

# ╔═╡ 24d5d19e-c9fd-4b4d-b90d-b2cc144586ce
levels = Trajectories.lod(trace)

# ╔═╡ 9c4e17ab-2f60-4d83-b915-6e0a7d3c81f4
trajectories = Trajectories.select(levels, selected_path)

# ╔═╡ 68f213e4-a64b-4aa1-bf77-9b131e657193
trajectory_view = if run_simulation && !ismissing(dark_mode)
    Bonito.App() do
        with_theme(dark_mode ? theme_dark() : Theme()) do
            DOM.div(
                WGLMakie.WithConfig(
                    Trajectories.render(
                        trajectories;
                        tracks=selected_tracks,
                        selected_genes=visible_genes,
                        group_colors,
                        aggregate_mode,
                    );
                    resize_to = (:parent, nothing),
                );
                style = "width: 100%; height: 100%;",
            )
        end
    end
end;

# ╔═╡ 2f83b1d6-9c07-4e58-a3f1-64d0b8ea7c19
analysis_control = @bind analysis PlutoUI.Select(
    ["none" => "none", "phase" => "phase", "distribution" => "distribution"];
    default="phase",
);

# ╔═╡ b81a8c99-9653-4288-846a-f56c873698cc
trajectory_header = isnothing(trajectory_view) ? nothing : @htl("""
<div class="dashboard-header stacked">
    <div style="display: flex; flex-direction: column; gap: 0.4rem;">
        <label class="dashboard-option">
            <span>path:</span>
            $(trajectory_path_control)
        </label>

        <label class="dashboard-option">
            <span>display:</span>
            $(trajectory_display_control)
        </label>

        <label class="dashboard-option">
            <span>analysis:</span>
            $(analysis_control)
        </label>
    </div>
    <label class="dashboard-option stacked">
        <span>genes:</span>
        $(trajectory_genes_control)
    </label>
    <label class="dashboard-option stacked">
        <span>tracks:</span>
        $(trajectory_tracks_control)
    </label>
</div>
""");

# ╔═╡ 5e913274-b08c-4366-b7ba-4c28e63b953a
dashboard_area(
    "trajectory",
    trajectory_header,
    dashboard_panel(trajectory_view; height=180 * max(1, length(selected_tracks))),
)

# ╔═╡ 4a1c8e77-2d95-4f3a-b0e6-1c7d9f2a4b58
analysis_track_control = @bind analysis_track PlutoUI.Select(
    string.(track_options);
    default=in("proteins", string.(track_options)) ? "proteins" : first(string.(track_options)),
);

# ╔═╡ 5a2d84f7-1c39-4b68-8e70-9f3a52c1d704
distribution_samples = if run_simulation && analysis == "distribution"
    Trajectories.counts(
        Trajectories.select(trace, selected_path),
        visible_genes,
        analysis_track,
    )
end

# ╔═╡ 8c40e91b-6f25-4d3a-97b1-4e08d5a2f637
distribution_reference =
    if !isnothing(distribution_samples) && !isnothing(secondary_path)
        Trajectories.counts(
            Trajectories.select(trace, secondary_path),
            visible_genes,
            analysis_track,
        )
    end

# ╔═╡ 9b1c47e3-5a82-4f07-b3e6-08d27f4a5169
distribution_view =
    if !isnothing(distribution_samples) && !ismissing(dark_mode)
        Bonito.App() do
            with_theme(dark_mode ? theme_dark() : Theme()) do
                figure = isnothing(distribution_reference) ?
                    Trajectories.render(distribution_samples; group_colors) :
                    Trajectories.render(
                        distribution_samples,
                        distribution_reference;
                        group_colors,
                    )
                DOM.div(
                    WGLMakie.WithConfig(figure; resize_to = (:parent, nothing));
                    style = "width: 100%; height: 100%;",
                )
            end
        end
    end

# ╔═╡ 2e79f5a8-b013-4c96-85d2-7a1fc6039e48
distribution_header = isnothing(distribution_samples) ? nothing : @htl("""
<div class="dashboard-header">
    <label class="dashboard-option">
        <span>track:</span>
        $(analysis_track_control)
    </label>

    <label class="dashboard-option">
        <span>compare:</span>
        $(distribution_path_control)
    </label>

    <span style="opacity: 0.5; font-size: 0.75rem;">
        $(size(distribution_samples.X, 2)) samples
    </span>
</div>
""");

# ╔═╡ be5ca395-d977-4082-a700-60730ba278df
phase_components_control = @bind phase_components PlutoUI.Select(
    [2 => "2D", 3 => "3D"];
    default=2,
);

# ╔═╡ 6d1f40a8-5b72-4c93-8e14-2a07c9f3bd65
phase_coloring_control = @bind phase_coloring PlutoUI.Select(
    [:genes => "genes", :time => "time"];
    default=:genes,
);

# ╔═╡ 9e8abb2f-4531-4c57-a15f-547938cbd6e4
phase_snapshot = if run_simulation && analysis == "phase"
    Trajectories.snapshots(
        Trajectories.select(trace, selected_path),
        visible_genes,
        analysis_track,
    )
end

# ╔═╡ ffa07386-4884-4f9d-9e2d-29da604aec87
phase_projection = if !isnothing(phase_snapshot)
    Trajectories.project(
        phase_snapshot;
        components=phase_components,
        group_colors,
        coloring=phase_coloring,
    )
end

# ╔═╡ 7f2b6c14-8a03-4d51-9e27-3b5c0a8f61d9
phase_header = !(run_simulation && analysis == "phase") ? nothing : @htl("""
<div class="dashboard-header">
    <label class="dashboard-option">
        <span>track:</span>
        $(analysis_track_control)
    </label>

    <label class="dashboard-option">
        $(phase_components_control)
    </label>

    <label class="dashboard-option">
        <span>colour:</span>
        $(phase_coloring_control)
    </label>

    $(isnothing(phase_projection) ? nothing : @htl("""
    <span style="opacity: 0.5; font-size: 0.75rem;">
        $(string(phase_projection.method))
    </span>
    """))
</div>
""");

# ╔═╡ b8a9bd86-5abc-48ed-87f5-534fe604bfb3
phase_view = if !isnothing(phase_projection) && !ismissing(dark_mode)
    Bonito.App() do
        with_theme(dark_mode ? theme_dark() : Theme()) do
            DOM.div(
                WGLMakie.WithConfig(
                    Trajectories.render(phase_projection);
                    resize_to = (:parent, nothing),
                );
                style = "width: 100%; height: 100%;",
            )
        end
    end
end

# ╔═╡ 3c7f1052-9e6a-4144-b5f8-2a06c419731e
dashboard_area(
    "left",
    schedule_panel,
    phase_header,
    dashboard_panel(phase_view; height=450),
    distribution_header,
    dashboard_panel(distribution_view; height=450),
)

# ╔═╡ a4bcea34-2ee9-48ee-bf1a-6d1fd6256c91
begin
    docs_url = "http://localhost:8001"
    app_header = @htl("""
<header class="grs-header">
    <div class="logo-wrapper">
        <a href="/" title="Pluto home">
            $(PlutoUI.LocalResource(joinpath(@__DIR__, "assets", "logo.png")))
        </a>
    </div>

    <div>
        <div class="grs-title">Gene Regulatory Systems</div>
    </div>
    <span style="opacity:0.6; font-size:0.85em">
        $(run_simulation ? "ran in $(round(simulation.elapsed; digits=2))s" : "")
    </span>

    <a class="docs-link" href="$(docs_url)" target="_blank">
        docs <span class="external-link">↗</span>
    </a>
    $(reload_control)
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
    .logo-wrapper a {
        display: block;
        line-height: 0;
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

    .logo-wrapper:hover img {
        animation: grs-logo-spin 0.6s ease-in-out;
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
""");
end;

# ╔═╡ 2b6e0f41-8d59-4033-a4e7-19f5b308620d
dashboard_area("header", app_header)

# ╔═╡ Cell order:
# ╟─29eb1509-3cea-4c50-9d8e-24bcd8a54bbb
# ╠═032084a8-991a-4eeb-a5e9-78110b2b7cbe
# ╠═116c3c09-d5b1-4c86-8c08-c470df26099b
# ╠═2b6e0f41-8d59-4033-a4e7-19f5b308620d
# ╠═1f6b30ca-8d47-4e52-9a03-6c8e1d47b259
# ╠═5a2d84f7-1c39-4b68-8e70-9f3a52c1d704
# ╠═8c40e91b-6f25-4d3a-97b1-4e08d5a2f637
# ╠═2e79f5a8-b013-4c96-85d2-7a1fc6039e48
# ╠═9b1c47e3-5a82-4f07-b3e6-08d27f4a5169
# ╠═3c7f1052-9e6a-4144-b5f8-2a06c419731e
# ╠═b7e4c1a0-5d92-4f3e-8a61-2c9f0d7b4e83
# ╠═c8f5d2b1-6ea3-4045-9b72-3da10e8c5f94
# ╠═4d802163-af7b-4255-a6a9-3b17d52a842f
# ╠═769e4fd2-71c3-418b-8838-0c8658dad724
# ╠═5e913274-b08c-4366-b7ba-4c28e63b953a
# ╠═8eb4f0d7-ef76-4ca7-a117-a48685c12667
# ╟─34a8deb2-0db9-456c-ba38-eabc3254ab50
# ╠═40d24905-5708-4580-9ae5-71e1c1de01a9
# ╠═12e75652-a3fd-11f1-9335-a986ee44237f
# ╟─76de62b9-be88-4b96-b801-b40872dfa0be
# ╠═c81d5f60-2a47-4e93-b6d1-5f39c284ae71
# ╠═a2866674-5e4e-4738-b18e-90d122d2d161
# ╠═d4a916c8-7b52-49f6-a8e3-1c06d5b74f29
# ╠═e57c02d3-9f18-4a6b-b30c-6d2e81f45c93
# ╠═b4f7c2d9-1e58-4a63-9c07-8d5e21f4a3b6
# ╠═7abbbe81-3d0b-4d94-a69f-32c66514687f
# ╠═8392f056-369c-45fb-9c00-719172e82616
# ╠═cdbda6eb-b5b9-4e9e-bac0-f334c19ae28e
# ╠═50fb4e60-7911-4f68-bd4a-166be8d4d44b
# ╠═69993fb8-daca-4b61-a711-907009dbee24
# ╠═e2b7f4a3-8c15-4d67-9e02-5fa31c7b8d46
# ╠═ea610eb8-caa3-4639-8d25-58904e600e52
# ╟─16d523e4-ebec-4d6e-b5e1-7aaf949ec5fc
# ╠═8bff584c-2987-4584-aaed-b1ce4104d891
# ╠═5e2c8a17-4b39-4f6d-9c81-7d0a3b5e62f4
# ╠═6c1e4b90-2a77-4d3e-9f58-31b0c7ea52d4
# ╠═d9a6e3c2-7fb4-4156-ac83-4eb21f9d6a05
# ╠═0d5a9f31-8e46-4c02-b7d1-9a2f6c48e713
# ╠═4d8d901f-d3c2-4b0a-b2c7-4ed8b2a5280d
# ╠═27655a09-73a0-4370-933b-2e4fe2e2ecf3
# ╠═aa943167-f20e-4349-a14b-d512c8005ab0
# ╟─8d1f0c53-6a24-4e91-b7d5-2c9a04f61b38
# ╠═5f0c8a71-93de-4b02-a6c4-1e78d3520fb9
# ╠═d1907f4c-3b26-4e85-9a70-52c8f31be6d4
# ╠═6b2e94f1-8d05-4c73-b1a8-97f2e6c04a3d
# ╠═c5a70e39-4b82-4d16-9f38-e07a2c61b845
# ╠═2b8e4ad5-6f93-4c71-9d02-af57b318ce03
# ╟─c7cd57d6-464d-4dbc-b955-699802e60be7
# ╠═a1d5e8f2-4b73-4c96-8e05-71fc32a9b6d8
# ╠═a8d4dcb7-33b0-44ef-8736-8cce40bd6eb7
# ╟─525c8855-9758-4d61-a8da-a616396be848
# ╠═180bf479-2269-4d94-a9f6-b3c060c17ab9
# ╠═5f81c7aa-273c-49e2-ae61-f473f810d6c3
# ╠═fc7d5237-870d-4549-9b95-d6eb7d508203
# ╠═1a7f39c4-5d82-4e60-b3a1-8c46f207de92
# ╠═7afe83d0-f91c-4d6d-80bb-82d634af59ed
# ╠═3c9f5be6-70a4-4d82-8e13-b0681c429df1
# ╠═5b1e9c47-3a82-4d0f-9e61-7c2f8a4d6b30
# ╠═24d5d19e-c9fd-4b4d-b90d-b2cc144586ce
# ╠═9c4e17ab-2f60-4d83-b915-6e0a7d3c81f4
# ╠═68f213e4-a64b-4aa1-bf77-9b131e657193
# ╠═b81a8c99-9653-4288-846a-f56c873698cc
# ╟─76d9a00b-10fa-4248-b460-6c5abb2424dd
# ╠═2f83b1d6-9c07-4e58-a3f1-64d0b8ea7c19
# ╠═4a1c8e77-2d95-4f3a-b0e6-1c7d9f2a4b58
# ╠═be5ca395-d977-4082-a700-60730ba278df
# ╠═6d1f40a8-5b72-4c93-8e14-2a07c9f3bd65
# ╠═9e8abb2f-4531-4c57-a15f-547938cbd6e4
# ╠═ffa07386-4884-4f9d-9e2d-29da604aec87
# ╠═b8a9bd86-5abc-48ed-87f5-534fe604bfb3
# ╠═7f2b6c14-8a03-4d51-9e27-3b5c0a8f61d9
# ╟─35e4e15a-e1f7-42c1-8ded-5714234811c6
# ╠═4a29ccba-5cdc-47e1-953c-5eee6d1edbe7
# ╠═6a3d1f28-04bc-4e7a-9c15-8b2f70d6ae41
# ╠═a4bcea34-2ee9-48ee-bf1a-6d1fd6256c91
# ╠═375349c2-33bb-45fb-b412-672c678b93b3
# ╠═b3f27a51-6c84-4d19-9f52-8e1c47a0d6b2
# ╠═0187c445-f7a0-4eae-b55d-2bd955df2c2c
