module CytoscapeJSExt

using GeneRegulatorySystems
using CytoscapeJS
using Colors

const Vis = GeneRegulatorySystems.Visualisation
const JS = CytoscapeJS.Bonito

const GENE_WIDTH = 144
const GENE_HEIGHT = 80
const GENE_PADDING = 6
const SPECIES_SIZE = 8
const GENE_CHILD_WIDTH = GENE_WIDTH - SPECIES_SIZE - 2 * GENE_PADDING
const GENE_CHILD_HEIGHT = GENE_HEIGHT - SPECIES_SIZE - 2 * GENE_PADDING

const DEFAULT_LAYOUT = (;
    name="fcose",
    quality="proof",
    randomize=true,
    animate=true,
    animationDuration=600,
    fit=true,
    padding=50,
    nodeDimensionsIncludeLabels=true,
    uniformNodeDimensions=false,
    packComponents=true,
    nodeRepulsion=50_000,
    idealEdgeLength=130,
    edgeElasticity=0.45,
    nestingFactor=0.1,
    gravity=32.8,
    numIter=1_000,
    tile=true,
    tilingPaddingVertical=30,
    tilingPaddingHorizontal=30,
    gravityRangeCompound=3.5,
    gravityCompound=1.0,
    gravityRange=3.8,
    initialEnergyOnIncremental=1,
)

const DEFAULT_TOOLTIP_ATTRIBUTES = (;
    class="cytoscapejs-tooltip",
    style="""
        max-width: 280px;
        padding: 8px 12px;
        color: #1f2937;
        background: rgba(255, 255, 255, 0.96);
        border: 1px solid #d1d5db;
        border-radius: 3px;
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.14);
        font-family: Montserrat, sans-serif;
        font-size: 10px;
        font-weight: 500;
        line-height: 1.4;
        white-space: pre-line;
    """
)

const DEFAULT_ATTRIBUTES = (;
    style="""
        width: 100%;
        height: calc(100vh);
        min-height: 560px;
        overflow: hidden;
        background-color: #ffffff;
        background-image: radial-gradient(circle, #d4d4d8 1px, transparent 1px);
        background-size: 30px 30px;
    """
)

function CytoscapeJS.Cytoscape(
    network::Vis.Network;
    layout=DEFAULT_LAYOUT,
    group_colors=Vis.group_colors(network.groups),
    stylesheet=stylesheet(network, group_colors),
    tooltip_attributes = DEFAULT_TOOLTIP_ATTRIBUTES,
    attributes = DEFAULT_ATTRIBUTES,
    wheelSensitivity=0.1,
    kwargs...
)
    strength_reference = get_strength_reference(network)
    gene_elements = elements(Vis.gene_view(network), network, Val(:gene), group_colors, strength_reference)
    species_elements = elements(Vis.species_view(network), network, Val(:species), group_colors, strength_reference)

    CytoscapeJS.Cytoscape(
        gene_elements;
        layout,
        stylesheet,
        setup=network_setup(gene_elements, species_elements),
        wheelSensitivity,
        tooltip_attributes,
        attributes,
        kwargs...
    )
end

function elements(view::Vis.Network, network::Vis.Network, mode::Val, group_colors, strength_reference)
    vcat(
        node_element.(view.nodes, Ref(network), mode, Ref(group_colors)),
        [link_element(link, network, mode, strength_reference) for link in view.links]
    )
end

css_color(color::Color) = "#$(hex(color))"
css_color(color::AbstractString) = color

node_color(node, group_colors) = css_color(
    get(group_colors, something(node.parent, node.name), "#808080")
)

function node_element(node::Vis.Node, network::Vis.Network, ::Val{:gene}, group_colors)
    (; data=(;
        id=string(node.name),
        label=Vis.node_label(node),
        kind=string(node.kind),
        colour=node_color(node, group_colors),
        tooltip=Vis.node_tooltip(node, network),
        parameters=Vis.parameter_lines(node, network),
    ), classes = string(node.kind))
end

function node_element(node::Vis.Node, network::Vis.Network, ::Val{:species}, group_colors)
    data=(;
        id=string(node.name),
        label=Vis.node_label(node),
        kind=string(node.kind),
        colour=node_color(node, group_colors),
        tooltip=Vis.node_tooltip(node, network),
        parameters=Vis.parameter_lines(node, network)
    )
    node.parent === nothing ||
        (data = merge(data, (; parent=string(node.parent))))
    orphan = node.kind === :species && node.parent === nothing ? " orphan-species" : ""
    (; data, classes="$(node.kind)$orphan")
end

function link_element(link::Vis.Link, network::Vis.Network, ::Val{V}, strength_reference) where V
    loop = link.from === link.to ? " loop" : ""
    data=(;
        id="$V:$(link.kind):$(link.from):$(link.to)",
        source=string(link.from),
        target=string(link.to),
        kind=string(link.kind),
        scope=string(link.scope),
        view=string(V),
        label=Vis.link_label(link),
        tooltip=Vis.link_tooltip(link, network),
    )
    if haskey(link.properties, :at)
        at = Float64(get(link.properties, :at, 1))
        data = merge(data, (; at, strengthNorm=inv(1 + sqrt(at / strength_reference))))
    end
    (; data, classes="$(link.kind)$loop")
end

function stylesheet(network, group_colors; fontfamily="Montserrat")
    rules = [
        # nodes
        (; selector="node", style=Dict(
            "label" => "data(label)",
            "font-family" => fontfamily,
            "text-halign" => "center",
            "text-valign" => "center",
            "background-color" => "data(colour)",
            "border-width" => 0
        )),
        (; selector="node.gene", style=Dict(
            "shape" => "round-rectangle",
            "width" => GENE_WIDTH,
            "height" => GENE_HEIGHT,
            "min-width" => GENE_WIDTH,
            "min-height" => GENE_HEIGHT,
            "font-size" => 20,
            "padding" => GENE_PADDING
        )),
        (; selector="node.compound-parent", style=Dict(
            "text-valign" => "top",
            "text-margin-y" => -8,
            "background-opacity" => 0.2,
            "text-wrap" => "none",
            "compound-sizing-wrt-labels" => "exclude"
        )),
        (; selector="node.species", style=Dict(
            "shape" => "ellipse",
            "width" => SPECIES_SIZE,
            "height" => SPECIES_SIZE,
            "font-size" => 2.4,
            "text-valign" => "bottom",
            "text-margin-y" => 1
        )),
        (; selector="node.reaction", style=Dict(
            "shape" => "ellipse",
            "width" => 12,
            "height" => 9,
            "font-size" => 2.0,
            "text-margin-y" => -2,
            "background-opacity" => 0.0,
        )),
        (; selector="node.dimmed", style=Dict(
            "opacity" => 0.3,
        )),
        (; selector="node.highlighted", style=Dict(
            "z-index" => 10,
        )),
        # edges
        (; selector="edge", style=Dict(
            "label" => "data(label)",
            "width" => 1,
            "curve-style" => "bezier",
            "target-arrow-shape" => "triangle",
            "font-family" => fontfamily,
            "font-size" => 7,
            "text-rotation" => "autorotate",
            "text-margin-y" => -8,
            "text-background-color" => "#ffffff",
            "text-background-opacity" => 0.7,
            "text-background-padding" => 2
        )),
        (; selector="edge[view =\"species\"]", style=Dict(
            "font-size" => 3,
            "arrow-scale" => 0.8,
            "z-index" => -100
        )),
        (; selector="edge.repression", style=Dict(
            "target-arrow-shape" => "tee"
        )),
        (; selector="edge.proteolysis", style=Dict(
            "width" => 1.5,
            "target-arrow-shape" => "diamond"
        )),
        (; selector="edge.catalyses", style=Dict(
            "width" => 1.2,
            "line-style" => "dashed",
            "line-dash-pattern" => [4,2],
            "target-arrow-shape" => "circle",
            "arrow-scale" => 0.7,
            "opacity" => 0.7
        )),
        (; selector="edge.substrate, edge.product", style=Dict(
            "width" => 0.5,
            "font-size" => 1.3,
            "arrow-scale" => 0.2,
            "curve-style" => "unbundled-bezier",
            "control-point-step-size" => 4,
            "text-margin-y" => -1,
            "text-background-opacity" => 0
        )),
        (; selector="edge.loop", style=Dict(
            "curve-style" => "unbundled-bezier",
            "control-point-step-size" => 100,
            "loop-sweep" => "60deg"
        )),
        (; selector="edge.dimmed", style=Dict(
            "opacity" => 0.3,
        )),
        (; selector="edge.activation[strengthNorm], edge.repression[strengthNorm]", style=Dict(
            "width" => "mapData(strengthNorm, 0, 1, 1, 5)",
        )),

        (; selector="""
            edge[view = "species"][strengthNorm]
            """,
        style=Dict(
            "width" => "mapData(strengthNorm, 0, 1, 0.75, 2.5)",
        )),
        (; selector="edge.promotes, edge.inhibits", style=Dict(
            "width" => 0.5,
            "font-size" => 1.5,
            "arrow-scale" => 0.2,
            "text-margin-y" => -1,
            "text-background-opacity" => 0.7,
            "text-background-padding" => 0.0,
        )),
        (; selector="edge.inhibits", style=Dict(
            "target-arrow-shape" => "tee",
        )),
    ]
    for (kind, color) in Vis.LINK_COLORS
        push!(rules, (;
            selector="edge.$kind",
            style=Dict(
                "line-color" => color,
                "target-arrow-color" => color
            )
        ))
    end
    rules
end

function selection_view()
    JS.js"""
    cy => {
        const selected = new Set();
        let pendingFrame = null;

        function selectionKey(node) {
            if (!node || node.empty()) return null;
            if (node.hasClass("gene") || node.hasClass("orphan-species")) {
                return node.id();
            }
            return node.data("parent") ?? null;
        }

        function updateSelection() {
            pendingFrame = null;

            cy.batch(() => {
                if (selected.size === 0) {
                    cy.elements().removeClass("dimmed highlighted");
                    return;
                }

                cy.nodes().forEach(node => {
                    const key = selectionKey(node);
                    const isSelected = key !== null && selected.has(key);
                    const isHighlightable =
                        node.hasClass("gene") || node.hasClass("orphan-species");

                    node.toggleClass("highlighted", isSelected && isHighlightable);
                    node.toggleClass("dimmed", !isSelected);
                });

                cy.edges().forEach(edge => {
                    const sourceSelected = selected.has(selectionKey(edge.source()));
                    const targetSelected = selected.has(selectionKey(edge.target()));
                    const kind = edge.data("kind");
                    const reactionEdge = kind === "substrate" || kind === "product";
                    const connected = reactionEdge
                        ? sourceSelected || targetSelected
                        : sourceSelected && targetSelected;

                    edge.toggleClass("dimmed", !connected);
                });
            });
        }

        function scheduleUpdate() {
            if (pendingFrame !== null) return;
            pendingFrame = requestAnimationFrame(updateSelection);
        }

        cy.on("tap", "node.gene, node.orphan-species", event => {
            const id = event.target.id();
            const original = event.originalEvent;
            const toggle = original?.ctrlKey === true || original?.metaKey === true;

            if (toggle) {
                selected.has(id) ? selected.delete(id) : selected.add(id);
            } else if (selected.size === 1 && selected.has(id)) {
                selected.clear();
            } else {
                selected.clear();
                selected.add(id);
            }

            updateSelection();
        });

        cy.on("add remove", scheduleUpdate);
        cy.on("destroy", () => {
            if (pendingFrame !== null) cancelAnimationFrame(pendingFrame);
            selected.clear();
        });
    }
    """
end

function adaptive_view(gene_elements, species_elements)
JS.js"""
cy => {
    const threshold = 2.0;
    const geneElements = $(gene_elements);
    const speciesElements = $(species_elements);
    const geneEdges = geneElements.filter(
        element => element.data.source !== undefined
    );

    let detailVisible = false;
    let timeout = null;
    const positions = new Map();

    function recenterGene(gene, target) {
        const current = gene.position();
        const offset = {
            x: target.x - current.x,
            y: target.y - current.y,
        };

        gene.children().shift(offset);
    }

    function fitGeneChildren(gene) {
        const children = gene.children();
        if (children.length < 2) return;

        const xs = children.map(node => node.position("x"));
        const ys = children.map(node => node.position("y"));
        const minX = Math.min(...xs);
        const maxX = Math.max(...xs);
        const minY = Math.min(...ys);
        const maxY = Math.max(...ys);
        const centreX = (minX + maxX) / 2;
        const centreY = (minY + maxY) / 2;
        const scaleX = maxX === minX ? 0 : $(GENE_CHILD_WIDTH) / (maxX - minX);
        const scaleY = maxY === minY ? 0 : $(GENE_CHILD_HEIGHT) / (maxY - minY);

        children.positions(node => {
            const position = node.position();
            return {
                x: centreX + (position.x - centreX) * scaleX,
                y: centreY + (position.y - centreY) * scaleY,
            };
        });
    }

    function layoutDetail(genePositions) {
        cy.nodes(".gene").forEach(gene => {
            const children = gene.children();
            if (children.empty()) return;

            const target = genePositions.get(gene.id());

            // Restored children already have a layout.
            const cached = children.every(node => positions.has(node.id()));
            if (cached) {
                if (target) recenterGene(gene, target);
                return;
            }

            const childIds = new Set(children.map(node => node.id()));
            const internalEdges = children.connectedEdges().filter(edge =>
                childIds.has(edge.source().id()) &&
                childIds.has(edge.target().id())
            );

            const layout = children.merge(internalEdges).layout({
                name: "fcose",
                quality: "default",
                randomize: true,
                animate: false,
                fit: false,
                nodeRepulsion: 1000,
                idealEdgeLength: 12,
                edgeElasticity: 0.9,
                numIter: 100,
                gravity: 1.2,
                gravityRange: 1.0,
                tile: false,
                packComponents: false,
            });

            layout.one("layoutstop", () => {
                if (target && !gene.removed() && gene.children().nonempty()) {
                    fitGeneChildren(gene);
                    recenterGene(gene, target);
                }
            });

            layout.run();
        });
    }
    function showDetail() {
        if (detailVisible) return;
        detailVisible = true;
        const genePositions = new Map(
            cy.nodes(".gene").map(gene => [
                gene.id(),
                { ...gene.position() },
            ])
        )
        cy.batch(() => {
            cy.remove(cy.edges());
            for (const element of speciesElements) {
                if (cy.getElementById(element.data.id).empty()) {
                    const added = cy.add(element);
                    const position = positions.get(element.data.id);
                    if (position) added.position(position);
                }
            }
            cy.nodes(".gene").addClass("compound-parent");
        })
        layoutDetail(genePositions);
    }
    function hideDetail() {
        if (!detailVisible) return;
        detailVisible = false;
        cy.nodes(":child").forEach(node => {
            positions.set(node.id(), node.position());
        })
        cy.batch(() => {
            cy.remove(cy.nodes(":child"));
            cy.remove(cy.edges());
            cy.add(geneEdges);
            cy.nodes(".gene").removeClass("compound-parent");
        })
    }
    function update() {
        cy.zoom() > threshold ? showDetail() : hideDetail()
    }
    cy.on("zoom", () => {
        clearTimeout(timeout);
        timeout = setTimeout(update, 50)
    })
    setTimeout(update, 500);
}
"""
end

function inline_parameters()
    JS.js"""
    cy => {
        const host = cy.container();

        if (getComputedStyle(host).position === "static") {
            host.style.position = "relative";
        }

        const anchors = new Map();
        let pendingFrame = null;

        function chipStyle(chip) {
            Object.assign(chip.style, {
                padding: "1px 3px",
                color: "#1f2937",
                background: "rgba(255, 255, 255, 0.7)",
                border: "none",
                borderRadius: "3px",
                boxShadow: "0 1px 2px rgba(0, 0, 0, 0.08)",
                fontFamily: "Montserrat, sans-serif",
                fontWeight: "400",
                lineHeight: "1.3",
                whiteSpace: "nowrap",
            });
        }

        function addAnchor(node) {
            if (!node.isNode() || !node.hasClass("reaction")) return;
            if (anchors.has(node.id())) return;

            const parameters = node.data("parameters") ?? [];
            if (parameters.length === 0) return;

            const container = document.createElement("div");

            Object.assign(container.style, {
                position: "absolute",
                left: "0",
                top: "0",
                display: "flex",
                gap: "2px",
                pointerEvents: "none",
                zIndex: "100",
                transformOrigin: "top center",
            });

            for (const parameter of parameters) {
                const chip = document.createElement("span");
                chip.textContent = parameter;
                chipStyle(chip);
                container.appendChild(chip);
            }

            host.appendChild(container);
            anchors.set(node.id(), { node, container });
            scheduleUpdate();
        }

        function removeAnchor(node) {
            const anchor = anchors.get(node.id());
            if (!anchor) return;

            anchor.container.remove();
            anchors.delete(node.id());
        }

        function updateAnchors() {
            pendingFrame = null;

            for (const { node, container } of anchors.values()) {
                if (node.removed() || !node.visible()) {
                    container.style.display = "none";
                    continue;
                }

                container.style.display = "flex";
                container.style.opacity = node.hasClass("dimmed") ? "0.3" : "1";

                const position = node.renderedPosition();

                container.style.fontSize =
                    `${1.8 * cy.zoom()}px`;

                container.style.transform =
                    `translate3d(${position.x}px, ${position.y - (2 / 3) * cy.zoom()}px, 0) ` +
                    `translate(-50%, 20%)`;
            }
        }

        function scheduleUpdate() {
            if (pendingFrame !== null) return;

            pendingFrame = requestAnimationFrame(updateAnchors);
        }

        cy.nodes().forEach(addAnchor);

        cy.on("add", "node.reaction", event => {
            addAnchor(event.target);
        });

        cy.on("remove", "node.reaction", event => {
            removeAnchor(event.target);
        });

        cy.on("pan zoom resize", scheduleUpdate);
        cy.on("style", "node.reaction", scheduleUpdate);
        cy.on("position", "node", scheduleUpdate);

        cy.on("destroy", () => {
            if (pendingFrame !== null) {
                cancelAnimationFrame(pendingFrame);
            }

            for (const { container } of anchors.values()) {
                container.remove();
            }

            anchors.clear();
        });

        scheduleUpdate();
    }
    """
end

function network_setup(gene_elements, species_elements)
    adaptive = adaptive_view(gene_elements, species_elements)
    selection = selection_view()
    parameters = inline_parameters()

    JS.js"""
    cy => {
        const setupAdaptive = $(adaptive);
        const setupSelection = $(selection);
        const setupParameters = $(parameters);

        setupAdaptive(cy);
        setupSelection(cy);
        setupParameters(cy);
    }
    """
end

function get_strength_reference(network)
    values = [
        Float64(link.properties[:at])
        for link in network.links
        if link.kind in (:activation, :repression) &&
           get(link.properties, :at, 0) > 0
    ]

    isempty(values) ? 1.0 : exp(sum(log, values) / length(values))
end

end
