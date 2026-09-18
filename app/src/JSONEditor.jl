
using HypertextLiteral: @htl

struct JSONEditor
    contents::String
    height::String
    on_save::Bool
end

Base.get(editor::JSONEditor) = editor.contents

JSONEditor(contents; height="400px", on_save=true) = JSONEditor(string(contents), string(height), on_save)

Base.show(io::IO, ::MIME"text/html", editor::JSONEditor) =
    show(io, MIME"text/html"(), @htl("""
<div class="json-editor" style="--editor-height: $(editor.height)"
    data-on-save="$(editor.on_save)">
    <textarea hidden>$(editor.contents)</textarea>
    <div class="editor"></div>

    <style>
        .json-editor,
        .json-editor .editor,
        .json-editor .cm-editor {
            width: 100%;
            min-width: 0;
            box-sizing: border-box;
        }

        .json-editor .cm-editor {
            height: var(--editor-height);
            border: 1px solid #d4d4d8;
            border-radius: 6px;
            overflow: hidden;
            font-size: 12px;
        }

        .json-editor {
            position: relative;
        }

        .json-editor::after {
            content: "";
            position: absolute;
            top: 8px;
            right: 14px;
            width: 7px;
            height: 7px;
            border-radius: 50%;
            background: currentColor;
            color: #a1a1aa;
            opacity: 0;
            transition: opacity 120ms ease;
            pointer-events: none;
        }

        .json-editor.dirty::after {
            opacity: 0.9;
        }

        .json-editor .color-swatch {
            display: inline-block;
            width: 0.75em;
            height: 0.75em;
            margin-left: 0.35em;
            border-radius: 2px;
            border: 1px solid rgba(128, 128, 128, 0.5);
            vertical-align: -0.05em;
        }

        .json-editor .cm-scroller {
            overflow: auto;
        }

        .json-editor .cm-gutters {
            background-color: #ffffff !important;
        }

        .json-editor .cm-lineNumbers .cm-gutterElement {
            color: #71717a !important;
            opacity: 1 !important;
        }

        @media (prefers-color-scheme: dark) {
            .json-editor .cm-editor {
                border-color: #6b7280;
            }

            .json-editor::after {
                color: #71717a;
            }

            .json-editor .cm-gutters {
                background-color: #282c34 !important;
            }

            .json-editor .cm-lineNumbers .cm-gutterElement {
                color: #a1a1aa !important;
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
        const COMMIT_ON_SAVE = root.dataset.onSave === "true"
        parent.addEventListener("input", event => event.stopPropagation())
        const initialValue = root.querySelector("textarea").value
        root.value = initialValue

        const { basicSetup, EditorView } =
            await import("https://esm.sh/codemirror@6.0.2")
        const { indentWithTab } =
            await import("https://esm.sh/@codemirror/commands@^6.0.0?target=es2022")
        const { keymap, Decoration, ViewPlugin, WidgetType } =
            await import("https://esm.sh/@codemirror/view@^6.0.0?target=es2022")
        const { json } =
            await import("https://esm.sh/@codemirror/lang-json@6.0.2")
        const { oneDark } =
            await import("https://esm.sh/@codemirror/theme-one-dark@6.1.3")

        const colorScheme =
            window.matchMedia("(prefers-color-scheme: dark)")

        const MIN_DELAY = 600
        const BURST_STEP = 150
        const MAX_DELAY = 1600
        const INVALID_DELAY = 1500
        const BURST_WINDOW = 1200
        const MAX_WAIT = 15000

        let timeout
        let view
        let burst = 0
        let lastChange = 0
        let pending = 0

        let contents = initialValue
        let dirty = false

        const setDirty = value => {
            dirty = value
            root.classList.toggle("dirty", value)
        }

        const commit = () => {
            pending = 0
            burst = 0
            setDirty(false)
            root.value = contents
            root.dispatchEvent(new CustomEvent("input"))
        }

        class ColorSwatch extends WidgetType {
            constructor(color) { super(); this.color = color }
            eq(other) { return other.color === this.color }
            toDOM() {
                const box = document.createElement("span")
                box.className = "color-swatch"
                box.style.backgroundColor = this.color
                return box
            }
        }

        const SWATCH_PATTERN =
            /#(?:[0-9a-fA-F]{8}|[0-9a-fA-F]{6}|[0-9a-fA-F]{3})/g

        const buildSwatches = view => {
            const widgets = []
            for (const { from, to } of view.visibleRanges) {
                const text = view.state.doc.sliceString(from, to)
                SWATCH_PATTERN.lastIndex = 0
                let match
                while ((match = SWATCH_PATTERN.exec(text)) !== null) {
                    widgets.push(
                        Decoration.widget({
                            widget: new ColorSwatch(match[0]),
                            side: 1,
                        }).range(from + match.index + match[0].length),
                    )
                }
            }
            return Decoration.set(widgets, true)
        }

        const colorSwatches = ViewPlugin.fromClass(
            class {
                constructor(view) {
                    this.decorations = buildSwatches(view)
                }
                update(update) {
                    if (update.docChanged || update.viewportChanged) {
                        this.decorations = buildSwatches(update.view)
                    }
                }
            },
            { decorations: plugin => plugin.decorations },
        )

        const createEditor = doc => new EditorView({
            doc,
            parent,
            extensions: [
                keymap.of([{
                    key: "Mod-s",
                    preventDefault: true,
                    run: () => { dirty && commit(); return true },
                }]),
                basicSetup,
                keymap.of([indentWithTab]),
                json(),
                colorSwatches,
                colorScheme.matches ? oneDark : [],
                EditorView.updateListener.of(update => {
                    if (!update.docChanged) return

                    const now = Date.now()
                    burst = now - lastChange < BURST_WINDOW ? burst + 1 : 1
                    lastChange = now
                    pending || (pending = now)

                    contents = update.state.doc.toString()
                    let valid = true
                    try {
                        JSON.parse(contents)
                    } catch {
                        valid = false
                    }

                    if (COMMIT_ON_SAVE) {
                        setDirty(true)
                        return
                    }

                    clearTimeout(timeout)

                    if (valid && now - pending >= MAX_WAIT) {
                        commit()
                        return
                    }

                    timeout = setTimeout(
                        commit,
                        valid ?
                            Math.min(
                                MIN_DELAY + BURST_STEP * (burst - 1),
                                MAX_DELAY,
                            ) :
                            INVALID_DELAY,
                    )
                }),
            ],
        })

        const updateTheme = () => {
            const current = view.state.doc.toString()
            view.destroy()
            view = createEditor(current)
        }

        view = createEditor(initialValue)

        colorScheme.addEventListener("change", updateTheme)

        invalidation.then(() => {
            clearTimeout(timeout)
            colorScheme.removeEventListener("change", updateTheme)
            view.destroy()
        })
    </script>
</div>
"""))
