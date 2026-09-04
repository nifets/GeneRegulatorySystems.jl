
using HypertextLiteral: @htl

struct JSONEditor
    contents::String
    height::String
end

Base.get(editor::JSONEditor) = editor.contents

JSONEditor(contents; height="400px") = JSONEditor(string(contents), string(height))

Base.show(io::IO, ::MIME"text/html", editor::JSONEditor) =
    show(io, MIME"text/html"(), @htl("""
<div class="json-editor" style="--editor-height: $(editor.height)">
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
        parent.addEventListener("input", event => event.stopPropagation())
        const initialValue = root.querySelector("textarea").value
        root.value = initialValue

        const { basicSetup, EditorView } =
            await import("https://esm.sh/codemirror@6.0.2")
        const { indentWithTab } =
            await import("https://esm.sh/@codemirror/commands@^6.0.0?target=es2022")
        const { keymap } =
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

        const createEditor = contents => new EditorView({
            doc: contents,
            parent,
            extensions: [
                basicSetup,
                keymap.of([indentWithTab]),
                json(),
                colorScheme.matches ? oneDark : [],
                EditorView.updateListener.of(update => {
                    if (!update.docChanged) return

                    const now = Date.now()
                    burst = now - lastChange < BURST_WINDOW ? burst + 1 : 1
                    lastChange = now
                    pending || (pending = now)

                    const contents = update.state.doc.toString()
                    let valid = true
                    try {
                        JSON.parse(contents)
                    } catch {
                        valid = false
                    }

                    const commit = () => {
                        pending = 0
                        burst = 0
                        root.value = contents
                        root.dispatchEvent(new CustomEvent("input"))
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
            const contents = view.state.doc.toString()
            view.destroy()
            view = createEditor(contents)
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
