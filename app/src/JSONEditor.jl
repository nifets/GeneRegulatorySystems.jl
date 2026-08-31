
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
        parent.addEventListener("input", event => event.stopPropagation())
        const initialValue = root.querySelector("textarea").value
        root.value = initialValue

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
                    }, 750)
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
