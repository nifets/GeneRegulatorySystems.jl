using Bonito: App, DOM, @js_str, onload
using HypertextLiteral: @htl

struct ValueMultiSelect
    options::Vector{String}
    default::Vector{String}
    size::Int
end

function ValueMultiSelect(options; default=String[], size=nothing)
    options = string.(collect(options))
    default = string.(collect(default))
    size = something(size, min(4, length(options)))
    ValueMultiSelect(options, default, size)
end

Base.get(select::ValueMultiSelect) = select.default

Base.show(io::IO, ::MIME"text/html", select::ValueMultiSelect) =
    show(io, MIME"text/html"(), @htl("""
    <span>
        <select
            multiple
            size=$(select.size)
            title="Cmd+Click or Ctrl+Click to select multiple items."
        >
            $([
                @htl("""
                <option
                    value=$(option)
                    selected=$(option in select.default)
                >$(option)</option>
                """)
                for option in select.options
            ])
        </select>

        <script>
            const root = currentScript.parentElement
            const select = root.querySelector("select")

            const update = event => {
                event?.stopPropagation()
                root.value = [...select.selectedOptions].map(
                    option => option.value
                )
                root.dispatchEvent(new CustomEvent("input"))
            }

            root.value = [...select.selectedOptions].map(
                option => option.value
            )
            select.addEventListener("input", update)
            invalidation.then(() =>
                select.removeEventListener("input", update)
            )
        </script>
    </span>
    """))

function observable_multiselect(selection, options; size=nothing)
    options = string.(collect(options))
    size = something(size, min(4, length(options)))

    App() do session
        select = DOM.select(
            (DOM.option(option; value=option) for option in options)...;
            multiple=true,
            size,
            title="Cmd+Click or Ctrl+Click to select multiple items.",
        )

        onload(session, select, js"""
            select => {
                if (!select.isConnected) return;
                const selection = $(selection);
                if (!selection) return;
                const applySelection = values => {
                    const selected = new Set(values);
                    [...select.options].forEach(option => {
                        option.selected = selected.has(option.value);
                    });
                };

                applySelection($(selection[]));
                selection.on(values => {
                    if (!select.isConnected) return false;
                    applySelection(values);
                });
                select.addEventListener("input", () => {
                    selection.notify(
                        [...select.selectedOptions].map(option => option.value)
                    );
                });
            }
        """)

        select
    end
end
