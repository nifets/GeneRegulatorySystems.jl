using Bonito: App, DOM, @js_str, onload

function value_multiselect(selection, options; size=nothing)
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
                const applySelection = values => {
                    const selected = new Set(values);
                    [...select.options].forEach(option => {
                        option.selected = selected.has(option.value);
                    });
                };

                applySelection($(selection[]));
                $(selection).on(values => {
                    if (!select.isConnected) return false;
                    applySelection(values);
                });
                select.addEventListener("input", () => {
                    $(selection).notify(
                        [...select.selectedOptions].map(option => option.value)
                    );
                });
            }
        """)

        select
    end
end
