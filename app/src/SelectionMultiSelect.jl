using Bonito: DOM, Session, @js_str, onload

struct SelectionMultiSelect{O}
    selection::O
    options::Vector{String}
    size::Int
end

function SelectionMultiSelect(selection, options; size=nothing)
    options = string.(collect(options))
    SelectionMultiSelect(selection, options, something(size, min(4, length(options))))
end

Base.get(select::SelectionMultiSelect) = copy(select.selection[])

function Base.show(io::IO, m::MIME"text/html", widget::SelectionMultiSelect)
    app = Bonito.App() do session
        return widget
    end
    show(io, m, app)
end

function Bonito.jsrender(session::Session, widget::SelectionMultiSelect)
    selection = widget.selection
    selected = Set(selection[])
    select = DOM.select(
        (
            DOM.option(option; value=option, selected=option in selected)
            for option in widget.options
        )...;
        multiple=true,
        size=widget.size,
        title="Cmd+Click or Ctrl+Click to select multiple items.",
    )

    onload(session, select, js"""
        select => {
            const selection = $(selection);
            if (!selection) return;

            let syncing = false;
            const values = () =>
                [...select.selectedOptions].map(option => option.value);
            const applySelection = values => {
                const selected = new Set(values);
                [...select.options].forEach(option => {
                    option.selected = selected.has(option.value);
                });
            };

            selection.on(values => {
                if (!select.isConnected) return false;
                syncing = true;
                applySelection(values);
                select.dispatchEvent(new Event("input", { bubbles: true }));
                syncing = false;
            });
            select.addEventListener("input", () => {
                if (!syncing) selection.notify(values());
            });
        }
    """)

    select
end
