module PageState

using HypertextLiteral

struct Shared
    key::String
end

const RUNTIME = @htl("""
<span>
<script>
    if (!window.PageState) {
        const values = new Map()
        const listeners = new Map()

        const same = (a, b) =>
            a.length === b.length && a.every((value, i) => value === b[i])

        const store = {
            get: key => values.get(key) ?? [],

            set(key, next, source) {
                const previous = store.get(key)
                const incoming = new Set(next)
                const ordered = [
                    ...previous.filter(value => incoming.has(value)),
                    ...next.filter(value => !previous.includes(value)),
                ]
                if (same(previous, ordered)) return

                values.set(key, ordered)
                for (const listener of listeners.get(key) ?? []) {
                    if (listener.source !== source) listener.fn(ordered)
                }
            },

            connect(key, source) {
                const entries = new Set()
                return {
                    get: () => store.get(key),
                    set: next => store.set(key, next, source),
                    changed(fn) {
                        const entry = { source, fn }
                        if (!listeners.has(key)) listeners.set(key, new Set())
                        listeners.get(key).add(entry)
                        entries.add(entry)
                        fn(store.get(key))
                    },
                    dispose() {
                        for (const entry of entries) listeners.get(key).delete(entry)
                    },
                }
            },
        }

        window.PageState = store
    }
</script>
</span>
""")

bridge(shared::Shared) = @htl("""
<span>
$(RUNTIME)
<script>
    const root = currentScript.parentElement
    const state = window.PageState.connect($(shared.key), "pluto")
    let timer
    let sent = null

    const same = (a, b) =>
        a !== null && a.length === b.length && a.every((value, i) => value === b[i])

    const busy = () =>
        document.querySelector("pluto-cell.running, pluto-cell.queued") !== null

    const emit = () => {
        if (busy()) {
            timer = setTimeout(emit, 200)
            return
        }
        if (same(sent, state.get())) return

        sent = [...state.get()]
        root.value = sent
        root.dispatchEvent(new CustomEvent("input"))
    }

    state.changed(() => {
        clearTimeout(timer)
        timer = setTimeout(emit, 200)
    })

    sent = [...state.get()]
    root.value = sent

    invalidation.then(() => {
        clearTimeout(timer)
        state.dispose()
    })
</script>
</span>
""")

picker(shared::Shared, options; size=4) = @htl("""
<span>
$(RUNTIME)
<select multiple size=$(size) title="Cmd+Click or Ctrl+Click to select multiple items.">
    $((@htl("<option value=$(option)>$(option)</option>") for option in options))
</select>
<script>
    const select = currentScript.parentElement.querySelector("select")
    const state = window.PageState.connect($(shared.key), "picker")

    select.addEventListener("input", () =>
        state.set([...select.selectedOptions].map(option => option.value))
    )

    state.changed(values => {
        const wanted = new Set(values)
        for (const option of select.options) option.selected = wanted.has(option.value)
    })

    invalidation.then(state.dispose)
</script>
</span>
""")

sync(shared::Shared, view) = @htl("""
<span>
$(RUNTIME)
$(view)
<script>
    const root = currentScript.parentElement
    const state = window.PageState.connect($(shared.key), "cytoscape")

    const api = () =>
        [...root.querySelectorAll("div")]
            .find(node => node.cytoscape?.selection)?.cytoscape.selection

    const adopt = selection => selection.set(state.get())

    root.addEventListener("cytoscape:ready", event => adopt(event.detail.selection))
    root.addEventListener("cytoscape:selection", event =>
        state.set(event.detail.selection)
    )

    state.changed(values => api()?.set(values))

    const existing = api()
    if (existing) adopt(existing)

    invalidation.then(state.dispose)
</script>
</span>
""")

end
