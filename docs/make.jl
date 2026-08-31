import ArgParse
using Documenter
import Markdown
using GeneRegulatorySystems

include_tool_script(m, name) =
    Base.include(m, "$(@__DIR__)/../tools/$name/script.jl")

function tool_help(settings)
    result = IOBuffer()
    ArgParse.show_help(result, settings, exit_when_done = false)
    Markdown.MD(Markdown.Code("", String(take!(result))))
end

makedocs(
    sitename = "GeneRegulatorySystems.jl",
    checkdocs = :exports,
    pages = [
        "User guide" => [
            "Overview" => "index.md"
            "guides/getting-started.md"
            #= "guides/v1-models.md" =#
            #= "guides/writing-schedules.md" =#
            "guides/usage-as-library.md"
        ]
        "Reference" => [
            "models/models.md"
            "models/scheduling.md"
            "models/plumbing.md"
            "models/regulation.md"
            "models/extraction.md"
        ]
        "CLI" => [
            "tools/experiment.md"
            "tools/inspect.md"
            "tools/reify.md"
            "tools/export.md"
        ]
    ],
    modules = [
        Models
        Specifications
    ]
)
