app = @__DIR__
root = normpath(joinpath(app, ".."))
julia = Base.julia_cmd()

docs_project = joinpath(root, "docs")
docs_code = "using LiveServer; cd($(repr(root))) do; servedocs(port=8001); end"
docs_process = run(`$julia --project=$docs_project -e $docs_code`; wait=false)

notebook = joinpath(app, "app_notebook.jl")
app_code = "using Pluto; Pluto.run(notebook=$(repr(notebook)))"

try
    run(`$julia --project=$app -e $app_code`)
finally
    Base.process_running(docs_process) && kill(docs_process)
end
