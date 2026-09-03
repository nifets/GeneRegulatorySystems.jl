app = @__DIR__
root = normpath(joinpath(app, ".."))
julia = Base.julia_cmd()

docs = joinpath(root, "docs")
build = joinpath(docs, "build")
docs_code = if isdir(build)
    "using LiveServer; serve(dir=$(repr(build)), port=8001)"
else
    "using LiveServer; cd($(repr(root))) do; servedocs(port=8001); end"
end

docs_command = `$julia --project=$docs -e $docs_code`
docs_process = run(Cmd(docs_command; detach=true); wait=false)

using Pluto

try
    Pluto.run(notebook=joinpath(app, "gene_regulatory_systems.jl"))
finally
    Base.process_running(docs_process) && kill(docs_process)
end
