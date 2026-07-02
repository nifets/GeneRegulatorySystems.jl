module GeneRegulatorySystems

using PrecompileTools

import Random
import SHA

const SPECIFICATION_EXAMPLES = "$(@__DIR__)/../examples/specification"

σ(x) = inv(one(x) + exp(-x))
logit(p) = log(p / (one(p) - p))

include("specifications.jl")
include("models/models.jl")

using .Specifications: Specification
using .Models: Model
using .Models.Scheduling: Scheduling, Schedule

export Specifications
export Specification
export Models
export Model
export Scheduling
export Schedule

@compile_workload begin
    dryrun(_primitive!, _x, _Δt; _...) = nothing
    trace(_state; _...) = nothing

    # Load and dry-run all specification examples:
    for filename in readdir(SPECIFICATION_EXAMPLES)
        schedule! = Models.load("$SPECIFICATION_EXAMPLES/$filename")
        schedule!(; dryrun)
    end

    # Load and additionally describe and simulate some examples:
    examples = [
        (filename = "kronecker.schedule.json", path = "+-2"),
        (filename = "templating.schedule.json", path = "+-1.do"),
        (filename = "differentiation.schedule.json", path = "+.do"),
    ]
    for (; filename, path) in examples
        schedule! = Models.load("$SPECIFICATION_EXAMPLES/$filename")
        Models.describe(Scheduling.reify(schedule!, path))
        schedule!(; trace)
    end
end

end
