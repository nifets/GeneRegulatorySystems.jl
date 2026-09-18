module Trajectories

using WGLMakie
using GeneRegulatorySystems: Models, Scheduling
using GeneRegulatorySystems.Visualisation:
    Catenation,
    Series,
    CountSeries,
    FractionSeries,
    Dimension,
    branch, by_kind, cut, place!, seriestype
import GeneRegulatorySystems.Visualisation: describe
import Colors
using Statistics: mean, std, var
using MultivariateStats: PCA, fit, predict
import UMAP

include("sink.jl")
include("lod.jl")
include("render.jl")
include("phase.jl")
include("distributions.jl")

end
