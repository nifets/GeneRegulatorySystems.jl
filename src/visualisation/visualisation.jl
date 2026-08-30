module Visualisation

using Printf
using ..Models
import Colors: Colors, Color, @colorant_str

struct GroupColors
    colors::Dict{String, Color}
end

GroupColors(::Nothing; kwargs...) = GroupColors(Dict{String, Color}())

function GroupColors(
    groups::AbstractVector;
    reserved=[colorant"white", colorant"black", colorant"crimson"],
    fixed=[
        colorant"#4196FF"
        colorant"#3CBC0B"
        colorant"#E38600"
        colorant"#D864F1"
        colorant"#A07CFF"
        colorant"#00C1D3"
        colorant"#A9A600"
        colorant"#FB53B1"
        colorant"#004D9D"
        colorant"#335D00"
        colorant"#8E2400"
        colorant"#6C00AB"
        colorant"#3D30AC"
        colorant"#006A67"
        colorant"#6B4900"
        colorant"#8D007C"
    ],
    seed=vcat(reserved, fixed),
    drop=length(reserved),
)
    groups = sort!(unique(string.(groups)))
    colors = Colors.distinguishable_colors(length(groups) + drop, seed)
    GroupColors(Dict(zip(groups, colors[(drop + 1):end])))
end

group_colors(groups) = GroupColors(groups)

Base.getindex(colors::GroupColors, group::Symbol) = colors[string(group)]
Base.getindex(colors::GroupColors, group::String) =
    get(colors.colors, group, colorant"gray")
Base.get(colors::GroupColors, group, default) =
    get(colors.colors, string(group), default)

const LINK_COLORS = Dict(
    :activation  => "#7ad9ac",
    :repression  => "#ea7e7e",
    :proteolysis => "#ff7f00",
    :promotes    => "#c4c4cb",
    :inhibits    => "#c4c4cb",
    :affects     => "#3aa6b9",
    :substrate   => "#c4c4cb",
    :product     => "#c4c4cb",
)

include("networks.jl")

end
