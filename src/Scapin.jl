module Scapin

using Base.Iterators
using LinearAlgebra

"""
    dimensionality(op)

Return the number of dimensions of the physical space that `op` operates on.
"""
function dimensionality end

export dimensionality

include("Elasticity.jl")
include("Grid.jl")
include("Brick.jl")
include("Bri17.jl")

end
