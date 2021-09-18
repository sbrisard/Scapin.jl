module Scapin

using Base.Iterators
using LinearAlgebra

include("Elasticity.jl")
include("Grid.jl")
include("Brick.jl")

include("bri17.jl")

export modal_strain_displacement!, modal_strain_displacement, modal_stiffness!, modal_stiffness

end
