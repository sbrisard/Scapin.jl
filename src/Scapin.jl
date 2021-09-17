module Scapin

using Base.Iterators
using LinearAlgebra

include("Grid.jl")
include("hooke.jl")
include("bri17.jl")

export Hooke, bulk_modulus, block_apply!, block_matrix, modal_strain_displacement!, modal_strain_displacement, modal_stiffness!, modal_stiffness, integrate, shape, gradient_operator, strain_displacement_operator, avg_strain_displacement_operator, stiffness_operator, cell_vertices, global_stiffness_operator, element_nodes

end
