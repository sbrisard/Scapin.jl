module Scapin

using Base.Iterators
using LinearAlgebra

include("hooke.jl")
include("bri17.jl")

export Hooke, bulk_modulus, block_apply!, block_matrix, modal_strain_displacement!, modal_strain_displacement, modal_stiffness!, modal_stiffness, element_nodes, integrate, shape, gradient_operator, strain_displacement_operator, avg_strain_displacement_operator, stiffness_operator, cell_vertices

end
