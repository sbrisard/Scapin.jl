module Scapin

include("hooke.jl")
include("bri17.jl")

export Hooke, bulk_modulus, block_apply!, block_matrix, modal_strain_displacement!, modal_strain_displacement, modal_stiffness!, modal_stiffness

end
