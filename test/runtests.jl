using LinearAlgebra
using Scapin
using Test

function isotropic_tensors(d)
    s = (d * (d + 1)) ÷ 2
    I₄ = Matrix{Float64}(I, s, s)
    I₂ = zeros(Float64, s)
    I₂[begin:begin+d-1] .= 1.0
    J₄ = I₂ * transpose(I₂) / d
    K₄ = I₄ - J₄
    return I₄, J₄, K₄
end


function asmatrix(C)
    (nrows, ncols) = size(C)
    C_mat = zeros(nrows, ncols)
    ε = zeros(Float64, ncols)
    for j = 1:ncols
        ε[j] = 1.0
        C_mat[:, j] = C * ε
        ε[j] = 0.0
    end
    return C_mat
end


@testset "Hooke 2D" begin
    d = 2
    (I₄, J₄, K₄) = isotropic_tensors(d)
    μ = 1.0
    ν = 0.3
    κ = μ / (1 - 2ν)
    C_exp = d * κ * J₄ + 2μ * K₄
    C_act = asmatrix(Hooke{d}(μ, ν))
    @show size(C_exp)
    @show size(C_act)
    @assert isapprox(C_act, C_exp, rtol=1e-15)
end
