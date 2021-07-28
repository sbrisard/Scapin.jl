using LinearAlgebra
using Scapin
using Test

function isotropic_tensors(T, d)
    s = (d * (d + 1)) ÷ 2
    I₄ = Matrix{T}(I, s, s)
    I₂ = zeros(T, s)
    I₂[begin:begin+d-1] .= one(T)
    J₄ = I₂ * transpose(I₂) / d
    K₄ = I₄ - J₄
    return I₄, J₄, K₄
end

for d = 2:3
    @testset "Hooke $(d)D" begin
        T = Float64
        (I₄, J₄, K₄) = isotropic_tensors(T, d)
        C = Hooke{T, d}(1.0, 0.3)
        C_exp = d * bulk_modulus(C) * J₄ + 2C.μ * K₄
        C_act = convert(Array, C)
        @test isapprox(C_act, C_exp, rtol = 1e-15)
    end
end
