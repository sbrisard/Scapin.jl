using Base.Iterators
using FFTW
using LinearAlgebra
using Scapin
using Scapin.Brick
using Scapin.Elasticity
using Scapin.Grid
using Test


@testset "Modal strain-displacement vector" begin
    for d = 2:3
        @testset "Modal strain-displacement vector, $(d)d" begin
            N = (3, 4, 5)[1:d]
            h = (1.1, 1.2, 1.3)[1:d]
            C = Hooke{Float64,d}(5.6, 0.3)
            B = global_strain_displacement_operator(N, h)
            B̂_exp = fft(B[fill(:, d + 3)..., fill(1, d)...], 3:(d+2))
            B̂_act = zeros(Complex{Float64}, d, d, N..., d)
            for n ∈ CartesianIndices(N)
                b̂ = modal_strain_displacement(n, N, h)
                for i = 1:d, j = 1:d
                    B̂_act[i, j, n, i] += b̂[j] / 2
                    B̂_act[i, j, n, j] += b̂[i] / 2
                end
            end
            @test B̂_act ≈ B̂_exp rtol = 1e-15 atol = 1e-15
        end
    end
end


#end

# @testset "Global stiffness operator" begin
#     N = [3, 4]
#     h = [1.1, 1.2]
#     μ = 5.6
#     ν = 0.3
#     K = global_stiffness_operator(N, h, μ, ν)
# end
