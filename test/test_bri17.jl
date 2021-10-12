using Base.Iterators
using FFTW
using LinearAlgebra
using Scapin
using Scapin.Bri17
using Scapin.Brick
using Scapin.Elasticity
using Scapin.Grid

using Test

@testset "Bri17 module" begin
    @testset "Modal strain-displacement vector" begin
        for d = 2:3
            @testset "Modal strain-displacement vector, $(d)d" begin
                N = (3, 4, 5)[1:d]
                h = (1.1, 1.2, 1.3)[1:d]
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

    @testset "Modal stiffness matrix" begin
        for d = 2:3
            @testset "Modal stiffness matrix, $(d)d" begin
                N = (3, 4, 5)[1:d]
                h = (1.1, 1.2, 1.3)[1:d]
                C = Hooke{d,Float64}(5.6, 0.3)
                K = global_stiffness_operator(N, h, C)
                K̂_exp = fft(K[fill(:, d + 2)..., fill(1, d)...], 2:(d+1))
                K̂_act = zeros(Complex{Float64}, d, N..., d)
                for n ∈ CartesianIndices(N)
                    K̂_act[:, n, :] = modal_stiffness(n, N, h, C)
                end
                @test K̂_act ≈ K̂_exp rtol = 1e-15 atol = 1e-15
            end
        end
    end
end
