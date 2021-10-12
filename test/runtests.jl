using LinearAlgebra
using Scapin
using Scapin.Elasticity
using StaticArrays
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
        C = Hooke{d,T}(1.0, 0.3)
        C_exp = d * bulk_modulus(C) * J₄ + 2C.μ * K₄
        C_act = convert(Array, C)
        @test isapprox(C_act, C_exp, rtol = 1e-15)
    end
end

function block_matrix_ref(hooke::Hooke{d,T}, k::SVector{d,T}) where {d,T}
    if d == 2
        sym = 3
        ij2i = [1, 2, 1]
        ij2j = [1, 2, 2]
    elseif d == 3
        sym = 6
        ij2i = [1, 2, 3, 2, 3, 1]
        ij2j = [1, 2, 3, 3, 1, 2]
    else
        throw(ArgumentError("d must be 2 or 3 (was $d)"))
    end

    mat = zeros(T, sym, sym)
    n = k / norm(k)
    for ij = 1:sym
        i = ij2i[ij]
        j = ij2j[ij]
        w_ij = ij <= d ? one(T) : sqrt(2 * one(T))
        for kl = 1:sym
            k = ij2i[kl]
            l = ij2j[kl]
            w_kl = kl <= d ? one(T) : sqrt(2 * one(T))
            δik_nj_nl = i == k ? n[j] * n[l] : zero(T)
            δil_nj_nk = i == l ? n[j] * n[k] : zero(T)
            δjk_ni_nl = j == k ? n[i] * n[l] : zero(T)
            δjl_ni_nk = j == l ? n[i] * n[k] : zero(T)
            aux1 = (δik_nj_nl + δil_nj_nk + δjk_ni_nl + δjl_ni_nk) / 4
            aux2 = n[i] * n[j] * n[k] * n[l] / (2 * (1 - hooke.ν))
            mat[ij, kl] = w_ij * w_kl * (aux1 - aux2) / hooke.μ
        end
    end
    mat
end

@testset "Green operator, Hooke 2D" begin
    hooke = Hooke{2,Float64}(5.6, 0.3)
    for k_norm ∈ [0.12, 2.3, 14.5]
        for θ ∈ LinRange(0.0, 2 * π, 21)[1:end-1]
            k = @SVector [k_norm * cos(θ), k_norm * sin(θ)]
            act = block_matrix(hooke, k)
            exp = block_matrix_ref(hooke, k)

            @test all(isapprox.(act, exp, atol = 1e-15))
        end
    end
end

@testset "Green operator, Hooke 3D" begin
    hooke = Hooke{3,Float64}(5.6, 0.3)
    for k_norm ∈ [0.12, 2.3, 14.5]
        for φ ∈ LinRange(0.0, 2 * π, 21)[1:end-1]
            for θ ∈ LinRange(0.0, π, 11)
                k = k_norm * (@SVector [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)])
                act = block_matrix(hooke, k)
                exp = block_matrix_ref(hooke, k)

                @test all(isapprox.(act, exp, atol = 1e-15))
            end
        end
    end
end

# @testset "Discrete Green operator [MS94], Hooke 2D" begin
#     T = Float64
#     DIM = 2
#     Γ = Hooke{T,DIM}(1.2, 0.3)
#     N = @SVector [5, 6]
#     L = @SVector [3.4, 5.6]
#     Γ_h = TruncatedGreenOperator{T,DIM}(Γ, N, L)
#     k = [(2π / L[1]) * [zero(T), 1, 2, -2, -1], (2π / L[2]) * [zero(T), 1, 2, -3, -2, -1]]
#     @test Γ_h.N == N
#     @test Γ_h.L == L
#     @test Γ_h.k == k
# end

# @testset "Discrete Green operator [MS94], Hooke 3D" begin
#     T = Float64
#     DIM = 3
#     Γ = Hooke{T,DIM}(1.2, 0.3)
#     N = @SVector [5, 6, 7]
#     L = @SVector [3.4, 5.6, 7.8]
#     Γ_h = TruncatedGreenOperator{T,DIM}(Γ, N, L)
#     k = [
#         (2π / L[1]) * [zero(T), 1, 2, -2, -1],
#         (2π / L[2]) * [zero(T), 1, 2, -3, -2, -1],
#         (2π / L[3]) * [zero(T), 1, 2, 3, -3, -2, -1],
#     ]
#     @test Γ_h.N == N
#     @test Γ_h.L == L
#     @test Γ_h.k == k
# end

include("test_grid.jl")
include("test_brick.jl")
include("test_bri17.jl")
