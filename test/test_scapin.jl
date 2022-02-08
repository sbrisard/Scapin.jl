using Base.Iterators
using FFTW
using Scapin
using Scapin.Bri17
using Scapin.Elasticity
using Test

@testset "Scapin module" begin
    @testset "Real vs. complex FFTs" begin
        Random.seed!(202202080909)
        T = Float64
        𝒩 = [8, 9, 15, 16]

        for d = 2:3
            C = Hooke{d,T}(1.0, 0.3)
            for N ∈ product(repeated(𝒩, d)...)
                τ = rand(Float64, div(d * (d + 1), 2), N...)
                Γ = DiscreteGreenOperatorBri17{d,T}(C, N, 1.0 ./ N)

                ε_c = apply(Γ, τ, plan_fft(τ, 2:(d+1)))
                ε_r = apply(Γ, τ, plan_rfft(τ, 2:(d+1)))

                @test real.(ε_c) ≈ ε_r rtol = 1e-15 atol = 1e-15
            end
        end
    end
end
