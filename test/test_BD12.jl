using Scapin
using Scapin.Elasticity
using Scapin.BD12
using Test

@testset "BD12 module" begin
    @testset "Creation of discrete Green operator for elasticity" begin
        T = Float64
	for d = 2:3
            s = d * (d + 1)
            @testset "d = $d" begin
                C = Hooke{d,T}(1.0, 0.3)
                N = Tuple(1:d)
                Γ = GreenOperatorHooke{d,T}(C)
                Γ_N = BrisardDormieux2012(Γ, N)

                @test eltype(typeof(Γ_N)) == T
                @test eltype(Γ_N) == T

                @test ndims(typeof(Γ_N)) == 2
                @test ndims(Γ_N) == 2

                @test size(Γ_N) == size(Γ)
                @test size(Γ_N, 1) == size(Γ, 1)
                @test size(Γ_N, 2) == size(Γ,2)

                @test dimensionality(typeof(Γ_N)) == d
                @test dimensionality(Γ_N) == d
            end
        end
    end
end
