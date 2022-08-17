using Scapin
using Scapin.Elasticity
using Scapin.BD12

const T = Float64
const d = 2

C = Hooke{d,T}(1.0, 0.3)
N = (16, 32)

Γ = GreenOperatorHooke{d,T}(C)
Γ_N = BrisardDormieux2012(Γ, N)

@show typeof(Γ_N)
@show eltype(typeof(Γ_N))
@show eltype(Γ_N)
@show ndims(typeof(Γ_N))
@show ndims(Γ_N)
@show size(Γ_N)
@show dimensionality(typeof(Γ_N))
@show dimensionality(Γ_N)
