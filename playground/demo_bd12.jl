using LinearAlgebra
using Plots
using Scapin
using Scapin.Elasticity
using Scapin.BD12

const T = Float64
const d = 2

C = Hooke{d,T}(1.0, 0.3)
N = (256, 256)

Γ = GreenOperatorHooke{d,T}(C)
Γ_N = BrisardDormieux2012(Γ, N)

τ = zeros(T, 3, N...)
τ[3, 1:(N[1]÷4), 1:(N[2]÷4)] .= 1
ε = apply(Γ_N, τ)

heatmap(real.(ε[3, :, :]), c=:viridis)
