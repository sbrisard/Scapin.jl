using ImageTransformations
using Interpolations
using LinearAlgebra
using Plots
using Scapin
using Scapin.Elasticity
using Scapin.BD12
using Scapin.Bri17

const T = Float64
const d = 2

C = Hooke{d,T}(1.0, 0.3)
Γ = GreenOperatorHooke{d,T}(C)

N_coarse = (4, 4)
r_max = 9
τ_coarse = zeros(T, size(Γ, 2), N_coarse...)
τ_coarse[end, 1, 1] = 1.0

N_ref = (2^r_max) .* N_coarse
ε_ref = Array{T}(undef, size(Γ, 2), N_ref...)

results = Dict()

for r ∈ 0:r_max
    N = (2^r) .* N_coarse
    𝒩 = CartesianIndices(N)
    h = 1.0 ./ N

    τ = zeros(T, size(Γ, 2), N...)
    τ[end, fill(1:2^r, d)...] .= one(T)

    Γ_N = BrisardDormieux2012(Γ, N)
    ε = apply(Γ_N, τ)
    # display(heatmap(real.(ε[3, :, :]), c=:viridis))
    results[r] = ε

    if r == r_max
        # Compute reference solution
        global ε_ref = apply(DiscreteGreenOperatorBri17{d,T}(C, N, h), τ)
    end
end

for r ∈ 0:(r_max-1)
    ε = imresize(results[r], size(ε_ref)..., method=BSpline(Constant()))
    println(norm(ε - ε_ref))
end
