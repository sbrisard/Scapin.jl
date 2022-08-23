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

function refine_and_apply(Γ, τ, r, ctor)
    T = eltype(Γ)
    N = (2^r) .* (size(τ)[2:end])
    Γ_N = ctor(Γ, N)
    τ_N = imresize(τ, size(Γ, 2), N..., method = BSpline(Constant()))
    return apply(Γ_N, τ_N)
end


N_coarse = (4, 4)
r_max = 8
τ_coarse = zeros(T, size(Γ, 2), N_coarse...)
τ_coarse[end, 1, 1] = 1.0

N_ref = (2^r_max) .* N_coarse
ε_ref = Array{T}(undef, size(Γ, 2), N_ref...)

results = Dict()
for r ∈ 0:r_max
    results[r] = refine_and_apply(Γ, τ_coarse, r, BrisardDormieux2012)
end

ε_ref = refine_and_apply(Γ, τ_coarse, r_max + 1, (Γ, N) -> DiscreteGreenOperatorBri17{d,T}(Γ.C, N, Tuple(fill(one(T), d))))

# for r ∈ 0:r_max+1
#     N = (2^r) .* N_coarse
#     𝒩 = CartesianIndices(N)
#     h = 1.0 ./ N

#     τ = zeros(T, size(Γ, 2), N...)
#     τ[end, fill(1:2^r, d)...] .= one(T)

#     if r <= r_max
#         Γ_N = BrisardDormieux2012(Γ, N)
#         ε = apply(Γ_N, τ)
#         # display(heatmap(real.(ε[3, :, :]), c=:viridis))
#         results[r] = ε
#     else
#         # Compute reference solution
#         global ε_ref = apply(DiscreteGreenOperatorBri17{d,T}(C, N, h), τ)
#     end
# end

x = [(2^r) .* N_coarse[1] for r ∈ 0:r_max]
err = zeros(Float64, r_max + 1)
for r ∈ 0:r_max
    ε = imresize(results[r], size(ε_ref)..., method = BSpline(Constant()))
    err[r+1] = norm(ε - ε_ref) / norm(ε_ref)
end

plot(
    x,
    err,
    xscale = :log10,
    yscale = :log10,
    markershape = :circle,
    xlabel = "Grid size, N",
    ylabel = "Relative L² error",
    label=""
)

savefig("convergence.png")
