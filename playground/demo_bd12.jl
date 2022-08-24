using ImageTransformations
using Interpolations
using LinearAlgebra
using Plots
using Scapin
using Scapin.Elasticity
using Scapin.BD12
using Scapin.Bri17


function refine_and_apply(Γ, τ, r, ctor)
    T = eltype(Γ)
    N = (2^r) .* (size(τ)[2:end])
    Γ_N = ctor(Γ, N)
    τ_N = imresize(τ, size(Γ, 2), N..., method = BSpline(Constant()))
    return apply(Γ_N, τ_N)
end

function convergence_analysis(Γ, τ_coarse, r_max, ctor)
    N_ref = (2^r_max) .* size(τ_coarse)
    ε_ref = Array{T}(undef, size(Γ, 2), N_ref...)
    ε_ref = refine_and_apply(
        Γ,
        τ_coarse,
        r_max + 1,
        (Γ, N) -> DiscreteGreenOperatorBri17{d,T}(Γ.C, N, Tuple(fill(one(T), d))),
    )

    [refine_and_apply(Γ, τ_coarse, r, ctor) for r ∈ 0:r_max], ε_ref
end

const T = Float64
const d = 2

C = Hooke{d,T}(1.0, 0.3)
Γ = GreenOperatorHooke{d,T}(C)
ctor = BrisardDormieux2012

N_coarse = (4, 4)
r_max = 8
τ_coarse = zeros(T, size(Γ, 2), N_coarse...)
τ_coarse[end, 1, 1] = 1.0



(ε, ε_ref) = convergence_analysis(Γ, τ_coarse, r_max, ctor)

x = [(2^r) .* N_coarse[1] for r ∈ 0:r_max]
y = [
    norm(imresize(ε_r, size(ε_ref)..., method = BSpline(Constant())) - ε_ref) / norm(ε_ref) for ε_r ∈ ε]

plot(
    x,
    y,
    xscale = :log10,
    yscale = :log10,
    markershape = :circle,
    markerstrokecolor = :white,
    xlabel = "Grid size, N",
    ylabel = "Relative L² error",
    label = "BD12",
)
plot!(x, 10 .* x .^ (-1), label="∝ N⁻¹")

savefig("convergence.pdf")
