using LinearAlgebra
using Plots
using Scapin
using Scapin.Elasticity
using Scapin.BD12
using Scapin.Bri17
using Scapin.ConvergenceAnalysis

function convergence_analysis(Γ, τ_coarse, r_max, ctor)
    # TODO Move to signature of function
    d = dimensionality(Γ)
    T = eltype(Γ)
    N_ref = (2^(r_max + 1)) .* size(τ_coarse)[2:end]
    h = Tuple(fill(one(T), d))

    Γ_ref = DiscreteGreenOperatorBri17{d,T}(Γ.C, N_ref, h)

    ε_ref = refine_and_apply(Γ_ref, τ_coarse)

    ε = []
    for r ∈ 0:r_max
        N = (2^r) .* (size(τ_coarse)[2:end])
        Γ_N = ctor(Γ, N)
        push!(ε, refine_and_apply(Γ_N, τ_coarse))
    end

    return ε, ε_ref
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
    norm(imresize(ε_r, size(ε_ref)..., method = BSpline(Constant())) - ε_ref) / norm(ε_ref) for ε_r ∈ ε
]

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
plot!(x, 10 .* x .^ (-1), label = "∝ N⁻¹")

savefig("convergence.pdf")
