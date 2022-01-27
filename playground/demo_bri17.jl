using CairoMakie
using ElectronDisplay
using FFTW
using LinearAlgebra
using Scapin
using Scapin.Elasticity
using Scapin.Bri17

const T = Float64
const d = 2

C = Hooke{d,T}(1.0, 0.3)
α = (0.25, 0.25) # Fraction of the domained that is polarized
N_coarse = (4, 4)
r_max = 9
N_fine = (2^r_max) .* N_coarse

results = Dict()

for r ∈ 0:r_max
    N = (2^r) .* N_coarse
    𝒩 = CartesianIndices(N)
    h = 1.0 ./ N
    Γ̂ = DiscreteGreenOperatorBri17{d,T}(C, N, h)

    τ = zeros(T, 3, N...)
    τ[3, fill(1:2^r, d)...] .= 1
    τ̂ = fft(τ, 2:(d+1))
    ε̂ = Array{eltype(τ̂)}(undef, size(τ̂)...)

    for n ∈ 𝒩
        apply_fourier!(view(ε̂, :, n), Γ̂, n, τ̂[:, n])
    end

    ε = real.(ifft(ε̂, 2:(d+1)))

    ε_fine = zeros(T, 3, N_fine...)
    s = 2^(r_max - r)
    for n ∈ map(Tuple, 𝒩)
        n₁ = CartesianIndex(s .* (n .- 1) .+ 1)
        n₂ = CartesianIndex(s .* n)
        ε_fine[:, n₁:n₂] .= ε[:, n...]
    end

    # fig, ax, hm = heatmap(ε_fine[3, :, :])
    # save("eps_xy-$(N[1])x$(N[2]).png", fig)
    results[r] = ε_fine
end

x = Array{T}(undef, r_max)
y = Array{T}(undef, r_max)
for r ∈ 0:(r_max-1)
    ε₁ = results[r]
    ε₂ = results[r_max]
    x[r+1] = 2^r * N_coarse[1]
    y[r+1] = norm(ε₂-ε₁)
end

C = x[end] * y[end]
fig = scatter(x, y,  axis = (xscale=log10, yscale = log10))
lines!(x, C ./ x)
#save("convergence.png", fig)
electrondisplay(fig)
