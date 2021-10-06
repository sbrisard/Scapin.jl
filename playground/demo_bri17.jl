using CairoMakie
using FFTW
using Scapin.Elasticity
using Scapin.Bri17

const T = Float64
const d = 2

C = Hooke{T,d}(1.0, 0.3)
N = (32, 32)
h = 1.0 ./ N

τ = zeros(T, 3, N...)
τ[3, 1:4, 1:4] .= 1
τ̂ = fft(τ, 2:(d+1))
ε̂ = Array{eltype(τ̂)}(undef, size(τ̂)...)

𝒩 = CartesianIndices(N)

for n ∈ 𝒩
    apply_discrete_green_operator!(view(ε̂, :, n), τ̂[:, n], n, N, h, C)
end

ε = real.(ifft(ε̂, 2:(d+1)))

fig, ax, hm = heatmap(ε[3, :, :])
save("eps_xy.png", fig)
