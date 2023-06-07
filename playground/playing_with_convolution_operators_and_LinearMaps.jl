using BenchmarkTools
using FFTW
using LinearAlgebra
using LinearMaps

struct LaplaceOperator <: LinearMaps.LinearMap{Float64}
    Lx::Float64
    Ly::Float64
    Nx::Int
    Ny::Int
    size::Dims{2}
    Δx²::Float64
    Δy²::Float64
    function LaplaceOperator(Lx::Float64, Ly::Float64, Nx::Int, Ny::Int)
        ncells = Nx * Ny
        return new(Lx, Ly, Nx, Ny, Dims([ncells, ncells]), (Lx / Nx)^2, (Ly / Ny)^2)
    end
end

Base.size(Δ::LaplaceOperator) = Δ.size

function LinearMaps._unsafe_mul!(v, Δ::LaplaceOperator, u::AbstractVector)
    u_arr = reshape(u, Δ.Nx, Δ.Ny)
    v_arr = reshape(v, Δ.Nx, Δ.Ny)
    for j₀ ∈ 1:Δ.Ny
        j₋₁ = j₀ == 1 ? Δ.Ny : j₀ - 1
        j₊₁ = j₀ == Δ.Ny ? 1 : j₀ + 1
        for i₀ ∈ 1:Δ.Nx
            i₋₁ = i₀ == 1 ? Δ.Nx : i₀ - 1
            i₊₁ = i₀ == Δ.Nx ? 1 : i₀ + 1
            v_arr[i₀, j₀] = (
                (u_arr[i₊₁, j₀] - 2 * u_arr[i₀, j₀] + u_arr[i₋₁, j₀]) / Δ.Δx² +
                    (u_arr[i₀, j₊₁] - 2 * u_arr[i₀, j₀] + u_arr[i₀, j₋₁]) / Δ.Δy²
            )
        end
    end
    return v
end

mutable struct ModalLaplaceOperator{T} <: LinearMaps.LinearMap{T}
    const Lx::Float64
    const Ly::Float64
    const Nx::Int
    const Ny::Int
    const kx²::Vector{Float64}
    const ky²::Vector{Float64}
    nx::Int
    ny::Int
    function ModalLaplaceOperator{T}(Lx, Ly, Nx, Ny) where {T}
        kx₀ = 2Nx / Lx
        ky₀ = 2Ny / Ly
        kx² = [(kx₀ * sin(π * n / Nx))^2 for n = 0:(Nx-1)]
        ky² = [(ky₀ * sin(π * n / Ny))^2 for n = 0:(Ny-1)]
        new{T}(Lx, Ly, Nx, Ny, kx², ky²)
    end
end

Base.size(Δ::ModalLaplaceOperator) = (1, 1)

function update_frequency!(Δ_hat:: ModalLaplaceOperator, nx, ny)
    Δ_hat.nx = nx
    Δ_hat.ny = ny
end

function LinearMaps._unsafe_mul!(v, Δ_hat::ModalLaplaceOperator, u)
    v .= -(Δ_hat.kx²[Δ_hat.nx] + Δ_hat.ky²[Δ_hat.ny]) .* u
    return v
end

Lx = 2.5
Ly = 5.0
Nx = 4
Ny = 8

Δ_ref = LaplaceOperator(Lx, Ly, Nx, Ny)

ux = [sin(π * (n + 0.5) / Nx) for n = 0:(Nx-1)]
uy = [sin(π * (n + 0.5) / Ny) for n = 0:(Ny-1)]
u = ux .* uy'
u_vec = reshape(u, :)
v_ref = reshape(Δ_ref * u_vec, Nx, Ny)

Δ_hat = ModalLaplaceOperator{ComplexF64}(Lx, Ly, Nx, Ny)
u_hat = reshape(fft(u), Nx, Ny, 1)
v_hat = similar(u_hat)

for nx = 1:Nx
    for ny = 1:Ny
        update_frequency!(Δ_hat, nx, ny)
        mul!(view(v_hat, nx, ny, :), Δ_hat, u_hat[nx, ny, :])
    end
end
v = ifft(v_hat)

@assert all(isapprox.(v, v_ref, atol = 1e-15, rtol = 1e-12))
