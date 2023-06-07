using BenchmarkTools
using FFTW
using LinearAlgebra
using LinearMaps
using SciMLOperators

abstract type AbstractGridOperator{T, N} <: SciMLOperators.AbstractSciMLOperator{T} end

grid_size(L::AbstractGridOperator{T, N}, d) where {T, N} = d::Integer <= N ? grid_size(L)[d] : 1
grid_step(L::AbstractGridOperator{T, N}, d) where {T, N} = d::Integer <= 2 ? grid_step(L)[d] : zero(T)

struct LaplaceOperator{T, DIM} <: AbstractGridOperator{T, DIM}
    size::NTuple{2, Int}
    grid_size::NTuple{DIM, Int}
    grid_step::NTuple{DIM, T}
    function LaplaceOperator{T, DIM}(N::NTuple{DIM, Integer}, h::NTuple{DIM, T}) where {T, DIM}
        num_cells = prod(N)
        return new{T, DIM}((num_cells, num_cells), N, h)
    end
end

Base.size(L::LaplaceOperator) = L.size
grid_size(L::LaplaceOperator) = L.grid_size
grid_step(L::LaplaceOperator) = L.grid_step

function LinearAlgebra.mul!(v::AbstractVecOrMat, Δ::LaplaceOperator{T, 2}, u::AbstractVecOrMat) where T
    (Nx, Ny) = grid_size(Δ)
    (hx², hy²) = grid_step(Δ) .^ 2
    u_arr = reshape(u, Nx, Ny)
    v_arr = reshape(v, Nx, Ny)

    for j₀ ∈ 1:Ny
        j₋₁ = j₀ == 1 ? Ny : j₀ - 1
        j₊₁ = j₀ == Ny ? 1 : j₀ + 1
        for i₀ ∈ 1:Nx
            i₋₁ = i₀ == 1 ? Nx : i₀ - 1
            i₊₁ = i₀ == Nx ? 1 : i₀ + 1
            v_arr[i₀, j₀] = (
                (u_arr[i₊₁, j₀] - 2 * u_arr[i₀, j₀] + u_arr[i₋₁, j₀]) / hx² +
                    (u_arr[i₀, j₊₁] - 2 * u_arr[i₀, j₀] + u_arr[i₀, j₋₁]) / hy²
            )
        end
    end
    nothing
end

function Base.:*(Δ::LaplaceOperator{T, 2}, u::AbstractVecOrMat) where T
    v = similar(u)
    mul!(v, Δ, u)
    return v
end

# function laplace_operator_SciML(Lx::Float64, Ly::Float64, Nx::Int, Ny::Int)
#     Δ = LaplaceOperatorSciML(Lx, Ly, Nx, Ny)
#     u = zeros(Float64, Nx * Ny)
#     v = similar(u)
#     return FunctionOperator(Δ, u, v; isinplace=true)
# end

# Δ_SciML = laplace_operator_SciML(Lx, Ly, Nx, Ny)
# v_SciML = reshape(Δ_SciML * u_vec, Nx, Ny)

# @assert all(isapprox.(v_SciML, v_ref, atol = 1e-15, rtol=1e-12))

Lx = 2.5
Ly = 5.0
Nx = 4
Ny = 8

Δ = LaplaceOperator{Float64, 2}((Nx, Ny), (Lx / Nx, Ly/Ny))

@assert isconstant(Δ) == true

ux = [sin(π * (n + 0.5) / Nx) for n = 0:(Nx-1)]
uy = [sin(π * (n + 0.5) / Ny) for n = 0:(Ny-1)]
u = ux .* uy'
u_vec = reshape(u, :)
v_vec = Δ * u_vec
v = reshape(v_vec, grid_size(Δ))
