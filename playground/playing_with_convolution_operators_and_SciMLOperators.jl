using BenchmarkTools
using FFTW
using LinearAlgebra
using SciMLOperators

struct RectilinearGrid{T,DIM}
    size::NTuple{DIM,Int}
    step::NTuple{DIM,T}
end

Base.size(g::RectilinearGrid{T,N}) where {T,N} = g.size
Base.step(g::RectilinearGrid{T,N}) where {T,N} = g.step
Base.size(g::RectilinearGrid{T,N}, k) where {T,N} = k::Integer <= N ? size(g)[k] : 1
Base.step(g::RectilinearGrid{T,N}, k) where {T,N} = k::Integer <= N ? step(g)[k] : zero(T)

Lx = 2.5
Ly = 5.0
Nx = 40
Ny = 80
g = RectilinearGrid{Float64,2}((Nx, Ny), (Lx / Nx, Ly / Ny))

# abstract type AbstractGridOperator{T, N} <: LinearMap{T} end

# grid_size(L::AbstractGridOperator{T, N}, d) where {T, N} = d::Integer <= N ? grid_size(L)[d] : 1
# grid_step(L::AbstractGridOperator{T, N}, d) where {T, N} = d::Integer <= N ? grid_step(L)[d] : zero(T)

# struct MinusLaplaceOperator{T, DIM} <: AbstractGridOperator{T, DIM}
#     size::NTuple{2, Int}
#     grid_size::NTuple{DIM, Int}
#     grid_step::NTuple{DIM, T}
#     function MinusLaplaceOperator{T, DIM}(N::NTuple{DIM, Integer}, h::NTuple{DIM, T}) where {T, DIM}
#         num_cells = prod(N)
#         return new{T, DIM}((num_cells, num_cells), N, h)
#     end
# end

# Base.size(L::MinusLaplaceOperator) = L.size
# grid_size(L::MinusLaplaceOperator) = L.grid_size
# grid_step(L::MinusLaplaceOperator) = L.grid_step

# function LinearMaps._unsafe_mul!(v::AbstractVecOrMat, L::MinusLaplaceOperator{T, 2}, u::AbstractVecOrMat) where T
#     (Nx, Ny) = grid_size(L)
#     (hx², hy²) = grid_step(L) .^ 2
#     u_arr = reshape(u, Nx, Ny)
#     v_arr = reshape(v, Nx, Ny)

#     for j₀ ∈ 1:Ny
#         j₋₁ = j₀ == 1 ? Ny : j₀ - 1
#         j₊₁ = j₀ == Ny ? 1 : j₀ + 1
#         for i₀ ∈ 1:Nx
#             i₋₁ = i₀ == 1 ? Nx : i₀ - 1
#             i₊₁ = i₀ == Nx ? 1 : i₀ + 1
#             v_arr[i₀, j₀] = -(
#                 (u_arr[i₊₁, j₀] - 2 * u_arr[i₀, j₀] + u_arr[i₋₁, j₀]) / hx² +
#                     (u_arr[i₀, j₊₁] - 2 * u_arr[i₀, j₀] + u_arr[i₀, j₋₁]) / hy²
#             )
#         end
#     end
#     return v
# end

# LinearMaps.MulStyle(L::MinusLaplaceOperator) = FiveArg()

# function LinearMaps._unsafe_mul!(v::AbstractVector, L::MinusLaplaceOperator{T, 2}, u::AbstractVector, α::Number, β::Number) where T
#     (Nx, Ny) = grid_size(L)
#     (hx², hy²) = grid_step(L) .^ 2
#     u_arr = reshape(u, Nx, Ny)
#     v_arr = reshape(v, Nx, Ny)

#     for j₀ ∈ 1:Ny
#         j₋₁ = j₀ == 1 ? Ny : j₀ - 1
#         j₊₁ = j₀ == Ny ? 1 : j₀ + 1
#         for i₀ ∈ 1:Nx
#             i₋₁ = i₀ == 1 ? Nx : i₀ - 1
#             i₊₁ = i₀ == Nx ? 1 : i₀ + 1
#             δ²_yy_u = (u_arr[i₀, j₊₁] - 2 * u_arr[i₀, j₀] + u_arr[i₀, j₋₁]) / hy²
#             δ²_xx_u = (u_arr[i₊₁, j₀] - 2 * u_arr[i₀, j₀] + u_arr[i₋₁, j₀]) / hx²
#             Δu = δ²_xx_u + δ²_yy_u
#             v_arr[i₀, j₀] = α * (-Δu) + β * v_arr[i₀, j₀]
#         end
#     end
#     return v
# end

# LinearAlgebra.issymmetric(L::MinusLaplaceOperator) = true
# LinearAlgebra.isposdef(L::MinusLaplaceOperator) = true

# Lx = 2.5
# Ly = 5.0
# Nx = 40
# Ny = 80

# L = MinusLaplaceOperator{Float64, 2}((Nx, Ny), (Lx / Nx, Ly/Ny))

# ux = [sin(π * (n + 0.5) / Nx) for n = 0:(Nx-1)]
# uy = [sin(π * (n + 0.5) / Ny) for n = 0:(Ny-1)]
# u = ux .* uy'
# u_vec = reshape(u, :)
# v_vec = L * u_vec
# v = reshape(v_vec, grid_size(L))

# α = 2.0
# β = 3.0

# @benchmark mul!(v_vec, L, u_vec, α, β)
