using BenchmarkTools
using FFTW
using Krylov
using LinearAlgebra
using LinearMaps
using Statistics

import FourierTools

abstract type AbstractGridOperator{T,N} <: LinearMap{T} end

grid_size(L::AbstractGridOperator{T,N}, d) where {T,N} =
    d::Integer <= N ? grid_size(L)[d] : 1
grid_step(L::AbstractGridOperator{T,N}, d) where {T,N} =
    d::Integer <= N ? grid_step(L)[d] : zero(T)

struct MinusLaplaceOperator{T,DIM} <: AbstractGridOperator{T,DIM}
    size::NTuple{2,Int}
    grid_size::NTuple{DIM,Int}
    grid_step::NTuple{DIM,T}
    function MinusLaplaceOperator{T,DIM}(
        N::NTuple{DIM,Integer},
        h::NTuple{DIM,T},
    ) where {T,DIM}
        num_cells = prod(N)
        return new{T,DIM}((num_cells, num_cells), N, h)
    end
end

MinusLaplaceOperator(N::NTuple{DIM,Integer}, h::NTuple{DIM,T}) where {T,DIM} =
    MinusLaplaceOperator{T,DIM}(N, h)

Base.size(L::MinusLaplaceOperator) = L.size
grid_size(L::MinusLaplaceOperator) = L.grid_size
grid_step(L::MinusLaplaceOperator) = L.grid_step
depth(L::MinusLaplaceOperator) = 1

function LinearMaps._unsafe_mul!(
    v::AbstractVecOrMat,
    L::MinusLaplaceOperator{T,2},
    u::AbstractVecOrMat,
) where {T}
    (Nx, Ny) = grid_size(L)
    (hx², hy²) = grid_step(L) .^ 2
    u_arr = reshape(u, depth(L), Nx, Ny)
    v_arr = reshape(v, size(u_arr))

    for j₀ ∈ 1:Ny
        j₋₁ = j₀ == 1 ? Ny : j₀ - 1
        j₊₁ = j₀ == Ny ? 1 : j₀ + 1
        for i₀ ∈ 1:Nx
            i₋₁ = i₀ == 1 ? Nx : i₀ - 1
            i₊₁ = i₀ == Nx ? 1 : i₀ + 1
            v_arr[1, i₀, j₀] = -(
                (u_arr[1, i₊₁, j₀] - 2 * u_arr[1, i₀, j₀] + u_arr[1, i₋₁, j₀]) / hx² +
                (u_arr[1, i₀, j₊₁] - 2 * u_arr[1, i₀, j₀] + u_arr[1, i₀, j₋₁]) / hy²
            )
        end
    end
    return v
end

LinearMaps.MulStyle(L::MinusLaplaceOperator) = FiveArg()

function LinearMaps._unsafe_mul!(
    v::AbstractVector,
    L::MinusLaplaceOperator{T,2},
    u::AbstractVector,
    α::Number,
    β::Number,
) where {T}
    (Nx, Ny) = grid_size(L)
    (hx², hy²) = grid_step(L) .^ 2
    u_arr = reshape(u, depth(L), Nx, Ny)
    v_arr = reshape(v, size(u_arr))

    for j₀ ∈ 1:Ny
        j₋₁ = j₀ == 1 ? Ny : j₀ - 1
        j₊₁ = j₀ == Ny ? 1 : j₀ + 1
        for i₀ ∈ 1:Nx
            i₋₁ = i₀ == 1 ? Nx : i₀ - 1
            i₊₁ = i₀ == Nx ? 1 : i₀ + 1
            δ²_yy_u = (u_arr[1, i₀, j₊₁] - 2 * u_arr[1, i₀, j₀] + u_arr[1, i₀, j₋₁]) / hy²
            δ²_xx_u = (u_arr[1, i₊₁, j₀] - 2 * u_arr[1, i₀, j₀] + u_arr[1, i₋₁, j₀]) / hx²
            Δu = δ²_xx_u + δ²_yy_u
            v_arr[1, i₀, j₀] = α * (-Δu) + β * v_arr[1, i₀, j₀]
        end
    end
    return v
end

LinearAlgebra.issymmetric(L::MinusLaplaceOperator) = true

struct MinusLaplaceOperatorFourier{T,DIM} <: AbstractGridOperator{T,DIM}
    size::NTuple{2,Int}
    grid_size::NTuple{DIM,Int}
    grid_step::NTuple{DIM,T}
    k²::NTuple{DIM,Vector{T}}
    cache::Array{Complex{T},DIM}
    function MinusLaplaceOperatorFourier{T,DIM}(
        N::NTuple{DIM,Integer},
        h::NTuple{DIM,T},
    ) where {T,DIM}
        num_cells = prod(N)
        f(N_, h_) = [(2 / h_ * sin(π * n / N_))^2 for n = 0:(N_-1)]
        k² = f.(N, h)
        cache = Array{complex(T)}(undef, FourierTools.rfft_size(N, 1:DIM))
        return new{T,DIM}((num_cells, num_cells), N, h, tuple(k²...), cache)
    end
end

MinusLaplaceOperatorFourier(N::NTuple{DIM,Integer}, h::NTuple{DIM,T}) where {T,DIM} =
    MinusLaplaceOperatorFourier{T,DIM}(N, h)

Base.size(L::MinusLaplaceOperatorFourier) = L.size
grid_size(L::MinusLaplaceOperatorFourier) = L.grid_size
grid_step(L::MinusLaplaceOperatorFourier) = L.grid_step
depth(::MinusLaplaceOperatorFourier) = 1

cache_size(L::MinusLaplaceOperatorFourier{T,DIM}) where {T,DIM} =
    FourierTools.rfft_size(grid_size(L), 1:DIM)
create_cache(L::MinusLaplaceOperatorFourier{T,DIM}) where {T,DIM} =
    Array{complex(T)}(undef, cache_size(L))

function mul_fourier!(v̂, L::MinusLaplaceOperatorFourier{T,DIM}, n, û) where {T,DIM}
    k² = zero(T)
    for (i, n_i) ∈ enumerate(n)
        k² += L.k²[i][n_i]
    end
    v̂ .= k² .* û
end

function LinearMaps._unsafe_mul!(v, L::MinusLaplaceOperatorFourier{T,2}, u) where {T}
    cache = isnothing(L.cache) ? create_cache(L) : L.cache
    N = grid_size(L)
    u_grid = reshape(u, depth(L), N...)
    v_grid = reshape(v, size(u_grid))
    ℱ = plan_rfft(u_grid[1, :, :])
    û_n = zeros(complex(T), 1)
    v̂_n = zeros(complex(T), 1)
    mul!(cache, ℱ, u_grid[1, :, :])
    for n in eachindex(IndexCartesian(), cache)
        û_n[1] = cache[n]
        mul_fourier!(v̂_n, L, Tuple(n), û_n)
        cache[n] = v̂_n[1]
    end
    ldiv!(view(v_grid, 1, :, :), ℱ, cache)
    return v
end

Lx = 2.5
Ly = 5.0
Nx = 80
Ny = 40

ℒ₁ = MinusLaplaceOperator{Float64,2}((Nx, Ny), (Lx / Nx, Ly / Ny))
ℒ₂ = MinusLaplaceOperatorFourier{Float64,2}((Nx, Ny), (Lx / Nx, Ly / Ny))

# ux = [sin(π * (n + 0.5) / Nx) for n = 0:(Nx-1)]
# uy = [sin(π * (n + 0.5) / Ny) for n = 0:(Ny-1)]
# u = ux .* uy'
# u_vec = reshape(u, :)
# v_vec = L * u_vec
# v = reshape(v_vec, grid_size(L))

# α = 2.0
# β = 3.0

# @benchmark mul!(v_vec, L, u_vec, α, β)

# u_true = rand(Float64, size(L, 2))
# f = L * u_true

# rtol = 1e-10
# atol = 1e-10
# (u, stats) = cg(L, f, atol=atol, rtol=rtol)

# v = u - u_true
# r = abs.(v .- mean(v))
# @assert all(isapprox.(r, 0.0, rtol=10rtol, atol=10atol))

rtol = 1e-12
atol = 1e-12

u = rand(Float64, size(ℒ₁, 2))
v0 = ℒ₁ * u

v1 = ℒ₂ * u
@assert all(isapprox.(v0, v1, rtol = rtol, atol = atol))

# cache = create_cache(ℒ₂)
# @benchmark __mul!(v1, ℒ₂, u, cache)
@benchmark ℒ₂ * u
