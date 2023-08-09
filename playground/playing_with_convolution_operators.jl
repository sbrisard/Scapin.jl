using BenchmarkTools
using FFTW
using Krylov
using LinearAlgebra
using LinearMaps
using Statistics

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

function LinearMaps._unsafe_mul!(
    v::AbstractVecOrMat,
    L::MinusLaplaceOperator{T,2},
    u::AbstractVecOrMat,
) where {T}
    (Nx, Ny) = grid_size(L)
    (hx², hy²) = grid_step(L) .^ 2
    u_arr = reshape(u, Nx, Ny)
    v_arr = reshape(v, Nx, Ny)

    for j₀ ∈ 1:Ny
        j₋₁ = j₀ == 1 ? Ny : j₀ - 1
        j₊₁ = j₀ == Ny ? 1 : j₀ + 1
        for i₀ ∈ 1:Nx
            i₋₁ = i₀ == 1 ? Nx : i₀ - 1
            i₊₁ = i₀ == Nx ? 1 : i₀ + 1
            v_arr[i₀, j₀] = -(
                (u_arr[i₊₁, j₀] - 2 * u_arr[i₀, j₀] + u_arr[i₋₁, j₀]) / hx² +
                (u_arr[i₀, j₊₁] - 2 * u_arr[i₀, j₀] + u_arr[i₀, j₋₁]) / hy²
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
    u_arr = reshape(u, Nx, Ny)
    v_arr = reshape(v, Nx, Ny)

    for j₀ ∈ 1:Ny
        j₋₁ = j₀ == 1 ? Ny : j₀ - 1
        j₊₁ = j₀ == Ny ? 1 : j₀ + 1
        for i₀ ∈ 1:Nx
            i₋₁ = i₀ == 1 ? Nx : i₀ - 1
            i₊₁ = i₀ == Nx ? 1 : i₀ + 1
            δ²_yy_u = (u_arr[i₀, j₊₁] - 2 * u_arr[i₀, j₀] + u_arr[i₀, j₋₁]) / hy²
            δ²_xx_u = (u_arr[i₊₁, j₀] - 2 * u_arr[i₀, j₀] + u_arr[i₋₁, j₀]) / hx²
            Δu = δ²_xx_u + δ²_yy_u
            v_arr[i₀, j₀] = α * (-Δu) + β * v_arr[i₀, j₀]
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
    function MinusLaplaceOperatorFourier{T,DIM}(
        N::NTuple{DIM,Integer},
        h::NTuple{DIM,T},
    ) where {T,DIM}
        num_cells = prod(N)
        f(N_, h_) = [(2 / h_ * sin(π * n / N_))^2 for n = 0:(N_-1)]
        k² = f.(N, h)
        return new{T,DIM}((num_cells, num_cells), N, h, tuple(k²...))
    end
end

MinusLaplaceOperatorFourier(N::NTuple{DIM,Integer}, h::NTuple{DIM,T}) where {T,DIM} =
    MinusLaplaceOperatorFourier{T,DIM}(N, h)

Base.size(L::MinusLaplaceOperatorFourier) = L.size
grid_size(L::MinusLaplaceOperatorFourier) = L.grid_size
grid_step(L::MinusLaplaceOperatorFourier) = L.grid_step

function mul_fourier!(v̂, L::MinusLaplaceOperatorFourier{T,DIM}, n, û) where {T,DIM}
    k² = zero(T)
    for (i, n_i) ∈ enumerate(n)
        k² += L.k²[i][n_i]
    end
    v̂ .= k² .* û
end

function LinearMaps._unsafe_mul!(
    v::AbstractVecOrMat,
    L::MinusLaplaceOperatorFourier{T,2},
    u::AbstractVecOrMat,
) where {T}
    (Nx, Ny) = grid_size(L)
    u_grid = reshape(u, Nx, Ny)
    v_grid = reshape(v, Nx, Ny)
    û_grid = fft(u_grid)
    v̂_grid = similar(û_grid)
    for ny ∈ 1:Ny
        for nx ∈ 1:Nx
            mul_fourier!(view(v̂_grid, nx, ny, :), L, (nx, ny), û_grid[nx, ny, :])
        end
    end
    v .= real.(reshape(ifft(v̂_grid), :))
    return v
end

function mul_rfft!(
    v::AbstractVecOrMat,
    L::MinusLaplaceOperatorFourier{T,2},
    u::AbstractVecOrMat,
) where {T}
    (Nx, Ny) = grid_size(L)
    u_grid = reshape(u, Nx, Ny)
    v_grid = reshape(v, Nx, Ny)
    û_grid = rfft(u_grid)
    v̂_grid = similar(û_grid)
    for i in eachindex(IndexCartesian(), û_grid)
        nx = i[1]
        ny = i[2]
        mul_fourier!(view(v̂_grid, nx, ny, :), L, (nx, ny), û_grid[nx, ny, :])
    end
    v .= reshape(irfft(v̂_grid, size(u_grid, 1)), :)
    return v
end

Lx = 2.5
Ly = 5.0
Nx = 40
Ny = 80

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

u = rand(Float64, size(ℒ₁, 2))
v_exp = ℒ₁ * u
v_act = ℒ₂ * u
v_act2 = similar(v_act)
mul_rfft!(v_act2, ℒ₂, u)
