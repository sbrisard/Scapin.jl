using CairoMakie
using FFTW
using LinearAlgebra
using SciMLOperators

struct LaplaceOperatorAsFunction
    Lx::Float64
    Ly::Float64
    Nx::Int
    Ny::Int
    dx::Float64
    dy::Float64
    ncells::Int
    function LaplaceOperatorAsFunction(Lx, Ly, Nx, Ny)
        dx = Lx / Nx
        dy = Ly / Ny
        ncells = Nx * Ny
        new(Lx, Ly, Nx, Ny, dx, dy, ncells)
    end
end

function apply!(y, Δ::LaplaceOperatorAsFunction, x)
    x_arr = reshape(x, Δ.Nx, Δ.Ny)
    y = similar(x)
    y_arr = reshape(y, Δ.Nx, Δ.Ny)
    for j₀ ∈ 1:Δ.Ny
        j₋₁ = j₀ == 1 ? Δ.Ny : j₀ - 1
        j₊₁ = j₀ == Δ.Ny ? 1 : j₀ + 1
        for i₀ ∈ 1:Δ.Nx
            i₋₁ = i₀ == 1 ? Δ.Nx : i₀ - 1
            i₊₁ = i₀ == Δ.Nx ? 1 : i₀ + 1
            y_arr[i₀, j₀] = (
                (x_arr[i₊₁, j₀] - 2 * x_arr[i₀, j₀] + x_arr[i₋₁, j₀]) / Δ.dx^2 +
                (x_arr[i₀, j₊₁] - 2 * x_arr[i₀, j₀] + x_arr[i₀, j₋₁]) / Δ.dy^2
            )
        end
    end
    y
end

struct LaplaceOperatorAsFunctionFourier
    Lx::Float64
    Ly::Float64
    Nx::Int
    Ny::Int
    dx::Float64
    dy::Float64
    ncells::Int
    kx²::Vector{Float64}
    ky²::Vector{Float64}
    function LaplaceOperatorAsFunctionFourier(Lx, Ly, Nx, Ny)
        dx = Lx / Nx
        dy = Ly / Ny
        ncells = Nx * Ny
        kx² = [(2 * sin(π * n / Nx) / dx)^2 for n = 0:(Nx-1)]
        ky² = [(2 * sin(π * n / Ny) / dy)^2 for n = 0:(Ny-1)]
        new(Lx, Ly, Nx, Ny, dx, dy, ncells, kx², ky²)
    end
end

function apply!(y, Δ::LaplaceOperatorAsFunctionFourier, x)
    x_arr = fft(reshape(x, Δ.Nx, Δ.Ny))
    for nx = 1:Δ.Nx
        for ny = 1:Δ.Ny
            x_arr[nx, ny] *= -(Δ.kx²[nx] + Δ.ky²[ny])
        end
    end
    y_arr = reshape(y, Δ.Nx, Δ.Ny)
    y_arr .= real.(ifft(x_arr))
    reshape(y_arr, Δ.ncells)
end

laplace = LaplaceOperatorAsFunction(1.0, 2.0, 4, 8)
m = n = laplace.Nx * laplace.Ny

ux = [sin(π * (n + 0.5) / laplace.Nx) for n = 0:(laplace.Nx-1)]
uy = [sin(π * (n + 0.5) / laplace.Ny) for n = 0:(laplace.Ny-1)]
u = ux .* uy'
u_vec = reshape(u, laplace.ncells)

x = zeros(Float64, n)
y = zeros(Float64, m)
Δ = FunctionOperator(
    (v, u, p, t) -> apply!(v, laplace, u),
    u_vec,
    u_vec;
    isinplace = true,
    issymmetric = true,
    isposdef = true,
)

laplace_fft =
    LaplaceOperatorAsFunctionFourier(laplace.Lx, laplace.Ly, laplace.Nx, laplace.Ny)
Δ = FunctionOperator(
    (v, u, p, t) -> apply!(v, laplace, u),
    u_vec,
    u_vec;
    isinplace = true,
    issymmetric = true,
    isposdef = true,
)

Δ_fft = FunctionOperator(
    (v, u, p, t) -> apply!(v, laplace_fft, u),
    u_vec,
    u_vec;
    isinplace = true,
    issymmetric = true,
    isposdef = true,
)

Δu = reshape(Δ * u_vec, laplace.Nx, laplace.Ny)
Δu_fft = reshape(Δ_fft * u_vec, laplace.Nx, laplace.Ny)

@assert all(isapprox.(Δu, Δu_fft, atol = 1e-15, rtol = 1e-12))
