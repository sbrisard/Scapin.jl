module Scapin

using Base.Iterators
using FFTW
using LinearAlgebra

"""
    dimensionality(op)

Return the number of dimensions of the physical space that `op` operates on.
"""
function dimensionality end


"""
    grid_size(ℱ)

Return the size of the grid the discrete operator `ℱ` operates on.
"""
function grid_size(ℱ) end


"""
    apply_fourier!(ŷ, ℱ, k, x̂) -> ŷ

Apply in-place operator `ℱ` to `x̂` in Fourier space and return the modified array `ŷ`.
More specifically,

    ŷ[:] = ℱ̂(k)⋅x̂,

where `k`, `x̂` and `ŷ` are column vectors: `k` is the vector of spatial
frequencies and `x̂` (resp. `ŷ`) is the Fourier mode of the input (resp. output),
for the spatial frequency `k`. The following must hold

    size(k, 1) == dimensionality(ℱ),
    size(x̂, 1) == size(ℱ, 2),
    size(ŷ, 1) == size(ℱ, 1).

The input vector `x̂` is an array of type `T` or `complex(T)`, where
`T == eltype(ℱ)`. The output vector `ŷ` can always be of type `complex(T)`. When
the Fourier coefficient `ℱ̂(k)` of the operator is real, it might make sense to
allow for `ŷ` to be of type `real(T)`. If this is disallowed, calling this
function should raise an exception.

!!! note "Continuous and discrete operators"

    For continuous operators, `k` would typically be a real-valued vector of wave
    numbers, while for discrete operators, `k` would be a integer-valued vector
    of indices.

!!! danger "Overriding the input vector"

    The present method *must* be implemented in such a way that aliasing `ŷ`
    with `x̂` is permitted. In other words, calling `apply_fourier!(x̂, ℱ, k, x̂)`
    *must* always deliver the correct answer, provided that this operation is
    allowed type-wise (when `x̂` is of type `real(T)`.

"""
function apply_fourier! end


"""
    apply_fourier(ℱ, k, x̂)

Apply operator `ℱ` to `x̂` in Fourier space and return the result. The returned
vector is promoted from `complex(T)`.

See [apply_fourier!](@ref)

"""
function apply_fourier(ℱ, k, x̂)
    T = promote_type(eltype(ℱ), eltype(x̂), eltype(k))
    ŷ = Array{T}(undef, size(ℱ, 1))
    return apply_fourier!(ŷ, ℱ, k, x̂)
end

"""
    fourier_matrix!(F̂, ℱ, k) -> F̂


Compute in-place the `k`-th mode of `ℱ`, `ℱ̂(k)`, as a matrix `F̂`.

In general, the coefficients of `F̂` should be s of type
`complex(eltype(ℱ))`. However, the Fourier coefficients of the operator ℱ might
be real, in which case, `eltype(F̂) == real(eltype(ℱ))` might be allowed.

!!! note "Performance of the default implementation"

    The default implementation relies on [apply_fourier!](@ref) to build the
    matrix column-by-column. It might be inefficient.

"""
function fourier_matrix!(F̂, ℱ, k)
    nrows, ncols = size(ℱ)
    T = eltype(F̂)
    x̂ = zeros(T, ncols)
    for i = 1:ncols
        x̂[i] = one(T)
        apply_fourier!(view(F̂, :, i), ℱ, k, x̂)
        x̂[i] = zero(T)
    end
    return F̂
end


"""
    fourier_matrix(ℱ, k)

Return the `k`-th mode of `ℱ`, `ℱ̂(k)`, as a matrix.

The coefficients of the returned matrix are of type `complex(eltype(ℱ))`. For a
different output type, use [fourier_matrix!](@ref).
"""
function fourier_matrix(ℱ, k)
    return fourier_matrix!(zeros(complex(eltype(ℱ)), size(ℱ)...), ℱ, k)
end


"""
    apply(ℱ, x)

Return the grid `y` that results from applying the convolution operator `ℱ` to
the input grid `x`.

The sizes of `x` and `y` are

    size(x) == (N..., ncols)
    size(y) == (N..., nrows)

where

    N = grid_size(ℱ)

and

    (nrows, ncols) == size(ℱ).

"""
function apply(ℱ, x)
    𝒩 = CartesianIndices(grid_size(ℱ))
    d = dimensionality(ℱ)
    y = fft(x, 2:(d+1))
    for n ∈ 𝒩
        apply_fourier!(view(y, :, n), ℱ, n, y[:, n])
    end
    # TODO: use in-place FFT
    return ifft(y, 2:(d+1))
end


function apply(ℱ, x, plan)
    y = plan * x
    for n ∈ CartesianIndices(size(y)[2:end])
        apply_fourier!(view(y, :, n), ℱ, n, y[:, n])
    end
    # TODO: Use in-place FFT
    return plan \ y
end


export dimensionality, grid_size, apply_fourier!, apply_fourier, fourier_matrix!, fourier_matrix, apply

include("Conductivity.jl")
include("Elasticity.jl")
include("Grid.jl")
include("Brick.jl")
include("Bri17.jl")

end
