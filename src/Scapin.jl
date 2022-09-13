module Scapin

using Base.Iterators
using FFTW
using LinearAlgebra


"""
   AbstractSpatialOperator

This is the root type of all operators defined in this library.

Such operators map vector fields `x : ℝᵈ → ℝⁿ` to vector fields
`y : ℝᵈ → ℝᵐ`.

The number of spatial dimensions, `d`, is referred to as the
*dimensionality* of the operator – see [dimensionality](@ref).

The pair `(m, n) == (dimension of the input vector field, dimension of
the output vector field)` is referred to as the *size* of the operator
– see [size](@ref). Consequently, `ndims(F) == 2` for all
`F::AbstractSpatialOperator`.

"""
abstract type AbstractSpatialOperator end


"""
   AbstractDiscreteSpatialOperator

This abstract type defines operators that are spatially discretized.

Such operators map vector fields `x : ℝᵈ → ℝⁿ` to vector fields
`y : ℝᵈ → ℝᵐ`.

This type extends [AbstractSpatialOperator](@ref). As such, it
implements [dimensionality](@ref), [size](@ref) and [ndims](@ref).

Discrete operators operate on discrete vector fields that are
represented as arrays. The size of the *input* arrays is
`(n, N[1], …, N[d])`, where `N[i]` is the number of grid cells in
direction `i = 1, …, d`. Note that the first (fast) index of the array
is the component of the vector field. Similarly, the size of the
*output* arrays is `(m, N[1], …, N[d])`.

The tuple `N` is referred to as the *grid size*, see
[grid_size](@ref).

"""
abstract type AbstractDiscreteSpatialOperator <: AbstractSpatialOperator end


"""
   ndims(F::AbstractSpatialOperator) -> 2

Return the number of dimensions of `F` – see
[AbstractSpatialOperator](@ref). Since this is a linear operator that
maps fields onto fields, the number of dimensions is always 2.

"""
Base.ndims(::AbstractSpatialOperator) = 2

"""
    dimensionality(F)

Return the number of dimensions of the physical space that `F`
operates on – see [AbstractSpatialOperator](@ref).

The default implementation assumes that `dimensionality(typeof(F))`
has been defined.

"""
dimensionality(F::AbstractSpatialOperator) = dimensionality(typeof(F))

"""
    grid_size(F)

Return the size of the grid the discrete operator `F` operates on.

- Size of *input* arrays of `F`: `(size(F, 2), grid_size(F)...)`,
- Size of *output* arrays of `F`: `(size(F, 1), grid_size(F)...)`.

"""
function grid_size(F::AbstractDiscreteSpatialOperator) end


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
    d = dimensionality(ℱ)
    y = fft(x, 2:(d+1))
    for n ∈ CartesianIndices(grid_size(ℱ))
        # NOTE: `y[:, n]` as the last argument is not a mistake, since `y`
        # temporarily stores the DFT of `x`. Since `y[:, n]` is not passed as a
        # view, a copy of `y[:, n]` is passed as an argument, which means that
        # `apply_fourier!` does not operate on the same first and last
        # arguments.
        apply_fourier!(view(y, :, n), ℱ, n, y[:, n])
    end
    # TODO: use in-place FFT
    return ifft(y, 2:(d+1))
end


function apply(ℱ, x, plan)
    y = plan * x
    for n ∈ CartesianIndices(size(y)[2:end])
        apply_fourier!(view(y, :, n), ℱ, n, x[:, n])
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
include("BD12.jl")
include("ConvergenceAnalysis.jl")

end
