module Scapin

using Base.Iterators
using LinearAlgebra

"""
    dimensionality(op)

Return the number of dimensions of the physical space that `op` operates on.
"""
function dimensionality end

"""
   eltype_real(type)

Return the type of the (scalar) elements the operator of given `type` operates
on in the *real* space (as opposed to the *Fourier* space). Note that this type
may be complex!

The definition `eltype_real(x) = eltype_real(typeof(x))` is provided for
convenience so that instances can be passed instead of types. However the form
that accepts a type argument should be defined for new types.
"""
function eltype_real end
eltype_real(op) = eltype_real(typeof(op))

"""
   eltype_fourier(type)

Return the type of the (scalar) elements the operator of given `type` operates
on in the *Fourier* space (as opposed to the *real* space). Note that this type
may be real!

The definition `eltype_fourier(x) = eltype_fourier(typeof(x))` is provided for
convenience so that instances can be passed instead of types. However the form
that accepts a type argument should be defined for new types.
"""
function eltype_fourier end
eltype_fourier(op) = eltype_fourier(typeof(op))

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

The input and output vectors `x̂` and `ŷ`  are arrays of type `T` or `Complex{T}`,
where `T = eltype(ℱ)`.

!!! note "Continuous and discrete operators"

    For continuous operators, `k` would typically be a real-valued vector of wave
    numbers, while for discrete operators, `k` would be a integer-valued vector
    of indices.

!!! danger "Overriding the input vector"

    The present method *must* be implemented in such a way that aliasing `ŷ` with
    `x̂` is permitted. In other words, calling `apply_fourier!(x̂, ℱ, x̂, k)` *must*
    always deliver the correct answer.

"""
function apply_fourier! end


"""
    apply_fourier(ℱ, x̂, k)

Apply in-place operator `ℱ` to `x̂` in Fourier space and return the result.

See [apply_fourier!](@ref)

"""
function apply_fourier(ℱ, x̂, k)
    T = promote_type(eltype(ℱ), eltype(x̂), eltype(k))
    ŷ = Array{T}(undef, size(ℱ, 1))
    return apply_fourier!(ŷ, ℱ, x̂, k)
end

"""
    fourier_matrix(ℱ, k)

Return the `k`-th mode of `ℱ`, `ℱ̂(k)`, as a matrix.

!!! note "Performance of the default implementation"

    The default implementation relies on [apply_fourier!](@ref) to build the
    matrix column-by-column. It might be inefficient.
"""
function fourier_matrix(ℱ, k)
    nrows, ncols = size(ℱ)
    T = eltype_fourier(ℱ)
    mat = zeros(T, nrows, ncols)
    x̂ = zeros(T, ncols)
    for i = 1:ncols
        x̂[i] = one(T)
        apply_fourier!(view(mat, :, i), ℱ, k, x̂)
        x̂[i] = zero(T)
    end
    return mat

end


export dimensionality, eltype_real, eltype_fourier, apply_fourier!, apply_fourier, fourier_matrix

include("Elasticity.jl")
include("Grid.jl")
include("Brick.jl")
include("Bri17.jl")

end
