module Scapin

using Base.Iterators
using LinearAlgebra

"""
    dimensionality(op)

Return the number of dimensions of the physical space that `op` operates on.
"""
function dimensionality end


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


export dimensionality, apply_fourier!

include("Elasticity.jl")
include("Grid.jl")
include("Brick.jl")
include("Bri17.jl")

end
