"""
Isotropic, linear elastic material.

    Hooke{DIM}(μ::Float64, ν::Float64)

Create a new instance with shear modulus `μ` and Poisson ratio `ν`.

!!! note "Material stability"

    Material stability requires that `μ > 0` and `-1 < ν < 1/2`; these conditions are *not*
    enforced here. In other words, *unstable* materials *can* be defined.

!!! tip "Plane stresses vs. plane strains"

    In the current implementation, `DIM = 2` refers to plane strain elasticity. For plane
    stresses, the *true* Poisson ratio `ν` should be replaced with the *fictitious* ratio
    `ν̃ = ν / (1 + ν)`.
"""
struct Hooke{DIM}
    "The shear modulus."
    μ::Float64
    "The Poisson ratio."
    ν::Float64
    """
    Lamé I coefficient, such that the constitutive law reads

        σ = λ⋅tr(ε) + 2⋅μ⋅ε,

    regardless of the dimensionality `DIM`.
    """
    λ::Float64
    Hooke{DIM}(μ, ν) where {DIM} = new(μ, ν, 2μ * ν / (1 - 2ν))
end

"""

    Base.ndims(mat::Hooke{DIM})

Return the number of dimensions of the space over which the constitutive law operates.
"""
Base.ndims(mat::Hooke{DIM}) where {DIM} = DIM

Base.size(C::Hooke{2}) = (3, 3)
Base.size(C::Hooke{2}, d) = 3
Base.size(C::Hooke{3}) = (6, 6)
Base.size(C::Hooke{3}, d) = 6

function Base.:*(C::Hooke{DIM}, ε::AbstractVector{Float64}) where {DIM}
    tr_ε = sum(ε[1:DIM])
    σ = 2C.μ * ε
    σ[1:DIM] .+= C.λ * tr_ε
    return σ
end

