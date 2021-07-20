module Scapin

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
end

"""

    Base.ndims(mat::Hooke{DIM})

Return the number of dimensions of the space over which the constitutive law operates.
"""
Base.ndims(mat::Hooke{DIM}) where {DIM} = DIM

end
