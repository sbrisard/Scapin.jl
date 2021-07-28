"""
Isotropic, linear elastic material.

    Hooke{T, DIM}(μ::T, ν::T)

Create a new instance with shear modulus `μ` and Poisson ratio `ν`.

!!! note "Material stability"

    Material stability requires that `μ > 0` and `-1 < ν < 1/2`; these conditions are *not*
    enforced here. In other words, *unstable* materials *can* be defined.

!!! tip "Plane stresses vs. plane strains"

    In the current implementation, `DIM = 2` refers to plane strain elasticity. For plane
    stresses, the *true* Poisson ratio `ν` should be replaced with the *fictitious* ratio
    `ν̃ = ν / (1 + ν)`.
"""
struct Hooke{T, DIM} <: AbstractMatrix{T}
    "The shear modulus."
    μ::T
    "The Poisson ratio."
    ν::T
    """
    Lamé I coefficient, such that the constitutive law reads

        σ = λ⋅tr(ε) + 2μ⋅ε,

    regardless of the dimensionality `DIM`.
    """
    λ::Float64
    Hooke{T, DIM}(μ::T, ν::T) where {T, DIM} = new(μ, ν, 2μ * ν / (1 - 2ν))
end

Base.size(::Hooke{T, DIM}) where {T, DIM} = ((DIM * (DIM + 1)) ÷ 2, (DIM * (DIM + 1)) ÷ 2)

function Base.getindex(C::Hooke{T, DIM}, i::Int, j::Int) where {T, DIM}
    value = zero(T)
    (i == j) && (value += 2C.μ)
    (i <= DIM) && (j <= DIM) && (value += C.λ)
    return value
end

function Base.:*(C::Hooke{T, DIM}, ε::AbstractVector{T}) where {T, DIM}
    tr_ε = sum(ε[1:DIM])
    σ = 2C.μ * ε
    σ[1:DIM] .+= C.λ * tr_ε
    return σ
end

"""
    bulk_modulus(C::Hooke)

Return the bulk modulus `κ` for the specified Hooke material.

For plane strain elasticity

```
κ = μ / (1 - 2ν),
```

and, for 3D elasticity

```
κ = 2/3 μ (1 + ν) / (1 - 2ν).
```
"""
function bulk_modulus(::Hooke) end

bulk_modulus(C::Hooke{T, 2}) where {T} = C.μ / (1-2C.ν)
bulk_modulus(C::Hooke{T, 3}) where {T} = 2C.μ * (1+C.ν)/3 / (1-2C.ν)
