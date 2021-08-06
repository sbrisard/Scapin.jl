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
struct Hooke{T,DIM} <: AbstractMatrix{T}
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
    Hooke{T,DIM}(μ::T, ν::T) where {T,DIM} = new(μ, ν, 2μ * ν / (1 - 2ν))
end

Base.size(::Hooke{T,DIM}) where {T,DIM} = ((DIM * (DIM + 1)) ÷ 2, (DIM * (DIM + 1)) ÷ 2)

function Base.getindex(C::Hooke{T,DIM}, i::Int, j::Int) where {T,DIM}
    value = zero(T)
    (i == j) && (value += 2C.μ)
    (i <= DIM) && (j <= DIM) && (value += C.λ)
    return value
end

function Base.:*(C::Hooke{T,DIM}, ε::AbstractVector{T}) where {T,DIM}
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

bulk_modulus(C::Hooke{T,2}) where {T} = C.μ / (1 - 2C.ν)
bulk_modulus(C::Hooke{T,3}) where {T} = 2C.μ * (1 + C.ν) / 3 / (1 - 2C.ν)

function block_apply!(out, hooke::Hooke{T,2}, k, τ) where {T}
    k² = sum(abs2, k)
    τk₁ = τ[1] * k[1] + τ[3] * k[2] / sqrt(2 * one(T))
    τk₂ = τ[2] * k[2] + τ[3] * k[1] / sqrt(2 * one(T))
    nτn = (k[1] * τk₁ + k[2] * τk₂) / k²
    const1 = nτn / (1 - hooke.ν)
    const2 = 1 / (2 * hooke.μ * k²)
    out[1] = const2 * (k[1] * (2 * τk₁ - const1 * k[1]))
    out[2] = const2 * (k[2] * (2 * τk₂ - const1 * k[2]))
    const3 = sqrt(2 * one(T)) * const2
    out[3] = const3 * (k[1] * τk₂ + k[2] * τk₁ - const1 * k[1] * k[2])
    return out
end

function block_apply!(out, hooke::Hooke{T,3}, k, τ) where {T}
    k² = sum(abs2, k)
    τk₁ = τ[1] * k[1] + (τ[6] * k[2] + τ[5] * k[3]) / sqrt(2 * one(T))
    τk₂ = τ[2] * k[2] + (τ[6] * k[1] + τ[4] * k[3]) / sqrt(2 * one(T))
    τk₃ = τ[3] * k[3] + (τ[5] * k[1] + τ[4] * k[2]) / sqrt(2 * one(T))
    nτn = (k[1] * τk₁ + k[2] * τk₂ + k[3] * τk₃) / k²
    const1 = nτn / (1 - hooke.ν)
    const2 = 1 / (2 * hooke.μ * k²)
    out[1] = const2 * (k[1] * (2 * τk₁ - const1 * k[1]))
    out[2] = const2 * (k[2] * (2 * τk₂ - const1 * k[2]))
    out[3] = const2 * (k[3] * (2 * τk₃ - const1 * k[3]))
    const3 = sqrt(2 * one(T)) * const2
    out[4] = const3 * (k[2] * τk₃ + k[3] * τk₂ - const1 * k[2] * k[3])
    out[5] = const3 * (k[3] * τk₁ + k[1] * τk₃ - const1 * k[3] * k[1])
    out[6] = const3 * (k[1] * τk₂ + k[2] * τk₁ - const1 * k[1] * k[2])
    return out
end

function block_matrix(op::Hooke{T,DIM}, k::AbstractVector{T}) where {T,DIM}
    nrows, ncols = size(op)
    mat = zeros(T, nrows, ncols)
    τ = zeros(T, ncols)
    for i = 1:ncols
        τ[i] = one(T)
        block_apply!(view(mat, :, i), op, k, τ)
        τ[i] = zero(T)
    end
    mat
end

# function ms94_frequencies(N, L) where {T}
#     # TODO Ensure that π and L are of the same type
#     [(2π / L) * (2n < N ? n : n - N) for n = 0:N-1]
# end

# struct TruncatedGreenOperator{T,DIM}
#     Γ::Hooke{T,DIM}
#     N::SVector{DIM,Int}
#     L::SVector{DIM,T}
#     # TODO: I would like k to be a SVector
#     k::AbstractArray{AbstractArray{T,1},1}
#     function TruncatedGreenOperator{T,DIM}(
#         Γ::Hooke{T,DIM},
#         N::SVector{DIM,Int},
#         L::SVector{DIM,T},
#     ) where {T,DIM}
#         k = [ms94_frequencies(N[i], L[i]) for i = 1:DIM]
#         new(Γ, N, L, k)
#     end
# end

# function block_apply!(out::AbstractArray{T, DIM}, Γ_h::TruncatedGreenOperator{T, DIM}, n::SVector{DIM, Int}, τ::AbstractArray{T, DIM}) where {T, DIM}
#     k = SVector [Γ_h.k[n[i]] for i=1:DIM]
#     block_apply!(out, Γ_h.Γ, k, τ)
#     return τ
# end

# function apply(out::AbstractArray{T, DIM+1}, Γ_h::TruncatedGreenOperator{T, DIM}, τ::AbstractArray{T, DIM+1}) where {T, DIM}
# end
