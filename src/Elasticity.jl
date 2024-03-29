module Elasticity

using LinearAlgebra
using Scapin

"""
    Hooke{d,T}

Isotropic, linear elastic material.

- `d` — number of spatial dimensions
- `T` — scalar type
- `μ::T` — shear modulus
- `ν::T` — Poisson ratio
- `λ::T` — first Lamé coefficient

Regardless of the number of space dimensions, `d`, the stress-strain relationship reads

    σ = λ⋅tr(ε)I + 2μ⋅ε.

!!! note "Material stability"

    Material stability requires that `μ > 0` and `-1 < ν < 1/2`; these conditions are *not*
    enforced here. In other words, *unstable* materials *can* be defined.

!!! tip "Plane stresses vs. plane strains"

    In the current implementation, `d = 2` refers to plane strain elasticity. For plane
    stresses, the *true* Poisson ratio `ν` should be replaced with the *fictitious* ratio
    `ν̃ = ν / (1 + ν)`.
"""
struct Hooke{d,T}
    μ::T
    ν::T
    λ::T
    Hooke{d,T}(μ::T, ν::T) where {d,T} = new(μ, ν, 2μ * ν / (1 - 2ν))
end

"""
    Hooke{d,T}(μ::T, ν::T)

Create a new instance of `Hooke{d,T}` with shear modulus `μ` and Poisson ratio `ν`.
"""
Hooke{d}(μ::T, ν::T) where {d,T} = Hooke{d,T}(μ, ν)

Base.eltype(::Type{Hooke{d,T}}) where {d,T} = T

Base.ndims(::Type{Hooke{d,T}}) where {d,T} = 2
Base.ndims(::Hooke{d,T}) where {d,T} = 2

Base.size(::Hooke{d,T}) where {d,T} = ((d * (d + 1)) ÷ 2, (d * (d + 1)) ÷ 2)

function Base.size(::Hooke{d,T}, n::Int) where {d,T}
    ((n <= 0) || (n > 2)) && throw(ErrorException("dimension must be 1 or 2 (got $n)"))
    (d * (d + 1)) ÷ 2
end

Scapin.dimensionality(::Type{Hooke{d, T}}) where {d, T} = d

function Base.convert(::Type{Array}, C::Hooke{d,T}) where {d,T}
    C₆₆ = 2C.μ
    C₁₁ = C.λ + C₆₆
    C₁₂ = C.λ
    if d == 2
        return [
            C₁₁ C₁₂ 0
            C₁₂ C₁₁ 0
            0 0 C₆₆
        ]
    elseif d == 3
        return [
            C₁₁ C₁₂ C₁₂ 0 0 0
            C₁₂ C₁₁ C₁₂ 0 0 0
            C₁₂ C₁₂ C₁₁ 0 0 0
            0 0 0 C₆₆ 0 0
            0 0 0 0 C₆₆ 0
            0 0 0 0 0 C₆₆
        ]
    end
end

function Base.getindex(C::Hooke{d,T}, i::Int, j::Int) where {d,T}
    value = zero(T)
    (i == j) && (value += 2C.μ)
    (i <= d) && (j <= d) && (value += C.λ)
    return value
end

function mul!(σ::AbstractVector{T}, C::Hooke{d,T}, ε::AbstractVector{T}) where {d,T}
    tr_ε = sum(ε[1:d])
    σ = 2C.μ * ε
    σ[1:d] .+= C.λ * tr_ε
    return σ
end

function Base.:*(C::Hooke{d,T}, ε::AbstractVector{T}) where {d,T}
    return mul!(AbstractVector{T}(undef, size(C, 1)), C, ε)
end

"""
    bulk_modulus(C::Hooke)

Return the bulk modulus `κ` for the specified Hooke material.

For plane strain elasticity, `κ = μ / (1 - 2ν)` and, for 3d elasticity
`κ = 2/3 μ (1 + ν) / (1 - 2ν)`.

"""
function bulk_modulus(::Hooke) end

bulk_modulus(C::Hooke{2,T}) where {T} = C.μ / (1 - 2C.ν)
bulk_modulus(C::Hooke{3,T}) where {T} = 2C.μ * (1 + C.ν) / 3 / (1 - 2C.ν)

"""
Continuous Green operator related to [Hooke](@ref) materials.

"""
struct GreenOperatorHooke{d,T}
    C::Hooke{d,T}
end

"""
   eltype(type)

Return the type of the (scalar) elements the operator of given `type` operates
on in the *real* space (as opposed to the *Fourier* space). Note that this type
may be complex!
"""
Base.eltype(::Type{GreenOperatorHooke{d,T}}) where {d,T} = T

Base.ndims(::Type{GreenOperatorHooke{d,T}}) where {d,T} = 2
Base.ndims(::GreenOperatorHooke{d,T}) where {d,T} = 2

Base.size(::GreenOperatorHooke{d,T}) where {d,T} = ((d * (d + 1)) ÷ 2, (d * (d + 1)) ÷ 2)

function Base.size(::GreenOperatorHooke{d,T}, n::Int) where {d,T}
    ((n <= 0) || (n > 2)) && throw(ErrorException("dimension must be 1 or 2 (got $n)"))
    (d * (d + 1)) ÷ 2
end

Scapin.dimensionality(::Type{GreenOperatorHooke{d, T}}) where {d, T} = d

function Scapin.apply_fourier!(ε̂, Γ::GreenOperatorHooke{2,T}, k, τ̂) where {T}
    k² = sum(abs2, k)
    if iszero(k²)
        # TODO: improve robustness of test
        ε̂[:] .= zero(eltype(ε̂))
    else
    τ̂k₁ = τ̂[1] * k[1] + τ̂[3] * k[2] / sqrt(2 * one(T))
    τ̂k₂ = τ̂[2] * k[2] + τ̂[3] * k[1] / sqrt(2 * one(T))
    nτ̂n = (k[1] * τ̂k₁ + k[2] * τ̂k₂) / k²
    const1 = nτ̂n / (1 - Γ.C.ν)
    const2 = 1 / (2 * Γ.C.μ * k²)
    ε̂[1] = const2 * (k[1] * (2 * τ̂k₁ - const1 * k[1]))
    ε̂[2] = const2 * (k[2] * (2 * τ̂k₂ - const1 * k[2]))
    const3 = sqrt(2 * one(T)) * const2
    ε̂[3] = const3 * (k[1] * τ̂k₂ + k[2] * τ̂k₁ - const1 * k[1] * k[2])
    end
    return ε̂
end

function Scapin.apply_fourier!(ε̂, Γ::GreenOperatorHooke{3,T}, k, τ̂) where {T}
    k² = sum(abs2, k)
    τ̂k₁ = τ̂[1] * k[1] + (τ̂[6] * k[2] + τ̂[5] * k[3]) / sqrt(2 * one(T))
    τ̂k₂ = τ̂[2] * k[2] + (τ̂[6] * k[1] + τ̂[4] * k[3]) / sqrt(2 * one(T))
    τ̂k₃ = τ̂[3] * k[3] + (τ̂[5] * k[1] + τ̂[4] * k[2]) / sqrt(2 * one(T))
    nτ̂n = (k[1] * τ̂k₁ + k[2] * τ̂k₂ + k[3] * τ̂k₃) / k²
    const1 = nτ̂n / (1 - Γ.C.ν)
    const2 = 1 / (2 * Γ.C.μ * k²)
    ε̂[1] = const2 * (k[1] * (2 * τ̂k₁ - const1 * k[1]))
    ε̂[2] = const2 * (k[2] * (2 * τ̂k₂ - const1 * k[2]))
    ε̂[3] = const2 * (k[3] * (2 * τ̂k₃ - const1 * k[3]))
    const3 = sqrt(2 * one(T)) * const2
    ε̂[4] = const3 * (k[2] * τ̂k₃ + k[3] * τ̂k₂ - const1 * k[2] * k[3])
    ε̂[5] = const3 * (k[3] * τ̂k₁ + k[1] * τ̂k₃ - const1 * k[3] * k[1])
    ε̂[6] = const3 * (k[1] * τ̂k₂ + k[2] * τ̂k₁ - const1 * k[1] * k[2])
    return ε̂
end


# function ms94_frequencies(N, L) where {T}
#     # TODO Ensure that π and L are of the same type
#     [(2π / L) * (2n < N ? n : n - N) for n = 0:N-1]
# end

# struct TruncatedGreenOperator{T,d}
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

export Hooke, bulk_modulus, GreenOperatorHooke

end
