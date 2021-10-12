module Bri17

using LinearAlgebra
using Scapin.Grid
using Scapin.Elasticity

"""
    modal_strain_displacement!(b, n, N, h)

Compute the modal strain-displacement vector `b` in place, for the spatial frequency `n`.

This function returns `b`.

The grid is defined by its size, `N` and spacing, `h`. The spatial frequency is
defined by `n ∈ CartesianIndices(1:N[1], …, 1:N[d])`.

"""
function modal_strain_displacement!(
    b::AbstractVector{Complex{T}},
    n::CartesianIndex{d},
    N::NTuple{d,Int},
    h::NTuple{d,T},
) where {d,T<:Number}
    ∑α = zero(T)
    c = Array{T}(undef, d)
    s = Array{T}(undef, d)
    for i = 1:d
        α = π * (n[i] - 1) / N[i]
        ∑α += α
        c[i] = cos(α)
        s[i] = sin(α) / h[i]
    end

    b₀ = -2sin(∑α) + 2im * cos(∑α)
    if d == 2
        b[1] = b₀ * s[1] * c[2]
        b[2] = b₀ * c[1] * s[2]
    elseif d == 3
        b[1] = b₀ * s[1] * c[2] * c[3]
        b[2] = b₀ * c[1] * s[2] * c[3]
        b[3] = b₀ * c[1] * c[2] * s[3]
    else
        throw(ArgumentError("not implemented"))
    end
    return b
end


"""
    modal_strain_displacement(n, N, h)

Return the modal strain-displacement vector `b`.

See [`modal_strain_displacement!`](@ref) for a description of the parameters.

"""
function modal_strain_displacement(
    n::CartesianIndex{d},
    N::NTuple{d,Int},
    h::NTuple{d,T},
) where {d,T<:Number}
    return modal_strain_displacement!(Array{Complex{T}}(undef, d), n, N, h)
end


"""
    modal_stiffness!(K, n, N, h, C)

Compute the modal stiffness matrix `K` in place, for the spatial frequency `n`.

This function returns `K`.

The grid is defined by its size, `N` and spacing, `h`. The spatial frequency is
defined by `n ∈ CartesianIndices(1:N[1], …, 1:N[d])`. The (homogeneous)
constitutive material is specified by `C`.

!!! danger "Scaling of the modal stiffness matrix"

    The modal stiffness matrix introduced above differs from the matrix initially
    introduced in Ref. [DOI:10.1002/nme.5263](https://doi.org/10.1002/nme.5263)
    by a factor `h[1] * … * h[d]`. More precisely,

    ```
    K̂{Scapin} = h[1] * … * h[d] * K̂{10.1002/nme.5263}
    ```

    Such rescaling makes the relation between modal and nodal stiffness operators
    more natural (the former is the DFT of the latter).

"""
function modal_stiffness!(
    K::AbstractMatrix{Complex{T}},
    n::CartesianIndex{d},
    N::NTuple{d,Int},
    h::NTuple{d,T},
    C::Hooke{T,d},
) where {d,T<:Number}
    # In the notation of [Bri17, see Eq. (B.17)]
    #
    # φ[i] = φ(zᵢ) * hᵢ
    # χ[i] = χ(zᵢ) * hᵢ²
    # ψ[i] = ψ(zᵢ) * hᵢ
    #
    # Which simplifies the expression of Hₖ and allows to express the
    # nodal stiffness matrix as the inverse DFT of the modal stiffness
    # matrix.

    φ = Array{T}(undef, d)
    χ = Array{T}(undef, d)
    ψ = Array{T}(undef, d)
    for i = 1:d
        β = 2π * (n[i] - 1) / N[i]
        φ[i] = 2 * (1 - cos(β)) / h[i]
        χ[i] = (2 + cos(β)) / 3 * h[i]
        ψ[i] = sin(β)
    end

    scaling = C.μ / (1 - 2 * C.ν)
    if d == 2
        H₁₁ = φ[1] * χ[2]
        H₂₂ = χ[1] * φ[2]
        K_diag = C.μ * (H₁₁ + H₂₂)
        K[1, 1] = scaling * H₁₁ + K_diag
        K[1, 2] = scaling * ψ[1] * ψ[2]
        K[2, 1] = K[1, 2]
        K[2, 2] = scaling * H₂₂ + K_diag
    elseif d == 3
        H₁₁ = φ[1] * χ[2] * χ[3]
        H₂₂ = χ[1] * φ[2] * χ[3]
        H₃₃ = χ[1] * χ[2] * φ[3]
        K_diag = C.μ * (H₁₁ + H₂₂ + H₃₃)
        K[1, 1] = scaling * H₁₁ + K_diag
        K[1, 2] = scaling * ψ[1] * ψ[2] * χ[3]
        K[1, 3] = scaling * ψ[1] * χ[2] * ψ[3]
        K[2, 1] = K[1, 2]
        K[2, 2] = scaling * H₂₂ + K_diag
        K[2, 3] = scaling * χ[1] * ψ[2] * ψ[3]
        K[3, 1] = K[1, 3]
        K[3, 2] = K[2, 3]
        K[3, 3] = scaling * H₃₃ + K_diag
    else
        throw(ArgumentError("not implemented"))
    end
    return K
end

"""
    modal_stiffness(n, N, h, C)

Return the modal stiffness matrix `K`.

See [`modal_stiffness!`](@ref) for a description of the parameters.

"""
function modal_stiffness(
    n::CartesianIndex{d},
    N::NTuple{d,Int},
    h::NTuple{d,T},
    C::Hooke{T,d},
) where {T<:Number,d}
    return modal_stiffness!(Array{Complex{T}}(undef, d, d), n, N, h, C)
end

struct DiscreteGreenOperatorBri17{T,d}
    C::Hooke{T,d}
    N::NTuple{d,Int}
    h::NTuple{d,T}
end

function apply!(
    ε̂::AbstractVector{Complex{T}},
    Γ̂::DiscreteGreenOperatorBri17{T, d},
    τ̂::AbstractVector{Complex{T}},
    n::CartesianIndex{d},
) where {T<:Number,d}
    if all(Tuple(n) .== 1)
        ε̂ .= zero(eltype(ε̂))
        return
    end
    s = one(T) / √(2 * one(T))
    b̂ = modal_strain_displacement(n, Γ̂.N, Γ̂.h)
    conj_b̂ = conj.(b̂)
    K̂ = modal_stiffness(n, Γ̂.N, Γ̂.h, Γ̂.C)
    # TODO: use static arrays here!
    τ̂_conj_b̂ = Vector{Complex{T}}(undef, d)
    if d == 2
        #     ┌               ┐
        # τ̂ = │ τ̂₁       τ̂₃/√2│
        #     │ τ̂₃/√2    τ̂₂   │
        #     └               ┘
        τ̂_conj_b̂[1] = τ̂[1] * conj_b̂[1] + s * τ̂[3] * conj_b̂[2]
        τ̂_conj_b̂[2] = τ̂[2] * conj_b̂[2] + s * τ̂[3] * conj_b̂[1]
    elseif d == 3
        #     ┌                        ┐
        #     │ τ̂₁       τ̂₆/√2    τ̂₅/√2│
        # τ̂ = │ τ̂₆/√2    τ̂₂       τ̂₄/√2│
        #     │ τ̂₅/√2    τ̂₄/√2    τ̂₃/√2│
        #     └                        ┘
        τ̂_conj_b̂[1] = τ̂[1] * conj_b̂[1] + s * (τ̂[6] * conj_b̂[2] + τ̂[5] * conj_b̂[3])
        τ̂_conj_b̂[2] = τ̂[2] * conj_b̂[2] + s * (τ̂[6] * conj_b̂[1] + τ̂[4] * conj_b̂[3])
        τ̂_conj_b̂[3] = τ̂[3] * conj_b̂[3] + s * (τ̂[5] * conj_b̂[1] + τ̂[4] * conj_b̂[2])
    end
    # Do not include the `-` sign!
    û = prod(Γ̂.h) .* (K̂ \ τ̂_conj_b̂)
    if d == 2
        ε̂[1] = û[1] * b̂[1]
        ε̂[2] = û[2] * b̂[2]
        ε̂[3] = s * (û[1] * b̂[2] + û[2] * b̂[1])
    elseif d == 3
        ε̂[1] = û[1] * b̂[1]
        ε̂[2] = û[2] * b̂[2]
        ε̂[3] = û[3] * b̂[3]
        ε̂[4] = s * (û[2] * b̂[3] + û[3] * b̂[2])
        ε̂[5] = s * (û[3] * b̂[1] + û[1] * b̂[3])
        ε̂[6] = s * (û[1] * b̂[2] + û[2] * b̂[1])
    end
    return ε̂
end

function apply(
    Γ̂::DiscreteGreenOperatorBri17{T, d},
    τ̂::AbstractVector{Complex{T}},
    n::CartesianIndex{d},
) where {T<:Number,d}
    return apply!(Array{Complex{T}}(undef, div(d * (d + 1), 2)), Γ̂, τ̂, n)
end

export modal_strain_displacement!,
    modal_strain_displacement,
    modal_stiffness!,
    modal_stiffness,
    DiscreteGreenOperatorBri17,
    apply!,
    apply
end
