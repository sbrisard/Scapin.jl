using Scapin.Grid
using Scapin.Elasticity

"""
    modal_strain_displacement!(B, n, N, h)

Compute the modal strain-displacement vector `B` in place, for the spatial frequency `n`.

The grid is defined by its size, `N` and spacing, `h`. The spatial frequency is
defined by `n ∈ CartesianIndices(1:N[1], …, 1:N[d])`.

See Sec. [The finite element discretization](@ref) for a definition of the modal
strain-displacement vector.

"""
function modal_strain_displacement!(
    B::AbstractVector{Complex{T}},
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

    b = -2sin(∑α) + 2im * cos(∑α)
    if d == 2
        B[1] = b * s[1] * c[2]
        B[2] = b * c[1] * s[2]
    elseif d == 3
        B[1] = b * s[1] * c[2] * c[3]
        B[2] = b * c[1] * s[2] * c[3]
        B[3] = b * c[1] * c[2] * s[3]
    else
        throw(ArgumentError("not implemented"))
    end
end


"""
    modal_strain_displacement(n, N, h)

Return the modal strain-displacement vector `B`.

See [`modal_strain_displacement!`](@ref) for a description of the parameters.
"""
function modal_strain_displacement(
    n::CartesianIndex{d},
    N::NTuple{d,Int},
    h::NTuple{d,T},
) where {d,T<:Number}
    B = Array{Complex{T}}(undef, d)
    modal_strain_displacement!(B, n, N, h)
    return B
end


"""
    modal_stiffness!(K, n, N, h, C)

Compute the modal stiffness matrix `K` in place, for the spatial frequency `n`.

The grid is defined by its size, `N` and spacing, `h`. The spatial frequency is
defined by `n ∈ CartesianIndices(1:N[1], …, 1:N[d])`. The (homogeneous)
constitutive material is specified by `C`.

See Sec. [The finite element discretization](@ref) for a definition of the modal
stiffness matrix.

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
    K = Array{Complex{T}}(undef, d, d)
    modal_stiffness!(K, n, N, h, C)
    return K
end
