using Scapin.Grid
using Scapin.Elasticity

"""
    modal_strain_displacement!(B, n, N, h)

Compute the modal strain-displacement vector `B` in place, for the spatial frequency `n`.

The grid is defined by its size, `N` and spacing, `h`. The spatial frequency is
defined by `n ∈ CartesianIndices(1:N[1], …, 1:N[d])`.

See “[The finite element discretization](@ref)” for a definition of the modal
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

Return the modal strain-displacement vector `B` for the spatial frequency `k`.

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
    modal_stiffness!(K, k, N, h, C)

Compute the modal stiffness matrix `K` in place, for the spatial frequency `k`.
The grid is defined by its size, `N` and spacing, `h`; `C` is the material.

"""
function modal_stiffness!(K, k, N, h, C)
    # In the notation of [Bri17, see Eq. (B.17)]
    #
    # φ[i] = φ(zᵢ) / hᵢ
    # χ[i] = χ(zᵢ) * hᵢ
    # ψ[i] = ψ(zᵢ)
    #
    # Which simplifies the expression of Hₖ (there are no hᵢs).
    dim = size(k, 1)

    β = 2π .* k ./ N
    φ = 2 .* (1 .- cos.(β)) ./ (h .^ 2)
    χ = (2 .+ cos.(β)) ./ 3
    ψ = sin.(β) ./ h

    scaling = C.μ / (1 - 2 * C.ν)
    if dim == 2
        H₁₁ = φ[1] * χ[2]
        H₂₂ = χ[1] * φ[2]
        K_diag = C.μ * (H₁₁ + H₂₂)
        K[1, 1] = scaling * H₁₁ + K_diag
        K[1, 2] = scaling * ψ[1] * ψ[2]
        K[2, 1] = K[1, 2]
        K[2, 2] = scaling * H₂₂ + K_diag
    elseif dim == 3
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
    modal_stiffness!(K, k, N, h, C)

Compute the modal stiffness matrix `K`, for the spatial frequency `k`.
The grid is defined by its size, `N` and spacing, `h`; `C` is the material.

"""
function modal_stiffness(k, N, h::AbstractVector{T}, C::Hooke{T,DIM}) where {T<:Number,DIM}
    K = Array{Complex{T}}(undef, DIM, DIM)
    modal_stiffness!(K, k, N, h, C)
    return K
end
