using Scapin.Grid

"""
    modal_strain_displacement!(B, k, N, h)

Compute the modal strain-displacement vector `B` in place, for the spatial frequency `k`.
The grid is defined by its size, `N` and spacing, `h`.
"""
function modal_strain_displacement!(B, k, N, h)
    dim = size(B, 1)
    ∑α = zero(T)

    for i in eachindex(k)
        α = π * k[i] / N[i]
        ∑α += α
        cos_α[i] = cos(α)
        sin_α[i] = sin(α) / h[i]
    end

    prefactor = -2sin(∑α) + 2im * cos(∑α)
    if dim == 2
        B[1] = prefactor * s[1] * c[2]
        B[2] = prefactor * c[1] * s[2]
    elseif dim == 3
        B[1] = prefactor * s[1] * c[2] * c[3]
        B[2] = prefactor * c[1] * s[2] * c[3]
        B[3] = prefactor * c[1] * c[2] * s[3]
    else
        throw(ArgumentError("not implemented"))
    end
end


"""
    modal_strain_displacement(k, N, h)

Compute the modal strain-displacement vector `B` for the spatial frequency `k`.
The grid is defined by its size, `N` and spacing, `h`.
"""
function modal_strain_displacement(k, N, h)
    B = similar(k, Complex{eltype(k)})
    modal_strain_displacement!(B, k, N, h)
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
    φ = 2 .* (1 .- cos(β)) ./ (h .^ 2)
    χ = (2 .+ cos(β)) ./ 3
    ψ = sin(β) ./ h

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
function modal_stiffness(k, N, h, C)
    K = similar(k, Complex{eltype(k)})
    modal_stiffness!(K, k, N, h, C)
    return K
end
