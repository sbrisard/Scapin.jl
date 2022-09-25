# Appendix

## On Fourier series

We consider the domain `Ω = (0, L[1]) × … × (0, L[d])` and a
`L`-periodic vector field `u: Ω → ℝᵐ`. The Fourier coefficients of `u`
are defined as follows, for all `k ∈ ℝᵈ`

    û(k) = |Ω|⁻¹ ∫ u(x) exp(-i k ⋅ x) dx.
                 Ω

The following *synthesis formula* then holds (under non-specified, but
non-restrictive, regularity conditions)

    u(x) =   ∑    û(kₚ) exp(i kₚ ⋅ x),
	       p ∈ ℤᵈ

where `kₚ` is the following wave-vector

            2π p[i]
    kₚ[i] = ───────,    for all    i ∈ {1, …, d}    and    p ∈ ℤᵈ.
	          L[i]
