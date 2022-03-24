module Conductivity

using Scapin

"""
    Ohm{d,T}

Isotropic, linear conductor that follows Ohm's law `𝐉 = σ𝐄`

- `𝐉` — current density
- `𝐄` — electric field
- `d` — number of spatial dimensions
- `T` — scalar type
- `σ::T` — conductivity (scalar)

!!! note "Material stability"

    Material stability requires that `σ > 0`; this condition is *not* enforced
    here. In other words, *unstable* conductors *can* be defined.
"""
struct Ohm{d,T}
    σ::T
end

"""
    Ohm{d}(σ::T)

Create a new instance of `Ohm{d, T}` with conductivity `σ`.
"""
Ohm{d}(σ::T) where {d,T} = Ohm{d,T}(σ)

Base.eltype(::Type{Ohm{d,T}}) where {d,T} = T

Base.ndims(::Ohm{d,T}) where {d,T} = 2

Base.size(::Ohm{d,T}) where {d,T} = (d, d)

function Base.size(::Ohm{d,T}, n::Int) where {d,T}
    ((n <= 0) || (n > 2)) && throw(ErrorException("dimension must be 1 or 2 (got $n)"))
    d
end

Scapin.dimensionality(::Ohm{d,T}) where {d,T} = d

Base.convert(::Type{Array}, conductor::Ohm{d,T}) where {d,T} = (conductor.σ * I)(d)

Base.getindex(conductor::Ohm{d,T}, i::Int, j::Int) where {d,T} = (i == j) ? conductor.σ : zero(T)

function mul!(J::AbstractVector{T}, conductor::Ohm{d,T}, E::AbstractVector{T}) where {d,T}
    J[1:d] = conductor.σ * E
    return J
end

function Base.:*(conductor::Ohm{d,T}, E::AbstractVector{T}) where {d,T}
    return mul!(AbstractVector{T}(undef, size(conductor, 1)), conductor, E)
end

"""
Continuous Green operator related to [Ohm](@ref) conductors.

"""
struct GreenOperatorOhm{d,T}
    conductor::Ohm{d,T}
end

"""
   eltype(type)

Return the type of the (scalar) elements the operator of given `type` operates
on in the *real* space (as opposed to the *Fourier* space). Note that this type
may be complex!
"""
Base.eltype(::Type{GreenOperatorOhm{d,T}}) where {d,T} = T

Base.ndims(::GreenOperatorOhm{d,T}) where {d,T} = 2

Base.size(::GreenOperatorOhm{d,T}) where {d,T} = (d, d)

function Base.size(::GreenOperatorOhm{d,T}, n::Int) where {d,T}
    ((n <= 0) || (n > 2)) && throw(ErrorException("dimension must be 1 or 2 (got $n)"))
    d
end

Scapin.dimensionality(::GreenOperatorOhm{d,T}) where {d,T} = d

function Scapin.apply_fourier!(Ê, Γ::GreenOperatorOhm{d,T}, k, P̂) where {d, T}
    (k ⋅ P̂) / (Γ.conductor.σ * sum(abs2, k)) .* k
end

export Ohm, GreenOperatorOhm

end
