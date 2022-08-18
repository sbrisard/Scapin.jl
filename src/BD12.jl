module BD12

using Base.Iterators
using Scapin

# CONT_OP: type of continuous operator
# d: number of space dimensions
# T: type of scalars for CONT_OP
# We must have: T == eltype(CONT_OP)

struct BrisardDormieux2012{CONT_OP,d,T}
    F::CONT_OP
    N::NTuple{d,Int}
    C::NTuple{d,Vector{T}}
    function BrisardDormieux2012{CONT_OP,d,T}(F::CONT_OP, N) where {CONT_OP,d, T}
        # TODO Check that d == dimensionality(F) == ndims(N)
        # TODO Check that T == eltype(CONT_OP)
        # TODO Check that ndims(F) == 2
        C = Tuple([cos(π * nᵢ / (2Nᵢ)) for nᵢ = 0:Nᵢ] for Nᵢ in N)
        new{CONT_OP,d, T}(F, N, C)
    end
end

BrisardDormieux2012(F::CONT_OP, N) where {CONT_OP} = BrisardDormieux2012{CONT_OP,dimensionality(CONT_OP),eltype(CONT_OP)}(F, N)

Base.eltype(::Type{BrisardDormieux2012{CONT_OP, d,T}}) where {CONT_OP,d,T} =
    eltype(CONT_OP)
Base.ndims(::Type{BrisardDormieux2012{CONT_OP,d,T}}) where {CONT_OP,d,T} = 2
Base.ndims(::BrisardDormieux2012{CONT_OP,d,T}) where {CONT_OP,d,T} = 2

Base.size(F_N::BrisardDormieux2012{CONT_OP, d, T}) where {CONT_OP,d, T} = size(F_N.F)
Base.size(F_N::BrisardDormieux2012{CONT_OP, d, T}, n::Int) where {CONT_OP, d, T} = size(F_N.F, n)

Scapin.dimensionality(::Type{BrisardDormieux2012{CONT_OP, d, T}}) where {CONT_OP, d, T} = d

Scapin.grid_size(F_N::BrisardDormieux2012{CONT_OP, d, T}) where {CONT_OP, d, T} = F_N.N

# function Scapin.apply_fourier!(ŷ, 𝔽::BrisardDormieux2012{d, CONT_OP}, n, x̂) where{d, CONT_OP}
#     # d = dimensionality(𝔽)
#     N = grid_size(𝔽)
#     ŷ[:] = zero(eltype(𝔽))
#     for i in 1:d
#     end

#     for m ∈ product(repeated(-1:0, d))
#         k = 2π .* ((n .- 1) .+ m .* N) ./ N
#         C = cos.(k ./ 4)
#         ŷ += C^2 * apply_fourier(𝔽.ℱ, k, x̂)
#     end
#     return ŷ
# end

export BrisardDormieux2012

end # of module
