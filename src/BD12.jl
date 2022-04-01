using Iterators
using Scapin

struct BD12{CONT_OP}
    ℱ::CONT_OP
    N::NTuple{ndims(CONT_OP),Int}
    C::NTuple{ndims(CONT_OP),Vector{eltype{CONT_OP}}} # TODO use static arrays
    function BD12{CONT_OP}(ℱ::CONT_OP, N)
        d = dimensionality(ℱ)
        C = Tuple([cos(0.5π * n / N) for n in 0:N] for i in 1:d)
        new(ℱ, N, C)
    end
end

BD12(ℱ::CONT_OP, N::NTuple{d,Int}) where {d,CONT_OP} = BD12{d,CONT_OP}(ℱ, N)

Base.eltype(::Type{BD12{T}}) where {T} = eltype(T)
Base.ndims(::BD12{T}) where {T} = ndims(T)
Base.size(::BD12{T}) where {T} = size(T)
Base.size(::BD12{T}, n::Int) where {d,T} = size(T, n)
Scapin.dimensionality(::BD12{T}) where {T} = dimensionality(T)

function Scapin.apply_fourier!(ŷ, 𝔽::BD12{CONT_OP}, n, x̂)
    d = dimensionality(𝔽)
    N = grid_size(𝔽)
    ŷ[:] = zero(eltype(T))
    for i in 1:d
    end

    for m ∈ product(repeated(-1:0, d))
        k = 2π .* ((n .- 1) .+ m .* N) ./ N
        C = cos.(k ./ 4)
        ŷ += C^2 * apply_fourier(𝔽.ℱ, k, x̂)
    end
    return ŷ
end
