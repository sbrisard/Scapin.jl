module Brick

using Base.Iterators
using LinearAlgebra
using Scapin.Grid

"""
    integrate(f, h)

Return the `N`-dimensional integral of `f` over `(0, h[1]) × (0, h[2]) × … × (0, h[N])`.

Uses 2-point Gauss-Legendre integration (tensorized over the `N` dimensions). `f` must
take a 1-dimensional array of size `N` as unique input. If `avg` is `true`, the
function returns the `N`-dimensional average.
"""
function integrate(f, h::AbstractArray{T,1}; avg = false) where {T<:Number}
    d = size(h, 1)
    ξ = [-1 / sqrt(3), 1 / sqrt(3)]
    weight = (avg ? one(T) : prod(h)) / 2^d
    x = map(collect, product((ξ .* h_ / 2 for h_ in h)...))
    return weight * sum(f, x)
end


"""
    shape(x, h)

Return the value of the shape functions of the element, at the specified point.

This function returns a `d`-dimensional array `N`, such that `N[n]` is the shape
function associated with node `n`, evaluated at `x`. In particular, `N[n]`
evaluated at node `m` is `δ[m, n]` (Kronecker).

See “[Shape functions](@ref _20210910114136)” in the docs.

"""
function shape(x::AbstractArray{T,1}, h::AbstractArray{T,1}) where {T<:Number}
    d = size(x, 1)
    # TODO — Check that x and h have same size
    ξ = 2 * x ./ h
    𝔑 = cell_vertices(d)
    return [prod((1 + (-1)^n[i] * ξ[i]) / 2 for i = 1:d) for n in 𝔑]
end

"""
    gradient_operator(x, h)

Return the gradient operator at the specified point.

This function returns a `(d+1)` dimensional array `D` of size `(d, 2, 2, …)`.
If `n` is the multi-index of the node, and `i` is the index of a component, then
`D[i, n]` is the partial derivative of `N[n]` w.r.t. `x[i]`, evaluated at `x`.

`h` is the size of the brick element.

See “[Geometry of the reference brick element](@ref _20210910120306)”.
"""
function gradient_operator(x::AbstractVector{T}, h::AbstractVector{T}) where {T<:Number}
    d = size(x, 1)
    # TODO — Check that x and h have same size
    ξ = 2 * x ./ h
    𝔑 = cell_vertices(d)
    return [
        prod(j == i ? (-1)^n[j] / h[j] : (1 + (-1)^n[j] * ξ[j]) / 2 for j = 1:d) for
        i = 1:d, n in 𝔑
    ]
end


"""
    avg_gradient_operator(h)

Return the cell average of the gradient operator.
"""
function avg_gradient_operator(h::AbstractArray{T,1}) where {T<:Number}
    return integrate(x -> gradient_operator(x, h), h, avg = true)
end


"""
    strain_displacement_operator(x, h)

Return the strain-displacement operator for the `d`-dimensional brick element of
size `h`, evaluated at point `x`.

This function returns a `(d+3)` dimensional array `B` of size `(d, d, 2, …, 2, d)`.
If `n` is the multi-index of the node, and `i`, `j`, `k` are component indices
then, the interpolated `(i, j)` component of the strain at `x` reads

```
ε[i, j] = Σₙ Σₖ B[i, j, n, k] * u[n, k].
```

See “[Gradient and strain-displacement operators](@ref _20210910114926)”.
"""
function strain_displacement_operator(
    x::AbstractVector{T},
    h::AbstractVector{T},
) where {T<:Number}
    d = size(x, 1)
    @assert size(h, 1) == d "x and h must have same size"

    𝔑 = cell_vertices(d)
    D = gradient_operator(x, h)
    B = zeros(T, d, d, fill(2, d)..., d)
    for i = 1:d, j = 1:d, n ∈ 𝔑
        B[i, j, n, i] += D[j, n] / 2
        B[i, j, n, j] += D[i, n] / 2
    end
    return B
end


"""
    avg_strain_displacement_operator(h)

Return the cell average of the strain-displacement operator.
"""
function avg_strain_displacement_operator(h::AbstractArray{T,1}) where {T<:Number}
    return integrate(x -> strain_displacement_operator(x, h), h, avg = true)
end

function stiffness_operator(h::AbstractArray{T,1}, μ::T, ν::T) where {T<:Number}
    # εₕₖ = Σₙⱼ B[h, k, n, j] * u[n, j]
    #
    # tr(ε) = Σₙⱼₖ B[k, k, n, j] * u[n, j]
    #       = Σₙⱼ tr_B[n, j] * u[n, j]
    #
    # where
    #
    # tr_B[n, j] = Σₖ B[k, k, n, j]
    #
    # σₕₖ = λ ⋅ tr(ε) ⋅ δₕₖ + 2μ ⋅ εₕₖ
    #     = Σₘᵢ {λ * tr_B[m, i] * δ[h, k] + 2 * μ * B[h, k, m, i]} * u[m, i]
    #
    # ½ σ : ε
    # = ½ Σₘᵢₙⱼₕₖ {λ * tr_B[m, i] * δ[h, k] + 2 * μ * B[h, k, m, i]} * B[h, k, n, j]
    #             * u[m, i] * u[n, j]
    # = ½ Σₘᵢₙⱼ {λ *  tr_B[m, i] * tr_B[n, j] + 2μ * Σₕₖ B[h, k, m, i] * B[h, k, n, j]}
    #           * u[m, i] * u[n, j]
    #
    # Kₘᵢₙⱼ = ∫ {λ *  tr_B[m, i] * tr_B[n, j] + 2μ * Σₕₖ B[h, k, m, i] * B[h, k, n, j]}

    d = size(h, 1)
    nodes = CartesianIndices(Tuple(fill(1:2, d)))
    λ = 2μ * ν / (1 - 2ν)
    function f(x)
        B = strain_displacement_operator(x, h)
        tr_B = [tr(B[:, :, n, i]) for n ∈ nodes, i = 1:d]
        σε = Array{T}(undef, fill(2, d)..., d, fill(2, d)..., d)
        for m ∈ nodes, i = 1:d, n ∈ nodes, j = 1:d
            σε[m, i, n, j] =
                λ * tr_B[m, i] * tr_B[n, j] +
                2μ * sum(B[h, k, m, i] * B[h, k, n, j] for h = 1:d, k = 1:d)
        end
        return σε
    end
    integrate(f, h)
end

"""
    global_stiffness_matrix(N, h, μ, ν)

Return the global stiffness matrix for periodic, homogeneous elasticity.

- `N`: grid size
- `h`: cell size
- `μ`: shear modulus
- `ν`: Poisson ratio

"""
# function global_stiffness_operator(
#     N::AbstractVector{Int},
#     h::AbstractVector{T},
#     μ::T,
#     ν::T,
# ) where {T<:Number}
#     d = size(N, 1)
#     cartesian = map(collect∘Tuple, CartesianIndices(Tuple(N)))
#     linear = LinearIndices(Tuple(N))
#     Ke = stiffness_operator(h, μ, ν)
#     K = zeros(T, N..., d, N..., d)
#     for e ∈ cartesian
#         nodes = element_nodes(e, N)
#         K[nodes, :, nodes, :] += Ke
#     end
#     return K
# end

export integrate, shape, gradient_operator, avg_gradient_operator, strain_displacement_operator, avg_strain_displacement_operator, stiffness_operator
end
