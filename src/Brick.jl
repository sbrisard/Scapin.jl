module Brick

using Base.Iterators
using LinearAlgebra
using Scapin.Grid

"""
    integrate(f, h)

Return the `d`-dimensional integral of `f` over `(0, h[1]) × (0, h[2]) × … × (0, h[d])`.

Uses 2-point Gauss-Legendre integration (tensorized over the `d` dimensions). `f` must
take a 1-dimensional array of size `d` as unique input. If `avg` is `true`, the
function returns the `N`-dimensional average.
"""
function integrate(f, h::NTuple{d, T}; avg = false) where {d, T<:Number}
    ξ = [-1 / √3, 1 / √3]
    weight = (avg ? one(T) : prod(h)) / 2^d
    return weight * sum(f, product((ξ .* h_ / 2 for h_ in h)...))
end


"""
    shape(x, h)

Return the value of the shape functions of the element, at the specified point.

This function returns a `d`-dimensional array `N`, such that `N[p]` is the shape
function associated with node `p`, evaluated at `x`. In particular, `N[p]`
evaluated at node `q` is `δ[p, q]` (Kronecker).

See “[Shape functions](@ref _20210910114136)” in the docs.

"""
function shape(x::NTuple{d, T}, h::NTuple{d, T}) where {d, T<:Number}
    ξ = 2 .* x ./ h
    ℒ = cell_vertices(d)
    return [prod((1 + (-1)^p[i] * ξ[i]) / 2 for i = 1:d) for p in ℒ]
end

"""
    gradient_operator(x, h)

Return the gradient operator at the specified point.

This function returns a `(d+1)` dimensional array `D` of size `(d, 2, 2, …)`.
If  `i` is the index of a component and `p` the `CartesianIndex` of the node, then
`D[i, p]` is the partial derivative of `N[p]` w.r.t. `x[i]`, evaluated at `x`.

`h` is the size of the brick element.

See “[Geometry of the reference brick element](@ref _20210910120306)”.
"""
function gradient_operator(x::NTuple{d, T}, h::NTuple{d, T}) where {d, T<:Number}
    ξ = 2 .* x ./ h
    ℒ = cell_vertices(d)
    return [
        prod(j == i ? (-1)^n[j] / h[j] : (1 + (-1)^n[j] * ξ[j]) / 2 for j = 1:d) for
        i = 1:d, n in ℒ
    ]
end


"""
    avg_gradient_operator(h)

Return the cell-average of the gradient operator.
"""
function avg_gradient_operator(h::NTuple{d, T}) where {d, T<:Number}
    return integrate(x -> gradient_operator(x, h), h, avg = true)
end


"""
    strain_displacement_operator(x, h)

Return the strain-displacement operator for the `d`-dimensional brick element of
size `h`, evaluated at point `x`.

This function returns a `(d+3)` dimensional array `B` of size `(d, d, 2, …, 2, d)`.
If `p` is the `CartesianIndex` of the node, and `i`, `j`, `k` are component
indices then, the interpolated `(i, j)` component of the strain at `x` reads

```
ε[i, j] = Σₖ Σₚ B[i, j, k, p] * u[k, p].
```

See “[Gradient and strain-displacement operators](@ref _20210910114926)”.
"""
function strain_displacement_operator(
    x::NTuple{d, T},
    h::NTuple{d, T},
) where {d, T<:Number}
    ℒ = cell_vertices(d)
    D = gradient_operator(x, h)
    B = zeros(T, d, d, d, fill(2, d)...)
    for i = 1:d, j = 1:d, p ∈ ℒ
        B[i, j, i, p] += D[j, p] / 2
        B[i, j, j, p] += D[i, p] / 2
    end
    return B
end


"""
    avg_strain_displacement_operator(h)

Return the cell average of the strain-displacement operator.
"""
function avg_strain_displacement_operator(h::NTuple{d, T}) where {d, T<:Number}
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
