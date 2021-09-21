module Brick

using Base.Iterators
using LinearAlgebra
using Scapin.Elasticity
using Scapin.Grid

"""
    integrate(f, h)

Return the `d`-dimensional integral of `f` over `(0, h[1]) × (0, h[2]) × … × (0, h[d])`.

Uses 2-point Gauss-Legendre integration (tensorized over the `d` dimensions). `f` must
take a 1-dimensional array of size `d` as unique input. If `avg` is `true`, the
function returns the `N`-dimensional average.
"""
function integrate(f, h::NTuple{d,T}; avg = false) where {d,T<:Number}
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
function shape(x::NTuple{d,T}, h::NTuple{d,T}) where {d,T<:Number}
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
function gradient_operator(x::NTuple{d,T}, h::NTuple{d,T}) where {d,T<:Number}
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
function avg_gradient_operator(h::NTuple{d,T}) where {d,T<:Number}
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
function strain_displacement_operator(x::NTuple{d,T}, h::NTuple{d,T}) where {d,T<:Number}
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
function avg_strain_displacement_operator(h::NTuple{d,T}) where {d,T<:Number}
    return integrate(x -> strain_displacement_operator(x, h), h, avg = true)
end


"""
    stiffness_operator(h, C)

Return the stifness operator for the brick element of size `h`, and Hooke material `C`.

The stiffness operator `K` delivers the strain energy associated to the nodal
displacements `u`

```
U = u[i, p] * K[i, p, j, q] * u[j, q] / 2,
```

where `i, j ∈ {1, …, d}` are component indices and `p, q ∈ CartesianIndices(1:2, …, 1:2)`.

"""
function stiffness_operator(h::NTuple{d,T}, C::Hooke{T,d}) where {d,T<:Number}
    ℒ = cell_vertices(d)
    function f(x)
        D = gradient_operator(x, h)
        B = strain_displacement_operator(x, h)
        σε = Array{T}(undef, d, fill(2, d)..., d, fill(2, d)...)
        for i = 1:d, p ∈ ℒ, j = 1:d, q ∈ ℒ
            σε[i, p, j, q] =
                C.λ * D[i, p] * D[j, q] +
                2C.μ * sum(B[h, k, i, p] * B[h, k, j, q] for h = 1:d, k = 1:d)
        end
        return σε
    end
    integrate(f, h)
end


"""
    global_strain_displacement_operator(N, h)

Return the global strain-displacement operator for periodic boundary conditions.

The grid is defined by its size `N` and its spacing `h`.

The global strain-displacement operator `B` is a `d+2`-dimensional array of size
`(d, d, N[1], …, N[d])`, such that the average strain within element `p` reads

```
ε[i, j, p] = B[i, j, p, k, q] * u[k, q],
```

where

- `i, j, k ∈ {1, …, d}`: component indices,
- `p, q ∈ CartesianIndices(1:N[1], …, 1:N[d])`: node indices,
- `ε[i, j, p]`: `(i, j)`-th component of the average strain in element `p`,
- `u[k, q]`: `i`-th component of the displacement of node `q`.

!!! note

    Assembly of the global strain-displacement operator is done under the
    assumption of periodicity.

"""
function global_strain_displacement_operator(
    N::NTuple{d,Int},
    h::NTuple{d,T},
) where {d,T<:Number}
    Be = avg_strain_displacement_operator(h)
    B = Array{T}(undef, d, d, N..., d, N...)
    𝓛 = cell_vertices(d) # Local node indices
    for p ∈ CartesianIndices(N)
        𝒢 = cell_vertices(p, N) # global node indices
        for q ∈ 𝓛
            B[:, :, p, :, 𝒢[q]] += Be[:, :, :, q]
        end
    end
    return B
end


"""
    global_stiffness_operator(N, h, μ, ν)

Return the global stiffness operator for periodic, homogeneous elasticity.

The grid size is `N`, the cell size is `h`. The constitutive material is
homogeneous, elastic linear and isotropic with stiffness `C`.

The global stiffness operator `K` is a `2d+2`-dimensional array of size
`(d, N[1], …, N[d], d, N[1], …, N[d])`, such that the strain energy of the
system reads

```
U = u[i, p] * K[i, p, j, q] * u[j, q] / 2,
```

where

- `i, j ∈ {1, …, d}`: component indices,
- `p, q ∈ CartesianIndices(1:N[1], …, 1:N[d])`: node indices,
- `u[i, p]`: `i`-th component of the displacement of node `p`.

!!! note

    Assembly of the global stiffness opertor is done under the assumption of
    periodicity.

"""
function global_stiffness_operator(
    N::NTuple{d,Int},
    h::NTuple{d,T},
    C::Hooke{T,d},
) where {d,T<:Number}
    Ke = stiffness_operator(h, C)
    K = zeros(T, d, N..., d, N...)
    𝓛 = cell_vertices(d) # Local node indices
    for e ∈ CartesianIndices(N)
        𝒢 = cell_vertices(e, N) # global node indices
        for m ∈ 𝓛, n ∈ 𝓛
            K[𝒢[m], :, 𝒢[n], :] += Ke[m, :, n, :]
        end
    end
    return K
end

export integrate,
    shape,
    gradient_operator,
    avg_gradient_operator,
    strain_displacement_operator,
    avg_strain_displacement_operator,
    stiffness_operator,
    global_strain_displacement_operator,
    global_stiffness_operator
end
