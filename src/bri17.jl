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


# Helper functions (for testing purposes)

"""
    integrate(f, h)

Return the `N`-dimensional integral of `f` over `(0, h[1]) × (0, h[2]) × … × (0, h[N])`.

Uses 2-point Gauss-Legendre integration (tensorized over the `N` dimensions). `f` must
take a 1-dimensional array of size `N` as unique input. If `avg` is `true`, the
function returns the `N`-dimensional average.
"""
function integrate(f, h::AbstractArray{T,1}; avg = false) where {T<:Number}
    ndims = size(h, 1)
    nvertices = 2^ndims
    ξ = [(1 - 1 / sqrt(3)) / 2, (1 + 1 / sqrt(3)) / 2]
    weight = (avg ? one(T) : prod(h)) / nvertices
    x = map(collect, product((ξ .* hᵢ for hᵢ in h)...))
    return weight * sum(f, x)
end


"""
    shape(ξ, k)

Return the value of the shape functions or their derivatives, at the specified point.

`ξ` is the `d`-dimensional array of the reduced coordinates (`d`: number of
spatial dimensions), such that `0 ≤ ξ[i] ≤ 1`.

Vertices of the element are referred to with their cartesian index `n`, a
`d`-dimensional array such that `n[1], …, n[d] ∈ {1, 2}`. The reduced coordinates
 of vertex `n` are `n .- 1`.

For `k == 0`, the function returns a `d`-dimensional array `N`, such that `N[n]`
is the shape function associated with vertex `n`. In other words,

```
size(N, 1) = … = size(N, d) = 2
```.

```
shape(n .- 1, 0)[m] = 0    (m ≠ n),
shape(n .- 1, 0)[n] = 1.
```

For `k == 1, …, d`, the function returns the array of the derivatives of the
shape functions at `ξ`, with respect to `ξ[k]`.

"""
function shape(ξ::AbstractArray{T,1}, k::Int) where {T<:Number}
    d = size(ξ, 1)
    N_1d = Array{T}(undef, d, 2)
    for i = 1:d
        N_1d[i, 1] = k == i ? -one(T) : 1 - ξ[i]
        N_1d[i, 2] = k == i ? one(T) : ξ[i]
    end
    cartesian = CartesianIndices(Tuple(fill(1:2, d)))
    return [prod(N_1d[i, n[i]] for i = 1:d) for n in cartesian]
end

"""
    strain_displacement_operator(x, h)

Return the strain-displacement matrix for the `d`-dimensional element of size `h`.

The strain-displacement operator maps the array of nodal displacements `u` onto
the array of strains at `x`. `u` is a `d+1` dimensional array (`d`: number of
spatial dimensions) such that

```
size(u, 1) = … = size(u, d) = 2``,
size(u, d+1) = d
```

and the `i`-th component of the displacement of vertex `n` is `u[n..., i]`. See
`shape_function` for the numbering of vertices. The strain at `x`, `ε` is a
2-dimensional array array such that

```
size(ε, 1) = size(ε, 2) = d
```

and

```
ε[i, j] = Σₙ Σₖ B[i, j, n..., k] * u[n..., k],
```

where `B` is the strain-displacement operator `B` returned by the present
function. This `d+3`-dimensional array is such that

```
size(B, 1) = size(B, 2) = d,
size(B, 3) = … = size(B, d+2) = 2,
size(B, d+3) = d.
```

"""
function strain_displacement_operator(
    x::AbstractVector{T},
    h::AbstractVector{T},
) where {T<:Number}
    # ∂[n, i] = ∂ᵢN[n] (partial derivative w.r.t x, not ξ!)
    #
    # uᵢ = Σₙ N[n] * u[n, i]
    #
    # εᵢⱼ = ½(∂ⱼuᵢ + ∂ᵢuⱼ)
    #     = ½Σₙ {∂ⱼN[n] * u[n, i] + ∂ᵢN[n] * u[n, j]}
    #     = ½Σₙ {∂[n, j] * u[n, i] + ∂[n, i] * u[n, j]}
    #     = ½ΣₙΣₖ{∂[n, j] * δ[i, k] + ∂[n, i] * δ[j, k]} * u[n, k]
    #
    # B[i, j, n, k] = ½{∂[n, j] * δ[i, k] + ∂[n, i] * δ[j, k]}
    d = size(x, 1)
    @assert size(h, 1) == d "x and h must have same size"

    ξ = x ./ h
    ∂ = cat((shape(ξ, i) ./ h[i] for i = 1:d)..., dims=d+1)
    cartesian = CartesianIndices(Tuple(fill(1:2, d)))

    B = zeros(T, d, d, fill(2, d)..., d)
    for i = 1:d, j = 1:d, n ∈ cartesian
        B[i, j, n, i] += ∂[n, j] / 2
        B[i, j, n, j] += ∂[n, i] / 2
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

function stiffness_matrix(h::AbstractArray{T,1}, μ::T, ν::T) where {T<:Number}
    d = size(h, 1)
    ndofs = d * 2^d
    λ = 2μ * ν / (1 - 2ν)
    function f(x)
        B = strain_displacement_matrix(x, h)
        tr_B = [tr(B[:, :, i]) for i = 1:ndofs]
        Ke = Array{T}(undef, ndofs, ndofs)
        for i = 1:ndofs, j = 1:ndofs
            Ke[i, j] =
                λ * tr_B[i] * tr_B[j] +
                2μ * sum(B[h, k, i] * B[h, k, j] for h = 1:d, k = 1:d)
        end
        return Ke
    end
    integrate(f, h)
end

"""
    cell_vertices(cell, linear)

Return the vertices of the `cell` (specified as a multi-index) as an
array of indices. The grid size is defined by `linear` (an instance of
`LinearIndices`).
"""
function cell_vertices(K, linear::LinearIndices)
    d = ndims(linear)
    bounds = [(K[i], (K[i] % axes(linear, i)[end] + 1)) for i = 1:d]
    return collect(flatten(linear[i...] for i in product(bounds...)))
end

"""
    global_stiffness_matrix(N, h, μ, ν)

Return the global stiffness matrix for periodic, homogeneous elasticity.

- `N`: grid size
- `h`: cell size
- `μ`: shear modulus
- `ν`: Poisson ratio

"""
function stiffness_matrix(
    N::AbstractVector{Int},
    h::AbstractVector{T},
    μ::T,
    ν::T,
) where {T<:Number}
    ndofs = prod(N) * size(N, 1)
    dims = Tuple(1:n for n in N)
    cartesian = CartesianIndices(dims)
    linear = LinearIndices(dims)
    Ke = stiffness_matrix(h, μ, ν)
    K = zeros(T, ndofs, ndofs)
    for J in CartesianIndices
        vertices = cell_vertices(J, linear)
        K[vertices, vertices] += Ke
    end
    return K
end
