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
    element_nodes(d)

Return the multi-indices of the `d`-dimensional brick element.

This function returns an object `𝔑` of type `CartesianIndices`. An
element `n ∈ 𝔑` represents the node with coordinates `(x[1], …, x[d])`

```
x[i] = (-1)^n[i] * h[i] / 2,    i = 1, …, d
```

See “[Geometry of the reference brick element](@ref _20210910120306)” in the docs.
"""
element_nodes(d::Int) = CartesianIndices(Tuple(fill(1:2, d)))


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
    𝔑 = element_nodes(d)
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
    𝔑 = element_nodes(d)
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

    𝔑 = element_nodes(d)
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
    element_nodes(e, N)

Return the nodes of element `e` within the grid of size `N`.

The element is specified as a `d`-dimensional multi-index, `e`. The function
returns an array of `2^d` multi-indices. The size of the grid is specified by the
array `N` (`size(N, 1) = d`).

Note that periodic boundary conditions are accounted for.

"""
function element_nodes(e, N::AbstractVector{Int})
    d = size(N, 1)
    e_ = collect(e)
    increments = map(collect, product(fill(0:1, d)...))
    reshape(map(n -> Tuple((e_ .+ n .- 1) .% N .+ 1), increments), 2^d)
end

"""
    global_stiffness_matrix(N, h, μ, ν)

Return the global stiffness matrix for periodic, homogeneous elasticity.

- `N`: grid size
- `h`: cell size
- `μ`: shear modulus
- `ν`: Poisson ratio

"""
function global_stiffness_operator(
    N::AbstractVector{Int},
    h::AbstractVector{T},
    μ::T,
    ν::T,
) where {T<:Number}
    d = size(N, 1)
    cartesian = CartesianIndices(Tuple(N))
    linear = LinearIndices(Tuple(N))
    Ke = stiffness_operator(h, μ, ν)
    K = zeros(T, N..., d, N..., d)
    for e ∈ cartesian
        nodes = cell_vertices(e, linear)
        K[nodes, :, nodes, :] += Ke
    end
    return K
end
