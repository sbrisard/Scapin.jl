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

The `ξ[i]` (`i = 1, …, d`) are the reduced coordinates, such that `0 ≤ ξ[i] ≤ 1`.

For `k == 0`, the function returns a vector `N`, such that `N[i]` is the value of the
`i`-th shape function at `ξ`. The shape functions are ordered using the column-major
convention. More precisely, if `ξ` are the reduced coordinates of a node (`ξ[i] ∈ {0, 1}`
for all `i = 1, …, d`), we can define the multi-index `I` such that `I[i] = ξ[i] + 1`.
Then, the index of the node under consideration can be found from

```
linear = LinearIndices(1:2, 1:2)      # d == 2
linear = LinearIndices(1:2, 1:2, 1:2) # d == 3
index = linear(I...)
```

For `k == 1, …, d`, the function returns the vector of the derivatives of the shape
functions at `ξ`, with respect to `ξ[k]`.

"""
function shape(ξ::AbstractArray{T,1}, k::Int) where {T<:Number}
    d = size(ξ, 1)
    N = Array{T}(undef, d, 2)
    for i = 1:d
        N[i, 1] = k == i ? -one(T) : 1 - ξ[i]
        N[i, 2] = k == i ? one(T) : ξ[i]
    end
    ranges = Tuple(fill(1:2, d))
    cartesian = CartesianIndices(ranges)
    [prod(N[i, cartesian[j][i]] for i = 1:d) for j = 1:length(cartesian)]
end

"""
    strain_displacement_matrix(x, h)

Return the strain-displacement matrix for the `d`-dimensional element of size `h`.

The strain-displacement matrix `B` is such that `B * q` is the strain at point `x`
(`d × d` matrix). Note that the degrees of freedom are ordered as follows

```
q = [u₁, v₁, u₂, v₂, u₃, v₃, u₄, v₄]           (2d)
q = [u₁, v₁, w₁, u₂, v₂, w₂, …, u₈, v₈, w₈]    (3d)
```

where `uₖ`, `vₖ` and `wₖ` are the components of the displacement of node `k` in the `x`,
`y` and `z` directions, respectively.
"""
function strain_displacement_matrix(
    x::AbstractVector{T},
    h::AbstractVector{T},
) where {T<:Number}
    d = size(x, 1)
    @assert size(h, 1) == d "x and h must have same size"
    nvertices = 2^d
    ndofs = d * nvertices

    ξ = x ./ h
    ∂ = hcat((shape(ξ, i) ./ h[i] for i = 1:d)...)
    u = zeros(T, d, nvertices, ndofs)
    for i = 1:d, j = 1:nvertices
        u[i, j, d*(j-1)+i] = one(T)
    end
    B = Array{T}(undef, d, d, ndofs)
    for i = 1:d, j = 1:d
        B[i, j, :] = (∂[:, i]' * u[j, :, :] + ∂[:, j]' * u[i, :, :]) / 2
    end
    return B
end


"""
    avg_strain_displacement_matrix(h)

Return the strain-displacement matrix, averaged over the whole element.
"""
function avg_strain_displacement_matrix(h::AbstractArray{T,1}) where {T<:Number}
    return integrate(x -> strain_displacement_matrix(x, h), h, avg = true)
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
