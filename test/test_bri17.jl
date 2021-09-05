using Base.Iterators
using LinearAlgebra

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


@testset "integrate" begin
    h = [1.1, 1.2, 1.3]
    n = [1.0, 2.0, 3.0]
    f(x) = prod(x .^ n)
    actual = integrate(f, h)
    expected = prod((h .^ (n .+ 1)) ./ (n .+ 1))
    @test isapprox(actual, expected, rtol = 1e-15)
end

for d = 1:3
    @testset "Node numbering, $(d)d" begin
        ranges = Tuple(fill(1:2, d))
        linear = LinearIndices(ranges)
        actual = zeros(Int, length(linear), length(linear))
        for ξ in map(collect, product(fill((0, 1), d)...))
            actual[:, linear[(ξ .+ 1)...]] = shape(ξ, 0)
        end
        @test actual == I
    end
end

@testset "Average strain-displacement matrix 2d" begin
    h = [1.1, 1.2]
    # Note: reference values where computed with maxima, using the
    # Kelvin–Mandel representation and a different numbering of dofs.
    nodes = [1, 5, 3, 7, 2, 6, 4, 8]
    B₁ = 0.4545454545454545
    B₂ = 0.4166666666666667
    B₃ = 0.2946278254943947
    B₄ = 0.3214121732666125
    B_exp = [
        -B₁ -B₁ B₁ B₁ 0 0 0 0
        0 0 0 0 -B₂ B₂ -B₂ B₂
        -B₃ B₃ -B₃ B₃ -B₄ -B₄ B₄ B₄
    ]

    B_act = avg_strain_displacement_matrix(h)

    @test isapprox(B_act[1, 1, nodes], B_exp[1, :], rtol = 1e-15)
    @test isapprox(B_act[2, 2, nodes], B_exp[2, :], rtol = 1e-15)
    @test isapprox(B_act[1, 2, nodes], B_exp[3, :] / sqrt(2), rtol = 1e-15)
end

@testset "Average strain-displacement matrix 3d" begin
    h = [1.1, 1.2, 1.3]
    # Note: reference values where computed with maxima, using the
    # Kelvin–Mandel representation and a different numbering of dofs.
    nodes = [1, 13, 7, 19, 4, 16, 10, 22,
             2, 14, 8, 20, 5, 17, 11, 23,
             3, 15, 9, 21, 6, 18, 12, 24]
    B₁ = 0.2272727272727273
    B₂ = 0.2083333333333333
    B₃ = 0.1923076923076923
    B₄ = 0.1359820733051053
    B₅ = 0.1473139127471974
    B₆ = 0.1607060866333062
    B_exp = [
        -B₁ -B₁ -B₁ -B₁ B₁ B₁ B₁ B₁ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 -B₂ -B₂ B₂ B₂ -B₂ -B₂ B₂ B₂ 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -B₃ B₃ -B₃ B₃ -B₃ B₃ -B₃ B₃;
        0 0 0 0 0 0 0 0 -B₄ B₄ -B₄ B₄ -B₄ B₄ -B₄ B₄ -B₅ -B₅ B₅ B₅ -B₅ -B₅ B₅ B₅;
        -B₄ B₄ -B₄ B₄ -B₄ B₄ -B₄ B₄ 0 0 0 0 0 0 0 0 -B₆ -B₆ -B₆ -B₆ B₆ B₆ B₆ B₆;
        -B₅ -B₅ B₅ B₅ -B₅ -B₅ B₅ B₅ -B₆ -B₆ -B₆ -B₆ B₆ B₆ B₆ B₆ 0 0 0 0 0 0 0 0
    ]

    B_act = avg_strain_displacement_matrix(h)

    @test isapprox(B_act[1, 1, nodes], B_exp[1, :], rtol = 1e-15)
    @test isapprox(B_act[2, 2, nodes], B_exp[2, :], rtol = 1e-15)
    @test isapprox(B_act[3, 3, nodes], B_exp[3, :], rtol = 1e-15)
    @test isapprox(B_act[2, 3, nodes], B_exp[4, :] / sqrt(2), rtol = 1e-15)
    @test isapprox(B_act[3, 1, nodes], B_exp[5, :] / sqrt(2), rtol = 1e-15)
    @test isapprox(B_act[1, 2, nodes], B_exp[6, :] / sqrt(2), rtol = 1e-15)
end



# @testset "bri17, 2D" begin
#     C = Hooke{Float64, 2}(5.6, 0.3)
#     N = [2, 4]
#     h = [1.1, 1.2]
# end
