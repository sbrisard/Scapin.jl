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

function global_stiffness_matrix(
    N::AbstractVector{Int},
    h::AbstractVector{T},
    μ::T,
    ν::T,
) where {T<:Number}
    d = size(N, 1)
    nnodes_per_cell = 2^d
    ndofs_per_cell = d * nnodes_per_cell
    ncells = prod(N)
    ndofs = ndofs_per_cell * ncells
    cartesian = CartesianIndices((1:Nᵢ for Nᵢ in N)...)
end


@testset "Gauss–Legendre integration" begin
    h = [1.1, 1.2, 1.3]
    n = [1.0, 2.0, 3.0]
    f(x) = prod(x .^ n)
    actual = integrate(f, h)
    expected = prod((h .^ (n .+ 1)) ./ (n .+ 1))
    @test isapprox(actual, expected, rtol = 1e-15)
end

@testset "Node numbering" begin
    for d = 1:3
        ranges = Tuple(fill(1:2, d))
        linear = LinearIndices(ranges)
        actual = zeros(Int, length(linear), length(linear))
        for ξ in map(collect, product(fill((0, 1), d)...))
            actual[:, linear[(ξ .+ 1)...]] = shape(ξ, 0)
        end
        @test actual == I
    end
end

@testset "2d quadrilateral element for linear elasticity" begin
    # Note: reference values where computed with maxima, using the
    # Kelvin–Mandel representation and a different numbering of dofs.
    h = [1.1, 1.2]
    μ = 5.6
    ν = 0.3
    nodes = [1, 5, 3, 7, 2, 6, 4, 8]

    @testset "Average strain-displacement matrix, 2d" begin
        B1 = 0.4545454545454545
        B2 = 0.4166666666666667
        B3 = 0.2946278254943947
        B4 = 0.3214121732666125
        B_exp = [
            -B1 -B1 B1 B1 0 0 0 0
            0 0 0 0 -B2 B2 -B2 B2
            -B3 B3 -B3 B3 -B4 -B4 B4 B4
        ]
        B_act = avg_strain_displacement_matrix(h)

        @test isapprox(B_act[1, 1, nodes], B_exp[1, :], rtol = 1e-15)
        @test isapprox(B_act[2, 2, nodes], B_exp[2, :], rtol = 1e-15)
        @test isapprox(B_act[1, 2, nodes], B_exp[3, :] / sqrt(2), rtol = 1e-15)
    end

    @testset "Stiffness matrix, 2d" begin
        # Note: reference values where computed with maxima, using the
        # Kelvin–Mandel representation and a different numbering of dofs.
        K1 = 8.83838383838384
        K2 = 1.852525252525252
        K3 = 6.271717171717172
        K4 = 4.41919191919192
        K5 = 3.5
        K6 = 0.7
        K7 = 8.025252525252526
        K8 = 4.97070707070707
        K9 = 0.9580808080808081
        K10 = 4.012626262626263
        K_exp = [
            K1 K2 -K3 -K4 K5 -K6 K6 -K5
            K2 K1 -K4 -K3 K6 -K5 K5 -K6
            -K3 -K4 K1 K2 -K6 K5 -K5 K6
            -K4 -K3 K2 K1 -K5 K6 -K6 K5
            K5 K6 -K6 -K5 K7 -K8 K9 -K10
            -K6 -K5 K5 K6 -K8 K7 -K10 K9
            K6 K5 -K5 -K6 K9 -K10 K7 -K8
            -K5 -K6 K6 K5 -K10 K9 -K8 K7
        ]

        K_act = stiffness_matrix(h, μ, ν)

        @test isapprox(K_act[nodes, nodes], K_exp, rtol = 1e-15)
    end
end

@testset "3d brick element for linear elasticity" begin
    # Note: reference values where computed with maxima, using the
    # Kelvin–Mandel representation and a different numbering of dofs.
    h = [1.1, 1.2, 1.3]
    μ = 5.6
    ν = 0.3
    nodes = vcat(
        [1, 13, 7, 19, 4, 16, 10, 22],
        [2, 14, 8, 20, 5, 17, 11, 23],
        [3, 15, 9, 21, 6, 18, 12, 24],
    )

    @testset "Average strain-displacement matrix" begin
        B1 = 0.2272727272727273
        B2 = 0.2083333333333333
        B3 = 0.1923076923076923
        B4 = 0.1359820733051053
        B5 = 0.1473139127471974
        B6 = 0.1607060866333062
        B_exp = [
            -B1 -B1 -B1 -B1 B1 B1 B1 B1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 -B2 -B2 B2 B2 -B2 -B2 B2 B2 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -B3 B3 -B3 B3 -B3 B3 -B3 B3
            0 0 0 0 0 0 0 0 -B4 B4 -B4 B4 -B4 B4 -B4 B4 -B5 -B5 B5 B5 -B5 -B5 B5 B5
            -B4 B4 -B4 B4 -B4 B4 -B4 B4 0 0 0 0 0 0 0 0 -B6 -B6 -B6 -B6 B6 B6 B6 B6
            -B5 -B5 B5 B5 -B5 -B5 B5 B5 -B6 -B6 -B6 -B6 B6 B6 B6 B6 0 0 0 0 0 0 0 0
        ]
        B_act = avg_strain_displacement_matrix(h)

        @test isapprox(B_act[1, 1, nodes], B_exp[1, :], rtol = 1e-15)
        @test isapprox(B_act[2, 2, nodes], B_exp[2, :], rtol = 1e-15)
        @test isapprox(B_act[3, 3, nodes], B_exp[3, :], rtol = 1e-15)
        @test isapprox(B_act[2, 3, nodes], B_exp[4, :] / sqrt(2), rtol = 1e-15)
        @test isapprox(B_act[3, 1, nodes], B_exp[5, :] / sqrt(2), rtol = 1e-15)
        @test isapprox(B_act[1, 2, nodes], B_exp[6, :] / sqrt(2), rtol = 1e-15)
    end

    @testset "Stiffness matrix" begin
        K1 = 4.461761201761202
        K2 = 1.283188293188293
        K3 = 1.118658378658379
        K4 = 0.08548303548303549
        K5 = 2.401846671846672
        K6 = 1.67476948976949
        K7 = 1.757034447034447
        K8 = 1.1154403004403
        K9 = 1.516666666666667
        K10 = 0.7583333333333333
        K11 = 0.3033333333333333
        K12 = 0.1516666666666667
        K13 = 1.4
        K14 = 0.28
        K15 = 0.7
        K16 = 0.14
        K17 = 4.109404299404299
        K18 = 1.107009842009842
        K19 = 1.838075628075628
        K20 = 1.392883967883968
        K21 = 0.7310657860657861
        K22 = 0.1083132608132608
        K23 = 1.580855995855996
        K24 = 1.027351074851075
        K25 = 1.283333333333333
        K26 = 0.2566666666666667
        K27 = 0.6416666666666667
        K28 = 0.1283333333333333
        K29 = 3.835187775187775
        K30 = 1.399329189329189
        K31 = 0.8053716653716654
        K32 = 1.255775705775706
        K33 = 0.593957523957524
        K34 = 1.361482776482776
        K35 = 0.2591323491323491
        K36 = 0.9587969437969438
        K_exp = [
            K1 K2 K3 K4 -K5 -K6 -K7 -K8 K9 K10 -K11 -K12 K11 K12 -K9 -K10 K13 -K14 K15 -K16 K14 -K13 K16 -K15
            K2 K1 K4 K3 -K6 -K5 -K8 -K7 K10 K9 -K12 -K11 K12 K11 -K10 -K9 K14 -K13 K16 -K15 K13 -K14 K15 -K16
            K3 K4 K1 K2 -K7 -K8 -K5 -K6 K11 K12 -K9 -K10 K9 K10 -K11 -K12 K15 -K16 K13 -K14 K16 -K15 K14 -K13
            K4 K3 K2 K1 -K8 -K7 -K6 -K5 K12 K11 -K10 -K9 K10 K9 -K12 -K11 K16 -K15 K14 -K13 K15 -K16 K13 -K14
            -K5 -K6 -K7 -K8 K1 K2 K3 K4 -K11 -K12 K9 K10 -K9 -K10 K11 K12 -K14 K13 -K16 K15 -K13 K14 -K15 K16
            -K6 -K5 -K8 -K7 K2 K1 K4 K3 -K12 -K11 K10 K9 -K10 -K9 K12 K11 -K13 K14 -K15 K16 -K14 K13 -K16 K15
            -K7 -K8 -K5 -K6 K3 K4 K1 K2 -K9 -K10 K11 K12 -K11 -K12 K9 K10 -K16 K15 -K14 K13 -K15 K16 -K13 K14
            -K8 -K7 -K6 -K5 K4 K3 K2 K1 -K10 -K9 K12 K11 -K12 -K11 K10 K9 -K15 K16 -K13 K14 -K16 K15 -K14 K13
            K9 K10 K11 K12 -K11 -K12 -K9 -K10 K17 K18 -K19 -K20 K21 -K22 -K23 -K24 K25 -K26 K26 -K25 K27 -K28 K28 -K27
            K10 K9 K12 K11 -K12 -K11 -K10 -K9 K18 K17 -K20 -K19 -K22 K21 -K24 -K23 K26 -K25 K25 -K26 K28 -K27 K27 -K28
            -K11 -K12 -K9 -K10 K9 K10 K11 K12 -K19 -K20 K17 K18 -K23 -K24 K21 -K22 -K26 K25 -K25 K26 -K28 K27 -K27 K28
            -K12 -K11 -K10 -K9 K10 K9 K12 K11 -K20 -K19 K18 K17 -K24 -K23 -K22 K21 -K25 K26 -K26 K25 -K27 K28 -K28 K27
            K11 K12 K9 K10 -K9 -K10 -K11 -K12 K21 -K22 -K23 -K24 K17 K18 -K19 -K20 K27 -K28 K28 -K27 K25 -K26 K26 -K25
            K12 K11 K10 K9 -K10 -K9 -K12 -K11 -K22 K21 -K24 -K23 K18 K17 -K20 -K19 K28 -K27 K27 -K28 K26 -K25 K25 -K26
            -K9 -K10 -K11 -K12 K11 K12 K9 K10 -K23 -K24 K21 -K22 -K19 -K20 K17 K18 -K28 K27 -K27 K28 -K26 K25 -K25 K26
            -K10 -K9 -K12 -K11 K12 K11 K10 K9 -K24 -K23 -K22 K21 -K20 -K19 K18 K17 -K27 K28 -K28 K27 -K25 K26 -K26 K25
            K13 K14 K15 K16 -K14 -K13 -K16 -K15 K25 K26 -K26 -K25 K27 K28 -K28 -K27 K29 -K30 K31 -K32 K33 -K34 -K35 -K36
            -K14 -K13 -K16 -K15 K13 K14 K15 K16 -K26 -K25 K25 K26 -K28 -K27 K27 K28 -K30 K29 -K32 K31 -K34 K33 -K36 -K35
            K15 K16 K13 K14 -K16 -K15 -K14 -K13 K26 K25 -K25 -K26 K28 K27 -K27 -K28 K31 -K32 K29 -K30 -K35 -K36 K33 -K34
            -K16 -K15 -K14 -K13 K15 K16 K13 K14 -K25 -K26 K26 K25 -K27 -K28 K28 K27 -K32 K31 -K30 K29 -K36 -K35 -K34 K33
            K14 K13 K16 K15 -K13 -K14 -K15 -K16 K27 K28 -K28 -K27 K25 K26 -K26 -K25 K33 -K34 -K35 -K36 K29 -K30 K31 -K32
            -K13 -K14 -K15 -K16 K14 K13 K16 K15 -K28 -K27 K27 K28 -K26 -K25 K25 K26 -K34 K33 -K36 -K35 -K30 K29 -K32 K31
            K16 K15 K14 K13 -K15 -K16 -K13 -K14 K28 K27 -K27 -K28 K26 K25 -K25 -K26 -K35 -K36 K33 -K34 K31 -K32 K29 -K30
            -K15 -K16 -K13 -K14 K16 K15 K14 K13 -K27 -K28 K28 K27 -K25 -K26 K26 K25 -K36 -K35 -K34 K33 -K32 K31 -K30 K29
        ]

        K_act = stiffness_matrix(h, μ, ν)

        @test isapprox(K_act[nodes, nodes], K_exp, rtol = 1e-15)
    end
end

@testset "Cell vertices" begin
    @testset "Cell vertices, 2d" begin
        """
        1────────2────────3────────4
        │        │        │        │
        │ (1, 4) │ (2, 4) │ (3, 4) │
        │        │        │        │
        10──────11────────12──────10
        │        │        │        │
        │ (1, 3) │ (2, 3) │ (3, 3) │
        │        │        │        │
        7────────8────────9────────7
        │        │        │        │
        │ (1, 2) │ (2, 2) │ (3, 2) │
        │        │        │        │
        4────────5────────6────────4
        │        │        │        │
        │ (1, 1) │ (2, 1) │ (3, 1) │
        │        │        │        │
        1────────2────────3────────1
        """
        N = [3, 4]
        d = size(N, 1)
        dims = Tuple(1:n for n in N)
        linear = LinearIndices(dims)
        @test cell_vertices((1, 1), linear) == [1, 2, 4, 5]
        @test cell_vertices((2, 1), linear) == [2, 3, 5, 6]
        @test cell_vertices((3, 1), linear) == [3, 1, 6, 4]
        @test cell_vertices((1, 2), linear) == [4, 5, 7, 8]
        @test cell_vertices((2, 2), linear) == [5, 6, 8, 9]
        @test cell_vertices((3, 2), linear) == [6, 4, 9, 7]
        @test cell_vertices((1, 3), linear) == [7, 8, 10, 11]
        @test cell_vertices((2, 3), linear) == [8, 9, 11, 12]
        @test cell_vertices((3, 3), linear) == [9, 7, 12, 10]
        @test cell_vertices((1, 4), linear) == [10, 11, 1, 2]
        @test cell_vertices((2, 4), linear) == [11, 12, 2, 3]
        @test cell_vertices((3, 4), linear) == [12, 10, 3, 1]
    end

    @testset "Cell vertices, 3d" begin
        """
        1───────────2───────────3───────────4
        │           │           │           │
        │ (1, 4, 1) │ (2, 4, 1) │ (3, 4, 1) │
        │           │           │           │
        10─────────11───────────12─────────10
        │           │           │           │
        │ (1, 3, 1) │ (2, 3, 1) │ (3, 3, 1) │
        │           │           │           │
        7───────────8───────────9───────────7
        │           │           │           │
        │ (1, 2, 1) │ (2, 2, 1) │ (3, 2, 1) │
        │           │           │           │
        4───────────5───────────6───────────4
        │           │           │           │
        │ (1, 1, 1) │ (2, 1, 1) │ (3, 1, 1) │
        │           │           │           │
        1───────────2───────────3───────────1


        13─────────14───────────15─────────13
        │           │           │           │
        │ (1, 4, 2) │ (2, 4, 2) │ (3, 4, 2) │
        │           │           │           │
        22─────────23───────────24─────────22
        │           │           │           │
        │ (1, 3, 2) │ (2, 3, 2) │ (3, 3, 2) │
        │           │           │           │
        19─────────20───────────21─────────19
        │           │           │           │
        │ (1, 2, 2) │ (2, 2, 2) │ (3, 2, 2) │
        │           │           │           │
        16─────────17───────────18─────────16
        │           │           │           │
        │ (1, 1, 2) │ (2, 1, 2) │ (3, 1, 2) │
        │           │           │           │
        13─────────14───────────15─────────13


        25─────────26───────────27─────────25
        │           │           │           │
        │ (1, 4, 3) │ (2, 4, 3) │ (3, 4, 3) │
        │           │           │           │
        34─────────35───────────36─────────34
        │           │           │           │
        │ (1, 3, 3) │ (2, 3, 3) │ (3, 3, 3) │
        │           │           │           │
        31─────────32───────────33─────────31
        │           │           │           │
        │ (1, 2, 3) │ (2, 2, 3) │ (3, 2, 3) │
        │           │           │           │
        28─────────29───────────30─────────28
        │           │           │           │
        │ (1, 1, 3) │ (2, 1, 3) │ (3, 1, 3) │
        │           │           │           │
        25─────────26───────────27─────────25


        37─────────38───────────39─────────37
        │           │           │           │
        │ (1, 4, 4) │ (2, 4, 4) │ (3, 4, 4) │
        │           │           │           │
        46─────────47───────────48─────────46
        │           │           │           │
        │ (1, 3, 4) │ (2, 3, 4) │ (3, 3, 4) │
        │           │           │           │
        43─────────44───────────45─────────43
        │           │           │           │
        │ (1, 2, 4) │ (2, 2, 4) │ (3, 2, 4) │
        │           │           │           │
        40─────────41───────────42─────────40
        │           │           │           │
        │ (1, 1, 4) │ (2, 1, 4) │ (3, 1, 4) │
        │           │           │           │
        37─────────38───────────39─────────37


        49─────────50───────────51─────────49
        │           │           │           │
        │ (1, 4, 5) │ (2, 4, 5) │ (3, 4, 5) │
        │           │           │           │
        58─────────59───────────60─────────58
        │           │           │           │
        │ (1, 3, 5) │ (2, 3, 5) │ (3, 3, 5) │
        │           │           │           │
        55─────────56───────────57─────────55
        │           │           │           │
        │ (1, 2, 5) │ (2, 2, 5) │ (3, 2, 5) │
        │           │           │           │
        52─────────53───────────54─────────52
        │           │           │           │
        │ (1, 1, 5) │ (2, 1, 5) │ (3, 1, 5) │
        │           │           │           │
        49─────────50───────────51─────────49


        1───────────2───────────3───────────4
        │           │           │           │
        │           │           │           │
        │           │           │           │
        10─────────11───────────12─────────10
        │           │           │           │
        │           │           │           │
        │           │           │           │
        7───────────8───────────9───────────7
        │           │           │           │
        │           │           │           │
        │           │           │           │
        4───────────5───────────6───────────4
        │           │           │           │
        │           │           │           │
        │           │           │           │
        1───────────2───────────3───────────1

        """
        N = [3, 4, 5]
        d = size(N, 1)
        dims = Tuple(1:n for n in N)
        linear = LinearIndices(dims)
        @test cell_vertices((1, 1, 1), linear) == [1, 2, 4, 5, 13, 14, 16, 17]
        @test cell_vertices((2, 1, 1), linear) == [2, 3, 5, 6, 14, 15, 17, 18]
        @test cell_vertices((3, 1, 1), linear) == [3, 1, 6, 4, 15, 13, 18, 16]
        @test cell_vertices((1, 2, 1), linear) == [4, 5, 7, 8, 16, 17, 19, 20]
        @test cell_vertices((2, 2, 1), linear) == [5, 6, 8, 9, 17, 18, 20, 21]
        @test cell_vertices((3, 2, 1), linear) == [6, 4, 9, 7, 18, 16, 21, 19]
        @test cell_vertices((1, 3, 1), linear) == [7, 8, 10, 11, 19, 20, 22, 23]
        @test cell_vertices((2, 3, 1), linear) == [8, 9, 11, 12, 20, 21, 23, 24]
        @test cell_vertices((3, 3, 1), linear) == [9, 7, 12, 10, 21, 19, 24, 22]
        @test cell_vertices((1, 4, 1), linear) == [10, 11, 1, 2, 22, 23, 13, 14]
        @test cell_vertices((2, 4, 1), linear) == [11, 12, 2, 3, 23, 24, 14, 15]
        @test cell_vertices((3, 4, 1), linear) == [12, 10, 3, 1, 24, 22, 15, 13]

        @test cell_vertices((1, 1, 2), linear) == [13, 14, 16, 17, 25, 26, 28, 29]
        @test cell_vertices((2, 1, 2), linear) == [14, 15, 17, 18, 26, 27, 29, 30]
        @test cell_vertices((3, 1, 2), linear) == [15, 13, 18, 16, 27, 25, 30, 28]
        @test cell_vertices((1, 2, 2), linear) == [16, 17, 19, 20, 28, 29, 31, 32]
        @test cell_vertices((2, 2, 2), linear) == [17, 18, 20, 21, 29, 30, 32, 33]
        @test cell_vertices((3, 2, 2), linear) == [18, 16, 21, 19, 30, 28, 33, 31]
        @test cell_vertices((1, 3, 2), linear) == [19, 20, 22, 23, 31, 32, 34, 35]
        @test cell_vertices((2, 3, 2), linear) == [20, 21, 23, 24, 32, 33, 35, 36]
        @test cell_vertices((3, 3, 2), linear) == [21, 19, 24, 22, 33, 31, 36, 34]
        @test cell_vertices((1, 4, 2), linear) == [22, 23, 13, 14, 34, 35, 25, 26]
        @test cell_vertices((2, 4, 2), linear) == [23, 24, 14, 15, 35, 36, 26, 27]
        @test cell_vertices((3, 4, 2), linear) == [24, 22, 15, 13, 36, 34, 27, 25]

        @test cell_vertices((1, 1, 3), linear) == [25, 26, 28, 29, 37, 38, 40, 41]
        @test cell_vertices((2, 1, 3), linear) == [26, 27, 29, 30, 38, 39, 41, 42]
        @test cell_vertices((3, 1, 3), linear) == [27, 25, 30, 28, 39, 37, 42, 40]
        @test cell_vertices((1, 2, 3), linear) == [28, 29, 31, 32, 40, 41, 43, 44]
        @test cell_vertices((2, 2, 3), linear) == [29, 30, 32, 33, 41, 42, 44, 45]
        @test cell_vertices((3, 2, 3), linear) == [30, 28, 33, 31, 42, 40, 45, 43]
        @test cell_vertices((1, 3, 3), linear) == [31, 32, 34, 35, 43, 44, 46, 47]
        @test cell_vertices((2, 3, 3), linear) == [32, 33, 35, 36, 44, 45, 47, 48]
        @test cell_vertices((3, 3, 3), linear) == [33, 31, 36, 34, 45, 43, 48, 46]
        @test cell_vertices((1, 4, 3), linear) == [34, 35, 25, 26, 46, 47, 37, 38]
        @test cell_vertices((2, 4, 3), linear) == [35, 36, 26, 27, 47, 48, 38, 39]
        @test cell_vertices((3, 4, 3), linear) == [36, 34, 27, 25, 48, 46, 39, 37]

        @test cell_vertices((1, 1, 4), linear) == [37, 38, 40, 41, 49, 50, 52, 53]
        @test cell_vertices((2, 1, 4), linear) == [38, 39, 41, 42, 50, 51, 53, 54]
        @test cell_vertices((3, 1, 4), linear) == [39, 37, 42, 40, 51, 49, 54, 52]
        @test cell_vertices((1, 2, 4), linear) == [40, 41, 43, 44, 52, 53, 55, 56]
        @test cell_vertices((2, 2, 4), linear) == [41, 42, 44, 45, 53, 54, 56, 57]
        @test cell_vertices((3, 2, 4), linear) == [42, 40, 45, 43, 54, 52, 57, 55]
        @test cell_vertices((1, 3, 4), linear) == [43, 44, 46, 47, 55, 56, 58, 59]
        @test cell_vertices((2, 3, 4), linear) == [44, 45, 47, 48, 56, 57, 59, 60]
        @test cell_vertices((3, 3, 4), linear) == [45, 43, 48, 46, 57, 55, 60, 58]
        @test cell_vertices((1, 4, 4), linear) == [46, 47, 37, 38, 58, 59, 49, 50]
        @test cell_vertices((2, 4, 4), linear) == [47, 48, 38, 39, 59, 60, 50, 51]
        @test cell_vertices((3, 4, 4), linear) == [48, 46, 39, 37, 60, 58, 51, 49]

        @test cell_vertices((1, 1, 5), linear) == [49, 50, 52, 53, 1, 2, 4, 5]
        @test cell_vertices((2, 1, 5), linear) == [50, 51, 53, 54, 2, 3, 5, 6]
        @test cell_vertices((3, 1, 5), linear) == [51, 49, 54, 52, 3, 1, 6, 4]
        @test cell_vertices((1, 2, 5), linear) == [52, 53, 55, 56, 4, 5, 7, 8]
        @test cell_vertices((2, 2, 5), linear) == [53, 54, 56, 57, 5, 6, 8, 9]
        @test cell_vertices((3, 2, 5), linear) == [54, 52, 57, 55, 6, 4, 9, 7]
        @test cell_vertices((1, 3, 5), linear) == [55, 56, 58, 59, 7, 8, 10, 11]
        @test cell_vertices((2, 3, 5), linear) == [56, 57, 59, 60, 8, 9, 11, 12]
        @test cell_vertices((3, 3, 5), linear) == [57, 55, 60, 58, 9, 7, 12, 10]
        @test cell_vertices((1, 4, 5), linear) == [58, 59, 49, 50, 10, 11, 1, 2]
        @test cell_vertices((2, 4, 5), linear) == [59, 60, 50, 51, 11, 12, 2, 3]
        @test cell_vertices((3, 4, 5), linear) == [60, 58, 51, 49, 12, 10, 3, 1]
    end
end

# @testset "bri17, 2D" begin
#     C = Hooke{Float64, 2}(5.6, 0.3)
#     N = [2, 4]
#     h = [1.1, 1.2]
# end
