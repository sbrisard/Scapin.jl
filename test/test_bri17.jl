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

# @testset "bri17, 2D" begin
#     C = Hooke{Float64, 2}(5.6, 0.3)
#     N = [2, 4]
#     h = [1.1, 1.2]
# end
