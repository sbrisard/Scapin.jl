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

@testset "integrate" begin
    h = [1.1, 1.2, 1.3]
    n = [1.0, 2.0, 3.0]
    f(x) = prod(x .^ n)
    actual = integrate(f, h)
    expected = prod((h .^ (n .+ 1)) ./ (n .+ 1))
    @test isapprox(actual, expected, rtol=1e-15)
end

# @testset "bri17, 2D" begin
#     C = Hooke{Float64, 2}(5.6, 0.3)
#     N = [2, 4]
#     h = [1.1, 1.2]
# end
