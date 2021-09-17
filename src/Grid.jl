module Grid

using Base.Iterators

"""
    cell_vertices(d)

Return the local multi-indices of vertices of the `d`-dimensional brick element.

This function returns an object `𝓛` of type `CartesianIndices`. An
element `v ∈ 𝓛` represents the vertex with coordinates `(x[1], …, x[d])`

```
x[i] = (-1)^v[i] * h[i] / 2,    i = 1, …, d
```

See “[Geometry of the reference brick element](@ref _20210910120306)” in the docs.
"""
cell_vertices(d::Int) = CartesianIndices(Tuple(fill(1:2, d)))


"""
    cell_vertices(p, N)

Return the global multi-indices of the vertices of a grid-cell.

The cell is specified through its `CartesianIndex`, `p`; the size of the grid is
defined by `N`.

This function returns a `d`-dimensional array `𝒢` of `CartesianIndex`, such that 𝒢[l] is
the global index of the vertex with local index `l ∈ CartesianIndices(1:2, …, 1:2)`.

Note that periodic boundary conditions are used.

"""
function cell_vertices(p::CartesianIndex{d}, N::NTuple{d, Int}) where d
    p_ = collect(Tuple(p))
    𝓛 = map(collect, product(fill(0:1, d)...))
    return [CartesianIndex(((p_ .+ q .- 1) .% N .+ 1)...) for q in 𝓛]
end

export cell_vertices

end
