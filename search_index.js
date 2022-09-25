var documenterSearchIndex = {"docs":
[{"location":"appendix/#Appendix","page":"Appendix","title":"Appendix","text":"","category":"section"},{"location":"appendix/#On-Fourier-series","page":"Appendix","title":"On Fourier series","text":"","category":"section"},{"location":"appendix/","page":"Appendix","title":"Appendix","text":"We consider the domain Ω = (0, L[1]) × … × (0, L[d]) and a L-periodic vector field u: Ω → ℝᵐ. The Fourier coefficients of u are defined as follows, for all k ∈ ℝᵈ","category":"page"},{"location":"appendix/","page":"Appendix","title":"Appendix","text":"û(k) = |Ω|⁻¹ ∫ u(x) exp(-i k ⋅ x) dx.\n             Ω","category":"page"},{"location":"appendix/","page":"Appendix","title":"Appendix","text":"The following synthesis formula then holds (under non-specified, but non-restrictive, regularity conditions)","category":"page"},{"location":"appendix/","page":"Appendix","title":"Appendix","text":"u(x) =   ∑    û(kₚ) exp(i kₚ ⋅ x),\n       p ∈ ℤᵈ","category":"page"},{"location":"appendix/","page":"Appendix","title":"Appendix","text":"where kₚ is the following wave-vector","category":"page"},{"location":"appendix/","page":"Appendix","title":"Appendix","text":"        2π p[i]\nkₚ[i] = ───────,    for all    i ∈ {1, …, d}    and    p ∈ ℤᵈ.\n          L[i]","category":"page"},{"location":"api/#API-of-the-Scapin-library","page":"Library","title":"API of the Scapin library","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [Scapin]","category":"page"},{"location":"api/#Scapin.AbstractContinuousSpatialOperator","page":"Library","title":"Scapin.AbstractContinuousSpatialOperator","text":"AbstractContinuousSpatialOperator\n\nThis is the root type of all operators that have not been discretized spatially. Such operators map continuous spatial fields x : ℝᵈ → ℝⁿ to vector fields y : ℝᵈ → ℝᵐ. In Fourier space, these operators are defined for an infinite set of wave-numbers.\n\n\n\n\n\n","category":"type"},{"location":"api/#Scapin.AbstractDiscreteSpatialOperator","page":"Library","title":"Scapin.AbstractDiscreteSpatialOperator","text":"AbstractDiscreteSpatialOperator\n\nThis abstract type defines operators that are spatially discretized.\n\nSuch operators map vector fields x : ℝᵈ → ℝⁿ to vector fields y : ℝᵈ → ℝᵐ.\n\nThis type extends AbstractSpatialOperator. As such, it implements dimensionality, size and ndims.\n\nDiscrete operators operate on discrete vector fields that are represented as arrays. The size of the input arrays is (n, N[1], …, N[d]), where N[i] is the number of grid cells in direction i = 1, …, d. Note that the first (fast) index of the array is the component of the vector field. Similarly, the size of the output arrays is (m, N[1], …, N[d]).\n\nThe tuple N is referred to as the grid size, see grid_size.\n\n\n\n\n\n","category":"type"},{"location":"api/#Scapin.AbstractSpatialOperator","page":"Library","title":"Scapin.AbstractSpatialOperator","text":"AbstractSpatialOperator\n\nThis is the root type of all operators defined in this library.\n\nSuch operators map vector fields x : ℝᵈ → ℝⁿ to vector fields y : ℝᵈ → ℝᵐ.\n\nThe number of spatial dimensions, d, is referred to as the dimensionality of the operator – see dimensionality.\n\nThe pair (m, n) == (dimension of the input vector field, dimension of the output vector field) is referred to as the size of the operator – see size. Consequently, ndims(F) == 2 for all F::AbstractSpatialOperator.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.ndims-Tuple{Scapin.AbstractSpatialOperator}","page":"Library","title":"Base.ndims","text":"ndims(F::AbstractSpatialOperator) -> 2\n\nReturn the number of dimensions of F – see AbstractSpatialOperator. Since this is a linear operator that maps fields onto fields, the number of dimensions is always 2.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.apply-Tuple{Any, Any}","page":"Library","title":"Scapin.apply","text":"apply(ℱ, x)\n\nReturn the grid y that results from applying the convolution operator ℱ to the input grid x.\n\nThe sizes of x and y are\n\nsize(x) == (N..., ncols)\nsize(y) == (N..., nrows)\n\nwhere\n\nN = grid_size(ℱ)\n\nand\n\n(nrows, ncols) == size(ℱ).\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.apply_fourier!","page":"Library","title":"Scapin.apply_fourier!","text":"apply_fourier!(ŷ, F, k, x̂) -> ŷ\n\nApply in-place operator F to x̂ in Fourier space and return the modified array ŷ. More specifically,\n\nŷ[:] = F̂(k) ⋅ x̂,\n\nwhere k, x̂ and ŷ are column vectors: k is the vector of spatial frequencies and x̂ (resp. ŷ) is the Fourier mode of the input (resp. output), for the spatial frequency k. The following must hold\n\nsize(k, 1) == dimensionality(F),\nsize(x̂, 1) == size(F, 2),\nsize(ŷ, 1) == size(F, 1).\n\nThe input vector x̂ is an array of type T or complex(T), where T == eltype(ℱ). The output vector ŷ can always be of type complex(T). When the Fourier coefficient ℱ̂(k) of the operator is real, it might make sense to allow for ŷ to be of type real(T). If this is disallowed, calling this function should raise an exception.\n\nnote: Continuous and discrete operators\nFor F::AbstractContinuousSpatialOperator,kwould typically be a real-valued vector of wave numbers, while forF::AbstractDiscreteSpatialOperator,k` would be a integer-valued vector of indices.\n\ndanger: Overriding the input vector\nThe present method must be implemented in such a way that aliasing ŷ with x̂ is permitted. In other words, calling apply_fourier!(x̂, ℱ, k, x̂) must always deliver the correct answer, provided that this operation is allowed type-wise [when x̂ is of type real(T)].\n\n\n\n\n\n","category":"function"},{"location":"api/#Scapin.apply_fourier-Tuple{Scapin.AbstractSpatialOperator, Any, Any}","page":"Library","title":"Scapin.apply_fourier","text":"apply_fourier(F, k, x̂)\n\nApply operator F to x̂ in Fourier space and return the result. The returned vector is promoted to a complex type if necessary.\n\nSee apply_fourier!\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.dimensionality-Tuple{Scapin.AbstractSpatialOperator}","page":"Library","title":"Scapin.dimensionality","text":"dimensionality(F)\n\nReturn the number of dimensions of the physical space that F operates on – see AbstractSpatialOperator.\n\nThe default implementation assumes that dimensionality(typeof(F)) has been defined.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.fourier_matrix!-Tuple{Any, Any, Any}","page":"Library","title":"Scapin.fourier_matrix!","text":"fourier_matrix!(F̂, ℱ, k) -> F̂\n\nCompute in-place the k-th mode of ℱ, ℱ̂(k), as a matrix F̂.\n\nIn general, the coefficients of F̂ should be s of type complex(eltype(ℱ)). However, the Fourier coefficients of the operator ℱ might be real, in which case, eltype(F̂) == real(eltype(ℱ)) might be allowed.\n\nnote: Performance of the default implementation\nThe default implementation relies on apply_fourier! to build the matrix column-by-column. It might be inefficient.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.fourier_matrix-Tuple{Any, Any}","page":"Library","title":"Scapin.fourier_matrix","text":"fourier_matrix(ℱ, k)\n\nReturn the k-th mode of ℱ, ℱ̂(k), as a matrix.\n\nThe coefficients of the returned matrix are of type complex(eltype(ℱ)). For a different output type, use fourier_matrix!.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.grid_size-Tuple{Scapin.AbstractDiscreteSpatialOperator}","page":"Library","title":"Scapin.grid_size","text":"grid_size(F)\n\nReturn the size of the grid the discrete operator F operates on.\n\nSize of input arrays of F: (size(F, 2), grid_size(F)...),\nSize of output arrays of F: (size(F, 1), grid_size(F)...).\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Bri17","page":"Library","title":"Scapin.Bri17","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [Scapin.Bri17]","category":"page"},{"location":"api/#Scapin.Bri17.DiscreteGreenOperatorBri17","page":"Library","title":"Scapin.Bri17.DiscreteGreenOperatorBri17","text":"Fourier representation of the discrete Green operator introduced in [Bri17].\n\n\n\n\n\n","category":"type"},{"location":"api/#Scapin.Bri17.modal_stiffness!-Union{Tuple{T}, Tuple{d}, Tuple{AbstractArray{Complex{T}, 2}, CartesianIndex{d}, Tuple{Vararg{Int64, d}}, Tuple{Vararg{T, d}}, Scapin.Elasticity.Hooke{d, T}}} where {d, T<:Number}","page":"Library","title":"Scapin.Bri17.modal_stiffness!","text":"modal_stiffness!(K, n, N, h, C)\n\nCompute the modal stiffness matrix K in place, for the spatial frequency n.\n\nThis function returns K.\n\nThe grid is defined by its size, N and spacing, h. The spatial frequency is defined by n ∈ CartesianIndices(1:N[1], …, 1:N[d]). The (homogeneous) constitutive material is specified by C.\n\ndanger: Scaling of the modal stiffness matrix\nThe modal stiffness matrix introduced above differs from the matrix initially introduced in Ref. DOI:10.1002/nme.5263 by a factor h[1] * … * h[d]. More precisely,K̂{Scapin} = h[1] * … * h[d] * K̂{10.1002/nme.5263}Such rescaling makes the relation between modal and nodal stiffness operators more natural (the former is the DFT of the latter).\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Bri17.modal_stiffness-Union{Tuple{d}, Tuple{T}, Tuple{CartesianIndex{d}, Tuple{Vararg{Int64, d}}, Tuple{Vararg{T, d}}, Scapin.Elasticity.Hooke{d, T}}} where {T<:Number, d}","page":"Library","title":"Scapin.Bri17.modal_stiffness","text":"modal_stiffness(n, N, h, C)\n\nReturn the modal stiffness matrix K.\n\nSee modal_stiffness! for a description of the parameters.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Bri17.modal_strain_displacement!-Union{Tuple{T}, Tuple{d}, Tuple{AbstractArray{Complex{T}, 1}, CartesianIndex{d}, Tuple{Vararg{Int64, d}}, Tuple{Vararg{T, d}}}} where {d, T<:Number}","page":"Library","title":"Scapin.Bri17.modal_strain_displacement!","text":"modal_strain_displacement!(b, n, N, h)\n\nCompute the modal strain-displacement vector b in place, for the spatial frequency n.\n\nThis function returns b.\n\nThe grid is defined by its size, N and spacing, h. The spatial frequency is defined by n ∈ CartesianIndices(1:N[1], …, 1:N[d]).\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Bri17.modal_strain_displacement-Union{Tuple{T}, Tuple{d}, Tuple{CartesianIndex{d}, Tuple{Vararg{Int64, d}}, Tuple{Vararg{T, d}}}} where {d, T<:Number}","page":"Library","title":"Scapin.Bri17.modal_strain_displacement","text":"modal_strain_displacement(n, N, h)\n\nReturn the modal strain-displacement vector b.\n\nSee modal_strain_displacement! for a description of the parameters.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Brick","page":"Library","title":"Scapin.Brick","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [Scapin.Brick]","category":"page"},{"location":"api/#Scapin.Brick.avg_gradient_operator-Union{Tuple{Tuple{Vararg{T, d}}}, Tuple{T}, Tuple{d}} where {d, T<:Number}","page":"Library","title":"Scapin.Brick.avg_gradient_operator","text":"avg_gradient_operator(h)\n\nReturn the cell-average of the gradient operator.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Brick.avg_strain_displacement_operator-Union{Tuple{Tuple{Vararg{T, d}}}, Tuple{T}, Tuple{d}} where {d, T<:Number}","page":"Library","title":"Scapin.Brick.avg_strain_displacement_operator","text":"avg_strain_displacement_operator(h)\n\nReturn the cell average of the strain-displacement operator.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Brick.global_stiffness_operator-Union{Tuple{T}, Tuple{d}, Tuple{Tuple{Vararg{Int64, d}}, Tuple{Vararg{T, d}}, Scapin.Elasticity.Hooke{d, T}}} where {d, T<:Number}","page":"Library","title":"Scapin.Brick.global_stiffness_operator","text":"global_stiffness_operator(N, h, μ, ν)\n\nReturn the global stiffness operator for periodic, homogeneous elasticity.\n\nThe grid size is N, the cell size is h. The constitutive material is homogeneous, elastic linear and isotropic with stiffness C.\n\nThe global stiffness operator K is a 2d+2-dimensional array of size (d, N[1], …, N[d], d, N[1], …, N[d]), such that the strain energy of the system reads\n\nU = u[i, p] * K[i, p, j, q] * u[j, q] / 2,\n\nwhere\n\ni, j ∈ {1, …, d}: component indices,\np, q ∈ CartesianIndices(1:N[1], …, 1:N[d]): node indices,\nu[i, p]: i-th component of the displacement of node p.\n\nnote: Note\nAssembly of the global stiffness opertor is done under the assumption of periodicity.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Brick.global_strain_displacement_operator-Union{Tuple{T}, Tuple{d}, Tuple{Tuple{Vararg{Int64, d}}, Tuple{Vararg{T, d}}}} where {d, T<:Number}","page":"Library","title":"Scapin.Brick.global_strain_displacement_operator","text":"global_strain_displacement_operator(N, h)\n\nReturn the global strain-displacement operator for periodic boundary conditions.\n\nThe grid is defined by its size N and its spacing h.\n\nThe global strain-displacement operator B is a d+2-dimensional array of size (d, d, N[1], …, N[d]), such that the average strain within element p reads\n\nε[i, j, p] = B[i, j, p, k, q] * u[k, q],\n\nwhere\n\ni, j, k ∈ {1, …, d}: component indices,\np, q ∈ CartesianIndices(1:N[1], …, 1:N[d]): node indices,\nε[i, j, p]: (i, j)-th component of the average strain in element p,\nu[k, q]: i-th component of the displacement of node q.\n\nnote: Note\nAssembly of the global strain-displacement operator is done under the assumption of periodicity.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Brick.gradient_operator-Union{Tuple{T}, Tuple{d}, Tuple{Tuple{Vararg{T, d}}, Tuple{Vararg{T, d}}}} where {d, T<:Number}","page":"Library","title":"Scapin.Brick.gradient_operator","text":"gradient_operator(x, h)\n\nReturn the gradient operator at the specified point.\n\nThis function returns a (d+1) dimensional array D of size (d, 2, 2, …). If  i is the index of a component and p the CartesianIndex of the node, then D[i, p] is the partial derivative of N[p] w.r.t. x[i], evaluated at x.\n\nh is the size of the brick element.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Brick.integrate-Union{Tuple{T}, Tuple{d}, Tuple{Any, Tuple{Vararg{T, d}}}} where {d, T<:Number}","page":"Library","title":"Scapin.Brick.integrate","text":"integrate(f, h)\n\nReturn the d-dimensional integral of f over (0, h[1]) × (0, h[2]) × … × (0, h[d]).\n\nUses 2-point Gauss-Legendre integration (tensorized over the d dimensions). f must take a 1-dimensional array of size d as unique input. If avg is true, the function returns the N-dimensional average.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Brick.shape-Union{Tuple{T}, Tuple{d}, Tuple{Tuple{Vararg{T, d}}, Tuple{Vararg{T, d}}}} where {d, T<:Number}","page":"Library","title":"Scapin.Brick.shape","text":"shape(x, h)\n\nReturn the value of the shape functions of the element, at the specified point.\n\nThis function returns a d-dimensional array N, such that N[p] is the shape function associated with node p, evaluated at x. In particular, N[p] evaluated at node q is δ[p, q] (Kronecker).\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Brick.stiffness_operator-Union{Tuple{T}, Tuple{d}, Tuple{Tuple{Vararg{T, d}}, Scapin.Elasticity.Hooke{d, T}}} where {d, T<:Number}","page":"Library","title":"Scapin.Brick.stiffness_operator","text":"stiffness_operator(h, C)\n\nReturn the stifness operator for the brick element of size h, and Hooke material C.\n\nThe stiffness operator K delivers the strain energy associated to the nodal displacements u\n\nU = u[i, p] * K[i, p, j, q] * u[j, q] / 2,\n\nwhere i, j ∈ {1, …, d} are component indices and p, q ∈ CartesianIndices(1:2, …, 1:2).\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Brick.strain_displacement_operator-Union{Tuple{T}, Tuple{d}, Tuple{Tuple{Vararg{T, d}}, Tuple{Vararg{T, d}}}} where {d, T<:Number}","page":"Library","title":"Scapin.Brick.strain_displacement_operator","text":"strain_displacement_operator(x, h)\n\nReturn the strain-displacement operator for the d-dimensional brick element of size h, evaluated at point x.\n\nThis function returns a (d+3) dimensional array B of size (d, d, 2, …, 2, d). If p is the CartesianIndex of the node, and i, j, k are component indices then, the interpolated (i, j) component of the strain at x reads\n\nε[i, j] = Σₖ Σₚ B[i, j, k, p] * u[k, p].\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Elasticity","page":"Library","title":"Scapin.Elasticity","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [Scapin.Elasticity]","category":"page"},{"location":"api/#Scapin.Elasticity.GreenOperatorHooke","page":"Library","title":"Scapin.Elasticity.GreenOperatorHooke","text":"Continuous Green operator related to Hooke materials.\n\n\n\n\n\n","category":"type"},{"location":"api/#Scapin.Elasticity.Hooke","page":"Library","title":"Scapin.Elasticity.Hooke","text":"Hooke{d,T}\n\nIsotropic, linear elastic material.\n\nd — number of spatial dimensions\nT — scalar type\nμ::T — shear modulus\nν::T — Poisson ratio\nλ::T — first Lamé coefficient\n\nRegardless of the number of space dimensions, d, the stress-strain relationship reads\n\nσ = λ⋅tr(ε)I + 2μ⋅ε.\n\nnote: Material stability\nMaterial stability requires that μ > 0 and -1 < ν < 1/2; these conditions are not enforced here. In other words, unstable materials can be defined.\n\ntip: Plane stresses vs. plane strains\nIn the current implementation, d = 2 refers to plane strain elasticity. For plane stresses, the true Poisson ratio ν should be replaced with the fictitious ratio ν̃ = ν / (1 + ν).\n\n\n\n\n\n","category":"type"},{"location":"api/#Scapin.Elasticity.Hooke-Union{Tuple{T}, Tuple{d}, Tuple{T, T}} where {d, T}","page":"Library","title":"Scapin.Elasticity.Hooke","text":"Hooke{d,T}(μ::T, ν::T)\n\nCreate a new instance of Hooke{d,T} with shear modulus μ and Poisson ratio ν.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.eltype-Union{Tuple{Type{Scapin.Elasticity.GreenOperatorHooke{d, T}}}, Tuple{T}, Tuple{d}} where {d, T}","page":"Library","title":"Base.eltype","text":"eltype(type)\n\nReturn the type of the (scalar) elements the operator of given type operates on in the real space (as opposed to the Fourier space). Note that this type may be complex!\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Elasticity.bulk_modulus-Tuple{Scapin.Elasticity.Hooke}","page":"Library","title":"Scapin.Elasticity.bulk_modulus","text":"bulk_modulus(C::Hooke)\n\nReturn the bulk modulus κ for the specified Hooke material.\n\nFor plane strain elasticity, κ = μ / (1 - 2ν) and, for 3d elasticity κ = 2/3 μ (1 + ν) / (1 - 2ν).\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Grid","page":"Library","title":"Scapin.Grid","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [Scapin.Grid]","category":"page"},{"location":"api/#Scapin.Grid.cell_vertices-Tuple{Int64}","page":"Library","title":"Scapin.Grid.cell_vertices","text":"cell_vertices(d)\n\nReturn the local multi-indices of vertices of the d-dimensional brick element.\n\nThis function returns an object 𝓛 of type CartesianIndices. An element v ∈ 𝓛 represents the vertex with coordinates (x[1], …, x[d])\n\nx[i] = (-1)^v[i] * h[i] / 2,    i = 1, …, d\n\nSee Geometry of the reference brick element.\n\n\n\n\n\n","category":"method"},{"location":"api/#Scapin.Grid.cell_vertices-Union{Tuple{d}, Tuple{CartesianIndex{d}, Tuple{Vararg{Int64, d}}}} where d","page":"Library","title":"Scapin.Grid.cell_vertices","text":"cell_vertices(p, N)\n\nReturn the global multi-indices of the vertices of a grid-cell.\n\nThe cell is specified through its CartesianIndex, p; the size of the grid is defined by N.\n\nThis function returns a d-dimensional array 𝒢 of CartesianIndex, such that 𝒢[l] is the global index of the vertex with local index l ∈ CartesianIndices(1:2, …, 1:2).\n\nNote that periodic boundary conditions are used.\n\n\n\n\n\n","category":"method"},{"location":"design/#On-the-design-of-the-Scapin-library","page":"Design","title":"On the design of the Scapin library","text":"","category":"section"},{"location":"design/#On-discrete-Green-operators","page":"Design","title":"On discrete Green operators","text":"","category":"section"},{"location":"design/#Real-vs.-complex-Fourier-transforms","page":"Design","title":"Real vs. complex Fourier transforms","text":"","category":"section"},{"location":"design/","page":"Design","title":"Design","text":"The interface for discrete Green operators should allow for switching between real and complex Fourier transforms.","category":"page"},{"location":"design/","page":"Design","title":"Design","text":"As of rev. 7bd9a4e, two functions are implemented, namely eltype_real and eltype_fourier. This might be overkill. Indeed, it is sufficient to define the scalar type in the Real space, T.","category":"page"},{"location":"convolution_operators/#On-convolution-operators","page":"Convolution operators","title":"On convolution operators","text":"","category":"section"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"In the present chapter, we define what we mean by convolution operators, and how they relate to the Julia library LinearOperators.jl.","category":"page"},{"location":"convolution_operators/#General-definitions:-vector-field,-spatial-operator","page":"Convolution operators","title":"General definitions: vector field, spatial operator","text":"","category":"section"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"We first define a vector field as a mapping from Ω to ℝᵐ, where Ω will be referred to as the real space or direct space (as opposed to the reciprocal space). The dimension m of the vector space is the size of the vector field.","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"tip: Tensor fields\nTensor fields can be considered, provided that they are first reduced to vector fields (see e.g. Voigt notation).","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"The field under consideration is discrete if Ω is a subset of ℕᵈ; typically","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"Ω = {1, 2, …, N[1]} × … × {1, 2, …, N[d]}.","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"Conversely, the field is continuous if Ω is a subset of ℝᵈ. Typically","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"Ω = (0, L[1]) × … × (0, L[d]).","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"In both cases, the number of spatial dimensions, d, will be referred to as the dimensionality of the field.","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"danger: Terminology\n“Continuous” is a very poor qualifier, since it does not relate to the regularity of the field itself.","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"An operator is then defined as a function that maps fields onto fields. The source and destination fields might have different vector dimensions; however, it is required that they share the same dimensionality d, which will be called by extension the dimensionality of the operator. The pair","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"(m, n) == (size of the output field, size of the input field)","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"is referred to as the size of the operator, which is consistent with the definition of the size of a linear mapping.","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"Discrete operators operate on rectangular grids of size N (d-tuple). The input field is stored in a (d + 1)-dimensional array of size (m, N...) in Julian notation. Similarly, the output field is stored in a (d + 1)-dimensional array of size (n, N...). Note that the components of both fields correspond to the first (fast) index.","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"We are now ready to define convolution operators in a periodic setting. Note that boundary conditions other than periodic are possible, but are irrelevant to this library.","category":"page"},{"location":"convolution_operators/#Convolution-operators","page":"Convolution operators","title":"Convolution operators","text":"","category":"section"},{"location":"convolution_operators/#Continuous-convolution-operators","page":"Convolution operators","title":"Continuous convolution operators","text":"","category":"section"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"We consider the operator F, that maps the field u: Ω → ℝᵐ onto v: Ω → ℝⁿ, where Ω = (0, L[1]) × … × (0, L[d]) and u and v are both L-periodic; F is a continuous convolution operator if v can be expressed as follows","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"v(x) = F(u)(x) = ∫ f(x - y) ⋅ u(y) dy,\n                 Ω","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"where f is the convolution kernel: for all x ∈ Ω, f(x) is a linear mapping from ℝᵐ to ℝⁿ. It is customary to identify F and its kernel. In other words, we often write","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"v(x) = ∫ F(x - y) ⋅ u(y) dy.","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"Whether we mean the operator or its kernel will be clear from the argument of F. The above expression is of course conveniently written in Fourier space as follows","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"v̂(k) = F̂(k) ⋅ û(k)    for all    k ∈ ℝᵈ,","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"where û(k), v̂(k) and F̂(k) are the Fourier coefficients of the fields u and v and the kernel F, respectively (see On Fourier series). Then, from the synthesis formula","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"v(x) = F(u)(x) =   ∑    F̂(kₚ) ⋅ û(kₚ) exp(i kₚ ⋅ x),\n                 p ∈ ℤᵈ","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"which shows that the operator F is fully defined by the Fourier coefficients of its kernel. In julia, we only need to define the method apply_fourier! to fully define a ContinuousConvolutionOperator. This method is called as follows","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"apply_fourier!(v, F, k, u)","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"and computes in-place, for all u ∈ ℂᵐ and k ∈ ℝᵈ, the vector v ∈ ℂⁿ such that v = F̂(k) ⋅ u. Note that, in general, both u and v are complex-valued vectors. More on that below.","category":"page"},{"location":"convolution_operators/#Discrete-convolution-operators","page":"Convolution operators","title":"Discrete convolution operators","text":"","category":"section"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"We now consider the operator F, that maps the field u: Ω → ℝᵐ onto v: Ω → ℝⁿ, where Ω = {1, 2, …, N[1]} × … × {1, 2, …, N[d]} and u and v are both N-periodic; F is a continuous convolution operator if v can be expressed as follows","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"v(x) = ∑ F(x - y) ⋅ u(y)    for all    y ∈ Ω.","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"Again, we have identified the operator and its kernel.","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"The above expression is of course conveniently written in Fourier space as follows","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"v̂(k) = F̂(k) ⋅ û(k)    for all    k ∈ ℝᵈ,","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"where û(k), v̂(k) and F̂(k) are the Fourier coefficients of the fields u and v and the kernel F, respectively (see On Fourier series). Then, from the synthesis formula","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"v(x) = F(u)(x) =   ∑    F̂(kₚ) ⋅ û(kₚ) exp(i kₚ ⋅ x),\n                 p ∈ ℤᵈ","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"which shows that the operator F is fully defined by the Fourier coefficients of its kernel. In julia, we only need to define the method apply_fourier! to fully define a ContinuousConvolutionOperator. This method is called as follows","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"apply_fourier!(v, F, k, u)","category":"page"},{"location":"convolution_operators/","page":"Convolution operators","title":"Convolution operators","text":"and computes in-place, for all u ∈ ℂᵐ and k ∈ ℝᵈ, the vector v ∈ ℂⁿ such that v = F̂(k) ⋅ u. Note that, in general, both u and v are complex-valued vectors. More on that below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Scapin","category":"page"},{"location":"#Scapin","page":"Home","title":"Scapin","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Scapin.","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nFor a theoretical overview of Lippmann–Schwinger solvers, please refer to the book An introduction to Lippmann–Scwhinger solvers. This book was written alongside the present code, which therefore follows closely the concepts and notations introduced in the book.","category":"page"}]
}
