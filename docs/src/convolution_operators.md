# On convolution operators

In the present chapter, we define what we mean by *convolution operators*, and
how they relate to the Julia library
[LinearOperators.jl](https://github.com/JuliaSmoothOptimizers/LinearOperators.jl).

We first define a *vector field* as a mapping from `Ω` to `ℝᵐ`, where `Ω` will
be referred to as the *real space* or *direct space* (as opposed to the
[reciprocal space](https://en.wikipedia.org/wiki/Reciprocal_lattice)). The
dimension `m` of the vector space is the *size* of the vector field.

!!! tip "Tensor fields"

    Tensor fields can be considered, provided that they are first
	reduced to vector fields
	(see e.g. [Voigt notation](https://en.wikipedia.org/wiki/Voigt_notation)).

The field under consideration is *discrete* if `Ω` is a subset of `ℕᵈ`;
typically

    Ω = {1, 2, …, N[1]} × … × {1, 2, …, N[d]}.

Conversely, the field is *continuous* if `Ω` is a subset of `ℝᵈ`. Typically

    Ω = (0, L[1]) × … × (0, L[d]).

In both cases, the number of spatial dimensions, `d`, will be referred to as the
*dimensionality* of the field.

!!! danger "Terminology"

    “Continuous” is a very poor qualifier, since it does not relate to the
	regularity of the field itself.

An *operator* is then defined as a function that maps fields onto fields. The
source and destination fields might have different vector dimensions; however,
it is required that they share the *same dimensionality* `d`, which will be
called by extension the dimensionality of the operator. The pair

    (m, n) == (size of the output field, size of the input field)

is referred to as the *size* of the operator, which is consistent with the
definition of the size of a linear mapping.

We consider the operator `F`, that maps the field `u: ℝᵈ → ℝᵐ` onto
`v: ℝᵈ → ℝⁿ`: `v = F(u)`.
