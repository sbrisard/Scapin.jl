# On convolution operators

In the present chapter, we define what we mean by *convolution operators*, and
how they relate to the Julia library
[LinearOperators.jl](https://github.com/JuliaSmoothOptimizers/LinearOperators.jl).

## General definitions: vector field, spatial operator

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

Discrete operators operate on rectangular grids of size `N` (`d`-tuple). The
input field is stored in a `(d + 1)`-dimensional array of size `(m, N...)` in
Julian notation. Similarly, the output field is stored in a
`(d + 1)`-dimensional array of size `(n, N...)`. Note that the *components* of
both fields correspond to the first (fast) index.

We are now ready to define convolution operators in a periodic setting. Note
that boundary conditions other than periodic are possible, but are irrelevant to
this library.


## Convolution operators

### Continuous convolution operators

We consider the operator `F`, that maps the field `u: Ω → ℝᵐ` onto `v: Ω → ℝⁿ`,
where `Ω = (0, L[1]) × … × (0, L[d])` and `u` and `v` are both `L`-periodic; `F`
is a continuous convolution operator if `v` can be expressed as follows

    v(x) = F(u)(x) = ∫ f(x - y) ⋅ u(y) dy,
                     Ω

where `f` is the convolution kernel: for all `x ∈ Ω`, `f(x)` is a linear mapping
from `ℝᵐ` to `ℝⁿ`. It is customary to identify `F` and its kernel. In other
words, we often write

    v(x) = ∫ F(x - y) ⋅ u(y) dy.

Whether we mean the operator or its kernel will be clear from the argument of
`F`. The above expression is of course conveniently written in Fourier space as
follows

    v̂(k) = F̂(k) ⋅ û(k)    for all    k ∈ ℝᵈ,

where `û(k)`, `v̂(k)` and `F̂(k)` are the Fourier coefficients of the fields `u`
and `v` and the kernel `F`, respectively (see [On Fourier series](@ref)). Then,
from the synthesis formula

    v(x) = F(u)(x) =   ∑    F̂(kₚ) ⋅ û(kₚ) exp(i kₚ ⋅ x),
	                 p ∈ ℤᵈ

which shows that the operator `F` is fully defined by the Fourier coefficients
of its kernel. In julia, we only need to define the method
[`apply_fourier!`](@ref) to fully define a `ContinuousConvolutionOperator`. This
method is called as follows

    apply_fourier!(v, F, k, u)

and computes *in-place*, for all `u ∈ ℂᵐ` and `k ∈ ℝᵈ`, the vector `v ∈ ℂⁿ` such
that `v = F̂(k) ⋅ u`. Note that, in general, both `u` and `v` are complex-valued
vectors. More on that below.


### Discrete convolution operators

We now consider the operator `F`, that maps the field `u: Ω → ℝᵐ` onto
`v: Ω → ℝⁿ`, where `Ω = {1, 2, …, N[1]} × … × {1, 2, …, N[d]}` and `u` and `v`
are both `N`-periodic; `F` is a continuous convolution operator if `v` can be
expressed as follows

    v(x) = ∑ F(x - y) ⋅ u(y)    for all    y ∈ Ω.

Again, we have identified the operator and its kernel.



The above expression is
of course conveniently written in Fourier space as follows

    v̂(k) = F̂(k) ⋅ û(k)    for all    k ∈ ℝᵈ,

where `û(k)`, `v̂(k)` and `F̂(k)` are the Fourier coefficients of the fields `u`
and `v` and the kernel `F`, respectively (see [On Fourier series](@ref)). Then,
from the synthesis formula

    v(x) = F(u)(x) =   ∑    F̂(kₚ) ⋅ û(kₚ) exp(i kₚ ⋅ x),
	                 p ∈ ℤᵈ

which shows that the operator `F` is fully defined by the Fourier coefficients
of its kernel. In julia, we only need to define the method
[`apply_fourier!`](@ref) to fully define a `ContinuousConvolutionOperator`. This
method is called as follows

    apply_fourier!(v, F, k, u)

and computes *in-place*, for all `u ∈ ℂᵐ` and `k ∈ ℝᵈ`, the vector `v ∈ ℂⁿ` such
that `v = F̂(k) ⋅ u`. Note that, in general, both `u` and `v` are complex-valued
vectors. More on that below.
