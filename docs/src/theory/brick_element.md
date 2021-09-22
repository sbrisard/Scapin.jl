# On the d-dimensional brick element

In this chapter, we formulate the linear brick finite element for conductivity
and linear elasticity in ``d`` dimensions (``d\in\{2, 3\}``).

!!! note

    In the remainder of this chapter, ``i``, ``j``, ``h`` and ``k`` denote scalar
	indices that refer to tensor components, and span ``\{1, 2, \cdots d\}``;
	``\tuple{p}``, ``\tuple{q}`` are ``d``-dimensional multi-indices that refer
	to nodes. Depending on the context (local, element operators or global,
	assembled operators) they span ``\{1, 2\}^d`` or the whole grid
	``\cellindices``. We adopt Einstein's convention for indices of both types.

## Geometry of the reference brick element

The reference element being centered at the origin, it occupies the following
domain

```math
\Omega^\element=\biggl(-\frac{h_1}2, \frac{h_1}2\biggr)\times\cdots
\times\biggl(-\frac{h_d}2, \frac{h_d}2\biggr),
```

where ``h_i`` is the size of the element along axis ``i`` (``i=1, \dots,
d``). The nodes of the element are indexed by ``\tuple{p}\in\{1, 2\}^d``, such
that the coordinates of node ``\tuple{p}`` read

```math
x_{i\tuple{p}}=(-1)^{p_i}\frac{h_i}2.
```

For integration purposes, it is useful to introduce the *reduced coordinates*
``\xi_i=2x_i/h_i``

```math
\int_{\Omega^\element}f(x_1, \ldots, x_d)\,\D x_1\cdots\D x_d
=\frac{h_1\cdots h_d}{2^d}\int_{(-1, 1)^d}f(h_1\xi_1/2, \ldots, h_d\xi_d/2)
\,\D\xi_1\cdots\D\xi_d
```

and it is observed that the reduced coordinates of the node ``\tuple{p}`` are
``\xi_i=(-1)^{p_i}``.

## Shape functions

The shape function associated with node ``\tuple{p}`` is
``N^\element_{\tuple{p}}(\vec x)``, wich is defined below as a function of the
reduced coordinates of the observation point

```math
N^\element_{\tuple{p}}(\vec x)=\frac1{2^d}\prod_{i=1}^d\bigl[1+(-1)^{p_i}\xi_i\bigr]
```

and we have as expected ``N^\element_{\tuple{p}}(\vec x_{\tuple q}) =
\delta_{\tuple{p}\tuple{q}}``. The nodal values of ``f`` being ``f_\tuple{p}``,
we have the interpolation formula: ``f(\vec x)=N^\element_\tuple{p}(\vec
x)f_\tuple{p}``.

## Gradient and strain-displacement operators

We consider a scalar interpolated field, ``f(\vec x)=N^\element_\tuple{p}(\vec
x)f_\tuple{p}``. The components of its gradient are given by the following
expression

```math
\bigl(f\otimes\nabla\bigr)_i=\partial_i f=D_{i\tuple{p}}^\element f_\tuple{p},
```

where the observation point ``\vec x`` has been dropped, and
``D_{i\tuple{p}}^\element`` denotes the element gradient operator

```math
D_{i\tuple{p}}^\element=\partial_iN_\tuple{p}
=\frac{(-1)^{p_i}}{2^{d-1}h_i}
\prod_{j\neq i}\bigl[1+(-1)^{p_j}\xi_j\bigr].
```

It will be useful to introduce the element average of this operator,
``\overline{D}_{i\tuple{n}}^\element``

```math
\overline{D}_{i\tuple{p}}^\element=\frac 1{h_1\cdots h_d}
\int_{\Omega^\element}D^\element_{i\tuple{p}}(\vec x)\,\D x_1\cdots\D x_d
```

and, observing that each factor in square brackets in the above expression of
``D^\element_{i\tuple{p}}`` varies linearly between 0 et 2 (and therefore
averages to 1), we find

```math
\overline{D}^\element_{i\tuple{p}}=\frac{(-1)^{p_i}}{2^{d-1}h_i}.
```

We now consider the interpolated displacement field, ``\vec u =
N^\element_\tuple{p} \vec u_\tuple{p}``. The components of the associated
strain, ``\tens\varepsilon=\vec u\symotimes\nabla`` are

```math
\varepsilon_{ij}=\tfrac12\bigl(\partial_iu_j+\partial_ju_i\bigr)
=\tfrac12\bigl(D^\element_{i\tuple{p}}u_{j\tuple{p}}
+D^\element_{j\tuple{p}}u_{i\tuple{p}}\bigr),
```

where ``u_{i\tuple{p}}`` denotes the ``i``-th component of the nodal value of
``\vec u`` at node ``\tuple{p}``, ``u_{i\tuple{p}}=\vec e_i\cdot\vec u(\vec
x_\tuple{p})``.

!!! note

    In the Julia implementation of the brick element, the nodal displacements
	are stored in a ``d+1``-dimensional array `u`, such that
	`u[i, p]```=u_{i\tuple{p}}``, where `i ∈ {1, ... d}` and
	`p ∈ CartesianIndices(1:N[1], …, 1:N[d])`. With this convention, the
	component index, `i`, is the *fast* index. Therefore, all dofs of the same
	node are contiguous, which is the usual convention in finite element.

Rearranging the above expression, we find
``\varepsilon_{ij}=B^\element_{ijk\tuple{p}}u_{k\tuple{p}}^\element``, where

```math
B_{ijk\tuple{p}}^\element
=\tfrac12\bigl(D^\element_{i\tuple{p}}\delta_{jk}
+D^\element_{j\tuple{p}}\delta_{ik}\bigr)
```

is the element strain-displacement operator, with volume average

```math
\overline{B}_{ijk\tuple{p}}^\element
=\tfrac12\bigl(\overline{D}^\element_{i\tuple{p}}\delta_{jk}
+\overline{D}^\element_{j\tuple{p}}\delta_{ik}\bigr)
=\frac1{2^d}\Bigl[\frac{(-1)^{p_i}\delta_{jk}}{h_i}
+\frac{(-1)^{p_j}\delta_{ik}}{h_j}\Bigr].
```

It is obverved that the *trace* of the strain tensor is given by the following
expression

```math
\varepsilon_{ii}=D_{i\tuple{p}}u_{i\tuple{p}}^\element
```

## Element stiffness operator

### Conductivity

### Linear elasticity

Assuming the material to be linearly elastic, with Lamé coefficients ``\lambda``
and ``\mu``, we have the following expression of the interpolated stresses as a
function of the nodal displacements ``u_{i\tuple{p}}``

```math
\sigma_{hk}=\bigl(\lambda\delta_{hk} D^\element_{i\tuple{p}}
+2\mu B^\element_{hki\tuple{p}}\bigr)u_{i\tuple{p}}
```

and the volume density of elastic energy reads

```math
\tfrac12\sigma_{hk}\varepsilon_{hk}
=\bigl(\lambda\delta_{hk} D^\element_{i\tuple{p}}
+2\mu B^\element_{hki\tuple{p}}\bigr)B_{hkj\tuple{q}}^\element
u_{i\tuple{p}}u_{j\tuple{q}}.
```

Upon expansion and integration, we get the expression of the elastic energy
within the element

```math
U^\element
=\int_{\Omega^\element}\tfrac12\sigma_{hk}\varepsilon_{hk}\,\D x_1\cdots\D x_d
=\tfrac12K^\element_{i\tuple{p}j\tuple{q}}u_{i\tuple{p}}u_{j\tuple{q}},
```

where the element stiffness operator reads

```math
K^\element_{i\tuple{p}j\tuple{q}}
=\lambda\,K^{\element, \lambda}_{i\tuple{p}j\tuple{q}}
+\mu\,K^{\element, \mu}_{i\tuple{p}j\tuple{q}},
```

with

```math
K^{\element, \lambda}_{i\tuple{p}j\tuple{q}}
=\int_{\Omega^\element}D^\element_{i\tuple{p}}D^\element_{j\tuple{q}}
\quad\text{and}\quad
K^{\element, \mu}_{i\tuple{p}j\tuple{q}}
=\int_{\Omega^\element} 2 B^\element_{hki\tuple{p}}B^\element_{hkj\tuple{q}}.
```
