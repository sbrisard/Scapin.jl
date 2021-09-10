# On the d-dimensional brick element

In this chapter, we formulate the brick element for conductivity and linear
elasticity in ``d`` dimensions (``d\in\{2, 3\}``).

!!! note

    In the remainder of this chapter, ``i``, ``j``, ``h``, ``k`` and ``l`` denote
	scalar indices that refer to tensor components, and span
	``\{1, 2, \cdots d\}``; ``\tuple{m}``, ``\tuple{n}`` are ``d``-dimensional
	multi-indices that refer to nodes. Depending on the context (local, element
	operators or global, assembled operators) they span ``\{1, 2\}^d`` or the
	whole grid ``\cellindices``. We adopt Einstein's convention for indices of
	both types.

## [Geometry of the reference brick element](@id 20210910120306)

The reference element being centered at the origin, it occupies the following
domain

```math
\Omega^\element=\biggl(-\frac{h_1}2, \frac{h_1}2\biggr)\times\cdots
\times\biggl(-\frac{h_d}2, \frac{h_d}2\biggr),
```

where ``h_i`` is the size of the element along axis ``i`` (``i=1, \dots,
d``). The nodes of the element are indexed by ``\tuple{n}\in\{1, 2\}^d``, such
that the coordinates of node ``\tuple{n}`` read

```math
x_{\tuple{n}i}=(-1)^{n_i}\frac{h_i}2.
```

For integration purposes, it is useful to introduce the *reduced coordinates*
``\xi_i=2x_i/h_i``

```math
\int_{\Omega^\element}f(x_1, \ldots, x_d)\,\D x_1\cdots\D x_d
=\frac{h_1\cdots h_d}{2^d}\int_{(-1, 1)^d}f(h_1\xi_1/2, \ldots, h_d\xi_d/2)
\,\D\xi_1\cdots\D\xi_d
```

and it is observed that the reduced coordinates of the node ``\tuple{n}`` are
``\xi_i=(-1)^{n_i}``.

## [Shape functions](@id sec:20210910114136)

The shape function associated with node ``\tuple{n}`` is
``N^\element_{\tuple{n}}(\vec x)``, wich is defined below as a function of the
reduced coordinates of the observation point

```math
N^\element_{\tuple{n}}(\vec x)=\frac1{2^d}\prod_{i=1}^d\bigl[1+(-1)^{n_i}\xi_i\bigr]
```

and we have as expected ``N^\element_{\tuple{m}}(\vec x_{\tuple n}) =
\delta_{\tuple{m}\tuple{n}}``. The nodal values of ``f`` being ``f_\tuple{n}``,
we have the interpolation formula: ``f(\vec x)=N^\element_\tuple{n}(\vec
x)f_\tuple{n}``.

## [Gradient and strain-displacement operators](@id 20210910114926)

We consider a scalar interpolated field, ``f(\vec x)=N^\element_\tuple{n}(\vec
x)f_\tuple{n}``. The components of its gradient are given by the following
expression

```math
\bigl(f\otimes\nabla\bigr)=\partial_i f=D_{i\tuple{n}}^\element f_\tuple{n},
```

where the observation point ``\vec x`` has been dropped, and
``D_{i\tuple{n}}^\element`` denotes the element gradient operator

```math
D_{i\tuple{n}}^\element=\partial_iN_\tuple{n}
=\frac{(-1)^{n_i}}{2^{d-1}h_i}
\prod_{j\neq i}\bigl[1+(-1)^{n_j}\xi_j\bigr].
```

It will be useful to introduce the element average of this operator,
``\overline{D}_{i\tuple{n}}^\element``

```math
\overline{D}_{i\tuple{n}}^\element=\frac 1{h_1\cdots h_d}
\int_{\Omega^\element}D^\element_{i\tuple{n}}(\vec x)\,\D x_1\cdots\D x_d
```

and, observing that each factor in square brackets in the above expression of
``D^\element_{i\tuple{n}}`` varies linearly between 0 et 2 (and therefore
averages to 1), we find

```math
\overline{D}^\element_{i\tuple{n}}=\frac{(-1)^{n_i}}{2^{d-1}h_i}.
```

We now consider the interpolated displacement field, ``\vec u =
N^\element_\tuple{n} \vec u_\tuple{n}``. The components of the associated
strain, ``\tens\varepsilon=\sym\bigl(\vec u\otimes\nabla\bigr)`` are

```math
\varepsilon_{ij}=\tfrac12\bigl(\partial_iu_j+\partial_ju_i\bigr)
=\tfrac12\bigl(D^\element_{i\tuple{n}}u_{\tuple{n}j}
+D^\element_{j\tuple{n}}u_{\tuple{n}i}\bigr),
```

where ``u_{\tuple{n}i}`` denotes the ``i``-th component of the nodal value of
``\vec u`` at node ``\tuple{n}``, ``u_{\tuple{n}i}=\vec u(\vec
x_\tuple{n})\cdot\vec e_i``. Rearranging the above expression, we find
``\varepsilon_{ij}=B^\element_{ij\tuple{n}k}u_{\tuple{n}k}^\element``, where

```math
B_{ij\tuple{n}k}^\element
=\tfrac12\bigl(D^\element_{i\tuple{n}}\delta_{jk}
+D^\element_{j\tuple{n}}\delta_{ik}\bigr)
```

is the element strain-displacement operator, with volume average

```math
\overline{B}_{ij\tuple{n}k}^\element
=\tfrac12\bigl(\overline{D}^\element_{i\tuple{n}}\delta_{jk}
+\overline{D}^\element_{j\tuple{n}}\delta_{ik}\bigr)
=\frac1{2^d}\Bigl[\frac{(-1)^{n_i}\delta_{jk}}{h_i}
+\frac{(-1)^{n_j}\delta_{ik}}{h_j}\Bigr].
```

It is obverved that the *trace* of the strain tensor is given by the following
expression

```math
\varepsilon_{ii}=D_{i\tuple{n}}u_{\tuple{n}i}^\element
```

## Element stiffness operator

### Conductivity

### Linear elasticity

Assuming the material to be linearly elastic, with Lamé coefficients ``\lambda``
and ``\mu``, we have the following expression of the interpolated stresses as a
function of the nodal displacements ``u_{\tuple{m}i}``

```math
\sigma_{hk}=\bigl(\lambda D^\element_{i\tuple{m}}\delta_{hk}
+2\mu B^\element_{hk\tuple{m}i}\bigr)u_{\tuple{m}i}
```

and the volume density of elastic energy reads

```math
\tfrac12\sigma_{hk}\varepsilon_{hk}
=\bigl(\lambda D^\element_{i\tuple{m}}\delta_{hk}
+2\mu B^\element_{hk\tuple{m}i}\bigr)B_{hk\tuple{n}j}^\element
u_{\tuple{m}i}u_{\tuple{n}j}.
```

Upon expansion and integration, we get the expression of the elastic energy
within the element

```math
U^\element
=\int_{\Omega^\element}\tfrac12\sigma_{hk}\varepsilon_{hk}\,\D x_1\cdots\D x_d
=\tfrac12K^\element_{\tuple{m}i\tuple{n}j}u_{\tuple{m}i}u_{\tuple{n}j},
```

where the element stiffness operator reads

```math
K^\element_{\tuple{m}i\tuple{n}j}
=\lambda\,K^{\element, \lambda}_{\tuple{m}i\tuple{n}j}
+\mu\,K^{\element, \mu}_{\tuple{m}i\tuple{n}j},
```

with

```math
K^{\element, \lambda}_{\tuple{m}i\tuple{n}j}
=\int_{\Omega^\element}D^\element_{i\tuple{m}}D^\element_{j\tuple{n}}
\quad\text{and}\quad
K^{\element, \mu}_{\tuple{m}i\tuple{n}j}
=\int_{\Omega^\element} 2 B^\element_{hk\tuple{m}i}B^\element_{hk\tuple{n}j}.
```
