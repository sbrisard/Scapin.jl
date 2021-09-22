# Continuous Green operators

In this chapter, we discuss various boundary-value problems in a periodic
setting. For each of these problems, we introduce the associated *Green
operator*.


## Conductivity

## Continuous Green operator for linear elasticity

We first define a few functional spaces; ``\tensors_2(\Omega)`` denotes the
space of second-order, symmetric, tensor fields, with square-integrable
components. Then, the space ``\stresses(\Omega)`` of periodic, self-equilibrated
stresses is defined as follows

```math
\tens\sigma\in\stresses(\Omega)\iff\left\{
\begin{gathered}
\tens\sigma\in\tensors_2(\Omega)\\
\tens\sigma\cdot\nabla=\vec 0\text{ a.e in }\Omega\\
\tens\sigma\cdot\vec e_i\text{ is }L_i\vec e_i
\text{-periodic for all }i=1, 2, \ldots, d\text{ (no summation),}
\end{gathered}
\right.
```

where the last condition expresses the periodicity of tractions in all
directions parallel to the sides of the unit-cell. The space
``\tens\strains(\Omega)`` of periodic, geometrically compatible strains is
defined as follows

```math
\tens\varepsilon\in\strains(\Omega)\iff\left\{
\begin{gathered}
\tens\varepsilon\in\tensors_2(\Omega)\\
\tens\varepsilon=\vec u\symotimes\nabla
\text{ a.e. in }\Omega\text{ for some vector field }\vec u\\
\vec u\text{ has square-integrable components}\\
\vec u\text{ is }\Omega\text{-periodic.}
\end{gathered}
\right.
```

Finally, we define the spaces of stresses and strains with zero average

```math
\stresses_0(\Omega)=\bigl\{\tens\sigma\in\stresses(\Omega),
\langle\tens\sigma\rangle=\tens0\bigr\}
\quad\text{and}\quad
\strains_0(\Omega)=\{\tens\varepsilon\in\strains(\Omega),
\langle\tens\varepsilon\rangle=\tens0\bigr\}.
```

We are now ready to define the periodic, fourth-order Green operator for strains
``\tens\Gamma``. Let ``\tens C`` be the homogeneous elastic stiffness of the
body ``\Omega``[^1]. Let ``\tens\tau\in\tensors_2(\Omega)`` be a prescribed
tensor field (*stress-polarization*). We want to find the equilibrium state of
the body ``\Omega``, subjected to the eigenstress ``\tens\tau`` and periodic
boundary conditions. In other words, we want to find the solution to the
following problem

```math
\text{Find }\tens\sigma\in\stresses_0(\Omega)
\text{ and }\tens\varepsilon\in\strains_0(\Omega)
\text{ such that }\tens\sigma=\tens C\dbldot\tens\varepsilon+\tens\tau
\text{ a.e. in }\Omega.
```

[^1]: In other words, ``\tens C`` is a constant, fourth-order tensor with major
      and minor symmetries; furthermore, ``\tens C`` is positive definite.

Owing to the periodic boundary conditions, we use [Fourier series](@ref)
expansions of ``\tens\tau``, ``\tens\sigma``, ``\tens\varepsilon`` and ``\vec
u``

```math
\begin{Bmatrix}
\tens\tau(\vec x)\\
\tens\sigma(\vec x)\\
\tens\varepsilon(\vec x)\\
\vec u(\vec x)
\end{Bmatrix}
=\sum_{\tuple{n}\in\integers^d}
\begin{Bmatrix}
\tilde{\tens\tau}_\tuple{n}\\
\tilde{\tens\sigma}_\tuple{n}\\
\tilde{\tens\varepsilon}_\tuple{n}\\
\tilde{\vec u}_\tuple{n}
\end{Bmatrix}
\exp(\I \vec k_\tuple{n}\cdot\vec x).
```

The Fourier modes ``\tilde{\tens\sigma}_\tuple{n}``,
``\tilde{\tens\varepsilon}_\tuple{n}`` and ``\tilde{\vec u}_\tuple{n}`` solve
the following equations (respectively: equilibrium, geometric compatibility,
constitutive relation)

```math
\begin{gathered}
\tilde{\tens\sigma}_\tuple{n}\cdot\vec k_\tuple{n}=\vec 0\\
\tilde{\tens\varepsilon}_\tuple{n}
=\frac{\I}{2}\bigl(\tilde{\vec u}_\tuple{n}\otimes\vec k_\tuple{n}
+\vec k_\tuple{n}\otimes\tilde{\vec u}_\tuple{n}\bigr)\\
\tilde{\tens\sigma}_\tuple{n}=\tens C\dbldot\tilde{\tens\varepsilon}_\tuple{n}
+\tilde{\tens\tau}_\tuple{n}.
\end{gathered}
```

Plugging the third equation into the second equation, and recalling that ``\tens
C`` has the minor symmetries, we find the following expression of
``\tilde{\tens\sigma}``

```math
\tilde{\tens\sigma}_\tuple{n}
=\I\bigl(\tens C\cdot\vec k_\tuple{n}\bigr)\cdot\tilde{\vec u}_\tuple{n}
+\tilde{\tens\tau}_\tuple{n}.
```

The Cauchy stress tensor being symmetric, the first of the above equations also
reads ``\vec k_\tuple{n}\cdot\tilde{\tens{\sigma}}_\tuple{n}=\vec 0`` and

```math
\tilde{\vec u}_\tuple{n}
=\I\bigl(\vec k_\tuple{n}\cdot\tens C\cdot\vec k_\tuple{n}\bigr)^{-1}
\cdot\tilde{\tens\tau}_\tuple{n}\cdot\vec k_\tuple{n}
```

which delivers the following expression for the Fourier modes of the strain
field

```math
\tilde{\tens\varepsilon}_\tuple{n}
=-\tfrac12\bigl[\bigl(\vec k_\tuple{n}\cdot\tens C\cdot\vec k_\tuple{n}\bigr)^{-1}
\cdot\tilde{\tens\tau}_\tuple{n}\cdot\vec k_\tuple{n}\bigr]\otimes\vec k_\tuple{n}
-\tfrac12\vec k_\tuple{n}
\otimes\bigl[\bigl(\vec k_\tuple{n}\cdot\tens C\cdot\vec k_\tuple{n}\bigr)^{-1}
\cdot\tilde{\tens\tau}_\tuple{n}\cdot\vec k_\tuple{n}\bigr].
```

The above relation defines a linear mapping between
``\tilde{\tens\tau}_\tuple{n}`` and ``\tilde{\tens\varepsilon}_\tuple{n}``. For
each Fourier mode ``\tuple{n}``, we therefore introduce the fourth-order tensor
``\tilde{\tens\Gamma}_\tuple{n}`` with major and minor symmetries, such that
``\tilde{\tens\varepsilon}_\tuple{n}=-\tilde{\tens\Gamma}_\tuple{n}\dbldot{\tilde{\tens\tau}}_\tuple{n}``. Our
analysis shows that ``\tilde{\tens\Gamma}_\tuple{n}=\hat{\tens\Gamma}(\vec
k_\tuple{n})`` where, for arbitrary wave-vector ``\vec k``,
``\hat{\tens\Gamma}(\vec k)`` is a fourth-order tensor with major and minor
symmetries, such that

```math
\hat{\tens\Gamma}(\vec k)\dbldot\tilde{\tens\tau}
=\tfrac12\bigl[\bigl(\vec n\cdot\tens C\cdot\vec n\bigr)^{-1}
\cdot\tilde{\tens\tau}\cdot\vec n\bigr]\otimes\vec n
+\tfrac12\vec n\otimes\bigl[\bigl(\vec n\cdot\tens C\cdot\vec n\bigr)^{-1}
\cdot\tilde{\tens\tau}\cdot\vec n\bigr],
```

where ``\vec n=\vec k/\lVert\vec k\rVert``. The above equation defines
``\hat{\tens\Gamma}(\vec k)`` by how it operates on second-order, symmetric
tensors. A closed-form expression of this tensor can be derived in the case of
an isotropic material, for which

```math
\tens C=\lambda\tens I_2\otimes\tens I_2+2\mu\tens I_4,
```

where ``\tens I_2`` (resp. ``\tens I_4``) is the second-order
(resp. fourth-order) identity tensor, and ``\lambda``, ``\mu`` are the Lamé
coefficients. Then

```math
\vec n\cdot\bigl(\tens I_2\otimes\tens I_2\bigr)\cdot\vec n=\vec n\otimes\vec n
```

then (recalling that ``\lVert\vec n\rVert=1``)

```math
\begin{aligned}
\vec n\cdot\tens I_4\cdot\vec n
&=\tfrac12 n_i\bigl(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}\bigr)n_l
\vec e_j\otimes\vec e_k
=\tfrac12\bigl(n_kn_j+n_in_i\delta_{jk}\bigr)\vec e_j\otimes\vec e_k\\
&=\tfrac12\bigl[\vec n\otimes\vec n+\bigl(\vec n\cdot\vec n\bigr)\tens I_2\bigr]
=\tfrac12\bigl(\vec n\otimes\vec n+\tens I_2\bigr)
=\vec n\otimes\vec n+\tfrac12\bigl(\tens I_2-\vec n\otimes\vec n\bigr)
\end{aligned}
```

and finally, we find the following expression of the “matrix-vector” product

```math
\begin{equation}
\label{_20210803063111}
\vec n\cdot\tens C\cdot\vec n
=\bigl(\lambda+2\mu\bigr)\vec n\otimes\vec n
+\mu\bigl(\tens I_2-\vec n\otimes\vec n\bigr)
=2\mu\frac{1-\nu}{1-2\nu}\vec n\otimes\vec n
+\mu\bigl(\tens I_2-\vec n\otimes\vec n\bigr),
\end{equation}
```

where ``\nu`` denotes the Poisson ratio. The above second-order tensor is easily
inverted, since ``\vec n\otimes\vec n`` and ``\tens I_2-\vec n\otimes\vec n``
are two orthogonal projectors (in the sense of the “``\dbldot``” product)

```math
2\mu\bigl(\vec n\cdot\tens C\cdot\vec n\bigr)^{-1}
=\frac{1-2\nu}{1-\nu}\vec n\otimes\vec n+2\bigl(\tens I_2-\vec n\otimes\vec n\bigr)
=2\tens I_2-\frac{\vec n\otimes\vec n}{1-\nu},
```

from which it results that

```math
2\mu\bigl(\vec n\cdot\tens C\cdot\vec n\bigr)^{-1}
\cdot\tilde{\tens\tau}\cdot\vec n
=2\tilde{\tens\tau}\cdot\vec n
-\frac{\vec n\cdot\tilde{\tens\tau}\cdot\vec n}{1-\nu}\vec n
```

and we finally get

```math
2\mu\hat{\tens \Gamma}(\vec k)\dbldot\tilde{\tens \tau}
=\bigl(\tilde{\tens \tau}\cdot\vec n\bigr)\otimes\vec n
+\vec n\otimes\bigl(\tilde{\tens \tau}\cdot\vec n\bigr)
-\frac{\vec n\cdot\tilde{\tens \tau}\cdot\vec n}{1-\nu}\vec n\otimes\vec n.
```

The components of the ``\hat{\tens \Gamma}`` tensor are then readily found

```math
\hat{\Gamma}_{ijkl}(\vec k)
=\frac{\delta_{ik}n_jn_l+\delta_{il}n_jn_k+\delta_{jk}n_in_l+\delta_{jl}n_in_k}{4\mu}
-\frac{n_in_jn_kn_l}{2\mu\bigl(1-\nu\bigr)},
```

which coincide with classical expressions (see e.g. [Suquet, 1990](@ref
suqu1990)). Implementation of the above equation is cumbersome; it is only used
for testing purposes. In Scapin, only the `matvec` product is required, and
Eq. \eqref{_20210803063111} was implemented.

## Hyperelasticity
