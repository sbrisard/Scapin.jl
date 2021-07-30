# Continuous Green operators

In this chapter, we discuss various boundary-value problems in a periodic
setting. For each of these problems, we introduce the associated *Green
operator*.

## On Fourier series

Owing to the periodic setting, the fields that are involved in the various BVPs
to be discussed in this chapter are expanded in Fourier series. ``\tens T``
being a ``Ω``-periodic tensor field

```math
\begin{equation}
  \tens T(\vec x)=\sum_{\tuple{n}∈\integers^d}\mathcal F(\tens T)(\vec
  k_{\tuple{n}})\exp(\I\vec k_{\tuple{n}}\cdot\vec x),
\end{equation}
```

where ``\tuple{n}`` denotes a ``d``-dimensional tuple of integers (see
[Nomenclature](@ref)). The wave vectors ``\vec k_{\tuple{n}}`` are given by

```math
\begin{equation}
  \vec k_{\tuple{n}}
  =\frac{2\pi n_1}{L_1}\vec e_1
  +\frac{2\pi n_2}{L_2}\vec e_2
  +\cdots
  +\frac{2\pi n_d}{L_d}\vec e_d,
\end{equation}
```

and the Fourier coefficients of ``\tens T`` are defined as follows

```math
\begin{equation}
  \mathcal F(\tens T)(\vec k)
  =\frac1V\int_{\vec x∈Ω}\tens T(\vec x)\exp(-\I\vec k\cdot\vec x)\,\D x_1\cdots\D x_d.
\end{equation}
```

It is recalled that the Fourier coefficients of the gradient and divergence of
``\tens T`` can readily be computed from the Fourier coefficients of ``\tens T``

```math
\begin{equation}
  \mathcal F(\tens T\otimes\nabla)(\vec k)=\mathcal F(\tens T)(\vec k)\otimes\I\vec k
	\quad\text{and}\quad
  \mathcal F(\tens T\cdot\nabla)(\vec k)=\mathcal F(\tens T)(\vec k)\cdot\I\vec k.
\end{equation}
```

When no confusion is possible, we will use the tilde to denote the Fourier
coefficients: ``\tilde{\tens T}_\tuple{n}=\mathcal F(\tens T)(\vec k_\tuple{n})``.

## Conductivity

## Elasticity

We first define a few functional spaces; ``\tensors_2(Ω)`` denotes the space of
second-order, symmetric, tensor fields, with square-integrable components. Then,
the space ``\tens\stresses(Ω)`` of periodic, self-equilibrated stresses is
defined as follows

```math
\begin{equation}
  \tens\sigma\in\stresses(Ω)\iff\left\{
  \begin{gathered}
    \tens\sigma\in\tensors_2(Ω)\\
    \tens\sigma\cdot\nabla=\vec 0\text{ a.e in }Ω\\
    \tens\sigma\cdot\vec e_i\text{ is }L_i\vec e_i
	\text{-periodic for all }i=1, 2, \ldots, d\text{ (no summation),}
  \end{gathered}
  \right.
\end{equation}
```

where the last condition expresses the periodicity of tractions in all
directions parallel to the sides of the unit-cell. The space
``\tens\strains(Ω)`` of periodic, geometrically compatible strains is defined as
follows

```math
\begin{equation}
  \tens\varepsilon\in\strains(Ω)\iff\left\{
  \begin{gathered}
    \tens\varepsilon\in\tensors_2(Ω)\\
    \tens\varepsilon=\sym\bigl(\vec u\otimes\nabla\bigr)
	\text{ a.e. in }Ω\text{ for some vector field }\vec u\\
    \vec u\text{ has square-integrable components}\\
    \vec u\text{ is }Ω\text{-periodic.}
  \end{gathered}
  \right.
\end{equation}
```

Finally, we define the spaces of stresses and strains with zero average

```math
\begin{equation}
  \stresses_0(Ω)=\bigl\{\tens\sigma\in\stresses(Ω),
  \langle\tens\sigma\rangle=\tens0\bigr\}
  \quad\text{and}\quad
  \strains_0(Ω)=\{\tens\varepsilon\in\strains(Ω),
  \langle\tens\varepsilon\rangle=\tens0\bigr\}.
\end{equation}
```

We are now ready to define the periodic, fourth-order Green operator for strains
``\tensΓ``. Let ``\tens C`` be the homogeneous elastic stiffness of the body
``Ω``[^1]. Let ``\tens\tau\in\tensors_2(Ω)`` be a prescribed tensor field
(*stress-polarization*). We want to find the equilibrium state of the body
``Ω``, subjected to the eigenstress ``\tens\tau`` and periodic boundary
conditions. In other words, we want to find the solution to the following
problem

```math
\begin{equation}
\text{Find }\tens\sigma\in\stresses_0(Ω)
\text{ and }\tens\varepsilon\in\strains_0(Ω)
\text{ such that }\tens\sigma=\tens C\dbldot\tens\varepsilon+\tens\tau
\text{ a.e. in }Ω.
\end{equation}
```

[^1]: In other words, ``\tens C`` is a constant, fourth-order tensor with major
      and minor symmetries; furthermore, ``\tens C`` is positive definite.

Owing to the periodic boundary conditions, we use Fourier series expansions of
``\tens\tau``, ``\tens\sigma``, ``\tens\varepsilon`` and ``\vec u``

```math
\begin{equation}
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
\end{equation}
```

The Fourier modes ``\tilde{\tens\sigma}_\tuple{n}``, ``\tilde{\tens\varepsilon}_\tuple{n}`` and
``\tilde{\vec u}_\tuple{n}`` solve the following equations (respectively: equilibrium,
geometric compatibility, constitutive relation)

```math
\begin{gather}
\label{eq:20210730094655}
\tilde{\tens\sigma}_\tuple{n}\cdot\vec k_\tuple{n}=\vec 0\\
\label{eq:20210730094504}
\tilde{\tens\varepsilon}_\tuple{n}=\frac{\I}{2}\bigl(\tilde{\vec u}_\tuple{n}\otimes\vec k_\tuple{n}
+\vec k_\tuple{n}\otimes\tilde{\vec u}_\tuple{n}\bigr)\\
\label{eq:20210730094514}
\tilde{\tens\sigma}_\tuple{n}=\tens C\dbldot\tilde{\tens\varepsilon}_\tuple{n}
+\tilde{\tens\tau}_\tuple{n}.
\end{gather}
```

Plugging Eq. \eqref{eq:20210730094514} into Eq. \eqref{eq:20210730094504}, and
recalling that ``\tens C`` has the minor symmetries, we find the following
expression of ``\tilde{\tens\sigma}``

```math
\begin{equation}
  \tilde{\tens\sigma}_\tuple{n}
  =\I\bigl(\tens C\cdot\vec k_\tuple{n}\bigr)\cdot\tilde{\vec u}_\tuple{n}+\tilde{\tens\tau}_\tuple{n}.
\end{equation}
```

The Cauchy stress tensor being symmetric, Eq. \eqref{eq:20210730094655} also
reads ``\vec k_\tuple{n}\cdot\tilde{\tens{\sigma}}_\tuple{n}=\vec 0`` and

```math
\begin{equation}
  \label{eq:16}
  \tilde{\vec u}_\tuple{n}=\I\bigl(\vec k_\tuple{n}\cdot\tens C\cdot\vec k_\tuple{n}\bigr)^{-1}
  \cdot\tilde{\tens\tau}_\tuple{n}\cdot\vec k_\tuple{n}
\end{equation}
```

which delivers the following expression for the Fourier modes of the strain
field

```math
\begin{equation}
\label{eq:20210730094915}
\tilde{\tens\varepsilon}_\tuple{n}=-\tfrac12\bigl[\bigl(\vec k_\tuple{n}\cdot\tens C\cdot
\vec k_\tuple{n}\bigr)^{-1}\cdot\tilde{\tens\tau}_\tuple{n}\cdot\vec k_\tuple{n}\bigr]\otimes\vec k_\tuple{n}
-\tfrac12\vec k_\tuple{n}\otimes\bigl[\bigl(\vec k_\tuple{n}\cdot\tens C\cdot\vec k_\tuple{n}
\bigr)^{-1}\cdot\tilde{\tens\tau}_\tuple{n}\cdot\vec k_\tuple{n}\bigr].
\end{equation}
```

The above relation defines a linear mapping between
``\tilde{\tens\tau}_\tuple{n}`` and ``\tilde{\tens\varepsilon}_\tuple{n}``. For
each Fourier mode ``\tuple{n}``, we therefore introduce the fourth-order tensor
``\tilde{\tens Γ}_\tuple{n}`` with major and minor symmetries, such that
``\tilde{\tens\varepsilon}_\tuple{n}=-\tilde{\tensΓ}_\tuple{n}\dbldot{\tilde{\tens\tau}}_\tuple{n}``. From
Eq. \eqref{eq:20210730094915}, it results that
``\tilde{\tensΓ}_\tuple{n}=\hat{\tensΓ(\vec k_\tuple{n})}`` where, for arbitrary
wave-vector ``\vec k``, ``\hat{\tensΓ}(\vec k)`` is a fourth-order tensor with
major and minor symmetries, such that

```math
\begin{equation}
  \label{eq:20210730095035}
  \hat{\tensΓ}(\vec k)\dbldot\tilde{\tens\tau}=\tfrac12\bigl[\bigl(\vec n
  \cdot\tens C\cdot\vec n\bigr)^{-1}\cdot\tilde{\tens\tau}\cdot\vec n\bigr]
  \otimes\vec n+\tfrac12\vec n\otimes\bigl[\bigl(\vec n\cdot\tens C\cdot\vec n
  \bigr)^{-1}\cdot\tilde{\tens\tau}\cdot\vec n\bigr],
\end{equation}
```

where ``\vec n=\vec k/\lVert\vec k\rVert``. Eq. \eqref{eq:20210730095035}
defines ``\hat{\tensΓ}(\vec k)`` by how it operates on second-order, symmetric
tensors. A closed-form expression of this tensor can be derived in the case of
an isotropic material, for which

```math
\begin{equation}
  \tens C=\lambda\tens I_2\otimes\tens I_2+2\mu\tens I_4,
\end{equation}
```

where ``\tens I_2`` (resp. ``\tens I_4``) is the second-order
(resp. fourth-order) identity tensor, and ``\lambda``, ``\mu`` are the Lamé
coefficients. Then

```math
\begin{equation}
  \vec n\cdot\bigl(\tens I_2\otimes\tens I_2\bigr)\cdot\vec n=\vec n\otimes\vec n
\end{equation}
```

then (recalling that ``\lVert\vec n\rVert=1``)

```math
\begin{equation}
  \begin{aligned}[b]
    \vec n\cdot\tens I_4\cdot\vec n&=\tfrac12 n_i\bigl(\delta_{ik}\delta_{jl}+
    \delta_{il}\delta_{jk}\bigr)n_l\vec e_j\otimes\vec e_k=\tfrac12\bigl(n_kn_j
    +n_in_i\delta_{jk}\bigr)\vec e_j\otimes\vec e_k\\
    &=\tfrac12\bigl[\vec n\otimes\vec n+\bigl(\vec n\cdot\vec n\bigr)\tens I_2
    \bigr]=\tfrac12\bigl(\vec n\otimes\vec n+\tens I_2\bigr)
    =\vec n\otimes\vec n+\tfrac12\bigl(\tens I_2-\vec n\otimes\vec n\bigr)
  \end{aligned}
\end{equation}
```

and finally

```math
\begin{equation}
  \vec n\cdot\tens C\cdot\vec n=\bigl(\lambda+2\mu\bigr)\vec n\otimes\vec n+\mu
  \bigl(\tens I_2-\vec n\otimes\vec n\bigr)=2\mu\frac{1-\nu}{1-2\nu}\vec n
  \otimes\vec n+\mu\bigl(\tens I_2-\vec n\otimes\vec n\bigr),
\end{equation}
```

where ``\nu`` denotes the Poisson ratio. The above second-order tensor is easily
inverted, since ``\vec n\otimes\vec n`` and ``\tens I_2-\vec n\otimes\vec n``
are two orthogonal projectors (in the sense of the “``\dbldot``” product)

```math
\begin{equation}
  2\mu\bigl(\vec n\cdot\tens C\cdot\vec n\bigr)^{-1}=\frac{1-2\nu}{1-\nu}\vec n
  \otimes\vec n+2\bigl(\tens I_2-\vec n\otimes\vec n\bigr)=2\tens I_2
  -\frac{\vec n\otimes\vec n}{1-\nu},
\end{equation}
```

from which it results that

```math
\begin{equation}
  2\mu\bigl(\vec n\cdot\tens C\cdot\vec n\bigr)^{-1}\cdot\tilde{\tens\tau}\cdot
  \vec n=2\tilde{\tens\tau}\cdot\vec n-\frac{\vec n\cdot\tilde{\tens\tau}\cdot
    \vec n}{1-\nu}\vec n
\end{equation}
```

and finally

```math
\begin{equation}
\label{eq:20210730095803}
2\mu\hat{\tensΓ}(\vec k)\dbldot\tilde{\tens\tau}=\bigl(\tilde{\tens\tau}
\cdot\vec n\bigr)\otimes\vec n+\vec n\otimes\bigl(\tilde{\tens\tau}\cdot\vec n
\bigr)-\frac{\vec n\cdot\tilde{\tens\tau}\cdot\vec n}{1-\nu}\vec n\otimes\vec n.
\end{equation}
```

The components of the ``\hat{\tensΓ}`` tensor are then readily found

```math
\begin{equation}
  \hat{Γ}_{ijkl}(\vec k)=\frac{\delta_{ik}n_jn_l+\delta_{il}n_jn_k
  +\delta_{jk}n_in_l+\delta_{jl}n_in_k}{4\mu}-\frac{n_in_jn_kn_l}{2\mu\bigl(1-
    \nu\bigr)},
\end{equation}
```

which coincide with classical expressions \parencite[see
e.g.][]{suqu1990}. Implementation of the above equation is cumbersome; it is
only used for testing purposes. In Scapin, only the `matvec` product is
required: Eq. \eqref{eq:20210730095803} was implemented.

## Hyperelasticity
