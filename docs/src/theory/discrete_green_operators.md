# Discrete Green operators

In this chapter, we introduce various discretizations of the Green operator; we
will adopt the vocabulary of linear elasticity, although the concepts apply to
all the various physical models presented in Chap. [Continuous Green
operators](@ref).


### The ``\fftfreq`` function

For ``n, N\in\naturals``, ``0\leq n<N``, we introduce ``\fftfreq(n, N)``

```math
\fftfreq(n, N)=
\begin{cases}
n & \text{if }2n<N,\\
n-N & \text{otherwise.}
\end{cases}
```

For ``n<0`` or ``n\geq N``, ``\fftfreq(n, N)`` is defined by
``N``-periodicity. ``\fftfreq`` is very similar to the NumPy
[fftfreq](https://numpy.org/doc/1.18/reference/generated/numpy.fft.fftfreq.html#numpy.fft.fftfreq)
function. We have the important result (see proof in Sec. [Properties of the
``\fftfreq`` function](@ref)).

```math
\fftfreq(N-n, N)=
\begin{cases}
\fftfreq(n) & \text{if }2n=N,\\
-\fftfreq(n) & \text{otherwise.}
\end{cases}
```

The ``\fftfreq`` function can be defined for ``d``-tuples as well

```math
\tuple \fftfreq(\tuple n, \tuple N)
=\bigl(\fftfreq(n_1, N_1), \ldots, \fftfreq(n_d, N_d)\bigr)
```

and we have again

```math
\tuple \fftfreq(\tuple N-\tuple n, \tuple N)=-\tuple \fftfreq(\tuple n)
```

if none of the ``n_i`` is such that ``2n_i=N_i``.

## The approximation space

In order to define a discrete Green operator, we need to introduce the
approximation space for the stress-polarizations. We will consider here
stress-polarizations that are constant over each cell of a regular grid of size
``N_1\times\cdots\times N_d``. The cells of this grid are

```math
\Omega_{\tuple{p}}^{\tuple{h}}=\{p_ih_i\leq x_i<\bigl(p_i+1\bigr)h_i, i=1,\ldots, d\},
```

where ``\tuple{p}=(p_1,\ldots,p_d)\in\cellindices`` denotes a ``d``-tuple of
integers, such that ``0\leq p_i<N_i``, ``i=1,\ldots, d`` and ``h_i=L_i/N_i`` is
the cell-size in the ``i``-th direction. Note that, in the above expression, no
summation is implied by the repeated index ``i``. The total number of cells is
``\lvert N\rvert=N_1\cdots N_d``.

We consider discrete stress-polarizations ``\tens\tau^{\tuple{h}}`` that are
constant over each cell of the grid: ``\tens\tau_{\tuple{p}}^{\tuple{h}}``
denotes the constant value of ``\tens\tau^{\tuple{h}}`` in cell
``\Omega_{\tuple{p}}^{\tuple{h}}``. The ``\tuple{h}`` superscript reminds that
``\tens\tau^{\tuple{h}}`` is a discrete approximation of the true
stress-polarization ``\tens\tau``. We will call this approximation subspace:
``\tensors_2^{\tuple{h}}(\Omega)``.

As discussed **TODO: xref**, the discrete Green operator is defined as the
restriction to this approximation space of the continuous Green operator, seen
as a bilinear form, or an approximation of it. In other words, we want to
propose an approximation of the quantity

```math
\langle\tens\varpi^{\tuple{h}}\dbldot\tens\Gamma(\tens\tau^{\tuple{h}})\rangle
\simeq\langle\tens\varpi^{\tuple{h}}\dbldot\tens\Gamma^{\tuple{h}}(\tens
\tau^{\tuple{h}})\rangle\quad\text{for all }\tens\tau^{\tuple{h}},
\tens\varpi^{\tuple{h}}\in\tensors_2^{\tuple{h}}(\Omega),
```

where ``\tens\Gamma^{\tuple{h}}`` is defined only over
``\tensors_2(\Omega)``. ``\tens\Gamma^{\tuple{h}}`` can therefore be seen as a
linear mapping between the cell values ``\tens\tau_{\tuple{p}}^{\tuple{h}}`` of
``\tens\tau^{\tuple{p}}`` and the cell values of
``\tens\Gamma^{\tuple{h}}(\tens\tau^{\tuple{h}})``; ``\tens\Gamma^{\tuple{h}}``
is therefore a *matrix*, and the above equation should be understood as

```math
\begin{equation}
\label{eq20210914103849}
\langle\tens\varpi^{\tuple{h}}\dbldot\tens\Gamma(\tens\tau^{\tuple{h}})\rangle
\simeq\frac1{\lvert\tuple{N}\rvert}\sum_{\tuple{p}, \tuple{q}\in\cellindices}
\tens\varpi_{\tuple{p}}^{\tuple{h}}\dbldot
\tens\Gamma_{\tuple{p}\tuple{q}}^{\tuple{h}}\dbldot\tens
\tau_{\tuple{q}}^{\tuple{h}}
\quad\text{for all}\quad
\tens\tau^{\tuple{h}},\tens\varpi^{\tuple{h}}\in\tensors_2^{\tuple{h}}(\Omega).
\end{equation}
```

The continuous Green operator is translation invariant, and this property will
of course be transferred to the “exact” discrete Green operator; we will in fact
require *all* dicretizations of the Green operator to have this property. In
other words, ``\tens\Gamma_{\tuple{p}\tuple{q}}^{\tuple{h}}=\tens\Gamma_{\tuple
p-\tuple q}^{\tuple{h}}`` and the above equation reads

```math
\langle\tens\varpi^{\tuple{h}}\dbldot\tens\Gamma(\tens\tau^{\tuple{h}})\rangle
\simeq\frac1{\lvert\tuple{N}\rvert}\sum_{\tuple{p}, \tuple{q}\in\cellindices}
\tens\varpi_{\tuple{p}}^{\tuple{h}}\dbldot
\tens\Gamma_{\tuple p-\tuple q}^{\tuple{h}}\dbldot\tens
\tau_{\tuple{q}}^{\tuple{h}}\quad\text{for all }\tens\tau^{\tuple{h}},\tens
\varpi^{\tuple{h}}\in\tensors_2^{\tuple{h}}(\Omega).
```

Note that ``\tens\Gamma_{\tuple{p}}^{\tuple{h}}`` is now indexed by *one* index
only and its DFT can be introduced unambiguously

```math
\begin{aligned}[b]
&\frac{1}{\lvert\tuple N\rvert}\sum_{\tuple{p}, \tuple{q}\in\cellindices}
\tens\varpi_{\tuple{p}}^{\tuple{h}}\dbldot\tens\Gamma_{\tuple p-\tuple q}^{\tuple h}
\dbldot\tens\tau_{\tuple{q}}^{\tuple{h}}\\
={}&\frac1{\lvert\tuple N\rvert^2}\sum_{\tuple p, \tuple q, \tuple n\in\cellindices}
\exp\Bigl[2\I\PI\sum_{j=1}^d\frac{n_j}{N_j}\bigl(p_j-q_j\bigr)\Bigr]
\tens\varpi_{\tuple{p}}^{\tuple{h}}\dbldot\hat{\tens\Gamma}_{\tuple n}^{\tuple h}
\dbldot\tens\tau_{\tuple q}^{\tuple h}\\
={}&\frac{1}{\lvert\tuple N\rvert^2}
\sum_{\tuple n\in\cellindices}\Bigl[\sum_{\tuple{p}\in\cellindices}
\exp\Bigl(2\I\PI\sum_{j=1}^d\frac{n_jp_j}{N_j}\Bigr)
\tens\varpi_{\tuple p}^{\tuple h}\Bigr]\dbldot\hat{\tens\Gamma}_{\tuple n}^{\tuple h}
\dbldot\Bigl[\sum_{\tuple q\in\cellindices}\exp\Bigl(-2\I\PI\sum_{j=1}^d
\frac{n_jq_j}{N_j}\tens\tau_{\tuple q}^{\tuple{h}}\Bigr)\Bigr].
\end{aligned}
```

Since ``\tens\varpi`` is real, we have
``\conj(\tens\varpi_{\tuple p}^h)=\tens\varpi_{\tuple p}^h`` and the first sum in
square brackets reads

```math
\sum_{\tuple{p}\in\cellindices}\exp\Bigl(2\I\PI\sum_{j=1}^d\frac{n_jp_j}{N_j}\Bigr)
\tens\varpi_{\tuple p}^{\tuple h}=\conj\Bigl[\sum_{\tuple{p}\in\cellindices}
\exp\Bigl(-2\I\PI\sum_{j=1}^d\frac{n_jp_j}{N_j}\Bigr)
\tens\varpi_{\tuple p}^{\tuple h}\Bigr]
=\conj(\hat{\tens\varpi}_{\tuple n}^{\tuple h}),
```

while the second sum in square brackets reduces to ``\hat{\tens\tau}_{\tuple
n}^{\tuple h}``. Gathering the above results, we find

```math
\frac{1}{\lvert\tuple N\rvert}\sum_{\tuple p, \tuple q\in\cellindices}
\tens\varpi_{\tuple p}^{\tuple h}\dbldot\tens\Gamma_{\tuple p-\tuple q}^{\tuple h}
\dbldot\tens\tau_{\tuple q}^{\tuple h}=\frac1{\lvert\tuple N\rvert^2}
\sum_{\tuple n\in\cellindices}\conj(\hat{\tens\varpi}_{\tuple n}^{\tuple h})
\dbldot\hat{\tens\Gamma}_{\tuple n}^{\tuple h}\dbldot
\hat{\tens\tau}_{\tuple n}^{\tuple h}.
```

The above equation can be understood as follows. ``\tens\Gamma^{\tuple h}`` is
a linear operator that maps the cell-wise constant field
``\tens\tau^{\tuple h}`` to the cell-wise constant field
``\tens\eta^{\tuple h}``, the cell-values of which are given by their DFT

```math
\tens{\eta}_{\tuple p}^{\tuple h}=
\dft^{-1}_{\tuple p}(\hat{\tens\eta}_{\bullet}^{\tuple h}),\quad\text{with}
\quad\hat{\tens\eta}_{\tuple n}^{\tuple h}
=\hat{\tens\Gamma}_{\tuple n}^{\tuple h}\dbldot
\hat{\tens\tau}_{\tuple n}^{\tuple h}.
```

Then, from the Plancherel theorem

```math
\frac1{\lvert\tuple N\rvert^2}
\sum_{\tuple n\in\cellindices}\conj(\hat{\tens\varpi}_{\tuple n}^{\tuple h})
\dbldot\hat{\tens\Gamma}_{\tuple n}^{\tuple h}\dbldot
\hat{\tens\tau}_{\tuple n}^{\tuple h}=\frac1{\lvert\tuple N\rvert^2}
\sum_{\tuple n\in\cellindices}\conj(\hat{\tens\varpi}_{\tuple n}^{\tuple h})
\dbldot\hat{\tens\eta}_{\tuple n}^{\tuple h}=\frac1{\lvert\tuple N\rvert}
\sum_{\tuple p\in\cellindices}\tens\varpi_{\tuple p}^{\tuple h}
\dbldot\tens\eta_{\tuple n}^{\tuple h}
```

and the last sum can be seen as the volume average ``\langle\tens\varpi^{\tuple
h}\dbldot\tens\eta^{\tuple h}\rangle``. Remembering that this expression was
proposed as an approximation of ``\langle\tens\varpi^{\tuple
h}\dbldot\tens\Gamma(\tens\tau^{\tuple h})\rangle``, we finally find

```math
\langle\tens\varpi^{\tuple h}\dbldot\tens\Gamma(\tens\tau^{\tuple h})\rangle
\simeq\langle\tens\varpi^{\tuple h}\dbldot\tens\eta^{\tuple h}\rangle
=\langle\tens\varpi^{\tuple h}\dbldot\tens\Gamma^{\tuple h}(\tens\tau^{\tuple h})\rangle
\quad\text{for all }\tens\varpi^{\tuple h}\in\tensors_2^{\tuple h}(\Omega),
```

from which we find

```math
\tens\Gamma(\tens\tau^{\tuple h})\simeq\tens\Gamma^{\tuple h}
(\tens\tau^{\tuple h}).
```

The discrete Green operator, which was first introduced as an approximation of
the continuous Green operator, seen as a bilinear form, can also be understood
as an approximation of the continuous Green operator, seen as a linear
mapping. This latter point of view will become extremely efficient when it comes
to discretizing the Lippmann–Schwinger equation.

It results from the above developments that an explicit expression of the
discrete Green operator as a (gigantic) matrix is never needed. Instead, the
matrix-vector product
``\tens\tau^{\tuple h}\mapsto\tens\Gamma^{\tuple h}(\tens\tau^{\tuple h})`` is
implemented in a matrix-free fashion as follows

1. Given ``\tens\tau^{\tuple h}\in\tensors_2^{\tuple h}(\Omega)``, compute the
   discrete Fourier transform ``\hat{\tens\tau}_{\tuple n}^{\tuple h}`` of its
   cell-values:``\hat{\tens\tau}_{\tuple n}^{\tuple h}=\dft_{\tuple
   n}(\tens\tau_\bullet^{\tuple h})``,
2. for each discrete frequency, compute ``\hat{\tens\eta}_{\tuple n}^{\tuple
   h}=\hat{\tens\Gamma}_{\tuple n}^{\tuple h}\dbldot\hat{\tens\tau}_{\tuple
   n}^{\tuple h}``,
3. compute the inverse discrete Fourier transform
   ``\tens\eta_{\tuple{p}}^{\tuple h}`` of ``\hat{\tens\eta}_{\tuple{n}}^{\tuple
   h}``,

discrete Fourier transforms being computed in steps 1 and 3 by means of the FFT.

### Condition for the discrete Green operator to map real fields onto real fields

In the remainder of this chapter, we propose various discretizations of the
Green operator. Before we proceed, though, it should be emphasized that the
discrete Green operator must map a *real* field onto a *real* field. In other
words, we must have ``\hat{\tens\eta}_{\tuple N-\tuple n}^{\tuple
h}=\conj(\hat{\tens\eta}_{\tuple n}^{\tuple h})`` for all ``\tuple n``. Since
``\hat{\tens\eta}_{\tuple n}^{\tuple h} =\hat{\tens\Gamma}_{\tuple n}^{\tuple
h}\dbldot\hat{\tens\tau}_{\tuple n}^{\tuple h}`` and ``\hat{\tens\tau}_{\tuple
n}`` already satisfies this condition (it is the DFT of a *real* field), we find
the following requirement.

Any discrete operator that will be considered below must ensure that

```math
\hat{\tens\Gamma}_{\tuple N-\tuple n}^{\tuple h}
=\conj(\hat{\tens\Gamma}_{\tuple n}^{\tuple h})
\quad\text{for all}\quad\tuple n\in\cellindices.
```

## Discretizations of the Green operator

### The finite element discretization

As discussed in Sec. [The approximation space](@ref) [see in particular
Eq. \eqref{eq20210914103849}], the discrete Green operator can be seen as an
operator that delivers the cell-averages of the strains induced by a cell-wise
constant eigenstress ``\tens\tau^\tuple{h}\in\tensors_2^\tuple{h}(\Omega)``. In
other words, let us consider the solution to the following problem (see
Sec. [Continuous Green operator for linear elasticity](@ref))

```math
\begin{equation}
\label{eq20210914140834}
\text{Find }\tens\sigma\in\stresses_0(\Omega)
\text{ and }\tens\varepsilon\in\strains_0(\Omega)
\text{ such that }\tens\sigma=\tens C\dbldot\tens\varepsilon+\tens\tau^\tuple{h}
\text{ a.e. in }\Omega.
\end{equation}
```

Note that, compared to the general problem that defines the continuous Green
operator (see Sec. [Continuous Green operator for linear elasticity](@ref)), the
above problem differs by the space where ``\tens\tau^\tuple{h}`` lives. An
approximate solution to this problem can be retrieved from a finite element
discretization (Brisard [2017](@ref bris2017))[^1], where the mesh coincides
with the grid introduced previously and each element is discretized with linear
[Shape functions](@ref). The resulting linear system can be solved efficiently
by means of a matrix-free approach that is outlined below.

[^1]: The preprint of this paper is freely available on the [HAL
      archive](https://hal-enpc.archives-ouvertes.fr/hal-01304603); it is not
      the final version, but the theoretical sections were unchanged through the
      review process.

#### The modal strain-displacement vector

Introducing the nodal displacement ``\vec u_\tuple{n}`` at node
``\tuple{n}\in\cellindices``, we define the inerpolated displacement ``\vec
u^\tuple{h}(\vec x)`` and strain ``\tens\strain^\tuple{h}(\vec x)``, as well as
the cell-average of the latter, ``\overline{\tens\strain}_\tuple{p}^\tuple{h}``

```math
\tens\varepsilon^\tuple{h}(\vec x)=\bigl(\vec
u^\tuple{h}\symotimes\nabla\bigr)(\vec x)
\quad\text{and}\quad
\overline{\tens\varepsilon}_\tuple{p}^\tuple{h}
=\frac1{\lvert\tuple{h}\rvert}\int_{\Omega_\tuple{p}}\tens\strain^\tuple{h}(\vec x)
\,\D x_1\cdots\D x_d,
```

where ``|\tuple{h}|=h_1\cdots h_d`` is the cell volume. In the present periodic
setting, the DFT ``\hat{\overline{\tens\varepsilon}}_\tuple{n}^\tuple{h}`` of
``\tens\varepsilon`` can be expressed as follows

```math
\begin{equation}
\label{eq20210914144114}
\hat{\overline{\tens\varepsilon}}_\tuple{n}^\tuple{h}
=\hat{\vec b}_\tuple{n}^\tuple{h}\symotimes\hat{\vec u}_\tuple{n}^\tuple{h},
\end{equation}
```

where ``\hat{\vec b}_\tuple{n}^\tuple{h}`` is the so-called *modal strain-displacement
vector* [see (Brisard [2017](@ref bris2017)) for its expression].

The modal strain-displacement vector introduced above is related to the *nodal*
strain-displacement operator introduced in Sec. [Gradient and
strain-displacement operators](@ref) of Chap. [On the d-dimensional brick
element](@ref). Indeed, in components, Eq. \eqref{eq20210914144114} reads

```math
\hat{\overline{\varepsilon}}_{ij\tuple{n}}^\tuple{h}
=\tfrac12\bigl(\hat{b}_{i\tuple{n}}^\tuple{h}\hat{u}_{j\tuple{n}}^\tuple{h}
+\hat{b}_{j\tuple{n}}^\tuple{h}\hat{u}_{i\tuple{n}}^\tuple{h}\bigr),
```

and, from the circular convolution theorem

```math
\overline{\varepsilon}_{ij\tuple{p}}^\tuple{h}
=\frac12\sum_{\tuple{q}\in\cellindices}\bigl(
b_{i,\tuple{p}-\tuple{q}+1}^\tuple{h}u_{j\tuple{q}}^\tuple{h}
+b_{j,\tuple{p}-\tuple{q}+1}^\tuple{h}u_{i\tuple{q}}^\tuple{h}\bigr)
=\sum_{k=1}^d\sum_{\tuple{q}\in\cellindices}B_{ij\tuple{p}k\tuple{q}}^\tuple{h}
u_{k\tuple{q}}^\tuple{h},
```

where

```math
B_{ij\tuple{p}k\tuple{q}}^\tuple{h}
=\frac12\bigl(\delta_{jk}b_{i,\tuple{p}-\tuple{q}+1}^\tuple{h}
+\delta_{ik}b_{j,\tuple{p}-\tuple{q}+1}^\tuple{h}\bigr)
```

is the global strain-displacement operator. The above identity reads, in Fourier
space

```math
\hat{B}_{ij\tuple{n}k\tuple{1}}^\tuple{h}
=\frac12\bigl(\delta_{jk}\hat{b}_{i\tuple{n}}^\tuple{h}
+\delta_{ik}\hat{b}_{j\tuple{n}}^\tuple{h}\bigr),
```

where ``\tuple{q}`` has been set to ``(1, \ldots, 1)``.

#### The modal stiffness matrix

It is recalled that the strain energy ``U`` is defined as the following integral
over the whole unit-cell ``\Omega`` (``\lambda``, ``\mu``: Lamé coefficients)

```math
U=\frac12\int_\Omega\bigl[\lambda\bigl(\tr\tens\varepsilon^\tuple{h}\bigr)^2
+2\mu\,\tens\varepsilon^\tuple{h}\dbldot\tens\varepsilon^\tuple{h}\bigr]
\D x_1\cdots\D x_d.
```

The strain energy is a quadratic form of the nodal displacements. Owing to the
homogeneity of the material (the Lamé coefficients are constant over the
unit-cell) and the periodic boundary conditions, the strain energy takes a
simple expression in Fourier space

```math
\begin{equation}
\label{eq20210914060318}
U=1{2\lvert\tuple{N}\rvert}
\sum_{\tuple{n}\in\cellindices}
\conj\bigl(\hat{\vec u}_\tuple{n}^\tuple{h}\bigr)
\cdot\hat{\tens K}_\tuple{n}^\tuple{h}
\cdot\hat{\vec u}_\tuple{n}^\tuple{h},
\end{equation}
```

where ``\hat{\tens K}_\tuple{n}^\tuple{h}`` is the *modal stiffness matrix*,
which is computed by the method XXX.

!!! danger

    The modal stiffness matrix introduced above differs from the matrix initially
	introduced by Brisard ([2017](@ref bris2017)) by a factor
	``\lvert\tuple{h}\rvert``. More precisely,

	```math
	\hat{\tens K}_\tuple{n}^\tuple{h}
	\bigr\rvert_\text{Scapin}
	=\lvert\tuple{h}\rvert\hat{\tens K}_\tuple{n}^\tuple{h}
	\bigr\rvert_\text{Brisard (2017)}.
	```

	Such rescaling makes the relation between modal and nodal stiffness operators
	more natural (the former is the DFT of the latter).

The modal stiffness matrix introduced above is related to the *nodal* (global)
stiffness operator introduced in Sec. [Element stiffness operator](@ref) of
Chap. [On the d-dimensional brick element](@ref). Indeed, plugging the
definition of the [Discrete Fourier transforms](@ref) into
Eq. \eqref{eq20210914060318} delivers the following expression

```math
\begin{aligned}
U={}&\frac1{2\lvert\tuple{N}\rvert}\sum_{\tuple{n}, \tuple{p}, \tuple{q}\in\cellindices}
\conj\bigl[\vec u_\tuple{p}^\tuple{h}\exp\bigl(-\I\phi_\tuple{np}\bigr)\bigr]
\cdot\hat{\tens K}_\tuple{n}^\tuple{h}
\cdot\vec u_\tuple{q}^\tuple{h}\exp\bigl(-\I\phi_\tuple{nq}\bigr)\\
={}&\frac1{2\lvert\tuple{N}\rvert}\sum_{\tuple{n}, \tuple{p}, \tuple{q}\in\cellindices}
\conj\bigl(\vec u_\tuple{p}^\tuple{h}\bigr)
\cdot\hat{\tens K}_\tuple{n}^\tuple{h}
\cdot\vec u_\tuple{q}^\tuple{h}
\exp\bigl[\I\bigl(\phi_\tuple{np}-\phi_\tuple{nq}\bigr)\bigr]\\
={}&\frac1{2\lvert\tuple{N}\rvert}\sum_{\tuple{n}, \tuple{p}, \tuple{q}\in\cellindices}
\conj\bigl(\vec u_\tuple{p}^\tuple{h}\bigr)
\cdot\hat{\tens K}_\tuple{n}^\tuple{h}
\cdot\vec u_\tuple{q}^\tuple{h}
\exp\bigl[\I\bigl(\phi_{\tuple{n}, \tuple{p}-\tuple{q}+1}\bigr)\bigr]\\
={}&\frac12\sum_{\tuple{p},\tuple{q}}
\conj\bigl(\vec u_\tuple{p}^\tuple{h}\bigr)\cdot\tens K_\tuple{pq}^\tuple{h}
\cdot\vec u_\tuple{q}^\tuple{h},
\end{aligned}
```

where

```math
\tens{K}_{\tuple{pq}}^\tuple{h}
=\frac1{\lvert\tuple{N}\rvert}\sum_{\tuple{n}}
\hat{\tens K}_\tuple{n}^\tuple{h}
\exp\bigl(\I\phi_{\tuple{n},\tuple{p}-\tuple{q}+1}\bigr)
=\dft_{\tuple{p}-\tuple{q}+1}^{-1}\bigl(\hat{\tens K}_\tuple{\bullet}^\tuple{h}\bigr)
```

is the *nodal stiffness matrix*, which appears as a block-circulant matrix. This
expresses the fact that the problem under consideration is
translation-invariant, owing to the homogeneity of the material and the periodic
boundary conditions (see also Sec. [On the d-dimensional brick element](@ref)).

In order to derive a FE-approximate solution to the problem described by
Eq. \eqref{eq20210914140834}, we need to account for the contribution
``U^\text{eigen}`` of the eigenstress ``\tens\tau^\tuple{h}`` to the total
potential energy ``\Pi``. This contribution reads, in the real space,

```math
U^\text{eigen}=\int_\Omega\bigl[\tens\tau^\tuple{h}(\vec x)
\dbldot\tens\varepsilon(\vec x)\bigr]\D x_1\cdots\D x_d
```

which can be expressed in Fourier space

```math
U^\text{eigen}=\frac{\lvert\tuple{h}\rvert}{\lvert\tuple{N}\rvert}
\sum_{\tuple{n}\in\cellindices}\conj\bigl(\hat{\vec u}_\tuple{n}\bigr)
\cdot\hat{\tens\tau}_\tuple{n}^\tuple{h}\cdot
\conj\hat{\tens B}_\tuple{n}^\tuple{h},
```

where ``\hat{\tens B}_\tuple{n}^\tuple{h}`` is the modal strain-displacement
operator introduced above, and ``\hat{\tens\tau}_\tuple{n}^\tuple{h}`` is the
DFT of the cell-values of the cell-wise constant eigenstrain.

Optimization of ``\Pi`` w.r.t. the nodal displacements delivers the following
equations

```math
\begin{equation}
\label{eq20210914145137}
\hat{\tens K}_\tuple{n}^\tuple{h}\cdot\hat{\vec u}_\tuple{n}^\tuple{h}
=-\hat{\tens\tau}_\tuple{n}^\tuple{h}\cdot\hat{\tens B}_\tuple{n}^\tuple{h},
\end{equation}
```

Note that these equations are *uncoupled*: they reduce to a ``d\times d`` linear
system *for each Fourier mode*. The solution to these equations delivers the
modal displacements

```math
\hat{\vec u}_\tuple{n}^\tuple{h}
=-\bigl(\hat{\tens K}_\tuple{n}^\tuple{h}\bigr)^{-1}
\cdot\hat{\tens\tau}_\tuple{n}^\tuple{h}
\cdot\hat{\tens B}_\tuple{n}^\tuple{h}.
```

Combining with Eq. \eqref{eq20210914144114}, we get the following expression of
the cell-average of the strains, in Fourier space

```math
\hat{\overline{\tens\varepsilon}}_\tuple{n}^\tuple{h}
=-\hat{\vec B}_\tuple{n}^\tuple{h}\symotimes
\bigl[\bigl(\hat{\tens K}_\tuple{n}^\tuple{h}\bigr)^{-1}
\cdot\hat{\tens\tau}_\tuple{n}^\tuple{h}\cdot\hat{\vec B}_\tuple{n}^\tuple{h}\bigr].
```

Upon appropriate symmetrization (``\tens\tau^\tuple{h}`` and
``\tens\varepsilon^\tuple{h}`` are both symmetric, second-rank tensors), we find

```math
\hat{\overline{\tens\varepsilon}}_\tuple{n}^\tuple{h}=-\hat{\tens\Gamma}_\tuple{n}^{\tuple{h}, \text{Bri17}}\dbldot\hat{\tens\tau}_\tuple{n}^\tuple{h},
```

with

```math
\hat{\tens\Gamma}_\tuple{n}^{\tuple{h}, \text{Bri17}}
=\hat{\vec B}_\tuple{n}^\tuple{h}
\symotimes\bigl(\tens K_\tuple{n}^\tuple{h}\bigr)^{-1}
\symotimes\hat{\vec B}_\tuple{n}^\tuple{h},
```

which defines the discrete Green operator.

!!! note

    Eq. \eqref{eq20210914145137} is singular for ``\vec k = \vec 0``. Indeed, in
	a periodic setting, the displacement is defined up to a constant translation.
	It is convenient to select the solution with zero average, that is
	``\hat{u}_\tuple{0}^\tuple{h}=\vec 0``.

### The discretization of Brisard and Dormieux

It was proved by Brisard and Dormieux [(2010)](@ref bris2010a) that, for all
``\tens\tau^{\tuple h}, \tens\varpi^{\tuple h}\in\tensors_2^{\tuple h}(\Omega)``

```math
\langle\tens\varpi^h\dbldot\tens\Gamma(\tens\tau^h)\rangle=\frac1{N^2}
\sum_{n_1=0}^{N_1-1}\cdots\sum_{n_d=0}^{N_d-1}\conj(\hat{\tens\varpi}_
{\tuple{n}}^h)\dbldot\hat{\tens\Gamma}_n^{h, \mathrm{BD10}}\dbldot
\hat{\tens\tau}_{\tuple{n}}^h,
```

where

```math
\hat{\tens\Gamma}_{\tuple n}^{\tuple h, \mathrm{BD10}}
=\sum_{\tuple m\in\integers^d}\bigl[F(\vec\alpha_{\tuple n+\tuple m\tuple N})
\bigr]^2\hat{\tens\Gamma}(\vec k_{\tuple n+\tuple m\tuple N}),
```

where ``\tuple{n+mN}`` denotes the ``d``-tuple: ``\tuple n+\tuple m\tuple
N=(n_1+m_1N_1, \ldots, n_d+m_dN_d)``, while ``\alpha_{\tuple{n}}`` is the
following dimensionless vector

```math
\vec\alpha_{\tuple{n}}
=\frac{2\pi h_1n_1}{L_1}\vec e_1+\cdots+\frac{2\pi h_dn_d}{L_d}\vec e_d,
```

finally, ``F`` is the tensor product of sine cardinal functions

```math
F(\vec\alpha)=\sinc\frac{\alpha_1}2\cdots\sinc\frac{\alpha_d}2.
```

Note that the above equations deliver the *exact* cell-averages of the Green
operator, applied exactly to any *cell-wise constant* polarization
field. Unfortunately, the series that defines the discrete Green operator can in
general not be evaluated, owing to very slow convergence. Therefore, this
discrete Green operator is unpractical, and is recalled here only for pedagocial
reasons.

### The discretization of Moulinec and Suquet

This is probably the most simple discretization, introduced first by Moulinec
and Suquet ([1994](@ref moul1994), [1998](@ref moul1998)). Only the lowest
(positive and negative) frequencies are kept

```math
\hat{\tens\Gamma}_{\tuple n}^{\tuple h, \mathrm{MS94}}=\hat{\tens\Gamma}
(\vec k_{\tuple Z(\tuple n, \tuple N)}).
```

We must check that the [Condition for the discrete Green operator to map real
fields onto real fields](@ref) is satisfied. Using the [Properties of the
``\fftfreq`` function](@ref) and assuming first that none of the ``n_i`` is such
that ``2n_i=N_i``

```math
\hat{\tens\Gamma}_{\tuple N-\tuple n}^{\tuple h, \mathrm{MS94}}
=\hat{\tens\Gamma}(\vec k_{\tuple Z(\tuple N-\tuple n, \tuple N)})
=\hat{\tens\Gamma}(\vec k_{-\tuple Z(\tuple n, \tuple N)})
=\hat{\tens\Gamma}(-\vec k_{\tuple Z(\tuple n, \tuple N)})
```

All Green operators presented in Chap. [Continuous Green operators](@ref) are
such that ``\hat{\tens\Gamma}(-\vec k)=\hat{\tens\Gamma}(\vec k)``, therefore

```math
\hat{\tens\Gamma}_{\tuple N-\tuple n}^{\tuple h, \mathrm{MS94}}
=\hat{\tens\Gamma}(\vec k_{\tuple Z(\tuple n, \tuple N)})
=\hat{\tens\Gamma}_{\tuple n}^{\tuple h, \mathrm{MS94}}
```

and the property is verified. Conversely, if all the ``n_i`` are such that
``2n_i=N_i``, then

```math
\hat{\tens\Gamma}_{\tuple N-\tuple n}^{\tuple h, \mathrm{MS94}}
=\hat{\tens\Gamma}(\vec k_{\tuple Z(\tuple N-\tuple n, \tuple N)})
=\hat{\tens\Gamma}(\vec k_{\tuple Z(\tuple n, \tuple N)})
=\hat{\tens\Gamma}_{\tuple n}^{\tuple h, \mathrm{MS94}}.
```

More problematic is the case when a few, but not all, ``n_i`` are such that
``2n_i=N_i``. Then the property does not hold for such frequencies. Moulinec and
Suquet ([1998](@ref moul1998)) use a specific treatment for such cases

```math
\hat{\tens\Gamma}(\vec k_{\tuple n})=\tens C^{-1},
```

if one of the ``n_i`` is such that ``2n_i=N_i``. This is implemented in
`Scapin`. Note that such cases occur only for even-sized grids.
