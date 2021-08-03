# Discrete Green operators

In this chapter, we introduce various discretizations of the Green operator; we
will adopt the vocabulary of linear elasticity, although the concepts apply to
all the various physical models presented in [Continuous Green operators](@ref).

## On the discrete Fourier transform

Let ``x=(x_{\tuple{p}})`` be a finite set of scalar values indexed by the
``d``-tuple ``\tuple{p}=(p_1, \ldots, p_d)`` where ``0\leq p_i<N_i`` (``N_i`` is
the number of data points in the ``i``-th direction). The discrete Fourier
transform is a discrete set of scalar values ``\dft_{\tuple{n}}(x)`` indexed by
the ``d``-tuple ``\tuple{n}\in\integers^d``, defined as follows

```math
\dft_{\tuple{n}}(x)=\sum_{p_1=0}^{N_1-1}\cdots\sum_{p_d=0}^{N_d-1}
\exp\Bigl[-2\I\PI\Bigl(\frac{n_1p_1}{N_1}
+\cdots+\frac{n_dp_d}{N_d}\Bigr)\Bigr] x_{\tuple{p}}.
```

Note that in the above definition, no restrictions are applied to the
multi-index ``\tuple{n}``. However, it can be verified that the above series of
tensors is in fact ``\tuple{N}``-periodic:
``\dft_{\tuple{n}+\tuple{N}}(x)=\dft_{\tuple{n}}(x)``, where
``\tuple{n}+\tuple{N}=(n_1+N_1, \ldots, n_d+N_d)``. Therefore, the
``\tuple{n}``-index is effectively restricted to ``0\leq n_i<N_i`` as well. The
most important results concerning the DFT are the *inversion formula*

```math
x_{\tuple{p}}=\frac1{\lvert\tuple{N}\rvert}\sum_{n_1=0}^{N_1-1}\cdots
\sum_{n_d=0}^{N_d-1}\exp\Bigl[2\I\PI\Bigl(\frac{n_1p_1}{N_1}+\cdots+
\frac{n_dp_d}{N_d}\Bigr)\Bigr]\dft_{\tuple{n}}(x),
```

the *Plancherel theorem*

```math
\sum_{p_1=1}^{N_1-1}\cdots\sum_{p_d=1}^{N_d-1}\conj(x_{\tuple p})y_{\tuple p}
=\frac1{\lvert\tuple N\rvert}\sum_{n_1=1}^{N_1-1}\cdots\sum_{n_d=1}^{N_d-1}
\conj[\dft_{\tuple n}(x)]\dft_{\tuple n}(y),
```

and the *circular convolution theorem*

```math
\dft_{\tuple{n}}(x\ast y)=\dft_{\tuple{n}}(x)\dft_{\tuple{n}}(y),
\quad\text{where}\quad
(x\ast y)_{\tuple p}=\sum_{q_1=0}^{N_1-1}\cdots\sum_{q_d=0}^{N_d-1}
x_{\tuple{q}}y_{\tuple{p}-\tuple{q}}.
```

The DFT is readily extended to tensor data points. In the absence of ambiguity,
the shorthand ``\hat{x}_{\tuple{n}}`` will be adopted for
``\dft_{\tuple{n}}(x)``.

To close this section, we observe that the DFT of a series of *real* data points
is a series of *complex* data points. However, these complex values have the
following property

```math
\dft_{\tuple{N}-\tuple{n}}(x)=\conj[\dft_{\tuple{n}}(x)].
```

The above condition is actually a *necessary and sufficient* condition for the
``x_{\tuple{p}}`` to be real.

### [The ``\fftfreq`` function](@id sec20210803055450)

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
function. We have the important result (see proof in
Sec. [Properties of the ``\fftfreq`` function](@ref) of the [Appendix](@ref))

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
\langle\tens\varpi^{\tuple{h}}\dbldot\tens\Gamma(\tens\tau^{\tuple{h}})\rangle
\simeq\frac1{\lvert\tuple{N}\rvert}\sum_{\tuple{p}, \tuple{q}\in\cellindices}
\tens\varpi_{\tuple{p}}^{\tuple{h}}\dbldot
\tens\Gamma_{\tuple{p}\tuple{q}}^{\tuple{h}}\dbldot\tens
\tau_{\tuple{q}}^{\tuple{h}}
\quad\text{for all}\quad
\tens\tau^{\tuple{h}},\tens\varpi^{\tuple{h}}\in\tensors_2^{\tuple{h}}(\Omega).
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

### [Condition for the discrete Green operator to map real fields onto real fields](@id sec20210803053100)

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

### The discretization of \textcite{bris2010a}

It was proved by \textcite{bris2010a} that, for all ``\tens\tau^{\tuple h},
\tens\varpi^{\tuple h}\in\tensors_2^{\tuple h}(\Omega)``

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

## The discretization of \textcite{moul1994, moul1998}

This is probably the most simple discretization, introduced first by
\textcite{moul1994}. Only the lowest (positive and negative) frequencies are
kept

```math
\hat{\tens\Gamma}_{\tuple n}^{\tuple h, \mathrm{MS94}}=\hat{\tens\Gamma}
(\vec k_{\tuple Z(\tuple n, \tuple N)}).
```

We must check that the [property ensuring realness of the discrete Green
operator](@ref sec20210803053100) is satisfied. Using the [fundamental
property](@ref sec20210803055450) of the ``\fftfreq`` function and assuming
first that none of the ``n_i`` is such that ``2n_i=N_i``

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
``2n_i=N_i``. Then the property does not hold for such
frequencies. \textcite{moul1998} use a specific treatment for such cases

```math
\hat{\tens\Gamma}(\vec k_{\tuple n})=\tens C^{-1},
```

if one of the ``n_i`` is such that ``2n_i=N_i``. This is implemented in
`Scapin`. Note that such cases occur only for even-sized grids.
