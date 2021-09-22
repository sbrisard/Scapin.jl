# Fourier transforms in a periodic setting

This chapter provides a brief overview of the two Fourier transforms that are
going to be used in the previous document, namely [Fourier series](@ref) and
[Discrete Fourier transforms](@ref).

## Fourier series

Owing to the periodic setting, the fields that are involved in the various BVPs
to be discussed in this document are expanded in Fourier series. ``\tens T``
being a ``\Omega``-periodic tensor field (with sufficient regularity), the
following decomposition holds

```math
\tens T(\vec x)
=\sum_{\tuple{n}∈\integers^d}\mathcal F(\tens T)(\vec k_{\tuple{n}})
\exp(\I\vec k_{\tuple{n}}\cdot\vec x),
```

where ``\tuple{n}`` denotes a ``d``-dimensional tuple of integers (see
[Nomenclature](@ref)). The wave vectors ``\vec k_{\tuple{n}}`` are given by

```math
\vec k_{\tuple{n}}=\frac{2\PI n_1}{L_1}\vec e_1+\frac{2\PI n_2}{L_2}\vec e_2
+\cdots+\frac{2\PI n_d}{L_d}\vec e_d,
```

and the Fourier coefficients of ``\tens T`` are defined as follows

```math
\mathcal F(\tens T)(\vec k)=\frac1{\lvert\Omega\rvert}\int_{\vec x\in\Omega}
\tens T(\vec x)\exp(-\I\vec k\cdot\vec x)\,\D x_1\cdots\D x_d.
```

It is recalled that the Fourier coefficients of the gradient and divergence of
``\tens T`` can readily be computed from the Fourier coefficients of ``\tens T``

```math
\mathcal F(\tens T\otimes\nabla)(\vec k)=\mathcal F(\tens T)(\vec k)\otimes\I\vec k
\quad\text{and}\quad
\mathcal F(\tens T\cdot\nabla)(\vec k)=\mathcal F(\tens T)(\vec k)\cdot\I\vec k.
```

When no confusion is possible, we will use a tilde to denote the Fourier
coefficients: ``\tilde{\tens T}_\tuple{n}=\mathcal F(\tens T)(\vec
k_\tuple{n})``.

## Discrete Fourier transforms

!!! note "Naming conventions for indices"

    In the present section, we use the following convention for naming indices:

	- italic, latin indices (``i``, ``j``, …) refer to components of a tensor in
	  the ``d``-dimensional space,
    - upright, sans-serif, indices are multi-indices over a cartesian grid,
    - ``\tuple{p}``, ``\tuple{q}`` are multi-indices in the real space (“pixels”),
    - ``\tuple{m}``, ``\tuple{n}`` are multi-indices in the frequency space.

We consider a ``d``-dimensional grid of size ``N_1\times\cdots\times N_d``. The
set of cell indices over this grid is

```math
\cellindices=\{1,\ldots, N_1\}\times\cdots\times\{1, \ldots, N_d\}.
```

For ``\tuple{n}, \tuple{p}\in\cellindices``, we define
$\Phi_{\tuple{n}\tuple{p}}$

```math
\phi_{\tuple{n}\tuple{p}}=2\PI\frac{\bigl(n_1-1\bigr)\bigl(p_1-1\bigr)}{N_1}
+\cdots+2\PI\frac{\bigl(n_d-1\bigr)\bigl(p_d-1\bigr)}{N_d}.
```

!!! warning "Indexing conventions"

    Owing to Julia being one-based, all multi-indices start at one, hence the
	“``-1``” correction in the above formula.

Let ``x=(x_{\tuple{p}})`` be a finite set of scalar values indexed by the
``d``-tuple ``\tuple{p}\in\cellindices``. The discrete Fourier transform is a
discrete set of scalar values ``\dft_{\tuple{n}}(x)`` indexed by the ``d``-tuple
``\tuple{n}\in\integers^d``, defined as follows

```math
\dft_{\tuple{n}}(x)=\sum_{\tuple{p}\in\cellindices}
\exp\bigl(-\I\phi_{\tuple{n}\tuple{p}}\bigr)x_{\tuple{p}}.
```

Note that in the above definition, no restrictions are applied to the
multi-index ``\tuple{n}``. However, it can be verified that the above series of
tensors is in fact ``\tuple{N}``-periodic:
``\dft_{\tuple{n}+\tuple{N}}(x)=\dft_{\tuple{n}}(x)``, where
``\tuple{n}+\tuple{N}=(n_1+N_1, \ldots, n_d+N_d)``. Therefore, the
``\tuple{n}``-index is effectively restricted to ``0\leq n_i<N_i`` as well. The
most important results concerning the DFT are the *inversion formula*

```math
x_{\tuple{p}}=\frac1{\lvert\tuple{N}\rvert}\sum_{\tuple{n}\in\cellindices}
\exp\bigl(\I\phi_{\tuple{n}\tuple{p}}\bigr)\dft_{\tuple{n}}(x),
```

the *Plancherel theorem*

```math
\sum_{\tuple{p}\in\cellindices}\conj(x_{\tuple p})y_{\tuple p}
=\frac1{\lvert\tuple N\rvert}\sum_{\tuple{n}\in\cellindices}
\conj[\dft_{\tuple n}(x)]\dft_{\tuple n}(y),
```

and the *circular convolution theorem*

```math
\dft_{\tuple{n}}(x\ast y)=\dft_{\tuple{n}}(x)\dft_{\tuple{n}}(y),
\quad\text{where}\quad
(x\ast y)_{\tuple p}=\sum_{\tuple{q}\in\cellindices}
x_{\tuple{q}}y_{\tuple{p}-\tuple{q}+1}.
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
