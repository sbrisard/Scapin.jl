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
  \mathcal F(\tgrad\tens T)(\vec k)=\mathcal F(\tens T)(\vec k)\otimes\I\vec k
  \quad\text{and}\quad
  \mathcal F(\tdiv\tens T)(\vec k)=\mathcal F(\tens T)(\vec k)\cdot\I\vec k.
\end{equation}
```

When no confusion is possible, we will use the tilde to denote the Fourier
coefficients: ``\tilde{\tens T}_n=\mathcal F(\tens T)(\vec k_n)``.
