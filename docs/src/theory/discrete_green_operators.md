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
\begin{equation}
\dft_{\tuple{n}}(x)=\sum_{p_1=0}^{N_1-1}\cdots\sum_{p_d=0}^{N_d-1}
\exp\Bigl[-2\I\PI\Bigl(\frac{n_1p_1}{N_1}
+\cdots+\frac{n_dp_d}{N_d}\Bigr)\Bigr] x_{\tuple{p}}.
\end{equation}
```

Note that in the above definition, no restrictions are applied to the
multi-index ``\tuple{n}``. However, it can be verified that the above series of
tensors is in fact ``\tuple{N}``-periodic:
``\dft_{\tuple{n}+\tuple{N}}(x)=\dft_{\tuple{n}}(x)``, where
``\tuple{n}+\tuple{N}=(n_1+N_1, \ldots, n_d+N_d)``. Therefore, the
``\tuple{n}``-index is effectively restricted to ``0\leq n_i<N_i`` as well. The
most important results concerning the DFT are the *inversion formula*

```math
\begin{equation}
  x_{\tuple{p}}=\frac1{\lvert\tuple{N}\rvert}\sum_{n_1=0}^{N_1-1}\cdots
  \sum_{n_d=0}^{N_d-1}\exp\Bigl[2\I\PI\Bigl(\frac{n_1p_1}{N_1}+\cdots+
  \frac{n_dp_d}{N_d}\Bigr)\Bigr]\dft_{\tuple{n}}(x),
\end{equation}
```

the *Plancherel theorem*

```math
\begin{equation}
  \label{eq:20203128093105}
  \sum_{p_1=1}^{N_1-1}\cdots\sum_{p_d=1}^{N_d-1}\conj(x_{\tuple p})y_{\tuple p}
  =\frac1{\lvert\tuple N\rvert}\sum_{n_1=1}^{N_1-1}\cdots\sum_{n_d=1}^{N_d-1}
  \conj[\dft_{\tuple n}(x)]\dft_{\tuple n}(y),
\end{equation}
```

and the *circular convolution theorem*

```math
\begin{equation}
  \dft_{\tuple{n}}(x\ast y)=\dft_{\tuple{n}}(x)\dft_{\tuple{n}}(y),
  \quad\text{where}\quad
  (x\ast y)_{\tuple p}=\sum_{q_1=0}^{N_1-1}\cdots\sum_{q_d=0}^{N_d-1}
  x_{\tuple{q}}y_{\tuple{p}-\tuple{q}}.
\end{equation}
```

The DFT is readily extended to tensor data points. In the absence of ambiguity,
the shorthand ``\hat{x}_{\tuple{n}}`` will be adopted for
``\dft_{\tuple{n}}(x)``.

To close this section, we observe that the DFT of a series of *real* data points
is a series of *complex* data points. However, these complex values have the
following property

```math
\begin{equation}
  \dft_{\tuple{N}-\tuple{n}}(x)=\conj[\dft_{\tuple{n}}(x)].
\end{equation}
```

The above condition is actually a *necessary and sufficient* condition for the
``x_{\tuple{p}}`` to be real.

### The ``\fftfreq`` function

For ``n, N\in\naturals``, ``0\leq n<N``, we introduce ``\fftfreq(n, N)``

```math
\begin{equation}
  \fftfreq(n, N)=
  \begin{cases}
    n & \text{if }2n<N,\\
    n-N & \text{otherwise.}
  \end{cases}
\end{equation}
```

For ``n<0`` or ``n\geq N``, ``\fftfreq(n, N)`` is defined by
``N``-periodicity. ``\fftfreq`` is very similar to the NumPy
[fftfreq](https://numpy.org/doc/1.18/reference/generated/numpy.fft.fftfreq.html#numpy.fft.fftfreq)
function. We have the important result (see proof in
Sec. [Properties of the ``\fftfreq`` function](@ref) of the [Appendix](@ref))

```math
\begin{equation}
  \label{eq:20202503052504}
  \tag{Prop. Z}
  \fftfreq(N-n, N)=
  \begin{cases}
    \fftfreq(n) & \text{if }2n=N,\\
    -\fftfreq(n) & \text{otherwise.}
  \end{cases}
\end{equation}
```

The ``\fftfreq`` function can be defined for ``d``-tuples as well

```math
\begin{equation}
  \tuple \fftfreq(\tuple n, \tuple N)
  =\bigl(\fftfreq(n_1, N_1), \ldots, \fftfreq(n_d, N_d)\bigr)
\end{equation}
```

and we have again

```math
\begin{equation}
  \tuple \fftfreq(\tuple N-\tuple n, \tuple N)=-\tuple \fftfreq(\tuple n)
\end{equation}
```

if none of the ``n_i`` is such that ``2n_i=N_i``.
