# Appendix

## Properties of the ``\fftfreq`` function

In this paragraph, we prove the following property of the ``\fftfreq``
function

```math
\fftfreq(N-n, N)=
\begin{cases}
\fftfreq(n) & \text{if }2n=N,\\
-\fftfreq(n) & \text{otherwise.}
\end{cases}
```

Several cases must be considered.

1. If ``n=0``
   ```math
   \fftfreq(N-n, N)=\fftfreq(N, N)=\fftfreq(0, N)=0=-\fftfreq(n, N)
   ```
2. If ``N`` is even, ``N=2M``
   1. If ``0<n<M``
      ```math
      \begin{gather*}
      2n<N\quad\Rightarrow\quad \fftfreq(n, N)=n,\\
      \begin{aligned}
      M<N-n&\quad\Rightarrow\quad N<2\bigl(N-n\bigr)\\
      &\quad\Rightarrow\quad \fftfreq(N-n, N)=N-n-N=-n.
      \end{aligned}
      \end{gather*}
      ```
   2. If ``n=M``
      ```math
      \begin{gather*}
      2n=N\quad\Rightarrow\quad \fftfreq(n, N)=n-N=-M,\\
      2(N-n)=2M=N\quad\Rightarrow\quad \fftfreq(N-n, N)=-M.
      \end{gather*}
      ```
   3. If ``M<n<N``
      ```math
      \begin{gather*}
      N<2n\quad\Rightarrow\quad \fftfreq(n, N)=n-N,\\
      \begin{aligned}
      N-n<M&\quad\Rightarrow\quad 2\bigl(N-n\bigr)<N\\
      &\quad\Rightarrow\quad \fftfreq(N-n, N)=N-n.
      \end{aligned}
      \end{gather*}
      ```
3. If ``N`` is odd, ``N=2M+1``
   1. If ``0<n\leq M``
      ```math
      \begin{gather*}
      2n\leq 2M<N\quad\Rightarrow\quad \fftfreq(n, N)=n,\\
      \begin{aligned}
      M+1\leq N-n&\quad\Rightarrow\quad N<2\bigl(N-n\bigr)\\
      &\quad\Rightarrow\quad \fftfreq(N-n, N)=N-n-N=-n.
      \end{aligned}
      \end{gather*}
      ```
   2. If ``M+1\leq n<N``
      ```math
      \begin{gather*}
      N+1\leq 2n\quad\Rightarrow\quad \fftfreq(n, N)=n-N\\
      \begin{aligned}
      N-n\leq N-M-1=M&\quad\Rightarrow\quad 2\bigl(N-n\bigr)\leq 2M<N\\
      &\quad\Rightarrow\quad \fftfreq(N-n, N)=N-n
      \end{aligned}
      \end{gather*}
      ```
