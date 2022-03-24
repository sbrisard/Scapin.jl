"""Various tests to create in-place FFT plans."""

using FFTW

N = (7, 8)
a = rand(complex(Float64), N...)
p = plan_fft!(a)

@assert a === p * a

b = rand(complex(Float64), N...)
b_ref = copy(b)

bb = p * b

bbb = p \ b

@assert bbb ≈ b_ref
