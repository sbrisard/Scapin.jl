# On the design of the Scapin library

## On discrete Green operators

### Real vs. complex Fourier transforms

The interface for discrete Green operators should allow for switching between
real and complex Fourier transforms.

As of
rev. [7bd9a4e](https://github.com/sbrisard/Scapin.jl/commit/7bd9a4ed36b3c9e8575dff384f34f9352e9e1845),
two functions are implemented, namely `eltype_real` and `eltype_fourier`. This
might be overkill. Indeed, it is sufficient to define the scalar type in the
Real space, `T`.
