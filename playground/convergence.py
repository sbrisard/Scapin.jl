import numpy as np
import matplotlib.pyplot as plt

import janus
import janus.material.elastic.linear.isotropic as material
import janus.operators as operators
import janus.fft.serial as fft
import janus.green as green


class ConvergenceTest:
    def __init__(self):
        self.dim = 2
        self.r_max = 9
        self.shape_coarse = self.dim * (4,)
        self.shape_fine = tuple((2 ** self.r_max) * n_ for n_ in self.shape_coarse)
        self.L = np.array(self.dim * (1.0,))
        self.patch_ratio = self.dim * (0.25,)
        self.C0 = material.create(1.0, 0.3, dim=self.dim)
        self.Γ0 = self.C0.green_operator()
        sym = (self.dim * (self.dim + 1)) // 2
        self.τ_in = np.zeros((sym,), dtype=np.float64)
        self.τ_in[-1] = 1.0
        self.τ_out = np.zeros_like(self.τ_in)

    def create_tau(self, r):
        shape = tuple(n_ >> (self.r_max - r) for n_ in self.shape_fine)
        patch_size = tuple(int(r_ * n_) for r_, n_ in zip(self.patch_ratio, shape))
        τ = np.empty(shape + self.τ_in.shape, dtype=np.float64)
        τ[...] = self.τ_out
        τ[tuple(slice(n_) for n_ in patch_size) + (...,)] = self.τ_in
        return τ

    def run(self, r):
        τ = self.create_tau(r)
        shape = τ.shape[:-1]
        Γ = green.filtered(self.Γ0, shape, 1.0, fft.create_real(shape))
        η = np.empty_like(τ)
        Γ.apply(τ, η)
        return η


def auxiliary_shapes(shape_coarse, shape_fine):
    aux_shape_coarse = []
    aux_shape_fine = []

    for i, (nc_i, nf_i) in enumerate(zip(shape_coarse, shape_fine)):
        aux_shape_coarse.append(nc_i)
        aux_shape_coarse.append(1)
        aux_shape_fine.append(nc_i)
        aux_shape_fine.append(nf_i // nc_i)

    return tuple(aux_shape_coarse[:-1]), tuple(aux_shape_fine[:-1])


if __name__ == "__main__":
    test = ConvergenceTest()
    results = []

    for r in range(test.r_max + 1):
        η = test.run(r)
        results.append(η)
        fig, ax = plt.subplots(1, η.shape[-1])
        for i, ax_i in enumerate(ax):
            ax_i.set_axis_off()
            ax_i.imshow(η[..., i].real)

    η_ref = results[-1]
    size = []
    err = []

    for η in results:
        new_shape_coarse, new_shape_fine = auxiliary_shapes(η.shape, η_ref.shape)
        η_zoom = np.empty_like(η_ref)
        η_zoom.reshape(new_shape_fine)[...] = η.reshape(new_shape_coarse)
        size.append(η.shape[0])
        err.append(np.linalg.norm(η_zoom - η_ref))

    size = np.asarray(size, dtype=np.float64)
    err = np.asarray(err, dtype=np.float64)
    C = size[-2] * err[-2]
    plt.figure()
    plt.loglog(size[:-1], C * size[:-1] ** -1)
    plt.loglog(size[:-1], err[:-1], "o-")
