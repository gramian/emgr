#!/usr/bin/env python3

"""
  project: emgr ( https://gramian.de )
  version: 5.9.py (2021-01-21)
  authors: Christian Himpe (0000-0003-2194-6754)
  license: BSD-2-Clause License (opensource.org/licenses/BSD-2-Clause)
  summary: Factorial empirical Gramian singular value decay testing
"""

import numpy as np
import matplotlib.pyplot as plt
from emgr import emgr


def emgrProbe():
    """Factorial test of emgr"""

    M = 1
    N = 16
    Q = M

    A = -2.0 * np.eye(N) + np.eye(N, N, 1) + np.eye(N, N, -1)
    A[0, 0] = -1.0

    B = np.zeros((N, 1))
    B[0] = 1.0

    C = B.T

    def F(x, u, p, t):  # vector field
        return A.dot(x) + B.dot(u) + p

    def H(x, u, p, t):  # adjoint vector field
        return A.T.dot(x) + C.T.dot(u)

    def G(x, u, p, t):  # output functional
        return C.dot(x)

    dt = 0.01
    Tf = 1.0

    p = np.zeros((N, 1))
    q = np.ones((N, 1)).dot([[0.5, 1.0]])

    gramian = ["c", "o", "x", "y", "s", "i", "j"]  # controllability, observability, minimality, linear minimality, sensitivity, identifiability, cross-identifiability

    kernels = [np.dot] #, quadratic, cubic, sigmoid]
    training = ["i", "s", "h", "a", "r"]  # impulse, step, havercosine-chirp, sinc, pseudo-random-binary
    weighting = [0] #, 1, 2, 3, 4]  # none, linear-time, quadratic-time, per-state, per-component
    centering = [0] #, 1, 2, 3, 4, 5]  # none, steady-state, final-state, arithmetic-mean, root-mean-square, midrange
    normalization = [0, 1, 2]  # none, steady-state, Jacobi-type
    stype = [0, 1]  # standard, (output-controllability, average-observability, non-symmetric-cross-gramian)
    extra = [0, 1]  # none, extra-input

    scales = [0]#, 1, 2, 3, 4]  # single, linear, geometric, logarithmic, sparse
    rotations = [0]#, 1]  # positive-negative, single

    pcentering = [0]#, 1, 2]  # none, linear, logarithmic
    ptype = [0, 1]  # standard (input-output-sensitivity, coarse-schur-complement)

    z = 0

    for w in gramian:

        W = []

        K = H if w == "y" else G

        for c in training:
            for d in centering:
                for e in normalization:
                    for f in stype:
                        for g in extra:
                            for h in weighting:
                                for i in kernels:

                                    if w in ["s", "i", "j"]:

                                        for j in ptype:
                                            for k in pcentering:
                                                W.append(np.linalg.svd(emgr(F, K, [M, N, Q], [dt, Tf], w, q, [d, 0, 0, 0, 0, e, f, g, 0, j, k, 0, h], c, 0.0, 0.0, 1.0, 1.0, i)[1], compute_uv=False))
                                    else:
                                        W.append(np.linalg.svd(emgr(F, K, [M, N, Q], [dt, Tf], w, p, [d, 0, 0, 0, 0, e, f, g, 0, 0, 0, 0, h], c, 0.0, 0.0, 1.0, 1.0, i), compute_uv=False))

        print(w)
        z += 1
        plotsingvals(W, (2, 4, z), w)

    plt.show()


def quadratic(x, y):
    return (x.dot(y))**2.0 + 1.0


def cubic(x, y):
    return (x.dot(y))**3.0 + 1.0


def sigmoid(x, y):
    return np.tanh(x.dot(y) - 1.0)


def plotsingvals(W, s, l):

    plt.subplot(s[0], s[1], s[2])

    for k in range(0, len(W)):

        plt.semilogy(W[k] / W[k].max())

    plt.xlim(0, len(W[0]) - 1)
    mi, ma = plt.ylim()
    plt.ylim(min(mi, 0.1), 1.1)
    plt.xlabel(l)

    return


emgrProbe()
