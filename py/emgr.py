"""
emgr - EMpirical GRamian Framework
==================================

  project: emgr ( https://gramian.de )
  version: 5.7-py ( 2019-02-26 )
  authors: Christian Himpe ( 0000-0003-2194-6754 )
  license: BSD-2-Clause License ( opensource.org/licenses/BSD-2-Clause )
  summary: Empirical system Gramians for (nonlinear) input-output systems.

USAGE:
------

  W = emgr(f,g,s,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm],[dp])

DESCRIPTION:
------------

  Empirical gramian matrix and empirical covariance matrix computation
  for model reduction, decentralized control, nonlinearity quantification,
  sensitivity analysis, parameter identification, uncertainty quantification &
  combined state and parameter reduction of large-scale input-output systems.
  Data-driven analysis of input-output coherence and system-gramian-based
  nonlinear model order reduction. Compatible with PYTHON2 and PYTHON3.

ALGORITHM:
----------

  C. Himpe (2018). emgr - The Empirical Gramian Framework. Algorithms 11(7):91
  `doi:10.3390/a11070091 <https://doi.org/10.3390/a11070091>`_

ARGUMENTS:
----------

   f {function} vector field handle: x' = f(x,u,p,t)
   g {function} output function handle: y = g(x,u,p,t)
   s {tuple} system dimensions: [inputs,states,outputs]
   t {tuple} time discretization: [time-step,time-horizon]
   w {string} single character encoding gramian type:
    * "c" empirical controllability gramian (Wc)
    * "o" empirical observability gramian (Wo)
    * "x" empirical cross gramian (Wx aka Wco or Xcg)
    * "y" empirical linear cross gramian (Wy)
    * "s" empirical sensitivity gramian (Ws)
    * "i" empirical identifiability gramian (Wi)
    * "j" empirical joint gramian (Wj)
  pr {matrix|0} parameters, each column is a set
  nf {tuple|0} option flags, twelve components, default zero:
    * center: no(0), steady(1), last(2), mean(3), rms(4), midr(5), geom(6)
    * input scales: single(0), linear(1), geom(2), log(3), sparse(4)
    * state scales: single(0), linear(1), geom(2), log(3), sparse(4)
    * input rotations: unit(0), single(1)
    * state rotations: unit(0), single(1)
    * normalization (only: Wc, Wo, Wx, Wy): none(0), Jacobi(1), steady(2)
    * state gramian variant:
      * cross gramian type (only: Wx, Wy, Wj): regular(0), non-symmetric(1)
      * observability gramian type (only: Wo, Wi): regular(0), averaged(1)
    * extra input (only: Wo, Wx, Ws, Wi, Wj): none(0), yes(1)
    * parameter centering (only: Ws, Wi, Wj): none(0), linear(1), log(2)
    * parameter gramian variant:
      * averaging type (only: Ws): input-state(0), input-output(1)
      * Schur-complement (only: Wi, Wj): detailed(0), approximate(1)
    * cross gramian partition size (only: Wx, Wj): full(0), partitioned(<N)
    * cross gramian partition index (only: Wx, Wj): partition(>=0)
  ut {function|"i"} input function handle: u_t = ut(t), or string
    * "i" delta impulse input
    * "s" step input / load vector / source term
    * "c" decaying exponential chirp input
    * "r" pseudo-random binary input
  us {vector|0} steady-state input
  xs {vector|0} steady-state and nominal initial state x_0
  um {matrix|1} input scales
  xm {matrix|1} initial-state scales
  dp {function|np.dot} custom inner product handle: xy = dp(x,y)

RETURNS:
--------

  W {matrix} Gramian Matrix (for: Wc, Wo, Wx, Wy)
  W {tuple}  [State-, Parameter-] Gramian (for: Ws, Wi, Wj)

CITE AS:
--------

  C. Himpe (2019). emgr - EMpirical GRamian Framework (Version 5.7)
  [Software]. Available from https://gramian.de . doi:10.5281/zenodo.2577980

KEYWORDS:
---------

  model reduction, system gramians, empirical gramians, cross gramian, MOR

SEE ALSO:
---------

  gram (Python Control Systems Library)
"""

import math
import numpy as np

__version__ = "5.7"
__copyright__ = "Copyright (C) 2019 Christian Himpe"
__author__ = "Christian Himpe"
__license__ = "BSD 2-Clause"


ODE = lambda f, g, t, x0, u, p: ssp2(f, g, t, x0, u, p)  # Integrator Handle


STAGES = 3  # Configurable number of stages for increased stability of ssp2


def emgr(f, g=None, s=None, t=None, w=None, pr=0, nf=0, ut="i", us=0.0, xs=0.0, um=1.0, xm=1.0, dp=np.dot):
    """ Compute empirical system Gramian matrix """

    # Version Info
    if f == "version":
        return __version__

    # Default Arguments
    if type(pr) in {int, float} or np.ndim(pr) == 1:
        pr = np.reshape(pr, (-1, 1))

    if nf == 0:
        nf = [0]

    """ SETUP """

    # System Dimensions
    M = int(s[0])                        # Number of inputs
    N = int(s[1])                        # Number of states
    Q = int(s[2])                        # Number of outputs
    A = int(s[3]) if len(s) == 4 else 0  # Number of augmented parameter-states
    P = pr.shape[0]                      # Dimension of parameter
    K = pr.shape[1]                      # Number of parameter-sets

    # Time Discretization
    dt = t[0]                            # Time-step width
    Tf = t[1]                            # Time horizon
    nt = int(math.floor(Tf / dt) + 1)    # Number of time-steps

    # Lazy Output Functional
    if type(g) == int and g == 1:
        g = ident
        Q = N

    # Pad Flag Vector
    if len(nf) < 12:
        nf = nf + [0] * (12 - len(nf))

    # Built-in input functions
    if type(ut) is str:
        if ut.lower() == "s":    # Step Input
            def ut(t):
                return 1

        elif ut.lower() == "c":  # Decaying Exponential Chirp Input
            a0 = (2.0 * math.pi) / (4.0 * dt) * Tf / math.log(4.0 * (dt / Tf))
            b0 = (4.0 * (dt / Tf)) ** (1.0 / Tf)

            def ut(t):
                return 0.5 * math.cos(a0 * (b0 ** t - 1)) + 0.5

        elif ut.lower() == "r":  # Pseudo-Random Binary Input
            def ut(t):
                return np.random.randint(0, 1, size=1)

        else:  # Delta Impulse Input

            def ut(t):
                return (t <= dt) / dt

    # Lazy Optional Arguments
    if type(us) in {int, float}: us = np.full(M, us)
    if type(xs) in {int, float}: xs = np.full(N, xs)
    if type(um) in {int, float}: um = np.full(M, um)
    if type(xm) in {int, float}: xm = np.full(N, xm)

    # Gramian Normalization
    if nf[5] and A == 0:

        if nf[5] == 1:    # Jacobi-type preconditioner
            NF = nf
            NF[5] = 0

            def DP(x, y):
                return np.sum(x * y.T, 1)  # Diagonal-only kernel
            WT = emgr(f, g, s, t, w, pr, NF, ut, us, xs, um, xm, DP)
            TX = np.sqrt(np.fabs(WT))

        elif nf[5] == 2:  # Steady-state preconditioner
            TX = xs

        TX[np.fabs(TX) < np.sqrt(np.spacing(1))] = 1

        def deco(f, g):
            def F(x, u, p, t):
                return f(TX * x, u, p, t) / TX

            def G(x, u, p, t):
                return g(TX * x, u, p, t)

            return F, G

        f, g = deco(f, g)

        xs = xs / TX

    # Non-symmetric cross Gramian or average observability Gramian
    if nf[6]:
        R = 1
    else:
        R = Q

    # Extra Input
    if nf[7]:
        def up(t):
            return us + ut(t)
    else:
        def up(t):
            return us

    # Scale Sampling
    if um.ndim == 1: um = np.outer(um, scales(nf[1], nf[3]))
    if xm.ndim == 1: vm = np.outer(xm[0:Q], scales(nf[1], nf[3]))
    if xm.ndim == 1: xm = np.outer(xm, scales(nf[2], nf[4]))

    C = um.shape[1]  # Number of input scales sets
    D = xm.shape[1]  # Number of state scales sets

    """ GRAMIAN COMPUTATION """

    W = 0.0  # Reserve gramian variable

    # General layout:
    #   Parameter gramians call state gramians
    #   For each {parameter, scale, input/state/parameter}:
    #     Perturb, simulate, center, normalize, accumulate
    #   Assemble, normalize, post-process

    if w == "c":    # Empirical Controllability Gramian

        for k in range(K):
            for c in range(C):
                for m in np.nditer(np.nonzero(um[:, c])):
                    em = np.zeros(M + P)
                    em[m + (A > 0) * M] = um[m, c]

                    def umc(t):
                        return up(t) + ut(t) * em[0:M]
                    pmc = pr[:, k] + em[M:M + P]
                    x = ODE(f, ident, t, xs, umc, pmc)
                    x -= avg(x, nf[0], xs)
                    x /= um[m, c]
                    if A > 0:
                        W += em[M:M + P] * dp(x, x.T)
                    else:
                        W += dp(x, x.T)
        W *= dt / (C * K)
        return W

    elif w == "o":  # Empirical Observability Gramian

            o = np.zeros((R * nt, N + A))  # Pre-allocate observability matrix
            for k in range(K):
                for d in range(D):
                    for n in np.nditer(np.nonzero(xm[:, d])):
                        en = np.zeros(N + P)
                        en[n] = xm[n, d]
                        xnd = xs + en[0:N]
                        pnd = pr[:, k] + en[N:N + P]
                        ys = g(xnd, us, pnd, 0)
                        y = ODE(f, g, t, xnd, up, pnd)
                        y -= avg(y, nf[0], ys)
                        y /= xm[n, d]
                        if nf[6]:  # Average observability gramian
                            o[:, n] = np.sum(y, 0)
                        else:      # Regular observability gramian
                            o[:, n] = y.flatten(1)
                    W += dp(o.T, o)
            W *= dt / (D * K)
            return W

    elif w == "x":  # Empirical Cross Gramian

        assert M == Q or nf[6], "emgr: non-square system!"

        i0 = 0
        i1 = N + A

        # Partitioned cross gramian
        if nf[10] > 0:
            sp = int(round(nf[10]))  # Partition size
            ip = int(round(nf[11]))  # Partition index
            i0 += ip * sp            # Start index
            i1 = min(i0 + sp, N)     # End index
            if i0 > N:
                i0 -= math.ceil(N / sp) * sp - N
                i1 = min(i0 + sp, N + A)

            if ip < 0 or i0 >= i1 or i0 < 0:
                return 0

        o = np.zeros((R, nt, i1 - i0))  # Pre-allocate observability 3-tensor
        for k in range(K):
            for d in range(D):
                for n in np.nditer(np.nonzero(xm[i0:i1, d])):
                    en = np.zeros(N + P)
                    en[i0 + n] = xm[i0 + n, d]
                    xnd = xs + en[0:N]
                    pnd = pr[:, k] + en[N:N + P]
                    ys = g(xs, us, pnd, 0)
                    y = ODE(f, g, t, xnd, up, pnd)
                    y -= avg(y, nf[0], ys)
                    y /= xm[i0 + n, d]
                    if nf[6]:  # Non-symmetric cross gramian
                        o[0, :, n] = np.sum(y, axis=0)
                    else:      # Regular cross gramian
                        o[:, :, n] = y
                for c in range(C):
                    for m in np.nditer(np.nonzero(um[:, c])):
                        em = np.zeros(M)
                        em[m] = um[m, c]

                        def umc(t):
                            return us + ut(t) * em
                        x = ODE(f, ident, t, xs, umc, pr[:, k])
                        x -= avg(x, nf[0], xs)
                        x /= um[m, c]
                        if nf[6]:  # Non-symmetric cross gramian
                            W += dp(x, o[0, :, :])
                        else:      # Regular cross gramian
                            W += dp(x, o[m, :, :])
        W *= dt / (C * D * K)
        return W

    elif w == "y":  # Empirical Linear cross Gramian

        assert M == Q or nf[6], "emgr: non-square system!"
        assert C == vm.shape[1], "emgr: scale count mismatch!"

        a = np.zeros((Q, nt, N))  # Pre-allocate adjoint 3-tensor
        for k in range(K):
            for c in range(C):
                for q in np.nditer(np.nonzero(vm[:, c])):
                    em = np.zeros(Q)
                    em[q] = vm[q, c]

                    def vqc(t):
                        return us + ut(t) * em

                    z = ODE(g, ident, t, xs, vqc, pr[:, k])
                    z -= avg(z, nf[0], xs)
                    z /= vm[q, c]
                    if nf[6]:  # Non-symmetric cross gramian
                        a[0, :, :] += z.T
                    else:      # Regular cross gramian
                        a[q, :, :] = z.T
                for m in np.nditer(np.nonzero(um[:, c])):
                    em = np.zeros(M)
                    em[m] = um[m, c]

                    def umc(t):
                        return us + ut(t) * em
                    x = ODE(f, ident, t, xs, umc, pr[:, k])
                    x -= avg(x, nf[0], xs)
                    x /= um[m, c]
                    if nf[6]:  # Non-symmetric cross gramian
                        W += dp(x, a[0, :, :])
                    else:      # Regular cross gramian
                        W += dp(x, a[m, :, :])
        W *= dt / (C * K)
        return W

    elif w == "s":  # Empirical Sensitivity Gramian
        pr, pm = pscales(pr, nf[8], C)
        WC = emgr(f, g, (M, N, Q), t, "c", pr, nf, ut, us, xs, um, xm, dp)
        if not nf[9]:  # Input-state sensitivity gramian
            def DP(x, y):
                return np.sum(x.dot(y))                # Trace pseudo-kernel
        else:          # Input-output sensitivity gramian
            def DP(x, y):
                return np.sum(np.reshape(y, (Q, -1)))  # Custom pseudo-kernel
            Y = emgr(f, g, (M, N, Q), t, "o", pr, nf, ut, us, xs, um, xm, DP)

            def DP(x, y):
                return np.fabs(np.sum(y * Y))          # Custom pseudo-kernel
        WS = emgr(f, g, (M, N, Q, P), t, "c", pr, nf, ut, us, xs, pm, xm, DP)
        return WC, WS

    elif w == "i":  # Empirical Augmented Observability Gramian

        pr, pm = pscales(pr, nf[8], D)
        V = emgr(f, g, (M, N, Q, P), t, "o", pr, nf, ut, us, xs, um, np.vstack((xm, pm)), dp)
        WO = V[0:N, 0:N]        # Observability gramian
        WM = V[0:N, N:N + P]
        WI = V[N:N + P, N:N + P]  # Identifiability gramian
        if not nf[9]:
            WI -= WM.T.dot(ainv(WO)).dot(WM)
        return WO, WI

    elif w == "j":  # Empirical Joint Gramian

        pr, pm = pscales(pr, nf[8], D)
        V = emgr(f, g, (M, N, Q, P), t, "x", pr, nf, ut, us, xs, um, np.vstack((xm, pm)), dp)
        if nf[10]:
            return V         # Joint gramian partition
        WX = V[0:N, 0:N]     # Cross gramian
        WM = V[0:N, N:N + P]
        if not nf[9]:        # Cross-identifiability gramian
            WI = -0.5 * WM.T.dot(ainv(WX + WX.T)).dot(WM)
        else:
            WI = -0.5 * WM.T.dot(WM)
        return WX, WI

    else:
        assert False, "emgr: unknown gramian type!"


def scales(nf1, nf2):
    """ Input and initial state perturbation scales """

    if nf1 == 1:    # Linear
        s = np.array([0.25, 0.50, 0.75, 1.0], ndmin=1)

    elif nf1 == 2:  # Geometric
        s = np.array([0.125, 0.25, 0.5, 1.0], ndmin=1)

    elif nf1 == 3:  # Logarithmic
        s = np.array([0.001, 0.01, 0.1, 1.0], ndmin=1)

    elif nf1 == 4:  # Sparse
        s = np.array([0.01, 0.50, 0.99, 1.0], ndmin=1)

    else:
        s = np.array([1.0], ndmin=1)

    if nf2 == 0:
        s = np.concatenate((-s, s))

    return s


def pscales(p, d, c):
    """ Parameter perturbation scales """

    assert p.shape[1] >= 2, "emgr: min and max parameter requires!"

    pmin = np.amin(p, axis=1)
    pmax = np.amax(p, axis=1)

    if d == 1:    # Linear centering and scales
        pr = 0.5 * (pmax + pmin)
        pm = np.outer(pmax - pmin, np.linspace(0, 1.0, c)) + (pmin - pr)[:, np.newaxis]

    elif d == 2:  # Logarithmic centering and scales
        lmin = np.log(pmin)
        lmax = np.log(pmax)
        pr = np.real(np.exp(0.5 * (lmax + lmin)))
        pm = np.real(np.exp(np.outer(lmax - lmin, np.linspace(0, 1.0, c)) + lmin[:, np.newaxis])) - pr[:, np.newaxis]

    else:         # No centering and linear scales
        pr = np.reshape(pmin, (pmin.size, 1))
        pm = np.outer(pmax - pmin, np.linspace(1.0 / c, 1.0, c))

    return pr, pm


def ident(x, u, p, t):
    """ Output identity function """

    return x


def avg(z, nf, zs):
    """ State and output trajectory centering """

    if nf == 1:    # Steady state / output
        a = zs

    elif nf == 2:  # Final state / output
        a = z[:, -1]

    elif nf == 3:  # Temporal mean state / output
        a = np.mean(z, 1)

    elif nf == 4:  # Temporal root-mean-square state / output
        a = np.sqrt(np.mean(z * z, 1))

    elif nf == 5:  # Midrange state / output
        a = 0.5 * (z.max(1) - z.min(1))

    elif nf == 6:  # Midrange state / output
        a = np.prod(np.sign(z), 1) * np.prod(np.fabs(z), 1) ** (1.0 / float(z.shape[1]))

    else:          # None
        a = np.zeros(z.shape[0])

    return a[:, np.newaxis]


def ainv(m):
    """ Quadratic complexity approximate inverse matrix """

    d = np.copy(np.diag(m))
    k = np.nonzero(np.fabs(d) > np.sqrt(np.spacing(1)))
    d[k] = 1.0 / d[k]
    x = m * (-d)
    x *= d.T
    x.flat[::np.size(d) + 1] = d
    return x


def ssp2(f, g, t, x0, u, p):
    """ Low-Storage Strong-Stability-Preserving Second-Order Runge-Kutta """

    dt = t[0]
    nt = int(math.floor(t[1] / dt) + 1)

    y0 = g(x0, u(0), p, 0)
    Q = y0.shape[0]        # Q = N when g = ident
    y = np.zeros((Q, nt))  # Pre-allocate trajectory
    y[:, 0] = y0

    xk1 = np.copy(x0)
    xk2 = np.copy(x0)
    for k in range(1, nt):
        tk = (k - 0.5) * dt
        uk = u(tk)
        for _ in range(STAGES - 1):
            xk1 += (dt / (STAGES - 1.0)) * f(xk1, uk, p, tk)

        xk2 += dt * f(xk1, uk, p, tk)
        xk2 /= STAGES
        xk2 += xk1 * ((STAGES - 1.0) / STAGES)
        xk1 = np.copy(xk2)
        y[:, k] = g(xk1, uk, p, tk).flatten(1)

    return y
