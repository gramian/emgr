"""
emgr - EMpirical GRamian Framework
==================================

  project: emgr ( https://gramian.de )
  version: 5.9.py (2021-01-21)
  authors: Christian Himpe (0000-0003-2194-6754)
  license: BSD-2-Clause License (opensource.org/licenses/BSD-2-Clause)
  summary: Empirical system Gramians for (nonlinear) input-output systems.

DESCRIPTION:
------------

  Empirical gramian matrix and empirical covariance matrix computation
  for model reduction, decentralized control, nonlinearity quantification,
  sensitivity analysis, parameter identification, uncertainty quantification &
  combined state and parameter reduction of large-scale input-output systems.
  Data-driven analysis of input-output coherence and system-gramian-based
  nonlinear model order reduction. Compatible with PYTHON3.

ALGORITHM:
----------

  C. Himpe (2018). emgr - The Empirical Gramian Framework. Algorithms 11(7):91
  doi:10.3390/a11070091

USAGE:
------

  W = emgr(f,g,s,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm],[dp])

MANDATORY ARGUMENTS:
--------------------

   f {function} vector field handle: x' = f(x,u,p,t)
   g {function} output function handle: y = g(x,u,p,t)
   s {tuple} system dimensions: [inputs, states, outputs]
   t {tuple} time discretization: [time-step, time-horizon]
   w {string} single character encoding gramian type:
    * "c" empirical controllability gramian (Wc)
    * "o" empirical observability gramian (Wo)
    * "x" empirical cross gramian (Wx aka Wco)
    * "y" empirical linear cross gramian (Wy)
    * "s" empirical sensitivity gramian (Ws)
    * "i" empirical identifiability gramian (Wi)
    * "j" empirical joint gramian (Wj)

OPTIONAL ARGUMENTS:
-------------------

  pr {matrix|0} parameter vector(s), each column is one parameter sample
  nf {vector|0} option flags, thirteen component vector, default all zero:
    * centering: none(0), steady(1), last(2), mean(3), rms(4), midrange(5)
    * input scales: single(0), linear(1), geometric(2), log(3), sparse(4)
    * state scales: single(0), linear(1), geometric(2), log(3), sparse(4)
    * input rotations: unit(0), single(1)
    * state rotations: unit(0), single(1)
    * normalization (only: Wc, Wo, Wx, Wy): none(0), steady(1), Jacobi(2)
    * state gramian variant:
      * controllability gramian type (only: Wc, Ws): regular(0), output(1)
      * observability gramian type (only: Wo, Wi): regular(0), averaged(1)
      * cross gramian type (only: Wx, Wy, Wj): regular(0), non-symmetric(1)
    * extra input (only: Wo, Wx, Ws, Wi, Wj): no(0), yes(1)
    * parameter centering (only: Ws, Wi, Wj): none(0), linear(1), log(2)
    * parameter gramian variant:
      * averaging type (only: Ws): input-state(0), input-output(1)
      * Schur-complement (only: Wi, Wj): approx(0), coarse(1)
    * cross gramian partition size (only: Wx, Wj): full(0), partitioned(<N)
    * cross gramian partition index (only: Wx, Wj): partition(>0)
    * weighting: none(0), linear(1), squared(2), state(3), scale(4)
  ut {handle|'i'} input function: u_t = ut(t) or single character string:
    * "i" delta impulse input
    * "s" step input / load vector / source term
    * "h" havercosine decaying exponential chirp input
    * "a" sinc (cardinal sine) input
    * "r" pseudo-random binary input
  us {vector|0} steady-state input (1 or M rows)
  xs {vector|0} steady-state and nominal initial state x_0 (1 or N rows)
  um {matrix|1} input scales (1 or M rows)
  xm {matrix|1} initial-state scales (1 or N rows)
  dp {handle|@mtimes} inner product or kernel: xy = dp(x,y)

RETURNS:
--------

  W {matrix} State-space system Gramian Matrix (for: Wc, Wo, Wx, Wy)
  W {tuple}  (State, Parameter)-space system Gramian (for: Ws, Wi, Wj)

CITE AS:
--------

  C. Himpe (2021). emgr - EMpirical GRamian Framework (Version 5.9)
  [Software]. Available from https://gramian.de . doi:10.5281/zenodo.4454679

KEYWORDS:
---------

  model reduction, system gramians, empirical gramians, cross gramian, MOR

SEE ALSO:
---------

  gram (Python Control Systems Library)

COPYRIGHT: Christian Himpe
---------

For more information, see: https://gramian.de
"""

import math
import numpy as np

__version__ = "5.9"
__date__ = "2021-01-21"
__copyright__ = "Copyright (C) 2021 Christian Himpe"
__author__ = "Christian Himpe"
__license__ = "BSD 2-Clause"


ODE = lambda f, g, t, x0, u, p: ssp2(f, g, t, x0, u, p)  # Preset default integrator


def emgr(f, g=None, s=None, t=None, w=None, pr=0, nf=0, ut="i", us=0.0, xs=0.0, um=1.0, xm=1.0, dp=np.dot):
    """ Compute empirical system Gramian matrix """

    # Version Info
    if f == "version":
        return __version__

    # Default Arguments
    if isinstance(pr, (int, float)) or np.ndim(pr) == 1:
        pr = np.reshape(pr, (-1, 1))

###############################################################################
# SETUP
###############################################################################

    # System Dimensions
    M = int(s[0])                        # Number of inputs
    N = int(s[1])                        # Number of states
    Q = int(s[2])                        # Number of outputs
    P = pr.shape[0]                      # Dimension of parameter
    K = pr.shape[1]                      # Number of parameter-sets

    # Time Discretization
    dt = t[0]                            # Time-step width
    Tf = t[1]                            # Time horizon
    nt = int(math.floor(Tf / dt) + 1)    # Number of time-steps

    # Force lower-case Gramian type
    w = w[0].lower()

    # Lazy Output Functional
    if isinstance(g, int) and g == 1:
        g = ident
        Q = N

    # Pad Flag Vector
    if nf == 0:
        nf = [0]

    if len(nf) < 13:
        nf = nf + [0] * (13 - len(nf))

    # Built-in input functions
    if isinstance(ut, str):
        if ut.lower() == "s":    # Step Input
            def ut(t):
                return 1

        elif ut.lower() == "h":  # Havercosine Chirp Input
            a0 = math.pi / (2.0 * dt) * Tf / math.log(4.0 * (dt / Tf))
            b0 = (4.0 * (dt / Tf)) ** (1.0 / Tf)

            def ut(t):
                return 0.5 * math.cos(a0 * (b0 ** t - 1.0)) + 0.5

        elif ut.lower() == "a":  # Sinc Input
            def ut(t):
                return math.sin(t / dt) / ((t / dt) + float(t == 0))

        elif ut.lower() == "r":  # Pseudo-Random Binary Input
            rt = np.random.randint(0, 2, size=nt)

            def ut(t):
                return rt[int(math.floor(t / dt))]

        else:                    # Delta Impulse Input
            def ut(t):
                return float(t <= dt) / dt

    # Lazy Optional Arguments
    if isinstance(us, (int, float)): us = np.full(M, us)
    if isinstance(xs, (int, float)): xs = np.full(N, xs)
    if isinstance(um, (int, float)): um = np.full(M, um)
    if isinstance(xm, (int, float)): xm = np.full(N, xm)

###############################################################################
# CONFIGURATION
###############################################################################

    # Trajectory Weighting
    if nf[12] == 1:    # Linear Time-Weighting
        def wei(m):
            return np.sqrt(np.linspace(0, Tf, nt))

    elif nf[12] == 2:  # Quadratic Time-Weighting
        def wei(m):
            return np.linspace(0, Tf, nt) * math.sqrt(0.5)

    elif nf[12] == 3:  # State-Weighting
        def wei(m):
            return 1.0 / np.maximum(math.sqrt(np.spacing(1)), np.linalg.norm(m, 2, axis=0))

    elif nf[12] == 4:  # Scale-Weighting
        def wei(m):
            return 1.0 / np.maximum(math.sqrt(np.spacing(1)), np.linalg.norm(m, np.inf, axis=1)[:, np.newaxis])

    else:              # None
        def wei(m):
            return 1.0

    # Trajectory Centering
    if nf[0] == 1:    # Steady-State / Output
        def avg(m, s):
            return s.reshape(-1, 1)

    elif nf[0] == 2:  # Final State / Output
        def avg(m, s):
            return m[:, -1].reshape(-1, 1)

    elif nf[0] == 3:  # Temporal Mean State / Output
        def avg(m, s):
            return np.mean(m, axis=1).reshape(-1, 1)

    elif nf[0] == 4:  # Temporal Root-Mean-Square / Output
        def avg(m, s):
            return np.sqrt(np.mean(m * m, axis=1)).reshape(-1, 1)

    elif nf[0] == 5:  # Temporal Mid-range of State / Output
        def avg(m, s):
            return 0.5 * (np.amax(m, axis=1) + np.amin(m, axis=1)).reshape(-1, 1)

    else:             # None
        def avg(m, s):
            return 0.0

    # Gramian Normalization
    if nf[5] and w in {"c", "o", "x", "y"}:

        nf[5] = 0
        TX = xs         # Steady-state preconditioner

        if nf[5] == 2:  # Jacobi-type preconditioner
            NF = nf
            NF[5] = 0
            if w == "c":
                NF[6] = 0
            PR = np.mean(pr, axis=1)

            def DP(x, y):
                return np.sum(x * y.T, 1)  # Diagonal-only kernel
            TX = np.sqrt(np.fabs(emgr(f, g, s, t, w, PR, NF, ut, us, xs, um, xm, DP)))

        TX[np.fabs(TX) < np.sqrt(np.spacing(1))] = 1.0

        tx = TX if w == "y" else 1.0

        def deco(f, g):
            def F(x, u, p, t):
                return f(TX * x, u, p, t) / TX

            def G(x, u, p, t):
                return g(TX * x, u, p, t) / tx

            return F, G

        f, g = deco(f, g)

        xs = xs / TX

    # State gramian variant setup
    if nf[6]:
        G = g                      # Output Controllability Gramian
        R = 1                      # Average Observability Cache Size

        def oavg(y):               # Average Observability Gramian
            return np.sum(y, axis=1)

        S = 0                      # Non-Symmetric (Linear) Cross Gramian
    else:
        G = ident                  # Regular Controllability Gramian
        R = Q                      # Regular Observability Cache Size

        def oavg(y):
            return y.flatten(0)  # Regular Observability Gramian

        S = 1                      # Regular (Linear) Cross Gramian

    # Extra input
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

    A = xm.shape[0]  # Number of total states (regular and augmented)
    C = um.shape[1]  # Number of input scales sets
    D = xm.shape[1]  # Number of state scales sets

###############################################################################
# EMPIRICAL SYSTEM GRAMIAN COMPUTATION
###############################################################################

    W = 0.0  # Initialize gramian variable

    # Common Layout:
    #   For each {parameter, scale, input/state/parameter component}:
    #     Perturb, simulate, weight, center, normalize, accumulate
    #   Parameter gramians call state gramians

###############################################################################
# EMPIRICAL CONTROLLABILITY GRAMIAN
###############################################################################

    if w == "c":

        for k in range(K):
            for c in range(C):
                for m in range(M):
                    if um[m, c] != 0:
                        em = np.zeros(M + P)
                        em[m] = um[m, c]

                        def umc(t):
                            return up(t) + ut(t) * em[0:M]
                        pmc = pr[:, k] + em[M:M + P]
                        x = ODE(f, G, t, xs, umc, pmc)
                        x *= wei(x)
                        x -= avg(x, G(xs, us, pmc, 0))
                        x /= um[m, c]
                        W += dp(x, x.T)
        W *= dt / (C * K)
        return W

###############################################################################
# EMPIRICAL OBSERVABILITY GRAMIAN
###############################################################################

    if w == "o":

        o = np.zeros((R * nt, A))  # Pre-allocate observability matrix
        for k in range(K):
            for d in range(D):
                for n in range(A):
                    if xm[n, d] != 0:
                        en = np.zeros(N + P)
                        en[n] = xm[n, d]
                        xnd = xs + en[0:N]
                        pnd = pr[:, k] + en[N:N + P]
                        y = ODE(f, g, t, xnd, up, pnd)
                        y *= wei(y)
                        y -= avg(y, g(xs, us, pnd, 0))
                        y /= xm[n, d]
                        o[:, n] = oavg(y.T)
                W += dp(o.T, o)
            W *= dt / (D * K)
            return W

###############################################################################
# EMPIRICAL CROSS GRAMIAN
###############################################################################

    if w == "x":

        assert M == Q or nf[6], "emgr: non-square system!"

        i0 = 0
        i1 = A

        # Partitioned cross gramian
        if nf[10] > 0:
            sp = int(round(nf[10]))  # Partition size
            ip = int(round(nf[11]))  # Partition index
            i0 += ip * sp            # Start index
            i1 = min(i0 + sp, N)     # End index
            if i0 > N:
                i0 -= math.ceil(N / sp) * sp - N
                i1 = min(i0 + sp, A)

            if ip < 0 or i0 >= i1 or i0 < 0:
                return 0

        o = np.zeros((R * nt, i1 - i0))  # Pre-allocate observability cache
        for k in range(K):
            for d in range(D):
                for n in range(i1 - i0):
                    if xm[n, d] != 0:
                        en = np.zeros(N + P)
                        en[i0 + n] = xm[i0 + n, d]
                        xnd = xs + en[0:N]
                        pnd = pr[:, k] + en[N:N + P]
                        y = ODE(f, g, t, xnd, up, pnd)
                        y *= wei(y)
                        y -= avg(y, g(xs, us, pnd, 0))
                        y /= xm[i0 + n, d]
                        o[:, n] = oavg(y.T)
                for c in range(C):
                    for m in range(M):
                        if um[m, c] != 0:
                            em = np.zeros(M)
                            em[m] = um[m, c]

                            def umc(t):
                                return us + ut(t) * em
                            x = ODE(f, ident, t, xs, umc, pr[:, k])
                            x *= wei(x)
                            x -= avg(x, xs)
                            x /= um[m, c]
                            W += dp(x, o[(nt * (S * (m - 1))):(nt * (S * m) + nt), :])
        W *= dt / (C * D * K)
        return W

###############################################################################
# EMPIRICAL LINEAR CROSS GRAMIAN
###############################################################################

    if w == "y":

        assert M == Q or nf[6], "emgr: non-square system!"
        assert C == vm.shape[1], "emgr: scale count mismatch!"

        a = np.zeros((N * nt, Q))  # Pre-allocate adjoint cache
        for k in range(K):
            for c in range(C):
                for q in range(Q):
                    if vm[q, c] != 0:
                        em = np.zeros(Q)
                        em[q] = vm[q, c]

                        def vqc(t):
                            return us + ut(t) * em
                        z = ODE(g, ident, t, xs, vqc, pr[:, k])
                        z *= wei(z)
                        z -= avg(z, xs)
                        z /= vm[q, c]
                        a[:, q] = z.flatten(0)
                if nf[6]:  # Non-symmetric cross gramian
                    a[:, 0] = np.sum(a, axis=1)
                for m in range(M):
                    if um[m, c] != 0:
                        em = np.zeros(M)
                        em[m] = um[m, c]

                        def umc(t):
                            return us + ut(t) * em
                        x = ODE(f, ident, t, xs, umc, pr[:, k])
                        x *= wei(x)
                        x -= avg(x, xs)
                        x /= um[m, c]
                        W += dp(x, np.reshape(a[:, S * m], (N, nt)).T)
        W *= dt / (C * K)
        return W

###############################################################################
# EMPIRICAL SENSITIVITY GRAMIAN
###############################################################################

    if w == "s":

        # Empirical Controllability Gramian
        pr, pm = pscales(pr, nf[8], C)
        WC = emgr(f, g, s, t, "c", pr, nf, ut, us, xs, um, xm, dp)

        if not nf[9]:  # Input-state sensitivity gramian
            def DP(x, y):
                return np.sum(x * y.T)                 # Trace pseudo-kernel
        else:          # Input-output sensitivity gramian
            def DP(x, y):
                return np.sum(np.reshape(y, (R, -1)))  # Custom pseudo-kernel

            Y = emgr(f, g, s, t, "o", pr, nf, ut, us, xs, um, xm, DP)

            def DP(x, y):
                return np.fabs(np.sum(y * Y))          # Custom pseudo-kernel

        # (Diagonal) Sensitivity Gramian
        WS = np.zeros((P, P))

        for p in range(P):
            pp = np.tile(pr, (1, pm.shape[1]))
            pp[p, :] = pp[p, :] + pm[p, :]
            WS[p, p] = emgr(f, g, s, t, "c", pp, nf, ut, us, xs, um, xm, DP)

        return WC, WS

###############################################################################
# EMPIRICAL IDENTIFIABILTY GRAMIAN
###############################################################################

    if w == "i":

        # Augmented Observability Gramian
        pr, pm = pscales(pr, nf[8], D)
        V = emgr(f, g, s, t, "o", pr, nf, ut, us, xs, um, np.vstack((xm, pm)), dp)

        WO = V[0:N, 0:N]      # Observability Gramian
        WM = V[0:N, N:N + P]  # Mixed Block

        # Identifiability Gramian
        if not nf[9]:         # Schur-complement via approximate inverse
            WI = V[N:N + P, N:N + P] - WM.T.dot(ainv(WO)).dot(WM)
        else:                 # Coarse Schur-complement via zero
            WI = V[N:N + P, N:N + P]

        return WO, WI

###############################################################################
# EMPIRICAL JOINT GRAMIAN
###############################################################################

    if w == "j":

        # Empirical Joint Gramian
        pr, pm = pscales(pr, nf[8], D)
        V = emgr(f, g, s, t, "x", pr, nf, ut, us, xs, um, np.vstack((xm, pm)), dp)

        if nf[10]:
            return V   # Joint gramian partition

        WX = V[0:N, 0:N]      # Cross gramian
        WM = V[0:N, N:N + P]  # Mixed Block

        if not nf[9]:         # Cross-identifiability gramian
            WI = 0.5 * WM.T.dot(ainv(WX + WX.T)).dot(WM)
        else:                 # Coarse Schur-complement via identity
            WI = 0.5 * WM.T.dot(WM)
        return WX, WI

    assert False, "emgr: unknown gramian type!"

###############################################################################
# LOCAL FUNCTION: scales
###############################################################################


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

###############################################################################
# LOCAL FUNCTION: pscales
###############################################################################


def pscales(p, nf, ns):
    """ Parameter perturbation scales """

    assert p.shape[1] >= 2, "emgr: min and max parameter requires!"

    pmin = np.amin(p, axis=1)
    pmax = np.amax(p, axis=1)

    if nf == 1:    # Linear centering and scales
        pr = 0.5 * (pmax + pmin)
        pm = np.outer(pmax - pmin, np.linspace(0, 1.0, ns)) + (pmin - pr)[:, np.newaxis]

    elif nf == 2:  # Logarithmic centering and scales
        lmin = np.log(pmin)
        lmax = np.log(pmax)
        pr = np.real(np.exp(0.5 * (lmax + lmin)))
        pm = np.real(np.exp(np.outer(lmax - lmin, np.linspace(0, 1.0, ns)) + lmin[:, np.newaxis])) - pr[:, np.newaxis]

    else:         # No centering and linear scales
        pr = np.reshape(pmin, (pmin.size, 1))
        pm = np.outer(pmax - pmin, np.linspace(1.0 / ns, 1.0, ns))

    return pr, pm

###############################################################################
# LOCAL FUNCTION: ident
###############################################################################


def ident(x, u, p, t):
    """ (Output) identity function """

    return x

###############################################################################
# LOCAL FUNCTION: ainv
###############################################################################


def ainv(m):
    """ Quadratic complexity approximate inverse matrix """

    # Based on truncated Neumann series: X = D^-1 - D^-1 (M - D) D^-1
    d = np.copy(np.diag(m))
    k = np.nonzero(np.fabs(d) > np.sqrt(np.spacing(1)))
    d[k] = 1.0 / d[k]
    x = m * (-d)
    x *= d.T
    x.flat[::np.size(d) + 1] = d
    return x

###############################################################################
# LOCAL FUNCTION: ssp2
###############################################################################


STAGES = 3  # Configurable number of stages for increased stability of ssp2


def ssp2(f, g, t, x0, u, p):
    """ Low-Storage Strong-Stability-Preserving Second-Order Runge-Kutta """

    dt = t[0]
    nt = int(math.floor(t[1] / dt) + 1)
    y0 = g(x0, u(0), p, 0)
    y = np.zeros((y0.shape[0], nt))  # Pre-allocate trajectory
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
        y[:, k] = g(xk1, uk, p, tk).flatten(0)

    return y
