import math
import numpy as np
import scipy as sp


def emgr(f, g=None, s=None, t=None, w=None, pr=0, nf=[0], ut=1, us=0.0, xs=0.0, um=1.0, xm=1.0, dp=np.dot):
## emgr - EMpirical GRamian Framework
#
#  project: emgr ( https://gramian.de )
#  version: 5.6-py ( 2019-01-02 )
#  authors: Christian Himpe ( 0000-0003-2194-6754 )
#  license: BSD-2-Clause License ( opensource.org/licenses/BSD-2-Clause )
#  summary: Empirical Gramians for (nonlinear) input-output systems.
#
# USAGE:
#
#  W = emgr(f,g,s,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm],[dp])
#
# DESCRIPTION:
#
#  Empirical gramian matrix and empirical covariance matrix computation
#  for model reduction, decentralized control, nonlinearity quantification,
#  sensitivity analysis, parameter identification, uncertainty quantification &
#  combined state and parameter reduction of large-scale input-output systems.
#  Data-driven analysis of input-output coherence and system-gramian-based
#  nonlinear model order reduction. Compatible with OCTAVE and MATLAB.
#
# ARGUMENTS:
#
#   f {function} vector field handle: x' = f(x,u,p,t)
#   g {function} output function handle: y = g(x,u,p,t)
#   s {tuple} system dimensions: [inputs,states,outputs]
#   t {tuple} time discretization: [time-step,time-horizon]
#   w {string} single character encoding gramian type:
#    * "c" empirical controllability gramian (Wc)
#    * "o" empirical observability gramian (Wo)
#    * "x" empirical cross gramian (Wx aka Wco or Xcg)
#    * "y" empirical linear cross gramian (Wy)
#    * "s" empirical sensitivity gramian (Ws)
#    * "i" empirical identifiability gramian (Wi)
#    * "j" empirical joint gramian (Wj)
#  pr {matrix|0} parameters, each column is a set
#  nf {tuple|0} option flags, twelve components, default zero:
#    * center: no(0), steady(1), last(2), mean(3), rms(4), midrange(5)
#    * input scales: single(0), linear(1), geometric(2), log(3), sparse(4)
#    * state scales: single(0), linear(1), geometric(2), log(3), sparse(4)
#    * input rotations: unit(0), single(1)
#    * state rotations: unit(0), single(1)
#    * normalization (only: Wc, Wo, Wx, Wy): none(0), Jacobi(1), steady(2)
#    * cross gramian type (only: Wx, Wy, Wj): regular(0), non-symmetric(1)
#    * extra input (only: Wo, Wx, Ws, Wi, Wj): none(0), yes(1)
#    * parameter centering (only: Ws, Wi, Wj): none(0), linear(1), log(2)
#    * parameter gramian variant:
#      * Averaging type (only: Ws): input-state(0), input-output(1)
#      * Schur-complement (only: Wi, Wj): detailed(0), approximate(1)
#    * cross gramian partition size (only: Wx, Wj): full(0), partitioned(<N)
#    * cross gramian partition index (only: Wx, Wj): partition(>0)
#  ut {function|1} input function handle: u_t = ut(t), default: impulse(1)
#  us {vector|0} steady-state input
#  xs {vector|0} steady-state and initial state x_0
#  um {matrix|1} input scales
#  xm {matrix|1} initial-state scales
#  dp {function|np.dot} custom inner product handle: z = dp(x,y)
#
# RETURNS:
#
#  W {matrix} Gramian Matrix (for: Wc, Wo, Wx, Wy)
#  W {tuple}  [State-, Parameter-] Gramian (for: Ws, Wi, Wj)
#
# CITATION:
#
#  C. Himpe (2019). emgr - EMpirical GRamian Framework (Version 5.6)
#  [Software]. Available from https://gramian.de . doi:10.5281/zenodo.2530021
#
# SEE ALSO:
#
#  gram
#
# KEYWORDS:
#
#  model reduction, system gramians, empirical gramians, cross gramian, MOR
#
# Further information: https://gramian.de

    # Integrator Handle
    global ODE

    if "ODE" not in globals():
        ODE = ssp2

    # Version Info
    if f is "version":
        return 5.6

    # Default Arguments
    if type(pr) in {int, float} or np.ndim(pr) == 1:
        pr = np.reshape(pr, (-1, 1))

## GENERAL SETUP

    # System Dimensions
    M = int(s[0])                        # Number of inputs
    N = int(s[1])                        # Number of states
    Q = int(s[2])                        # Number of outputs
    A = int(s[3]) if len(s) == 4 else 0  # Number of augmented parameter-states
    P = pr.shape[0]                      # Dimension of parameter
    K = pr.shape[1]                      # Number of parameter-sets
    h = t[0]                             # Time-step width
    L = int(math.floor(t[1]/h) + 1)      # Number of time-steps

    # Lazy Output Functional
    if type(g) == int and g == 1:
        g = ident
        Q = N

    # Ensure lower case gramian type
    w = w.lower()

    # Ensure flag vector length
    if len(nf) < 12:
        nf = nf + [0] * (12 - len(nf))

    # Built-in input functions
    if type(ut) in {int, float}:

        if ut == 0:  # Pseudo-Random Binary Input
            def ut(t):
                return np.random.randint(0, 1, size=M)

        elif ut == np.inf:  # Decaying Exponential Chirp Input
            mh = 0.5 * np.ones(M)
            gr = (10.0/L)**(1.0/(L*h))
            st = 2.0 * math.pi * (0.1/h) / math.log(gr)

            def ut(t):
                return mh * math.cos(st * (gr**t - 1.0)) + 0.5

        else:  # Delta Impulse Input
            mh = np.ones(M) / h

            def ut(t):
                return mh * (t <= h)

    # Lazy Optional Arguments
    if type(us) in {int, float}: us = np.full(M, us)
    if type(xs) in {int, float}: xs = np.full(N, xs)
    if type(um) in {int, float}: um = np.full(M, um)
    if type(xm) in {int, float}: xm = np.full(M if w == "y" else N, xm)

    if um.ndim == 1: um = scales(um, nf[1], nf[3])
    if xm.ndim == 1: xm = scales(xm, nf[2-(w == "y")], nf[4-(w == "y")])

    C = um.shape[1]  # Number of input scales sets
    D = xm.shape[1]  # Number of state scales sets

## GRAMIAN SETUP

    # Gramian Normalization
    if (w == "c" or "o" or "x" or "y") and nf[5] and A == 0:

        TX = np.ones(N)
        if nf[5] == 1:    # Jacobi-type preconditioner
            def DP(x, y):
                return np.sum(x * y.T, 1)  # Diagonal-only kernel
            WT = emgr(f, g, s, t, w, pr, nf[0:4]+[0]+nf[6:11], ut, us, xs, um, xm, DP)
            TX = np.sqrt(np.fabs(WT))

        elif nf[5] == 2:  # Steady-state preconditioner
            TX[xs != 0] = xs[xs != 0]

        def deco(f, g):
            def F(x, u, p, t):
                return f(TX * x, u, p, t) / TX

            def G(x, u, p, t):
                return g(TX * x, u, p, t)

            return F, G

        f, g = deco(f, g)

        xs = xs / TX

    # Extra Input
    if nf[7]:
        def up(t):
            return us + ut(t)
    else:
        def up(t):
            return us

## GRAMIAN COMPUTATION

    W = 0.0  # Reserve gramian variable

    # General layout:
    #   Parameter gramians call state gramians
    #   For each {parameter, scale, input/state/parameter}:
    #     Perturb, simulate, center, normalize, accumulate
    #   Assemble, normalize, post-process

    if w == "c":    ## Empirical Controllability Gramian

        for k in range(K):
            pk = pr[:, k]
            for c in range(C):
                for m in np.nditer(np.nonzero(um[:, c])):
                    em = np.zeros(M+P)
                    em[m+(A > 0)*M] = um[m, c]

                    def uu(t):
                        return up(t) + ut(t) * em[0:M]
                    pp = pk + em[M:M+P]
                    x = ODE(f, ident, t, xs, uu, pp)
                    x -= avg(x, nf[0], xs, t)[:, np.newaxis]
                    x /= um[m, c]
                    if A > 0:
                        W += em[M:M+P] * dp(x, x.T)
                    else:
                        W += dp(x, x.T)
        W *= h/(C*K)
        return W

    elif w == "o":  ## Empirical Observability Gramian

            o = np.zeros((Q*L, N+A))  # Pre-allocate observability matrix
            for k in range(K):
                pk = pr[:, k]
                for d in range(D):
                    for n in np.nditer(np.nonzero(xm[:, d])):
                        en = np.zeros(N+P)
                        en[n] = xm[n, d]
                        xx = xs + en[0:N]
                        pp = pk + en[N:N+P]
                        y = ODE(f, g, t, xx, up, pp)
                        y -= avg(y, nf[0], g(xs, us, pp, 0), t)[:, np.newaxis]
                        y /= xm[n, d]
                        o[:, n] = y.flatten(1)
                    W += dp(o.T, o)
            W *= h/(D*K)
            return W

    elif w == "x":  ## Empirical Cross Gramian

        assert M == Q or nf[6], "emgr: non-square system!"

        i0 = 0
        i1 = N + A

        # Partitioned cross gramian
        if nf[10] > 0:
            sp = int(round(nf[10]))    # Partition size
            ip = int(round(nf[11]))    # Partition index
            i0 += ip * sp         # Start index
            i1 = min(i0 + sp, N)  # End index
            if i0 > N:
                i0 -= math.ceil(N / sp) * sp - N
                i1 = min(i0 + sp, N+A)

            if ip < 0 or i0 >= i1 or i0 < 0:
                return 0

        o = np.zeros((Q, L, i1-i0))  # Pre-allocate observability 3-tensor
        for k in range(K):
            pk = pr[:, k]
            for d in range(D):
                for n in np.nditer(np.nonzero(xm[i0:i1, d])):
                    en = np.zeros(N+P)
                    en[i0+n] = xm[i0+n, d]
                    xx = xs + en[0:N]
                    pp = pk + en[N:N+P]
                    y = ODE(f, g, t, xx, up, pp)
                    y -= avg(y, nf[0], g(xs, us, pp, 0), t)[:, np.newaxis]
                    y /= xm[i0+n, d]
                    o[:, :, n] = y
                if nf[6]:  # Non-symmetric cross gramian: cache average
                    o[0, :, :] = np.sum(o, axis=0)
                for c in range(C):
                    for m in np.nditer(np.nonzero(um[:, c])):
                        em = np.zeros(M)
                        em[m] = um[m, c]

                        def uu(t):
                            return us + ut(t) * em
                        x = ODE(f, ident, t, xs, uu, pk)
                        x -= avg(x, nf[0], xs, t)[:, np.newaxis]
                        x /= um[m, c]
                        if nf[6]:  # Non-symmetric cross gramian
                            W += dp(x, o[0, :, :])
                        else:      # Regular cross gramian
                            W += dp(x, o[m, :, :])
        W *= h/(C*D*K)
        return W

    elif w == "y":  ## Empirical Linear cross Gramian

        assert M == Q or nf[6], "emgr: non-square system!"
        assert C == D, "emgr: scale count mismatch!"

        a = np.zeros((Q, L, N))  # Pre-allocate adjoint 3-tensor
        for k in range(K):
            pk = pr[:, k]
            for c in range(C):
                for q in np.nditer(np.nonzero(xm[:, c])):
                    em = np.zeros(Q)
                    em[q] = xm[q, c]

                    def uu(t):
                        return us + ut(t) * em
                    z = ODE(g, ident, t, xs, uu, pk)
                    z -= avg(z, nf[0], xs, t)[:, np.newaxis]
                    z /= xm[q, c]
                    a[q, :, :] = z.T
                if nf[6]:  # Non-symmetric cross gramian: cache average
                    a[0, :, :] = np.sum(a, axis=0)
                for m in np.nditer(np.nonzero(um[:, c])):
                    em = np.zeros(M)
                    em[m] = um[m, c]

                    def uu(t):
                        return us + ut(t) * em
                    x = ODE(f, ident, t, xs, uu, pk)
                    x -= avg(x, nf[0], xs, t)[:, np.newaxis]
                    x /= um[m, c]
                    if nf[6]:  # Non-symmetric cross gramian
                        W += dp(x, a[0, :, :])
                    else:      # Regular cross gramian
                        W += dp(x, a[m, :, :])
        W *= h/(C*K)
        return W

    elif w == "s":  ## Empirical Sensitivity Gramian
        pr, pm = pscales(pr, nf[8], C)
        WC = emgr(f, g, (M, N, Q), t, "c", pr, nf, ut, us, xs, um, xm, dp)
        if not nf[9]:  # Input-state sensitivity gramian
            def DP(x, y):
                return np.sum(x.dot(y))  # Trace pseudo-kernel
        else:          # Input-output sensitivity gramian
            av = np.kron(np.eye(L), np.zeros(1, Q))

            def DP(x, y):
                return av.dot(y)  # Custom pseudo-kernel
            V = emgr(f, g, (M, N, Q), t, "o", pr, nf, ut, us, xs, um, xm, DP).T

            def DP(x, y):
                return np.fabs(np.sum(x.dot(V)))  # Custom pseudo-kernel
        WS = emgr(f, g, (M, N, Q, P), t, "c", pr, nf, ut, us, xs, pm, xm, DP)
        WS /= np.max(WS)
        return WC, WS

    elif w == "i":  ## Empirical Augmented Observability Gramian

        pr, pm = pscales(pr, nf[8], D)
        V = emgr(f, g, (M, N, Q, P), t, "o", pr, nf, ut, us, xs, um, np.vstack((xm, pm)), dp)
        WO = V[0:N, 0:N]      # Observability gramian
        WM = V[0:N, N:N+P]
        WI = V[N:N+P, N:N+P]  # Identifiability gramian
        if not nf[9]:
            WI -= WM.T.dot(ainv(WO)).dot(WM)
        return WO, WI

    elif w == "j":  ## Empirical Joint Gramian

        pr, pm = pscales(pr, nf[8], D)
        V = emgr(f, g, (M, N, Q, P), t, "x", pr, nf, ut, us, xs, um, np.vstack((xm, pm)), dp)
        if nf[10]:
            return V         # Joint gramian partition
        WX = V[0:N, 0:N]     # Cross gramian
        WM = V[0:N, N:N+P]
        if not nf[9]:        # Cross-identifiability gramian
            WI = -0.5 * WM.T.dot(ainv(WX + WX.T)).dot(WM)
        else:
            WI = -0.5 * WM.T.dot(WM)
        return WX, WI

    else:
        assert False, "emgr: unknown gramian type!"


## LOCAL FUNCTION: scales
def scales(s, d, c):
#  summary: Input and initial state perturbation scales

    if d == 1:    # Linear
        sc = np.array([0.25, 0.50, 0.75, 1.0], ndmin=1)

    elif d == 2:  # Geometric
        sc = np.array([0.125, 0.25, 0.5, 1.0], ndmin=1)

    elif d == 3:  # Logarithmic
        sc = np.array([0.001, 0.01, 0.1, 1.0], ndmin=1)

    elif d == 4:  # Sparse
        sc = np.array([0.01, 0.50, 0.99, 1.0], ndmin=1)

    else:
        sc = np.array([1.0], ndmin=1)

    if c == 0:
        sc = np.concatenate((-sc, sc))

    return np.outer(s, sc)


## LOCAL FUNCTION: pscales
def pscales(p, d, c):
#  summary: Parameter perturbation scales

    assert p.shape[1] >= 2, "emgr: min & max parameter requires!"

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
        pm = np.outer(pmax - pmin, np.linspace(1.0/c, 1.0, c))

    return pr, pm


## LOCAL FUNCTION: id
def ident(x, u, p, t):
#  summary: Output identity function

    return x


## LOCAL FUNCTION: avg
def avg(x, d, c, t):
#  summary: State and output trajectory centering

    if d == 1:    # Steady state / output
        return c

    elif d == 2:  # Final state / output
        return x[:, -1]

    elif d == 3:  # Temporal mean state / output
        return np.sum(np.fabs(x), 1) * (t[0]/t[1])

    elif d == 4:  # Temporal root-mean-square state / output
        return np.sqrt(np.sum(x*x, 1)) * (t[0]/t[1])

    elif d == 5:  # Midrange state / output
        return 0.5 * (x.max(1) - x.min(1))

    else:         # None
        return np.zeros(x.shape[0])


## LOCAL FUNCTION: ainv
def ainv(m):
#  summary: Quadratic complexity approximate inverse matrix

    d = np.copy(np.diag(m))
    k = np.nonzero(np.fabs(d) > np.sqrt(np.spacing(1)))
    d[k] = 1.0/d[k]
    x = m * (-d)
    x *= d.T
    x.flat[::np.size(d)+1] = d
    return x


## LOCAL FUNCTION: ssp2
def ssp2(f, g, t, x0, u, p):
#  summary: Low-Storage Strong-Stability-Preserving Second-Order Runge-Kutta

    global STAGES  # Configurable number of stages for increased stability

    if "STAGES" not in globals():
        STAGES = 3

    h = t[0]
    K = int(math.floor(t[1]/h) + 1)

    y0 = g(x0, u(0), p, 0)
    y = np.pad(y0[:, np.newaxis], ((0, 0), (0, K-1)), mode="constant")  # Pre-allocate trajectory

    xk1 = np.copy(x0)
    xk2 = np.copy(x0)

    for k in range(1, K):
        tk = (k - 0.5) * h
        uk = u(tk)
        for s in range(STAGES-1):
            xk1 += (h/(STAGES-1.0)) * f(xk1, uk, p, tk)

        xk2 += h * f(xk1, uk, p, tk)
        xk2 /= STAGES
        xk2 += xk1 * ((STAGES-1.0)/STAGES)
        xk1 = np.copy(xk2)
        y[:, k] = g(xk1, uk, p, tk).flatten(1)

    return y
