"""
emgr - EMpirical GRamian Framework
==================================

  project: emgr ( https://gramian.de )
  version: 5.99.py (2022-04-13)
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

BRIEF:
------

  Unsupervised learning of I/O system properties for data-driven control.

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
  nf {tuple|0} option flags, thirteen component vector, default all zero:
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
    * parameter centering (only: Ws, Wi, Wj): none(0), lin(1), log(2), nom(3)
    * parameter gramian variant:
      * averaging type (only: Ws): input-state(0), input-output(1)
      * Schur-complement (only: Wi, Wj): approx(0), coarse(1)
    * cross gramian partition size (only: Wx, Wj): full(0), partitioned(<N)
    * cross gramian partition index (only: Wx, Wj): partition(>0)
    * weighting: none(0), linear(1), squared(2), state(3), scale(4), rsqrt(5)
  ut {handle|'i'} input function: u_t = ut(t) or single character string:
    * "i" delta impulse input
    * "s" step input / load vector / source term
    * "h" havercosine decaying exponential chirp input
    * "a" sinc (cardinal sine) input
    * "r" pseudo-random binary input
  us {vector|0} steady-state input (1 or #inputs rows)
  xs {vector|0} steady-state and nominal initial state x_0 (1 or #states rows)
  um {matrix|1} input scales (1 or #inputs rows)
  xm {matrix|1} initial-state scales (1 or #states rows)
  dp {handle|@mtimes} inner product or kernel: xy = dp(x,y)

RETURNS:
--------

  W {matrix} State-space system Gramian Matrix (for: Wc, Wo, Wx, Wy)
  W {tuple}  (State, Parameter)-space system Gramian (for: Ws, Wi, Wj)

CITE AS:
--------

  C. Himpe (2022). emgr - EMpirical GRamian Framework (Version 5.99)
  [Software]. Available from https://gramian.de . doi:10.5281/zenodo.6457616

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

__version__ = "5.99"
__date__ = "2022-04-13"
__copyright__ = "Copyright (C) 2022 Christian Himpe"
__author__ = "Christian Himpe"
__license__ = "BSD 2-Clause"


ODE = lambda f, g, t, x0, u, p: ssp2(f, g, t, x0, u, p)  # Preset default integrator


def emgr(f, g=None, s=None, t=None, w=None, pr=0, nf=0, ut="i", us=0.0, xs=0.0, um=1.0, xm=1.0, dp=np.dot):
    """ Compute empirical system Gramian matrix """

    # Version Info
    if f == "version":
        return __version__

    fState = f

    # Default Arguments
    if isinstance(pr, (int, float)) or np.ndim(pr) == 1:
        pr = np.reshape(pr, (-1, 1))

###############################################################################
# SETUP
###############################################################################

    # System Dimensions
    nInputs = int(s[0])   # Number of inputs
    nStates = int(s[1])   # Number of states
    nOutputs = int(s[2])  # Number of outputs

    # Parameter Dimensions
    nParams = pr.shape[0]        # Dimension of parameter
    nParamSamples = pr.shape[1]  # Number of parameter-sets

    # Time Discretization
    tStep = t[0]                           # Time-step width
    tFinal = t[1]                          # Time horizon
    nSteps = int(math.floor(tFinal / tStep) + 1)  # Number of time-steps

    # Gramian Type
    gramianType = w[0].lower()

    # Flag Vector
    if nf == 0:
        flags = [0]
    else:
        flags = nf

    if len(flags) < 13:
        flags = flags + [0] * (13 - len(flags))

    # Built-in Input Functions
    if isinstance(ut, str):

        a0 = math.pi / (2.0 * tStep) * tFinal / math.log(4.0 * (tStep / tFinal))
        b0 = (4.0 * (tStep / tFinal)) ** (1.0 / tFinal)

        if ut.lower() == "i":    # Delta Impulse Input
            def fExcite(t):
                return float(t <= tStep) / tStep

        elif ut.lower() == "s":  # Step Input
            def fExcite(t):
                return 1.0

        elif ut.lower() == "h":  # Havercosine Chirp Input
            def fExcite(t):
                return 0.5 * math.cos(a0 * (b0 ** t - 1.0)) + 0.5

        elif ut.lower() == "a":  # Sinc Input
            def fExcite(t):
                return math.sin(t / tStep) / ((t / tStep) + float(t == 0))

        elif ut.lower() == "r":  # Pseudo-Random Binary Input
            def fExcite(t):
                return np.random.randint(0, 1, 1)

        else:
            assert False, "emgr: unknown input ut!"
    else:
        fExcite = ut

###############################################################################
# CONFIGURATION
###############################################################################

    # Output Function
    if (isinstance(g, int) and g == 1) or ((gramianType == "c") and not flags[6]) or (gramianType == "y"):
        fOutput = ident
        fAdjoint = g
    else:
        fOutput = g

    # Trajectory Weighting
    tInstances = np.linspace(0, tFinal, nSteps)[np.newaxis, :]
    tInstances[0, 0] = 0.5 * tStep
    if flags[12] == 1:    # Linear Time-Weighting
        def fWeight(traj):
            return traj * np.sqrt(tInstances)

    elif flags[12] == 2:  # Quadratic Time-Weighting
        def fWeight(traj):
            return traj * (tInstances / math.sqrt(2.0))

    elif flags[12] == 3:  # State-Weighting
        def fWeight(traj):
            return traj / np.maximum(math.sqrt(np.spacing(1)), np.linalg.norm(traj, 2, axis=0))

    elif flags[12] == 4:  # Scale-Weighting
        def fWeight(traj):
            return traj / np.maximum(math.sqrt(np.spacing(1)), np.linalg.norm(traj, np.inf, axis=1)[:, np.newaxis])

    elif flags[12] == 5:  # Reciprocal Square-Root Time-Weighting
        def fWeight(traj):
            return traj / (math.pi * tInstances) ** 0.25

    else:              # None
        def fWeight(traj):
            return traj

    # Trajectory Centering
    if flags[0] == 1:    # Steady-State / Output
        def fCenter(traj, xs):
            return traj - xs.reshape(-1, 1)

    elif flags[0] == 2:  # Final State / Output
        def fCenter(traj, xs):
            return traj - traj[:, -1].reshape(-1, 1)

    elif flags[0] == 3:  # Temporal Mean State / Output
        def fCenter(traj, xs):
            return traj - np.mean(traj, axis=1).reshape(-1, 1)

    elif flags[0] == 4:  # Temporal Root-Mean-Square / Output
        def fCenter(traj, xs):
            return traj - np.sqrt(np.mean(traj * traj, axis=1)).reshape(-1, 1)

    elif flags[0] == 5:  # Temporal Mid-range of State / Output
        def fCenter(traj, xs):
            return traj - 0.5 * (np.amax(traj, axis=1) + np.amin(traj, axis=1)).reshape(-1, 1)

    else:             # None
        def fCenter(traj, xs):
            return traj

    # Steady State
    vSteadyInput = np.full((nInputs, 1), us) if np.isscalar(us) else us
    vSteadyState = np.full((nStates, 1), xs) if np.isscalar(xs) else xs

    # Gramian Normalization
    if flags[5] in {1, 2} and w in {"c", "o", "x", "y"}:

        if flags[5] == 2:  # Jacobi-type preconditioner
            NF = list(flags)
            NF[5] = 0
            if w == "c":
                NF[6] = 0

            def DP(x, y):
                return np.sum(x[:nStates, :] * y[:, :nStates].T, 1)  # Diagonal-only kernel

            TX = np.sqrt(np.abs(emgr(f, g, s, t, w, np.mean(pr, axis=1), NF, ut, us, xs, um, xm, DP)))[:, np.newaxis]

        if flags[5] == 1:  # Steady-state preconditioner
            TX = vSteadyState

        TX[np.fabs(TX) < np.sqrt(np.spacing(1))] = 1.0

        vSteadyState = vSteadyState / TX

        def fState(x, u, p, t):
            return f(TX * x, u, p, t) / TX

        def fAdjoint(x, u, p, t):
            return g(TX * x, u, p, t) / TX

        if fOutput == ident:

            def fOutput(x, u, p, t):
                return ident(TX * x, u, p, t)
        else:

            def fOutput(x, u, p, t):
                return g(TX * x, u, p, t)

    # Output Averaging
    nPages = 1 if flags[6] != 0 else nOutputs

    # Extra Input (for control explicit observability)
    if flags[7] != 0:
        def fSteady(t):
            return vSteadyInput + fExcite(t)

    else:
        def fSteady(t):
            return vSteadyInput

    # Perturbation Scales
    vInputMax = np.full((nInputs, 1), um) if np.isscalar(um) else um
    vStateMax = np.full((nStates, 1), xm) if np.isscalar(xm) else xm
    vOutputMax = np.full((nOutputs, 1), xm) if np.isscalar(xm) else xm

    mInputScales = (vInputMax * scales(flags[1], flags[3])) if vInputMax.shape[1] == 1 else vInputMax
    mStateScales = (vStateMax * scales(flags[2], flags[4])) if vStateMax.shape[1] == 1 else vStateMax
    mOutputScales = (vOutputMax * scales(flags[1], flags[3])) if vOutputMax.shape[1] == 1 else vOutputMax

    nTotalStates = mStateScales.shape[0]

    nInputScales = mInputScales.shape[1]
    nStateScales = mStateScales.shape[1]
    nOutputScales = mOutputScales.shape[1]

###############################################################################
# EMPIRICAL SYSTEM GRAMIAN COMPUTATION
###############################################################################

    W = 0.0  # Initialize gramian variable

    # Common Layout:
    #   For each {parameter, scale, input/state/parameter component}:
    #     Perturb, simulate, weight, center, normalize, accumulate
    #   Output and adjoint trajectories are cached to prevent recomputation
    #   Parameter gramians "s", "i", "j" call state gramians "c", "o", "x"

###############################################################################
# EMPIRICAL CONTROLLABILITY GRAMIAN
###############################################################################

    if w == "c":

        for k in range(nParamSamples):
            vParam = pr[:, [k]]
            vSteadyOutput = fOutput(vSteadyState, vSteadyInput, vParam, 0)
            for c in range(nInputScales):
                for m in range(nInputs):
                    sPerturb = mInputScales[m, c]
                    if sPerturb != 0.0:
                        vUnit = np.zeros(nInputs)
                        vUnit[m] = sPerturb

                        def fInput(t):
                            return vSteadyInput + vUnit * fExcite(t)

                        mTraj = fWeight(fCenter(ODE(fState, fOutput, t, vSteadyState, fInput, vParam), vSteadyOutput)) / sPerturb
                        W += dp(mTraj, mTraj.T)
        W *= tStep / (nInputScales * nParamSamples)
        return W

###############################################################################
# EMPIRICAL OBSERVABILITY GRAMIAN
###############################################################################

    if w == "o":

        obsCache = np.zeros((nPages * nSteps, nTotalStates))
        for k in range(nParamSamples):
            vParam = pr[:, [k]]
            for d in range(nStateScales):
                for n in range(nTotalStates):
                    sPerturb = mStateScales[n, d]
                    if sPerturb != 0.0:
                        vUnit = np.zeros((nTotalStates, 1))
                        vUnit[n] = sPerturb
                        vInit = vSteadyState + vUnit[:nStates]
                        vParamInit = np.copy(vParam)
                        if nTotalStates > nStates:
                            vParamInit += vUnit[nStates:]
                        vSteadyOutput = fOutput(vSteadyState, vSteadyInput, vParamInit, 0)
                        mTraj = fWeight(fCenter(ODE(fState, fOutput, t, vInit, fSteady, vParamInit), vSteadyOutput)) / sPerturb
                        if flags[6] != 0:
                            obsCache[:, n] = np.sum(mTraj, axis=0)
                        else:
                            obsCache[:, n] = mTraj.flatten(order='F')
                W += dp(obsCache.T, obsCache)
        W *= tStep / (nStateScales * nParamSamples)
        return W

###############################################################################
# EMPIRICAL CROSS GRAMIAN
###############################################################################

    if w == "x":

        assert nInputs == nOutputs or nf[6], "emgr: non-square system!"

        colFirst = 0            # Start partition column index
        colLast = nTotalStates  # Final partition column index

        # Partitioned cross gramian
        if flags[10] > 0:
            parSize = int(round(nf[10]))          # Partition size
            parIndex = int(round(nf[11]))         # Partition index
            colFirst += (parIndex - 1) * parSize  # Start index
            colLast = min(colFirst + (parSize - 1), nStates)
            if colFirst > nStates:
                colFirst -= (math.ceil(nStates / parSize) * parSize - nStates)
                colLast = min(colFirst + parSize - 1, nTotalStates)

            if parIndex < 0 or colFirst >= colLast or colFirst < 0:
                return 0

        obsCache = np.zeros((nSteps * nPages, colLast - colFirst))

        for k in range(nParamSamples):
            vParam = pr[:, [k]]
            for d in range(nStateScales):
                for n in range(colLast - colFirst):
                    sPerturb = mStateScales[colFirst + n, d]
                    if sPerturb != 0.0:
                        vUnit = np.zeros((nTotalStates, 1))
                        vUnit[colFirst + n] = sPerturb
                        vInit = vSteadyState + vUnit[:nStates]
                        vParamInit = np.copy(vParam)
                        if nTotalStates > nStates:
                            vParamInit += vUnit[nStates:]
                        vSteadyOutput = fOutput(vSteadyState, vSteadyInput, vParamInit, 0)
                        mTraj = fWeight(fCenter(ODE(fState, fOutput, t, vInit, fSteady, vParamInit), vSteadyInput)) / sPerturb
                        if flags[6] != 0:
                            obsCache[:, n] = np.sum(mTraj, axis=0).T
                        else:
                            obsCache[:, n] = mTraj.T.flatten(0)
                for c in range(nInputScales):
                    for m in range(nInputs):
                        sPerturb = mInputScales[m, c]
                        if sPerturb != 0.0:
                            vUnit = np.zeros((nInputs, 1))
                            vUnit[m] = sPerturb

                            def fInput(t):
                                return vSteadyInput + vUnit * fExcite(t)
                            mTraj = fWeight(fCenter(ODE(fState, ident, t, vSteadyState, fInput, vParam), vSteadyInput)) / sPerturb
                            nBlock = 0 if flags[6] else m * nSteps
                            W += dp(mTraj, obsCache[nBlock:nBlock + nSteps, :])
        W *= tStep / (nInputScales * nStateScales * nParamSamples)
        return W

###############################################################################
# EMPIRICAL LINEAR CROSS GRAMIAN
###############################################################################

    if w == "y":

        assert nInputs == nOutputs or nf[6], "emgr: non-square system!"
        assert nInputScales == nOutputScales, "emgr: scale count mismatch!"

        adjCache = np.zeros((nSteps, nStates, nPages))

        for k in range(nParamSamples):
            vParam = pr[:, [k]]
            for c in range(nInputScales):
                for q in range(nOutputs):
                    sPerturb = mOutputScales[q, c]
                    if sPerturb != 0.0:
                        vUnit = np.zeros((nOutputs, 1))
                        vUnit[q] = sPerturb

                        def fInput(t):
                            return vSteadyInput + vUnit * fExcite(t)
                        mTraj = fWeight(fCenter(ODE(fAdjoint, ident, t, vSteadyState, fInput, vParam), vSteadyInput)) / sPerturb
                        adjCache[:, :, q] = mTraj.T
                if flags[6] != 0:
                    adjCache[:, :, 0] = np.sum(adjCache, axis=2)
                for m in range(nInputs):
                    sPerturb = mInputScales[m, c]
                    if sPerturb != 0.0:
                        vUnit = np.zeros((nInputs, 1))
                        vUnit[m] = sPerturb

                        def fInput(t):
                            return vSteadyInput + vUnit * fExcite(t)
                        mTraj = fWeight(fCenter(ODE(fState, ident, t, vSteadyState, fInput, vParam), vSteadyInput)) / sPerturb
                        W += dp(mTraj, adjCache[:, :, 0 if flags[6] != 0 else m])
        W *= tStep / (nInputScales * nParamSamples)
        return W

###############################################################################
# EMPIRICAL SENSITIVITY GRAMIAN
###############################################################################

    if w == "s":

        # Controllability Gramian
        pr, mParamScales = paramScales(pr, flags[8], nInputScales)
        WC = emgr(f, g, s, t, "c", pr, flags, ut, us, xs, um, xm, dp)

        if not flags[9]:  # Input-state sensitivity gramian
            def DP(x, y):
                return np.sum(x * y.T)         # Trace pseudo-kernel
        else:             # Input-output sensitivity gramian
            def DP(x, y):
                return y  # Custom pseudo-kernel
            flags[6] = 1
            Y = emgr(f, g, s, t, "o", pr, flags, ut, us, xs, um, xm, DP)
            flags[6] = 0

            def DP(x, y):
                return np.abs(np.sum(y * Y))  # Custom pseudo-kernel

        # (Diagonal) Sensitivity Gramian
        WS = np.zeros((nParams, 1))

        for p in range(nParams):
            paramSamples = np.tile(pr, (1, mParamScales.shape[1]))
            paramSamples[p, :] += mParamScales[p, :]
            WS[p] = emgr(f, g, s, t, "c", paramSamples, flags, ut, us, xs, um, xm, DP)

        return WC, WS

###############################################################################
# EMPIRICAL IDENTIFIABILTY GRAMIAN
###############################################################################

    if w == "i":

        # Augmented Observability Gramian
        pr, mParamScales = paramScales(pr, flags[8], nStateScales)
        V = emgr(f, g, s, t, "o", pr, flags, ut, us, xs, um, np.vstack((mStateScales, mParamScales)), dp)

        # Return augmented observability gramian
        if flags[10] != 0:
            return V

        WO = V[:nStates, :nStates]  # Observability Gramian
        WM = V[:nStates, nStates:]  # Mixed Block
        WI = V[nStates:, nStates:]  # Parameter Gramian

        # Identifiability Gramian
        if flags[9] == 2:    # Exact Schur-complement via pseudo-inverse
            WI -= (WM.T).dot(np.linalg.pinv(WO)).dot(WM)
        elif flags[9] == 0:  # Approximate Schur-complement via approximate inverse
            WI -= (WM.T).dot(ainv(WO)).dot(WM)

        return WO, WI

###############################################################################
# EMPIRICAL JOINT GRAMIAN
###############################################################################

    if w == "j":

        # Empirical Joint Gramian
        pr, mParamScales = paramScales(pr, flags[8], nStateScales)
        V = emgr(f, g, s, t, "x", pr, flags, ut, us, xs, um, np.vstack((mStateScales, mParamScales)), dp)

        # Return joint gramian (partition)
        if flags[10] != 0:
            return V

        WX = V[:nStates, :nStates]  # Cross gramian
        WM = V[:nStates, nStates:]  # Mixed Block

        # Cross-identifiability Gramian via Schur Complement
        if flags[9] == 1:    # Coarse Schur-complement via identity
            WI = 0.5 * (WM.T).dot(WM)
        elif flags[9] == 2:  # Exact Schur-complement via pseudo-inverse
            WI = 0.5 * (WM.T).dot(np.linalg.pinv(WX + WX.T)).dot(WM)
        else:                # Approximate Schur-complement via approximate inverse
            WI = 0.5 * (WM.T).dot(ainv(WX + WX.T)).dot(WM)

        return WX, WI

    assert False, "emgr: unknown gramian type!"

###############################################################################
# LOCAL FUNCTION: ident
###############################################################################


def ident(x, u, p, t):
    """ (Output) identity function """

    return x

###############################################################################
# LOCAL FUNCTION: scales
###############################################################################


def scales(flScales, flRot):
    """ Input and initial state perturbation scales """

    if flScales == 1:    # Linear
        mScales = np.array([0.25, 0.50, 0.75, 1.0], ndmin=1)

    elif flScales == 2:  # Geometric
        mScales = np.array([0.125, 0.25, 0.5, 1.0], ndmin=1)

    elif flScales == 3:  # Logarithmic
        mScales = np.array([0.001, 0.01, 0.1, 1.0], ndmin=1)

    elif flScales == 4:  # Sparse
        mScales = np.array([0.01, 0.50, 0.99, 1.0], ndmin=1)

    else:
        mScales = np.array([1.0], ndmin=1)

    if flRot == 0:
        mScales = np.concatenate((-mScales, mScales))

    return mScales

###############################################################################
# LOCAL FUNCTION: pscales
###############################################################################


def paramScales(p, flScales, nParamScales):
    """ Parameter perturbation scales """

    vParamMin = np.amin(p, axis=1)
    vParamMax = np.amax(p, axis=1)

    if flScales == 1:    # Linear centering and scales
        assert p.shape[1] >= 2, "emgr: min and max parameter requires!"
        vParamSteady = 0.5 * (vParamMax + vParamMin)
        vScales = np.linspace(0.0, 1.0, nParamScales)

    elif flScales == 2:  # Logarithmic centering and scales
        assert p.shape[1] >= 2, "emgr: min and max parameter requires!"
        vParamSteady = np.sqrt(vParamMax * vParamMin)
        vParamMin = np.log(vParamMin)
        vParamMax = np.log(vParamMax)
        vScales = np.linspace(0.0, 1.0, nParamScales)

    elif flScales == 3:  # Nominal centering and scaling
        assert p.shape[1] >= 3, "emgr: min, nom, max parameter requires!"
        vParamSteady = p[:, 1]
        vParamMin = p[:, 0]
        vParamMax = p[:, 2]
        vScales = np.linspace(0.0, 1.0, nParamScales)

    else:                # No centering and linear scales
        assert p.shape[1] >= 2, "emgr: min and max parameter requires!"
        vParamSteady = np.copy(vParamMin)
        vParamMin = np.full(p.shape[0], 1.0 / nParamScales)
        vScales = np.linspace(1.0 / nParamScales, 1.0, nParamScales)

    mParamScales = np.outer(vParamMax - vParamMin, vScales.T) + np.expand_dims(vParamMin, -1)
    if flScales == 2:
        mParamScales = np.exp(mParamScales)
    vParamSteady = np.expand_dims(vParamSteady, -1)
    mParamScales -= vParamSteady

    return vParamSteady, mParamScales

###############################################################################
# LOCAL FUNCTION: ainv
###############################################################################


def ainv(m):
    """ Quadratic complexity approximate inverse matrix """

    # Based on truncated Neumann series: X = D^-1 - D^-1 (M - D) D^-1
    d = np.copy(np.diag(m))[:, np.newaxis]
    k = np.nonzero(np.abs(d) > np.sqrt(np.spacing(1)))
    d[k] = 1.0 / d[k]
    x = (m * (-d)) * d.T
    x.flat[::np.size(d) + 1] = d

    return x

###############################################################################
# LOCAL FUNCTION: ssp2
###############################################################################


STAGES = 3  # Configurable number of stages for enhanced stability


def ssp2(f, g, t, x0, u, p):
    """ Low-Storage Strong-Stability-Preserving Second-Order Runge-Kutta """

    nStages = STAGES if isinstance(STAGES, int) else 3

    tStep = t[0]
    nSteps = int(math.floor(t[1] / tStep) + 1)
    y0 = g(x0, u(0), p, 0)
    y = np.empty((y0.shape[0], nSteps))  # Pre-allocate trajectory
    y[:, 0] = y0.T

    xk1 = np.copy(x0)
    for k in range(1, nSteps):
        xk2 = np.copy(xk1)
        tCurr = (k - 0.5) * tStep
        uCurr = u(tCurr)
        for _ in range(1, nStages):
            xk1 += (tStep / (STAGES - 1.0)) * f(xk1, uCurr, p, tCurr)
        xk1 = (xk1 * (nStages - 1) + xk2 + tStep * f(xk1, uCurr, p, tCurr)) / nStages
        y[:, k] = g(xk1, uCurr, p, tCurr).flatten(0)

    return y
