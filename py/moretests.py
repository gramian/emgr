from __future__ import print_function
import numpy as np
from scipy.linalg import toeplitz
from emgr import emgr
import matplotlib.pyplot as plt


def moretests(ut=1):
### summary: moretests
### project: emgr - EMpirical GRamian Framework ( https://gramian.de )
### authors: Christian Himpe ( 0000-0003-2194-6754 )
### license: BSD-2-Clause (2019)

    global ODE

    try:
        print("emgr (version: {0})".format(emgr("version")))
    except:
        print("emgr not found! Get emgr at: https://gramian.de")
        exit()

## SYSTEM SETUP
    M = 4      # number of inputs
    N = M * M  # number of states
    Q = M      # number of outputs
    h = 0.01   # time step size
    T = 1.0    # time horizon

    on = np.arange(1,N+1)
    om = np.arange(1,M+1)
    X = np.linspace(0,1,N).T                       # initial state
    #def U(t): np.ones((M,1)) * ((t<=h)/h)          # impulse input function
    P = 0.5 + 0.5 * np.cos(on).T                   # parameter
    R = np.ones((N,1)).dot(np.array([[0.5,1.0]]))  # parameter range

    A = np.diag(N*(np.sin(on)-1.0))    # system matrix
    B = np.sin(toeplitz(on,om)) + 1.0  # input matrix
    C = np.cos(toeplitz(om,on)) + 1.0  # output matrix

    def LIN(x,u,p,t): return A.dot(x) + B.dot(u) + p  # vector field
    def ADJ(x,u,p,t): return A.T.dot(x) + C.T.dot(u)  # adjoint vector field
    def OUT(x,u,p,t): return C.dot(x)                 # output functional

## Controllability Gramian

    print("Empirical Controllability Gramian: ", end=""),
    W = []

    # Baseline
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Centering
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Input Scales
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[0,1,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[0,2,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[0,3,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="") # fast drop
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[0,4,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Input Rotations
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[0,0,0,1,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Normalization
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[0,0,0,0,0,1,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"c",P,[0,0,0,0,0,2,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    plotsingvals(W,(2,6,1),"Controllability Gramian",8)
    print("")

## Observability Gramian

    print("Empirical Observability Gramian:   ", end="")
    W = []

    # Baseline
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Centering
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # State scales
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[0,0,1,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[0,0,2,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[0,0,3,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")  # fast drop
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[0,0,4,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # State Rotations
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[0,0,0,0,1,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Normalization
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[0,0,0,0,0,1,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[0,0,0,0,0,2,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Extra Input
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"o",P,[0,0,0,0,0,0,0,1,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    plotsingvals(W,(2,6,2),"Observability Gramian",8)
    print("")

## Linear Cross Gramian

    print("Empirical Linear Cross Gramian:    ", end="")
    W = []

    # Baseline
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Centering
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Input scales
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,1,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,2,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,3,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")  # fast drop
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,4,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Adjoint scales
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,1,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,2,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,3,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,4,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Input Rotations
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,0,1,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Adjoint Rotations
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,0,0,1,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Normalization
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,0,0,0,1,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,0,0,0,2,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    plotsingvals(W,(2,6,3),"Linear Cross Gramian",8)
    W = []

    # Non-Symmetric Cross Gramian
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Centering
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[1,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[2,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[3,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[4,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[5,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Input Scales
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,1,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,2,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,3,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")  # fast drop
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,4,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Adjoint Scales
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,1,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,2,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,3,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,4,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Input Rotations
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,0,1,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Adjoint Rotations
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,0,0,1,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Normalization
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,0,0,0,1,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,ADJ,(M,N,Q),(h,T),"y",P,[0,0,0,0,0,2,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    plotsingvals(W,(2,6,4),"Linear Cross Gramian (Non-Symmetric)",8)
    print("")

## Cross Gramian

    print("Empirical Cross Gramian:           ", end="")
    W = []

    # Baseline
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Centering
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Input scales
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,1,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,2,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,3,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,4,0,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # State scales
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,1,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,2,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,3,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,4,0,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Input Rotations
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,1,0,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # State Rotations
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,0,1,0,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")  # fast drop

    # Normalization
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,0,0,1,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,0,0,2,0,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Extra Input
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,0,0,0,0,1,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    plotsingvals(W,(2,6,5),"Cross Gramian",15)
    W = []

    # Baseline
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Centering
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[1,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[2,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[3,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[4,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[5,0,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Input scales
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,1,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,2,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,3,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,4,0,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # State scales
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,1,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,2,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,3,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,4,0,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # Input Rotations
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,1,0,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    # State Rotations
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,0,1,0,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")  # fast drop

    # Normalization
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,0,0,1,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,0,0,2,1,0,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="") 

    # Extra Input
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"x",P,[0,0,0,0,0,0,1,1,0,0,0,0],ut,0,X),compute_uv=False)); print("#", end="")

    plotsingvals(W,(2,6,6),"Cross Gramian (Non-Symmetric)",15)
    print("")

## Sensitivity Gramian

    print("Empirical Sensitivity Gramian:     ", end="")
    W = []

    # Baseline
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # Centering
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # Input scales
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,1,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,2,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,3,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,4,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # Input Rotations
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,0,1,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # Extra Input
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,0,0,0,0,0,1,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # Parameter Scaling
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,0,0,0,0,0,0,1,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")  # fast drop
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,0,0,0,0,0,0,2,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    plotsingvals(W,(2,6,9),"Sensitivity Gramian",12)
    W = []

    # Baseline
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # Centering
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")  # fast drop
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # Input scales
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,1,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,2,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,3,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,4,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # Input scales
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,1,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,2,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,3,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,4,0,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # Input Rotations
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,0,1,0,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # State Rotations
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,0,0,1,0,0,0,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # Extra Input
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,0,0,0,0,0,1,0,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")

    # Parameter Scaling
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,0,0,0,0,0,0,1,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")
    W.append(np.sort(emgr(LIN,OUT,(M,N,Q),(h,T),"s",R,[0,0,0,0,0,0,0,0,2,0,0,0],ut,0,X)[1])[::-1]); print("#", end="")  # fast drop

    plotsingvals(W,(2,6,10),"Sensitivity Gramian (Input-Output)",17)
    print("")

## Identifiability Gramian

    print("Empirical Identifiability Gramian: ", end="")
    W = []

    # Baseline
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Centering
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # State scales
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[0,0,1,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[0,0,2,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[0,0,3,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[0,0,4,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # State Rotations
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[0,0,0,0,1,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Extra Input
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[0,0,0,0,0,0,0,1,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")  # fast drop

    # Parameter Scaling
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[0,0,0,0,0,0,0,0,1,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[0,0,0,0,0,0,0,0,2,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Schur Complement
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"i",R,[0,0,0,0,0,0,0,0,0,1,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    plotsingvals(W,(2,6,8),"Identifiability Gramian",11)
    print("")

## Joint Gramian

    print("Empirical Joint Gramian:           ", end="");
    W = []

    # Baseline
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Centering
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Input scales
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,1,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,2,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,3,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,4,0,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # State scales
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,1,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,2,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,3,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,4,0,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Input Rotations
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,1,0,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # State Rotations
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,1,0,0,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Extra Input
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,0,0,0,1,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")  # fast drop

    # Parameter Scaling
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,0,0,0,0,1,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,0,0,0,0,2,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Schur Complement
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,0,0,0,0,0,1,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    plotsingvals(W,(2,6,11),"Cross-Identifiability Gramian",16)
    W = []

    # Baseline
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Centering
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[1,0,0,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[2,0,0,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[3,0,0,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[4,0,0,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[5,0,0,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Input scales
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,1,0,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,2,0,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,3,0,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,4,0,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # State scales
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,1,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,2,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,3,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,4,0,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Input Rotations
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,1,0,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # State Rotations
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,1,0,1,0,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Extra Input
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,0,0,1,1,0,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")  # fast drop

    # Parameter Scaling
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,0,0,1,0,1,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,0,0,1,0,2,0,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    # Schur Complement
    W.append(np.linalg.svd(emgr(LIN,OUT,(M,N,Q),(h,T),"j",R,[0,0,0,0,0,0,1,0,0,1,0,0],ut,0,X)[1],compute_uv=False)); print("#", end="")

    plotsingvals(W,(2,6,12),"Cross-Identifiability Gramian (Non-Symmetric)",16)
    print("")

    plt.show()

    return


def plotsingvals(W,s,l,b):

    plt.subplot(s[0],s[1],s[2])

    for k in range(1,len(W)):

        plt.semilogy(W[k]/W[k].max())

    plt.semilogy(W[0]/W[0].max(),color="r")
    plt.semilogy(W[b]/W[b].max(),color="k")

    plt.xlim(1,W[0].size)
    mi, ma = plt.ylim()
    plt.ylim(min(mi,0.1),1)
    plt.xlabel(l)

    return

moretests()
