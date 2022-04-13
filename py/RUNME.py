#!/usr/bin/env python3

"""
  project: emgr ( https://gramian.de )
  version: 5.99.py (2022-04-13)
  authors: Christian Himpe (0000-0003-2194-6754)
  license: BSD-2-Clause License (opensource.org/licenses/BSD-2-Clause)
  summary: RUNME (Minimal emgr test script)
"""

import numpy as np
from emgr import emgr


print("emgr version: {0}".format(emgr("version")))

# Linear System
A = -np.eye(4)
B = np.array([[0.0,1.0,0.0,1.0]]).T
C = np.array([[0.0,0.0,1.0,1.0]])

P = np.zeros((4,1))
Q = np.array([[0.01,0.25],[0.01,0.5],[0.01,0.75],[0.01,1.0]])

def f(x,u,p,t): return A.dot(x) + B.dot(u) + p
def g(x,u,p,t): return C.dot(x)
def h(x,u,p,t): return A.T.dot(x) + C.T.dot(u)

s = (1,4,1)
t = (0.01,1.0)

# (Empircal) Controllability Gramian
WC = emgr(f,g,s,t,"c",P)
print(WC)

# (Empirical) Observability Gramian
WO = emgr(f,g,s,t,"o",P)
print(WO)

# (Empirical) Cross Gramian
WX = emgr(f,g,s,t,"x",P)
print(WX)

# (Empirical) Linear Cross Gramian
WY = emgr(f,h,s,t,"y",P)
print(WY)

# (Empircal) Controllability Gramian
WC,WS = emgr(f,g,s,t,"s",Q)
print(WC)
print(WS)

# (Empirical) Observability Gramian
WO,WI = emgr(f,g,s,t,"i",Q)
print(WO)
print(WI)

# (Empirical) Cross Gramian
WX,WJ = emgr(f,g,s,t,"j",Q)
print(WX)
print(WJ)

