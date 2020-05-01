"""
  project: emgr ( https://gramian.de )
  version: 5.8.py (2020-05-01)
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

def f(x,u,p,t): return A.dot(x) + B.dot(u)
def g(x,u,p,t): return C.dot(x)
def h(x,u,p,t): return A.T.dot(x) + C.T.dot(u)

# (Empircal) Controllability Gramian
WC = emgr(f,g,(1,4,1),(0.01,1.0),"c")
print(WC)

# (Empirical) Observability Gramian
WO = emgr(f,g,(1,4,1),(0.01,1.0),"o")
print(WO)

# (Empirical) Cross Gramian
WX = emgr(f,g,(1,4,1),(0.01,1.0),"x")
print(WX)

# (Empirical) Linear Cross Gramian
WY = emgr(f,h,(1,4,1),(0.01,1.0),"y")
print(WY)

Q = np.array([[0.0,0.25],[0.0,0.5],[0.0,0.75],[0.0,1.0]])

def f(x,u,p,t): return A.dot(x) + B.dot(u) + p

# (Empircal) Controllability Gramian
WC,WS = emgr(f,g,(1,4,1),(0.01,1.0),"s",Q)
print(WS)

# (Empirical) Observability Gramian
WO,WI = emgr(f,g,(1,4,1),(0.01,1.0),"i",Q)
print(WI)

# (Empirical) Cross Gramian
WX,WJ = emgr(f,g,(1,4,1),(0.01,1.0),"j",Q)
print(WJ)
