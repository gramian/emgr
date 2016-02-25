![emgr Logo](emgr.png) emgr
===========================
Empirical Gramian Framework (Version 3.9) ![DOI badge](https://zenodo.org/badge/doi/10.5281/zenodo.46523.png)

Empirical gramians can be computed for linear and nonlinear state-space control systems for purposes of model order reduction (MOR), system identification (SYSID) and uncertainty quantification (UQ).
Model reduction using empirical gramians can be applied to the state-space, to the parameter-space or to both through combined reduction.
For state reduction, balanced truncation of the empirical controllability gramian and the empirical observability gramian, or alternatively, direct truncation (approximate balancing) of the empirical cross gramian or the empirical linear cross gramian for large-scale linear systems, is available.
For parameter reduction, parameter identification and sensitivity analysis the empirical sensitivity gramian (controllability of parameters) or the empirical identifiability gramian (observability of parameters) are provided.
Combined state and parameter reduction is enabled by the empirical joint gramian, which computes controllability and observability of states (cross gramian) and observability of parameters (cross-identifiability gramian) concurrently.
The empirical gramian framework is a compact open-source toolbox for (empirical) GRAMIAN-based model reduction and compatible with OCTAVE and MATLAB.
emgr provides a common interface for the computation of empirical gramians and empirical covariance matrices.

More information at: http://gramian.de

Code Meta Information: [CODE](CODE)
