![code meta-data](https://img.shields.io/badge/code_meta--data-%E2%9C%93-brightgreen.svg) 
![zenodo listed](https://zenodo.org/badge/doi/10.5281/zenodo.3779889.png)
![matlab compatible](https://img.shields.io/badge/matlab-compatible-lightgrey.svg)
![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/gramian/emgr/)

![emgr logo](emgr.png) emgr - EMpirical GRamian Framework
=========================================================

[Website](https://gramian.de) |
[Twitter](https://twitter.com/modelreduction) |
[Feedback](&#x63;&#x68;&#x40;&#x67;&#x72;&#x61;&#x6D;&#x69;&#x61;&#x6E;&#x2E;&#x64;&#x65;)

* [emgr - EMpirical GRamian Framework](https://gramian.de)
* version: **5.8** (2020-05-01)
* by: [Christian Himpe](https://orcid.org/0000-0003-2194-6754)
* under: [BSD-2-Clause](https://opensource.org/licenses/BSD-2-Clause) License
* summary: Empirical system Gramians for (nonlinear) input-output systems.

## Scope

* Model Reduction / Model Order Reduction (MOR)
    - Nonlinear Model Order Reduction (nMOR)
    - Parametric Model Order Reduction (pMOR)
    - Robust Reduction
* Parameter Identification (Structural Identifiability)
* Parameter Reduction
* Combined Reduction (**Combined State and Parameter Reduction**)
* Decentralized Control
* State Sensitivity Analysis
* Parameter Sensitivity Analysis
* Nonlinearity Quantification
* Uncertainty Quantification
* Gramian Indices
    - Optimal Sensor Placement
    - Optimal Actuator Placement
* System Indices
* Tau Function

## Empirical Gramians (Empirical Gramian Matrix Types)

* Empirical Controllability Gramian
* Empirical Observability Gramian
* Empirical Cross Gramian
* Empirical Linear Cross Gramian
* Empirical Sensitivity Gramian
* Empirical Augmented Observability Gramian (Empirical Identifiability Gramian)
* Empirical Joint Gramian (Empirical Cross-Identifiability Gramian)

### Features

* Interfaces for:
    - custom solvers / integrators
    - custom inner products / dot products / kernels / pseudo-kernels
    - distributed / partitioned / low-rank cross Gramian
* Non-symmetric cross gramian option for all cross Gramians:
    - Empirical Cross Gramian
    - Empirical Linear Cross Gramian
    - Empirical Joint Gramian
* Program Design:
    - Implicit parallelization via vectorization of bulk matrix operations
    - Explicit parallelization via `parfor` hints
    - Functional design via closures

### Algorithm

For a mathematical summary and technical documentation of the empirical Gramian framework (5.4),
detailing most features and capabilities see:

* C. Himpe.
  [emgr -- the Empirical Gramian Framework](https://doi.org/10.3390/a11070091).
  Algorithms 11(7): 91, 2018.

and references therein.
See also the [reference lists](emgr-est.md) for further theoretical backgrounds on empirical Gramians.

## Compatibility

* Mathworks MATLAB >=2017b
* GNU Octave >=5.2.0 (in Ubuntu **20.04** LTS)
* Python >=3.6.9 (NumPy >=1.13.3)

## Citation

* C. Himpe (2020). emgr -- EMpirical GRamian Framework (Version 5.8) [Software].
  https://gramian.de doi:[doi:10.5281/zenodo.3779889](https://doi.org/10.5281/zenodo.3779889)

## Getting Started

Run a minimal example in a Matlab interpreter like OCTAVE or MATLAB:
```
RUNME
```

To run all tests use:
```
emgrTest
```

To run demos use:
```
emgrDemo(id) % with id one of 'hnm', 'isp', 'fss', 'nrc', 'rqo', 'lte', 'fbc', 'qso'
```

## Empirical System Theory

The Empirical System Theory (EST) emgr frontend wraps the empirical Gramian framework
to enable system-theoretic computations for linear and nonlinear, small- and medium-scale
input-output systems. The following tasks (and methods) are available:

* Model reduction ([References](emgr-est.md#L104))
    - Empirical poor man (aka POD)
    - Empirical dominant subspaces
    - Approximate balancing (aka Modified POD)
    - Empirical balanced truncation (via Balanced POD)
* Parameter identification ([References](emgr-est.md#L239))
* Parameter reduction
* Combined state and parameter reduction ([References](emgr-est.md#L281))
    - Observability-based
    - Minimality-based
* Decentralized control ([References](emgr-est.md#L308))
    - Relative gain array
    - Input-output coherence
    - Input-output pairing
    - Participation matrix
    - Hardy-2 measure
    - Hardy-∞ measure
    - Hankel interaction
    - Root-mean-square HSV
* State sensitivity ([References](emgr-est.md#L372))
    - Controllability-based
    - Observability-based
    - Minimality-based
* Parameter sensitivity ([References](emgr-est.md#L378))
    - Controllability-based
    - Output-Controllability-based
    - Minimality-based
* Nonlinearity quantification ([References](emgr-est.md#L388))
    - Controllability-based
    - Observability-based
    - Gain-based
    - Correlation-based
* Uncertainty quantification ([References](emgr-est.md#L407))
    - Controllability-based
    - Output-Controllability-based
* Gramian indices ([References](emgr-est.md#L418))
    - [Generalized means](https://en.wikipedia.org/wiki/Generalized_mean) of singular values
    - Storage efficiency
    - Performance index
* System indices ([References](emgr-est.md#L448))
    - Cauchy index
    - System entropy
    - System symmetry
    - Input-output coherence
    - System gain
    - Gramian distance
    - Network sensitivity
    - HSV geometric mean
    - RV-coefficient
* System norms
    - Approximate Hardy-2 norm
    - Approximate Hardy-∞ norm
    - Approximate Hilbert-Schmidt-Hankel norm
    - Approximate Hankel norm
* Tau function ([References](emgr-est.md#L599))

## Files and Folders

[`README.md`](README.md) Basic Information (this file)

[`AUTHORS`](AUTHORS) Author identity and contributions

[`CHANGELOG`](CHANGELOG) Version Information

[`CITATION`](CITATION) Citation Information

[`CODE`](CODE) Meta Information

[`LICENSE`](LICENSE) License Information

[`VERSION`](VERSION) Current Version Number

[`RUNME.m`](RUNME.m) Minimal Code Example

[`emgr.m`](emgr.m) Empirical Gramian Framework (main file, crc32:`08b4da00`)

[`emgrTest.m`](emgrTest.m) System Tests

[`est.m`](est.m) Empirical System Theory (EST) emgr frontend

[`estTest.m`](estTest.m) EST System Tests

[`estProbe.m`](estProbe.m) emgr factorial comparison of singular value decays

[`estDemo.m`](estDemo.m) Run demo (sample applications)

* `'hnm'` Hyperbolic Network Model
* `'isp'` Inverse Sylvester Procedure
* `'fss'` Flexible Space Structures
* `'nrc'` Nonlinear Resistor-Capacitor Cascade
* `'rqo'` Random Diagonal System with Quadratic Output
* `'lte'` Linear Transport Equation
* `'fbc'` Five Body Choreography
* `'qso'` Quasi-Stable Orbits Inside Black Holes

[`est.ogdl`](est.odl) EST tree-based documentation in [ordered graph data lanuage](https://ogdl.org)

[`emgr-est.md`](emgr-est.md) emgr / est reference list

[`emgr-ref.pdf`](emgr-ref.pdf) emgr reference cheat sheet

[`py/emgr.py`](py/emgr.py) Empirical Gramian Framework (Python variant)

[`py/RUNME.py`](RUNME.py) Minimal Code Example (Python variant)

[`py/emgrProbe.py`](emgrProbe.py) Factorial emgr.py tests

## Documentation Table-of-Contents

### [Summary](https://gramian.de/#summary)

### [Scope](https://gramian.de/#scope)

### [Download](https://gramian.de/#download)

### [License](https://gramian.de/#license)

### [Disclaimer](https://gramian.de/#disclaimer)

### [Algorithm](https://gramian.de/#algorithm)

### [Usage](https://gramian.de/#usage)

### [Arguments](https://gramian.de/#arguments)

### [Empirical Gramian Types](https://gramian.de/#gramians)

### [Option Flags](https://gramian.de/#options)

### [Custom Solver](https://gramian.de/#solver)

### [Frontend](https://gramian.de/#extra)

### [Minimal Example](https://gramian.de/#example)

### [Tests](https://gramian.de/#tests)

### [Demos](https://gramian.de/#demos)

### [About](https://gramian.de/#about)

### [References](https://gramian.de/#references)

### [Contact](https://gramian.de/#contact)

### [Links](https://gramian.de/#links)

### [Notes](https://gramian.de/#notes)

### [Troubleshooting](https://gramian.de/#trouble)

## More information at: https://gramian.de

Follow [@modelreduction](https://twitter.com/modelreduction) for `emgr` news.
