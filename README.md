![code meta-data](https://img.shields.io/badge/code_meta--data-%E2%9C%93-brightgreen.svg) 
![zenodo listed](https://zenodo.org/badge/doi/10.5281/zenodo.4454679.png)
![matlab compatible](https://img.shields.io/badge/matlab-compatible-lightgrey.svg)
![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/gramian/emgr/)

![emgr logo](emgr.png) emgr -- EMpirical GRamian Framework (5.9)
================================================================

[Website](https://gramian.de) |
[Twitter](https://twitter.com/modelreduction) |
[Feedback](&#x63;&#x68;&#x40;&#x67;&#x72;&#x61;&#x6D;&#x69;&#x61;&#x6E;&#x2E;&#x64;&#x65;)

* [emgr - EMpirical GRamian Framework](https://gramian.de)
* version: **5.9** (2021-01-21)
* by: [Christian Himpe](https://orcid.org/0000-0003-2194-6754)
* under: [BSD-2-Clause](https://opensource.org/licenses/BSD-2-Clause) License
* summary: Empirical system Gramians for (nonlinear) input-output systems.

## Scope

* Model Reduction / Model Order Reduction (MOR)
    - Nonlinear Model Order Reduction (nMOR)
    - Parametric Model Order Reduction (pMOR) | Robust Reduction
    - Parameter Reduction
    - Combined State and Parameter Reduction
* Decentralized Control
* Sensitivity Analysis
* Parameter Identification | Structural Identifiability | Input-Output Identifiability
* Nonlinearity Quantification
* Uncertainty Quantification
* System Norms | System Indices | System Invariants
* Optimal Sensor Placement | Optimal Actuator Placement
* Matrix Equations
* Tau Functions

## Empirical Gramians

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
* Configurable non-symmetric cross gramian for:
    - Empirical Cross Gramian
    - Empirical Linear Cross Gramian
    - Empirical Joint Gramian
* Program design:
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
See also the [reference lists](https://gramian.de/emgr-est.md) for further theoretical backgrounds on empirical Gramians.

## Compatibility

Successfully tested on:

* Mathworks MATLAB 2017b
* Mathworks MATLAB 2020b
* GNU Octave 5.2.0
* GNU Octave 6.1.0
* Python 3.8.5 (NumPy 1.17.4)

## Citation

* C. Himpe (2020). emgr -- EMpirical GRamian Framework (Version 5.9) [Software].
  https://gramian.de [doi:10.5281/zenodo.4454679](https://doi.org/10.5281/zenodo.4454679)

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
emgrDemo(id) % with id one of 'hnm', 'isp', 'fss', 'nrc', 'rqo', 'lte', 'aps', 'fbc', 'qso'
```

## Files and Folders

[`README.md`](README.md) Basic Information (this file)

[`AUTHORS`](AUTHORS) Author identity and contributions

[`CHANGELOG`](CHANGELOG) Version Information

[`CITATION`](CITATION) Citation Information

[`CODE`](CODE) Meta Information

[`LICENSE`](LICENSE) License Information

[`VERSION`](VERSION) Current Version Number

[`RUNME.m`](RUNME.m) Minimal Code Example

[`emgr.m`](emgr.m) Empirical Gramian Framework (main file, crc32:`c136507c`)

[`emgrTest.m`](emgrTest.m) Run all tests

[`est.m`](est.m) Empirical System Theory (EST) emgr frontend

[`estTest.m`](estTest.m) EST System Tests

[`estDemo.m`](estDemo.m) Run demo (sample applications)

* `'hnm'` Hyperbolic Network Model
* `'isp'` Inverse Sylvester Procedure
* `'fss'` Flexible Space Structures
* `'nrc'` Nonlinear Resistor-Capacitor Cascade
* `'rqo'` Random Diagonal System with Quadratic Output
* `'aps'` All Pass System
* `'lte'` Linear Transport Equation
* `'fbc'` Five Body Choreography
* `'qso'` Quasi-Stable Orbits Inside Black Holes

[`estProbe.m`](estProbe.m) emgr factorial comparison of singular value decays

[`est.scm`](est.scm) EST tree-based documentation in [Scheme](https://en.wikipedia.org/wiki/Scheme_(programming_language))

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

### [Model](https://gramian.de/#model)

### [Usage](https://gramian.de/#usage)

### [Arguments](https://gramian.de/#arguments)

### [Empirical Gramian Types](https://gramian.de/#gramians)

### [Option Flags](https://gramian.de/#options)

### [Return Values](https://gramian.de/#return)

### [Solver Interface](https://gramian.de/#solver)

### [EST Frontend](https://gramian.de/#extra)

### [Minimal Example](https://gramian.de/#example)

### [Tests](https://gramian.de/#tests)

### [Demos](https://gramian.de/#demos)

### [About](https://gramian.de/#about)

### [References](https://gramian.de/#references)

### [Contact](https://gramian.de/#contact)

### [Cite](https://gramian.de/#cite)

### [Links](https://gramian.de/#links)

### [References](https://gramian.de/emgr-est.md)

### [Notes](https://gramian.de/#notes)

### [Troubleshooting](https://gramian.de/#trouble)

### [Frequently Asked Questions](https://gramian.de/#faq)

## More information at: https://gramian.de

Follow [@modelreduction](https://twitter.com/modelreduction) for `emgr` news.
