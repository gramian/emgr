![emgr logo](emgr.png) emgr - EMpirical GRamian Framework
=========================================================

![code meta-data.](https://img.shields.io/badge/code_meta--data-%E2%9C%93-brightgreen.svg) 
![zenodo listed.](https://zenodo.org/badge/doi/10.5281/zenodo.1401500.png)
![matlab compatible](https://img.shields.io/badge/matlab-compatible-lightgrey.svg)

* emgr - EMpirical GRamian Framework ([gramian.de](https://gramian.de))
* version: **5.5** (2018-08-22)
* by: Christian Himpe ([0000-0003-2194-6754](https://orcid.org/0000-0003-2194-6754))
* under: [BSD-2-Clause](https://opensource.org/licenses/BSD-2-Clause) License
* summary: Empirical Gramians for (model reduction of) input-output systems.

## Scope

* Model Reduction / Model Order Reduction (MOR)
  * parametric Model Order Reduction (pMOR) / Robust Reduction
  * nonlinear Model Order Reduction (nMOR)
  * Parameter Identification / Parameter Reduction
  * **Combined State and Parameter Reduction** (Combined Reduction)
* Sensitivity Analysis
* Decentralized Control
* Optimal Sensor Placement / Optimal Actuator Placement
* Nonlinearity Quantification

## Empirical Gramians (Empirical Gramian Matrix Types)

* Empirical Controllability Gramian
* Empirical Observability Gramian
* Empirical Cross Gramian
* Empirical Linear Cross Gramian
* Empirical Sensitivity Gramian
* Empirical Identifiability Gramian
* Empirical Joint Gramian (Empirical Cross-Identifiability Gramian)

### Features

* Interfaces for:
  * custom solvers / integrators
  * custom inner products / dot products / kernels / pseudo-kernels
  * distributed / partitioned / low-rank cross Gramian
* Non-symmetric cross gramian option for:
  * Empirical Cross Gramian
  * Empirical Linear Cross Gramian
  * Empirical Joint Gramian
* Compatible with:
  * GNU Octave
  * Mathworks MATLAB
* Vectorized and parallelizable
* Functional design

### Details

For a mathematical summary and technical documentation of the empirical Gramian framework
detailing all features and capabilities see:

* C. Himpe. "[emgr -- the Empirical Gramian Framework](https://doi.org/10.3390/a11070091)". Algorithms 11(7): 91, 2018.

## Compatbility

* GNU Octave >=4.2.0 (available in Ubuntu **18.04** LTS repositories)
* Mathworks MATLAB >=2016b (use `emgr_lgc.m` for earlier version)

## Basic Usage

Run a minimal example in a Matlab interpreter like OCTAVE or MATLAB:
```
RUNME
```

To run all tests use:
```
emgrtest
```

To run demos use:
```
examples(id) % with id one of 'hnm', 'isp', 'fss', 'nrc', 'lte', 'fbc', 'qso'
```

## Files and Folders

[`README.md`](README.md) Basic Information

[`RUNME.m`](RUNME.m) Minimal Code Example

[`CODE`](CODE) Code Meta Information

[`CITATION`](CITATION) Citation Information

[`LICENSE`](LICENSE) License Information

[`CHANGELOG`](CHANGELOG) Version Change Information

[`emgr.m`](emgr.m) Empirical Gramian Framework (Main File)

[`emgr_oct.m`](emgr_oct.m) Empirical Gramian Framework (Optional Octave Variant) 

[`emgr_lgc.m`](emgr_lgc.m) Empirical Gramian Framework (Pre 2016b Matlab Variant)

[`curios.m`](curios.m) Clearing Up Reducibility of Input-Output Systems (Simple frontend)

[`emgrtest.m`](emgrtest.m) Run tests (system tests)

[`moretests.m`](moretests.m) More tests (functionality tests)

[`examples.m`](examples.m) Run demo
  * `'hnm'` Hyperbolic Network Model
  * `'isp'` Inverse Sylvester Procedure
  * `'fss'` Flexible Space Structures
  * `'nrc'` Nonlinear Resistor-Capacitor Cascade
  * `'lte'` Linear Transport Equation
  * `'fbc'` Five Body Choreography
  * `'qso'` Quasi-Stable Orbits Inside Black Holes

[`emgr-ref.pdf`](emgr-ref.pdf) emgr reference cheat sheet

## Documentation Table-of-Contents

### [Summary](https://gramian.de/#summary)

### [Scope](https://gramian.de/#scope)

### [Download](https://gramian.de/#download)

### [License](https://gramian.de/#license)

### [Disclaimer](https://gramian.de/#disclaimer)

### [Usage](https://gramian.de/#usage)

### [Arguments](https://gramian.de/#arguments)

### [Empirical Gramian Types](https://gramian.de/#gramians)

### [Option Flags](https://gramian.de/#options)

### [Custom Solver](https://gramian.de/#solver)

### [Extra Utilities](https://gramian.de/#extra)

### [Tests](https://gramian.de/#tests)

### [Demos](https://gramian.de/#demos)

### [About](https://gramian.de/#about)

### [References](https://gramian.de/#references)

### [Links](https://gramian.de/#links)

### [Notes](https://gramian.de/#notes)

### [Troubleshooting](https://gramian.de/#trouble)

## More information at: https://gramian.de
