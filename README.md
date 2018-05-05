![emgr logo](emgr.png) emgr - EMpirical GRamian Framework
=========================================================

![code meta-data.](https://img.shields.io/badge/code_meta--data-%E2%9C%93-brightgreen.svg) 
![zenodo listed.](https://zenodo.org/badge/doi/10.5281/zenodo.1241532.png)
![matlab compatible](https://img.shields.io/badge/matlab-compatible-lightgrey.svg)

* emgr - EMpirical GRamian Framework ( http://gramian.de )
* version: **5.4** ( 2018-05-05 )
* by: Christian Himpe ( 0000-0003-2194-6754 )
* under: BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
* summary: Empirical Gramians for model reduction of input-output systems.

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
  * distributed memory / column-wise computation of the cross Gramian
* Non-symmetric cross gramian option for:
  * Empirical Cross Gramian
  * Empirical Linear Cross Gramian
  * Empirical Joint Gramian
* Compatible with:
  * GNU Octave
  * Mathworks MATLAB
* Vectorized and parallelizable
* Functional design

## Basic Usage

Run a minimal example in a Matlab interpreter like OCTAVE or MATLAB:
```
RUNME
```

To use specific utilities, tests or demos include:
```
addpath('demos')
```

## Files and Folders

[`README.md`](README.md) Basic Information

[`RUNME.m`](RUNME.m) Minimal Code Example

[`CODE`](CODE) Code Meta Information

[`CITATION`](CITATION) Citation Information

[`LICENSE`](LICENSE) License Indormation

[`emgr.m`](emgr.m) Empirical Gramian Framework (Main File)

[`emgr_oct.m`](emgr_oct.m) Empirical Gramian Framework (Optional Octave Variant) 

[`emgr_lgc.m`](emgr_lgc.m) Empirical Gramian Framework (Pre 2016b Matlab)

[`emgrtest.m`](emgrtest.m) Run tests

[`emgr_ref.pdf`](emgr_ref.pdf) emgr reference cheat sheet

[`examples.m`](examples.m) Run demos

## Documentation ToC

### [Summary](http://gramian.de/#summary)

### [Scope](http://gramian.de/#scope)

### [Download](http://gramian.de/#download)

### [License](http://gramian.de/#license)

### [Disclaimer](http://gramian.de/#disclaimer)

### [Usage](http://gramian.de/#usage)

### [Arguments](http://gramian.de/#arguments)

### [Empirical Gramian Types](http://gramian.de/#gramians)

### [Option Flags](http://gramian.de/#options)

### [Custom Solver](http://gramian.de/#solver)

### [Extra Utilities](http://gramian.de/#utilities)

### [Tests](http://gramian.de/#tests)

### [Demos](http://gramian.de/#demos)

### [About](http://gramian.de/#about)

### [References](http://gramian.de/#references)

### [Links](http://gramian.de/#links)

### [Notes](http://gramian.de/#notes)

### [Troubleshooting](http://gramian.de/#troubleshooting)

## More information at: http://gramian.de
