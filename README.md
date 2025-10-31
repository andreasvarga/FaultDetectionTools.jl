# FaultDetectionTools.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5957950.svg)](https://doi.org/10.5281/zenodo.5957950)
[![DocBuild](https://github.com/andreasvarga/FaultDetectionTools.jl/workflows/CI/badge.svg)](https://github.com/andreasvarga/FaultDetectionTools.jl/actions) 
[![codecov.io](https://codecov.io/gh/andreasvarga/FaultDetectionTools.jl/coverage.svg?branch=master)](https://codecov.io/gh/andreasvarga/FaultDetectionTools.jl?branch=master)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://andreasvarga.github.io/FaultDetectionTools.jl/dev/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://andreasvarga.github.io/FaultDetectionTools.jl/stable/)
[![LinkedIn Badge](https://img.shields.io/badge/LinkedIn-Article-informational?style=flat&logo=linkedin&logoColor=white&color=0D76A8)](https://www.linkedin.com/pulse/faultdetectiontools-julia-package-model-based-fault-diagnosis-varga-osr6e/?trackingId=xcVaQcjNQ8uY2%2FqCFcxRTQ%3D%3D)
[![The MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square)](https://github.com/andreasvarga/FaultDetectionTools.jl/blob/main/LICENSE.md)

## Fault Detection and Isolation Tools in Julia

## Compatibility

Julia 1.10 and higher.

## How to install

````JULIA
pkg> add FaultDetectionTools
pkg> test FaultDetectionTools
````

<!-- For a short interactive demonstration of the main functions execute 

````JULIA
using FaultDetectionTools
cd(joinpath(pkgdir(FaultDetectionTools), "test"))
include("FDIToolsDemo.jl")
````
-->
To execute all examples in Chapters 5, 6 and 7, and all case studies in Chapter 8 in the fault diagnosis book (see below), execute 

````JULIA
using FaultDetectionTools
cd(joinpath(pkgdir(FaultDetectionTools), "Examples"))
include("runexamples.jl")
include("runcasestudies.jl")
````
To execute a particular example, say Example 5.4 and its compact variant 5.4c, execute 

````JULIA
using FaultDetectionTools
cd(joinpath(pkgdir(FaultDetectionTools), "Examples"))
include("Ex5_4.jl")
include("Ex5_4c.jl")
````
To execute a particular case study, say Case Study 2.1, execute 

````JULIA
using FaultDetectionTools
cd(joinpath(pkgdir(FaultDetectionTools), "Examples"))
include("CS2_1.jl")
````
_Note:_ For the execution of the test examples and case study examples, the packages **Polynomials**, **Measurements**, **GenericLinearAlgebra**, **Makie**, **CairoMakie**, **LaTeXStrings**, **JLD2** and **Optim** are also required and must be additionally installed. 

## About

`FaultDetectionTools` is a collection of Julia functions for the analysis and solution 
of fault detection and model detection problems. The functions of this collection rely on 
the [DescriptorSystems](https://github.com/andreasvarga/DescriptorSystems.jl) package. 

The implemented functions are based on the computational procedures described in Chapters 5, 6 and 7 of the book:

**Andreas Varga**, [Solving Fault Diagnosis Problems - Linear Synthesis Techniques with Julia Code Examples](https://link.springer.com/book/10.1007/978-3-031-35767-1), vol. 482 of Studies in Systems, Decision and Control, Springer International Publishing, 2024.

This book describes the mathematical background of solving synthesis problems of fault detection and model detection filters and gives detailed descriptions of the underlying synthesis procedures. 

The implemented functionality parallels the functionality of the MATLAB collection of tools [FDITOOLS](https://github.com/andreasvarga/FDITools), whose User's Guide of the version V1.0 is provided in the file [`fditoolsdoc.pdf`](https://github.com/andreasvarga/FDITools/blob/master/fditoolsdoc.pdf).  


## Supplementary information

Supplementary information on the targeted functionality is also available on [arXiv](https://arxiv.org/abs/1703.08480) in the documentation of the companion MATLAB [`FDITOOLS`](https://github.com/andreasvarga/FDITools) collection.
