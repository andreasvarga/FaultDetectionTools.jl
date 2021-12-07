# FaultDetectionTools.jl

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4568159.svg)](https://doi.org/10.5281/zenodo.4568159) -->
[![DocBuild](https://github.com/andreasvarga/FaultDetectionTools.jl/workflows/CI/badge.svg)](https://github.com/andreasvarga/FaultDetectionTools.jl/actions) 
[![codecov.io](https://codecov.io/gh/andreasvarga/FaultDetectionTools.jl/coverage.svg?branch=master)](https://codecov.io/gh/andreasvarga/FaultDetectionTools.jl?branch=master)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://andreasvarga.github.io/FaultDetectionTools.jl/dev/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://andreasvarga.github.io/FaultDetectionTools.jl/stable/)
[![The MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square)](https://github.com/andreasvarga/FaultDetectionTools.jl/blob/main/LICENSE.md)

## FDITools - Fault Detection and Isolation Tools in Julia

## Compatibility

Julia 1.6 and higher.

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
To execute all examples in Chapter 5 and case studies in Chapter 8 in the fault diagnosis book (see below), execute 

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

## About

`FDITools` is a collection of Julia functions for the analysis and solution 
of fault detection problems. The functions of this collection relies on 
the [DescriptorSystems](https://github.com/andreasvarga/DescriptorSystems.jl) package. 

The implemented functions are based on the computational procedures described in Chapters 5, 6 and 7 of the book:

**Andreas Varga**, [Solving Fault Diagnosis Problems, Linear Synthesis Techniques](https://www.springer.com/us/book/9783319515588), vol. 84 of Studies in Systems, Decision and Control, Springer International Publishing, xxviii+394, 2017.

This book describes the mathematical background of solving synthesis problems of fault detection and model detection filters and gives detailed descriptions of the underlying synthesis procedures. 

The targetted functionality parallels the functionality of the MATLAB collection of tools [FDITOOLS](https://github.com/andreasvarga/FDITools), whose User's Guide of the version V1.0 is provided in the file [`fditoolsdoc.pdf`](https://github.com/andreasvarga/FDITools/blob/master/fditoolsdoc.pdf).  


## Supplementary information

Supplementary information on the targetted functionality is also available on [arXiv](https://arxiv.org/abs/1703.08480) in the documentation of the companion MATLAB [`FDITOOLS`](https://github.com/andreasvarga/FDITools) collection.
