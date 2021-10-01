```@meta
CurrentModule = FaultDetectionTools
DocTestSetup = quote
    using FaultDetectionTools
end
```

# FaultDetectionTools.jl

[![DocBuild](https://github.com/andreasvarga/FaultDetectionTools.jl/workflows/CI/badge.svg)](https://github.com/andreasvarga/FaultDetectionTools.jl/actions)
[![Code on Github.](https://img.shields.io/badge/code%20on-github-blue.svg)](https://github.com/andreasvarga/FaultDetectionTools.jl)

The `FaultDetectionTools.jl` package (or shortly`FDITools`) is a collection of Julia functions for the analysis and solution 
of fault detection problems. The functions of this collection relies on 
the [DescriptorSystems](https://github.com/andreasvarga/DescriptorSystems.jl) package [2], whose underlying computational functions are based on the
[`MatrixPencils.jl`](https://github.com/andreasvarga/MatrixPencils.jl) [3] and
[`MatrixEquations.jl`](https://github.com/andreasvarga/MatrixEquations.jl) [4] packages. 

The implemented functions are based on the computational procedures described in Chapters 5, 6 and 7 of the book:

Andreas Varga, "[Solving Fault Diagnosis Problems, Linear Synthesis Techniques](http://www.springer.com/us/book/9783319515588)", vol. 84 of Studies in Systems, Decision and Control, Springer International Publishing, xxviii+394, 2017.

This book describes the mathematical background of solving synthesis problems of fault detection and model detection filters and gives detailed descriptions of the underlying synthesis procedures. 

The targetted functionality parallels the functionality of the MATLAB collection of tools [`FDITOOLS`](https://github.com/andreasvarga/FDITools), whose User's Guide of the version V1.0 is provided in the file [`fditoolsdoc.pdf`](https://github.com/andreasvarga/FDITools/blob/master/fditoolsdoc.pdf).  
Supplementary information on the targetted functionality is also available on [arXiv](https://arxiv.org/abs/1703.08480) in the documentation of the companion MATLAB [`FDITOOLS`](https://github.com/andreasvarga/FDITools) collection.

The available functions in the current version of the `FaultDetectionTools.jl` package are:


**Building fault detection and isolation related objects**

* **[`FDIModel`](@ref)**  Fault detection and isolation synthesis model.
* **[`FDFilter`](@ref)**  Fault detection filter object.
* **[`FDIFilter`](@ref)**  Fault detection and isolation filter object.
* **[`FDFilterIF`](@ref)**  Fault detection filter internal form object.
* **[`FDIFilterIF`](@ref)**  Fault detection and isolation filter internal form object.
* **[`fdimodset`](@ref)**  Setup of synthesis models for solving fault detection and isolation problems.
* **[`fdIFeval`](@ref)**  Evaluation of the internal forms of fault detection and isolation filters. 

**Performance evaluation of FDI filters**

* **[`fditspec`](@ref)**  Computation of the weak or strong structure matrix.
* **[`fdisspec`](@ref)**  Computation of the strong structure matrix.
* **[`fdscond`](@ref)**  Computation of the fault detection sensitivity condition.
* **[`fdif2ngap`](@ref)**  Computation of the fault-to-noise gap.

**Solving fault detection and isolation problems**

* **[`efdsyn`](@ref)**  Exact synthesis of fault detection filters.
* **[`efdisyn`](@ref)**  Exact synthesis of fault detection and isolation filters.

**FDI related computational utilities**

* **[`fdhinfminus`](@ref)**  Evaluation of the `H∞-` index of the transfer function matrix of a descriptor system model.  
* **[`fdhinfmax`](@ref)**  Evaluation of the maximum of column norm of the transfer function matrix of a descriptor system model.  
* **[`fditspec_`](@ref)**  Computation of the weak or strong structure matrix of a descriptor system model.
* **[`fdisspec_`](@ref)**  Computation of the strong structure matrix of a descriptor system model.
* **[`fdscond_`](@ref)**  Computation of the column-gains sensitivity condition of the transfer function matrix of a descriptor system model.

## [Release Notes](https://github.com/andreasvarga/FaultDetectionTools.jl/blob/main/ReleaseNotes.md)

## Main developer

[Andreas Varga](https://sites.google.com/view/andreasvarga/home)

License: MIT (expat)

## References

[1]   A. Varga, Solving Fault Diagnosis Problems – Linear Synthesis Techniques, Vol. 84 of Studies in Systems, Decision and Control, Springer International Publishing, 2017.

[2]  A. Varga, [DescriptorSystems.jl: Manipulation of generalized state-space (descriptor) system representations](https://github.com/andreasvarga/DescriptorSystems.jl). [Zenodo: https://doi.org/10.5281/zenodo.5148319](https://doi.org/10.5281/zenodo.5148319).

[3]  A. Varga, [MatrixPencils.jl: Matrix pencil manipulation using Julia](https://github.com/andreasvarga/MatrixPencils.jl).
[Zenodo: https://doi.org/10.5281/zenodo.3894503](https://doi.org/10.5281/zenodo.3894503).

[4]  A. Varga, [MatrixEquations.jl: Solution of Lyapunov, Sylvester and Riccati matrix equations using Julia](https://github.com/andreasvarga/MatrixEquations.jl). [Zenodo: https://doi.org/10.5281/zenodo.3556867](https://doi.org/10.5281/zenodo.3556867).


