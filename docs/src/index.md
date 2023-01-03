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
of _fault detection and isolation_ (FDI) problems and _model detection_ problems. The functions of this collection rely on 
the [DescriptorSystems](https://github.com/andreasvarga/DescriptorSystems.jl) package [2], whose underlying computational functions are based on the
[`MatrixPencils.jl`](https://github.com/andreasvarga/MatrixPencils.jl) [3] and
[`MatrixEquations.jl`](https://github.com/andreasvarga/MatrixEquations.jl) [4] packages. 

The implemented functions are based on the computational procedures described in Chapters 5, 6 and 7 of the book:

Andreas Varga, "[Solving Fault Diagnosis Problems, Linear Synthesis Techniques](http://www.springer.com/us/book/9783319515588)", vol. 84 of Studies in Systems, Decision and Control, Springer International Publishing, 2017.

This book describes the mathematical background of solving synthesis problems of fault detection and model detection filters and gives detailed descriptions of the underlying synthesis procedures. 

The targeted functionality parallels the functionality of the MATLAB collection of tools [`FDITOOLS`](https://github.com/andreasvarga/FDITools), whose User's Guide of the version V1.0 is provided in the file [`fditoolsdoc.pdf`](https://github.com/andreasvarga/FDITools/blob/master/fditoolsdoc.pdf).  
Supplementary information on the targeted functionality is also available on [arXiv](https://arxiv.org/abs/1703.08480) in the documentation of the companion MATLAB [`FDITOOLS`](https://github.com/andreasvarga/FDITools) collection.

The available functions in the current version of the `FaultDetectionTools.jl` package are:


**Building FDI related objects**

* **[`FDIModel`](@ref)**  Fault detection and isolation synthesis model.
* **[`FDFilter`](@ref)**  Fault detection filter object.
* **[`FDIFilter`](@ref)**  Fault detection and isolation filter object.
* **[`FDFilterIF`](@ref)**  Fault detection filter internal form object.
* **[`FDIFilterIF`](@ref)**  Fault detection and isolation filter internal form object.
* **[`fdimodset`](@ref)**  Setup of synthesis models for solving fault detection and isolation problems.
* **[`fdIFeval`](@ref)**  Evaluation of the internal forms of fault detection and isolation filters. 

**Analysis of FDI synthesis models**

* **[`fdigenspec`](@ref)**  Generation of achievable FDI specifications.
* **[`fdichkspec`](@ref)**  Feasibility analysis of a set of FDI specifications.

**Solving FDI problems**

* **[`efdsyn`](@ref)**  Exact synthesis of fault detection filters.
* **[`efdisyn`](@ref)**  Exact synthesis of fault detection and isolation filters.
* **[`afdsyn`](@ref)**  Approximate synthesis of fault detection filters.
* **[`afdisyn`](@ref)**  Approximate synthesis of fault detection and isolation filters.
* **[`emmsyn`](@ref)**  Exact model-matching based synthesis of fault detection filters.
* **[`ammsyn`](@ref)**  Approximate model-matching based synthesis of fault detection filters.

**Performance evaluation of FDI filters**

* **[`fditspec`](@ref)**  Computation of the weak or strong structure matrix.
* **[`fdisspec`](@ref)**  Computation of the strong structure matrix.
* **[`fdiscond`](@ref)**  Computation of the fault detection sensitivity condition.
* **[`fdif2ngap`](@ref)**  Computation of the fault-to-noise gap.
* **[`fdimmperf`](@ref)**  Computation of the model-matching performace.

**Building model detection related objects**

* **[`MDModel`](@ref)**  Model detection component synthesis model.
* **[`MDMModel`](@ref)**  Model detection multiple synthesis model.
* **[`mdmodset`](@ref)**  Setup of multiple synthesis models for solving model detection problems.
* **[`MDFilter`](@ref)**  Model detection filter object.
* **[`MDFilterIF`](@ref)**  Model detection filter internal form object.
* **[`mdIFeval`](@ref)**  Evaluation of the internal forms of model detection filters. 

**Analysis of model detection synthesis models**

* **[`mdgenspec`](@ref)**  Generation of achievable model detection specifications.
* **[`mddist`](@ref)**  Computation of distances between component models.
* **[`mddist2c`](@ref)**  Computation of pairwise distances between two sets of component models.

**Solving model detection problems**

* **[`emdsyn`](@ref)**  Exact synthesis of model detection filters.
* **[`amdsyn`](@ref)**  Approximate synthesis of model detection filters.

**Performance evaluation of model detection filters**

* **[`mdspec`](@ref)**  Computation of the weak structure matrix.
* **[`mdsspec`](@ref)**  Computation of the strong structure matrix.
* **[`mdperf`](@ref)**  Computation of the distance-matching performace.
* **[`mdmatch`](@ref)**  Computation of the distance-matching performace to a component model.
* **[`mdgap`](@ref)**  Computation of the noise gaps.

**Computational utilities**

* **[`fdhinfminus`](@ref)**  Evaluation of the `H∞-` index of the transfer function matrix of a descriptor system model.  
* **[`fdhinfmax`](@ref)**  Evaluation of the maximum of column norm of the transfer function matrix of a descriptor system model.  
* **[`fditspec_`](@ref)**  Computation of the weak or strong structure matrix of a descriptor system model.
* **[`fdisspec_`](@ref)**  Computation of the strong structure matrix of a descriptor system model.
* **[`fdiscond_`](@ref)**  Computation of the column-gains sensitivity condition of the transfer function matrix of a descriptor system model.

## [Release Notes](https://github.com/andreasvarga/FaultDetectionTools.jl/blob/master/ReleaseNotes.md)

## Main developer

[Andreas Varga](https://sites.google.com/view/andreasvarga/home)

License: MIT (expat)

## References

[1]   A. Varga, Solving Fault Diagnosis Problems – Linear Synthesis Techniques, Vol. 84 of Studies in Systems, Decision and Control, Springer International Publishing, 2017.

[2]  A. Varga, [DescriptorSystems.jl: Manipulation of generalized state-space (descriptor) system representations](https://github.com/andreasvarga/DescriptorSystems.jl). [Zenodo: https://doi.org/10.5281/zenodo.5148319](https://doi.org/10.5281/zenodo.5148319).

[3]  A. Varga, [MatrixPencils.jl: Matrix pencil manipulation using Julia](https://github.com/andreasvarga/MatrixPencils.jl).
[Zenodo: https://doi.org/10.5281/zenodo.3894503](https://doi.org/10.5281/zenodo.3894503).

[4]  A. Varga, [MatrixEquations.jl: Solution of Lyapunov, Sylvester and Riccati matrix equations using Julia](https://github.com/andreasvarga/MatrixEquations.jl). [Zenodo: https://doi.org/10.5281/zenodo.3556867](https://doi.org/10.5281/zenodo.3556867).


