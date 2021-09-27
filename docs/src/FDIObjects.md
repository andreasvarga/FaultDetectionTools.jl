# Building fault detection and isolation related objects

* **[`FDIModel`](@ref)**  Fault detection and isolation synthesis model.
* **[`FDFilter`](@ref)**  Fault detection filter object.
* **[`FDIFilter`](@ref)**  Fault detection and isolation filter object.
* **[`FDFilterIF`](@ref)**  Fault detection filter internal form object.
* **[`FDIFilterIF`](@ref)**  Fault detection and isolation filter internal form object.
* **[`fdimodset`](@ref)**  Setup of synthesis models for solving fault detection and isolation problems.
* **[`fdIFeval`](@ref)**  Evaluation of the internal forms of fault detection and isolation filters. 

```@docs
FDIModel
fdimodset
FDFilter
FDFilter(::DescriptorSystems.DescriptorStateSpace, ::Int, ::Int)
FDFilterIF
FDFilterIF(::DescriptorSystems.DescriptorStateSpace, ::Int, ::Int, ::Int)
FDIFilter
FDIFilter(::Array{DescriptorSystems.DescriptorStateSpace{T},1}, ::Int, ::Int) where T
FDIFilterIF
FDIFilterIF(::Array{DescriptorSystems.DescriptorStateSpace{T},1}, ::Int, ::Int, ::Int) where T
fdIFeval
```
