# Building fault detection and isolation related objects

* **[`FDIModel`](@ref)**  Falt detection and isolation synthesis model.
* **[`FDFilter`](@ref)**  Fault detection filter object.
* **[`FDFilterIF`](@ref)**  Fault detection filter internal form object.
* **[`fdimodset`](@ref)**  Setup of synthesis models for solving fault detection and isolation problems.

```@docs
FDIModel
fdimodset
FDFilter
FDFilter(::DescriptorStateSpace, ::Int, ::Int)
FDFilterIF
FDFilterIF(::DescriptorStateSpace, ::Int, ::Int, ::Int)
```
