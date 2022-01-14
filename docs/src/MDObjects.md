# Building model detection related objects

* **[`MDModel`](@ref)**  Model detection component synthesis model.
* **[`MDMModel`](@ref)**  Model detection multiple synthesis model.
* **[`mdmodset`](@ref)**  Setup of multiple synthesis models for solving model detection problems.
* **[`MDFilter`](@ref)**  Model detection filter object.
* **[`MDFilterIF`](@ref)**  Model detection filter internal form object.
* **[`mdIFeval`](@ref)**  Evaluation of the internal forms of model detection filters. 

```@docs
MDModel
MDModel(sys::DescriptorStateSpace; mu, md, mw, ma) 
MDMModel
MDMModel(sys::Vector{<:DescriptorStateSpace}; mu, md, mw, ma)
mdmodset
MDFilter
MDFilterIF
mdIFeval
```
