"""
    FDSystem <: AbstractFDDObject

Type for a FDD system for fault detection or fault isolation.
    
If `fddsys::FDSystem` is the FDD system object, then the matrices of the underlying standard state-space system model 
of the discrete-time residual generator can be obtained via `fddsys.A`, `fddsys.By`, `fddsys.Bu`, `fddsys.C`, `fddsys.Dy`, `fddsys.Du`, `fddsys.A`, 
the sampling time can be obtained via `fddsys.Ts`, the initial state and initial time can be obatined 
via `fddsys.x0` and `fddsys.t0`, while `fddsys.SFDI` contains the structure matrix underlying the synthesis. 
For decision making `fddsys.τ` contains a vector of detection thresholds. 
The discrete-time Narendra-type fault evaluation filter parameters are contained in the vectors `fddsys.α`, `fddsys.β` and `fddsys.γ`.
The current state of the FDD system is formed from the current residual state `fddsys.x`, current time `fddsys.t`,
current residual value `fddsys.r`, current residual evaluation filter state `fddsys.re`, 
and current evaluation signals `fddsys.θ`. For decision making purposes, the vector `fddsys.isig` contains 
the current fault detection status, such that the `i`-th component is set to `1` if the `i`-th evaluation signal
exceeds the `i`-th threshold value.  _Strong fault isolation_ is performed if the structure matrix `fddsys.SFDI = I`
to allow the detection of simultaneous faults. In this case, `fddsys.strongfdi = true` and `fddsys.indfault` contains the indices of currently 
detected faults. If `fddsys.strongfdi = false`, weak fault isolation is performed and `fddsys.indfault = k`, where 
`k` is the index of currently matched fault signature (i.e., the column number of `fddsys.SFDI`). 
`fddsys.indfault` is an empty vector if no signature match occurred. 

This FDD system type is relevant for residual generators determined by one of the functions  [`efdsyn`](@ref), [`afdsyn`](@ref),  
[`emmsyn(::FDIModel, ::FDFilterIF)`](@ref) or [`ammsyn(::FDIModel, ::FDFilterIF)`](@ref).
"""
mutable struct FDSystem{T,T1} <: AbstractFDDObject 
    # filter parameters
    A::Matrix{T}
    By::Matrix{T}
    Bu::Matrix{T}
    C::Matrix{T}
    Dy::Matrix{T}
    Du::Matrix{T}
    Ts::T # sampling time
    SFDI::BitMatrix # structure matrix
    const x0::Vector{T1} # initial state
    const t0::T # initial time
    # detection threshold
    τ::Vector{T}
    # evaluation filter parameters (α,β,γ)
    α::Vector{T}
    β::Vector{T}
    γ::Vector{T}
    x::Vector{T1} # current state
    t::T # current time
    r::Vector{T1} # residual value
    re::Vector{T1} # residual evaluation state 
    θ::Vector{T1} # current evaluation signal 
    isig::Vector{Int} # fault flag (0 - no fault, 1 - fault)
    strongfdi::Bool # indicates strong FDI
    indfault::Vector{Int} # indices of detected faults 
end
"""
    FDSystem(filter::FDFilter[,SFDI::BitMatrix = falses(1,1)]; Ne = size(SFDI,1), 
    x0::Union{Vector,Missing} = missing, t0 = 0, τ::Vector = ones(Ne), α::Vector = zeros(Ne), β::Vector = ones(Ne), γ::Vector = .9*ones(Ne), 
    ) -> fdsys:FDSystem

Build for a `filter::FDFilter` and, optionally, for a structure matrix `SFDI` (Default: `SFDI = BitMatrix([1;;]`) ), 
a FDD system of type [`FDSystem`](@ref) for fault detection or fault isolation. 
`Ne` is the number of residual evaluation signals (i.e., the number of rows of SFDI).
If `Ne = 1` a FDD system for fault detection is built, 
while for `Ne > 1` a FDD system for fault isolation is built.      
The initial state vector of the residual generator can be set using the keyword argument `x0` 
(default: `x0 = 0` if `x0 = missing`) and the initial time can be set via the keyword argument `t0` 
(default: `t0 = 0`).
The vectors `τ`, `α`, `β`, `γ` have the same number of components `Ne`, and contain the 
detection thresholds and the Narendra-type evaluation filter parameters `(α,β,γ)`, respectively.
"""
function FDSystem(f::FDFilter{T},SFDI::Union{BitMatrix,BitVector,Array{Bool,2},Array{Bool,1}} = falses(1,1); 
    Ne = size(SFDI,1), x0::Union{Vector,Missing} = missing, τ = ones(T,Ne), α::Vector = zeros(T,Ne), β::Vector = ones(T,Ne), γ::Vector = .9*ones(T,Ne), t0::T = zero(T),
    ) where {T} 
    f.sys.Ts > 0 || throw(ArgumentError("only discrete-time residual generators supported"))
    Ne == size(SFDI,1) || throw(ArgumentError("number of residuals evaluation signals must match the row size of SFDI"))   
    Ne > 1 && Ne != size(f.sys,1) && throw(ArgumentError("number of residuals evaluation signals must match the number of residuals")) 
    ismissing(x0) && (x0 = zeros(T,order(f.sys)))
    length(x0) == order(f.sys) || throw(DimensionMismatch("order of detector must match the initial condition vector"))
    T1 = eltype(x0)
    Ne == length(τ) || throw(ArgumentError("number of residuals evaluation signals must match the dimensions of vector τ"))
    (Ne == length(α) == length(β) == length(γ)) || throw(ArgumentError("number of residuals evaluation signals must match the dimensions of vectors α, β and γ"))    
    all(τ .>= 0) || throw(ArgumentError("only nonnegative threshold values allowed"))
    any(α .== 0 .&& β .== 0) && throw(ArgumentError("components α[i] and β[i] cannot be simultaneously zero"))
    all(0 .<= γ .&& γ .<= 1) || throw(ArgumentError("only values γ ∈ [0,1] are allowed")) 
    myu = f.ny+f.mu
    return FDSystem{T,T1}(
        f.sys.A, f.sys.B[:,1:f.ny], f.sys.B[:,f.ny+1:myu], f.sys.C, f.sys.D[:,1:f.ny], f.sys.D[:,f.ny+1:myu],
        f.sys.Ts, SFDI, x0, t0, τ, α, β, γ, copy(x0), t0, zeros(T,f.sys.ny), zeros(T,Ne), zeros(T,Ne), zeros(Int,Ne), SFDI == I, Int[])
end
"""
    FDISystem <: AbstractFDDObject

Type for a FDD system for fault isolation.
    
If `fddsys::FDISystem` is the FDD system object, then `fddsys.fdsys` is a vector of FDD systems for fault detection, where the `i`-th FDD system 
is setup based on the `i`-th residual generation filter of a bank of fault detection and isolation filters (see [`FDSystem`](@ref)).
The sampling time can be obtained via `fddsys.Ts`, `fddsys.SFDI` contains the structure matrix underlying the synthesis, the current time is 
contained in `fddsys.t` and the current evaluation signals are contained in the vector `fddsys.θ`. 
For decision making purposes, the vector `fddsys.isig` contains 
the current fault detection status, such that the `i`-th component is set to `1` if the `i`-th evaluation signal
exceeds the `i`-th threshold value.  _Strong fault isolation_ is performed if the structure matrix `fddsys.SFDI = I`
to allow the detection of simultaneous faults. In this case, `fddsys.strongfdi = true` and `fddsys.indfault` contains the indices of currently 
detected faults. If `fddsys.strongfdi = false`, weak fault isolation is performed and `fddsys.indfault = k`, where 
`k` is the index of currently matched fault signature (i.e., the column number of `fddsys.SFDI`). 
`fddsys.indfault` is an empty vector if no signature match occurred. 

This FDD system type is relevant for residual generators determined by one of the functions  [`efdisyn`](@ref), [`afdisyn`](@ref),  
[`emmsyn(::FDIModel, ::FDIFilterIF)`](@ref emmsyn(::FDIModel, ::FDIFilterIF)) or [`ammsyn(::FDIModel, ::FDIFilterIF)`](@ref).
"""
mutable struct FDISystem{T,T1} <: AbstractFDDObject 
    fdsys::Vector{FDSystem{T,T1}}
    Ts::T # sampling time
    SFDI::BitMatrix # structure matrix
    t::T # current time
    θ::Vector{T1} # current evaluation signals 
    isig::Vector{Int} # fault flags (0 - no fault, 1 - fault)
    strongfdi::Bool # indicates strong FDI
    indfault::Vector{Int} # indices of detected faults 
end
"""
    FDISystem(filter::FDIFilter,SFDI::BitMatrix; Ne = size(SFDI,1), 
    x0::Union{Vector,Missing} = missing, t0 = 0, τ::Vector = ones(Ne), α::Vector = zeros(Ne), β::Vector = ones(Ne), γ::Vector = .9*ones(Ne), 
    ) -> fdsys:FDISystem

Build for a `filter::FDIFilter` and a structure matrix `SFDI`, 
a FDD system of type [`FDISystem`](@ref) for fault detection and isolation. 
`Ne` is the number of component residual generators in `filter.sys` (also the number of evaluation signals and the number of rows of `SFDI`).   
The initial state vectors of the bank of `Ne` residual generators contained in `filter.sys` can be set using the keyword argument `x0` as a vector of initial conditions
of appropriate dimensions, where the `i`-th vector `x0[i]` is the initial condition of the `i`-th residual generator contained in `filter.sys[i]`
(default: `x0[i] = 0` for `i = 1, ..., Ne` if `x0 = missing`). The initial time can be set via the keyword argument `t0` 
(default: `t0 = 0`).
The vectors `τ`, `α`, `β`, `γ` have the same number of components `Ne`, and contain the 
detection thresholds and the Narendra-type evaluation filter parameters `(α,β,γ)`, respectively.
"""
function FDISystem(f::FDIFilter{T},SFDI::Union{BitMatrix,BitVector,Array{Bool,2},Array{Bool,1}};
    Ne = size(SFDI,1), x0::Union{Vector,Missing} = missing, 
    τ = ones(T,Ne), α::Vector = zeros(T,Ne), β::Vector = ones(T,Ne), γ::Vector = .9*ones(T,Ne), t0::T = zero(T)) where {T} 
    Ne == length(f.sys) || throw(ArgumentError("number of residuals evaluation signals must match the number of residuals")) 
    if ismissing(x0) 
       x0 = [zeros(T,order(f.sys[i])) for i in 1:Ne]
    else
       Ne == length(x0) || throw(DimensionMismatch("number of filters must be equal to the number of initial condition vectors"))
       all(length.(x0) == order.(f.sys)) || throw(DimensionMismatch("orders of filters must match the initial condition vectors"))
    end
    T1 = eltype(x0[1]) 
    f.sys[1].Ts > 0 || throw(ArgumentError("only discrete-time residual generators supported"))
    (Ne == length(α) == length(β) == length(γ)) || throw(ArgumentError("number of residuals evaluation signals must match the dimensions of vector τ"))
    Ne == length(τ) || throw(ArgumentError("number of residuals evaluation signals must match the dimensions of vectors α, β and γ"))
    all(τ .>= 0) || throw(ArgumentError("only nonnegative threshold values allowed"))
    any(α .== 0 .&& β .== 0) && throw(ArgumentError("α and β cannot be simulataneously zero"))
    all(0 .<= γ .&& γ .<= 1) || throw(ArgumentError("only values γ ∈ [0,1] are allowed")) 
    ny = f.ny; mu = f.mu
    return FDISystem{T,T1}(
    [FDSystem(FDFilter{T}(f.sys[i],ny,mu), SFDI[i:i,:]; x0 = x0[i], t0, τ = [τ[i]], α = [α[i]], β = [β[i]], γ = [γ[i]]) for i in 1:Ne], 
        f.sys[1].Ts, SFDI, t0, zeros(T1,Ne), zeros(Int,Ne), SFDI == I, Int[])
end
"""
    MDSystem <: AbstractFDDObject

Type for a FDD system for model detection.
    
If `mfddsys::MDSystem` is the FDD system object, then `mfddsys.mdsys` is a vector of FDD systems for fault detection, where the `i`-th FDD system 
is setup based on the `i`-th residual generation filter of a bank of model detection filters (see [`MDSystem`](@ref)).
The sampling time can be obtained via `mfddsys.Ts`, the current time is 
contained in `mfddsys.t` and the current evaluation signals are contained in the vector `mfddsys.θ`. 
For decision making purposes, the vector `mfddsys.isig` contains 
the current model detection status, such that the `i`-th component is set to `1` if the `i`-th evaluation signal
exceeds the `i`-th threshold value. `mfddsys.indmodel` contains the index of currently 
detected model and is zero if no signature match occurred. `mfddsys.indminim` provides the index of best matched model, 
corresponding to the least component of the evaluation vector.

This FDD system type is relevant for model detection filters determined by one of the functions  [`emdsyn`](@ref) or [`amdsyn`](@ref).
"""
mutable struct MDSystem{T,T1} <: AbstractFDDObject 
    mdsys::Vector{FDSystem{T,T1}}
    Ts::T # sampling time
    t::T # current time
    θ::Vector{T1} # current evaluation signals 
    isig::Vector{Int} # diagnosis flags (0 - no match, 1 - match)
    indmodel::Int # index of detected model 
    indminim::Int # index of best matched model 
end
"""
    MDSystem(filter::MDFilter; N = length(filter.sys), 
    x0::Union{Vector,Missing} = missing, t0 = 0, τ::Vector = ones(N), α::Vector = zeros(N), β::Vector = ones(N), γ::Vector = .9*ones(N), 
    ) -> mdsys:MDSystem

Build for a `filter::MDFilter` a FDD system of type [`MDSystem`](@ref) for model detection. 
`N` is the number of component model detection filters in `filter.sys` (also the number of evaluation signals).   
The initial state vectors of the bank of `N` model detection filters contained in `filter.sys` can be set using the keyword argument `x0` 
as a vector of initial conditions
of appropriate dimensions, where the `i`-th vector `x0[i]` is the initial condition of the `i`-th model detection filter contained in `filter.sys[i]`
(default: `x0[i] = 0` for `i = 1, ..., N` if `x0 = missing`). The initial time can be set via the keyword argument `t0` 
(default: `t0 = 0`).
The vectors `τ`, `α`, `β`, `γ` have the same number of components `N`, and contain the 
detection thresholds and the Narendra-type evaluation filter parameters `(α,β,γ)`, respectively.
"""
function MDSystem(f::MDFilter{T}; N = length(f.sys), x0::Union{Vector,Missing} = missing, 
    τ = ones(T,N), α::Vector = zeros(T,N), β::Vector = ones(T,N), γ::Vector = .9*ones(T,N), t0::T = zero(T)) where {T} 
    N == length(f.sys) || throw(ArgumentError("number of models must match the number of residuals")) 
    f.sys[1].Ts > 0 || throw(ArgumentError("only discrete-time residual generators supported"))
    if ismissing(x0) 
       x0 = [zeros(T,order(f.sys[i])) for i in 1:N]
    else
       N == length(x0) || throw(DimensionMismatch("number of filters must be equal to the number of initial condition vectors"))
       all(length.(x0) == order.(f.sys)) || throw(DimensionMismatch("orders of filters must match the initial condition vectors"))
    end
    T1 = eltype(x0[1]) 
    # x0::Vector{Vector{T1}} = [zeros(T1,order(f.sys[i])) for i in 1:length(f.sys)], 
    # N == length(x0) || throw(DimensionMismatch("number of filters must be equal to the number of initial condition vectors"))
    # all(length.(x0) == order.(f.sys)) || throw(DimensionMismatch("orders of filters must match the initial condition vectors"))
    (N == length(α) == length(β) == length(γ)) || throw(ArgumentError("number of residuals evaluation signals must match the dimensions of vectors α, β and γ"))
    N == length(τ) || throw(ArgumentError("number of residuals evaluation signals must match the dimensions of vector τ")) 
    all(τ .>= 0) || throw(ArgumentError("only nonnegative threshold values allowed"))
    any(α .== 0 .&& β .== 0) && throw(ArgumentError("α and β cannot be simulataneously zero"))
    all(0 .<= γ .&& γ .<= 1) || throw(ArgumentError("only values γ ∈ [0,1] are allowed")) 
    ny = f.ny; mu = f.mu
    return MDSystem{T,T1}(
    [FDSystem(FDFilter{T}(f.sys[i],ny,mu); x0 = x0[i], t0, τ = [τ[i]], α =[α[i]], β = [β[i]], γ=[γ[i]]) for i in 1:N], 
        f.sys[1].Ts, t0, zeros(T1,N), zeros(Int,N), 0, 0)
end
