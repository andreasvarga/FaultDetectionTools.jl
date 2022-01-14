"""
    MDModel <: AbstractFDDObject

Type for component synthesis models employed to solve model detection problems.
    
If `sysc::MDModel` is a component synthesis model object, the underlying descriptor system model
can be obtained via `sysc.sys` and the dimensions of the control, disturbance, noise and auxiliary
vectors are contained in the integers `sysm.mu`, `sysm.md`, `sysm.mw` and `sysm.ma`, respectively.
"""
struct MDModel <: AbstractFDDObject 
    sys::DescriptorStateSpace
    mu::Int
    md::Int
    mw::Int
    ma::Int
    function MDModel(sys::DescriptorStateSpace, mu::Int, md::Int, mw::Int, ma::Int)  
        mdmodel_validation(sys, mu, md, mw, ma)
        new(sys, mu, md, mw, ma)
    end
end
"""
    MDModel(sys; mu, md, mw, ma) -> sysc::MDModel

Build for a linear time-invariant descriptor system `sys = (A-λE,B,C,D)` 
a component synthesis model object `sysm::MDModel` 
to be used in conjunction with the analysis and synthesis functions of model detection filters. 

The resulting synthesis model object `sysc` contains the component model with partitioned inputs, 
`sysc.sys = (A-λE,[Bu Bd Bw Bv],C,[Du Dd Dw Dv])`, where 
`Bu`, `Bd`, `Bw` and `Bv` are formed from the successive columns of `B` and are
 the input matrices from the control inputs `u`, disturbance inputs `d`, 
noise inputs `w` and auxiliary inputs `v`, respectively, 
and `Du`, `Dd`, `Dw` and `Dv` are formed from the successive columns of `D` and 
are the feedthrough matrices from those inputs.
The dimensions of control, disturbance, noise and auxiliary input vectors are contained in 
`sysm.mu`, `sysm.md`, `sysm.mw` and `sysm.ma`, respectively.  

The information on the partition of the input components 
in control, disturbance, noise and auxiliary inputs can be specified using the following keyword arguments:

`mu = nu` specifies the dimension `nu` of the control input vector `u` (default: `nu = 0`)

`md = nd` specifies the dimension `nd` of the disturbance input vector `d` (default: `nd = 0`)

`mw = nw` specifies the dimension `nw` of the noise input vector `w` (default: `nw = 0`)

`ma = na` specifies the dimension `na` of the auxiliary input vector `v` (default: `na = 0`)
"""
function MDModel(sys::DescriptorStateSpace; mu::Int = 0, md::Int = 0, mw::Int = 0, ma::Int = 0)  
    MDModel(sys, mu, md, mw, ma)
end
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, sysc::MDModel)
    summary(io, sysc); println(io)
    display(sysc.sys)
    mu = sysc.mu
    md = sysc.md
    mw = sysc.mw
    maux = sysc.ma
    if mu+md+mw+maux > 0
       println(io, "Input groups:")
       println(io, "Name         Channels")
       mu > 0 &&   println(io,"controls     $(1:mu)")
       md > 0 &&   println(io,"disturbances $(mu+1:mu+md)")
       mw > 0 &&   println(io,"noise        $(mu+md+1:mu+md+mw)")
       maux > 0 && println(io,"aux          $(mu+md+mw+1:mu+md+mw+maux)")
    end      
end
"""
    MDMModel <: AbstractFDDObject

Type for multiple synthesis models employed to solve model detection problems.
    
If `sysm::MDMModel` is the multiple synthesis model object, the underlying vector of descriptor system models
can be obtained via `sysm.sys`, the common dimension of the control vectors is contained in the integer
`sysm.mu`, and the dimensions of the disturbance, noise 
and auxiliary input vectors for the `i`-th model `sysm.sys[i]` are contained in the `i`-th 
components of the integer vectors `sysm.md`, `sysm.mw` and `sysm.ma`, respectively.
"""
struct MDMModel <: AbstractFDDObject 
    sys::Vector{<:DescriptorStateSpace}
    mu::Int
    md::Vector{Int}
    mw::Vector{Int}
    ma::Vector{Int}
    function MDMModel(sys::Vector{<:DescriptorStateSpace}, mu::Int, md::Vector{Int}, 
                     mw::Vector{Int}, ma::Vector{Int})  
        mdmodel_validation(sys, mu, md, mw, ma)
        new(sys, mu, md, mw, ma)
    end
end
"""
    MDMModel(sys; mu, md, mw, ma) -> sysm::MDMModel

Build for a vector of linear time-invariant system models `sys`, with the `i`-th component model
`sys[i] = (Ai-λEi,Bi,Ci,Di)`, a multiple synthesis model object `sysm::MDMModel` 
to be used in conjunction with the analysis and synthesis functions of model detection filters. 

The resulting multiple synthesis model object `sysm` contains the vector `sysm.sys` 
of component models with partitioned inputs, 
with the `i`-th model `sysm.sys[i] = (Ai-λEi,[Bui Bdi Bwi Bvi],Ci,[Dui Ddi Dwi Dvi])`, where 
`Bui`, `Bdi`, `Bwi` and `Bvi` are formed from the successive columns of `Bi` and are
 the input matrices from the control inputs `u`, disturbance inputs `di`, 
noise inputs `wi` and auxiliary inputs `vi`, respectively, 
and `Dui`, `Ddi`, `Dwi` and `Dvi` are formed from the successive columns of `Di` and 
are the feedthrough matrices from those inputs.
The dimensions of control, disturbance, noise and auxiliary input vectors are contained in 
`sysm.mu`, `sysm.md[i]`, `sysm.mw[i]` and `sysm.ma[i]`, respectively.  

If `N` is the number of component models, the information on the partition of the input components 
in control, disturbance, noise and auxiliary inputs can be specified using the following keyword arguments:

`mu = nu` specifies the dimension `nu` of control input vector `u` (default: `nu = 0`)

`md = nd` specifies the vector `nd` containing the `N` dimensions of disturbance input vectors, such that
          the `i`-th disturbance vector `di` has dimension `nd[i]` (default: `nd = zeros(Int,N)`)

`mw = nw` specifies the vector `nw` containing the `N` dimensions of noise input vectors, such that
          the `i`-th noise vector `wi` has dimension `nw[i]` (default: `nw = zeros(Int,N)`)

`ma = na` specifies the vector `na` containing the `N` dimensions of auxiliary input vectors, such that
          the `i`-th auxiliary vector `vi` has dimension `na[i]` (default: `na = zeros(Int,N)`)
"""
function MDMModel(sys::Vector{<:DescriptorStateSpace}; mu::Int = 0, 
                  md::Vector{Int} = fill(0,length(sys)), 
                  mw::Vector{Int} = fill(0,length(sys)), 
                  ma::Vector{Int} = fill(0,length(sys)))  
    MDMModel(sys, mu, md, mw, ma)
end
function mdmodel_validation(sys::DST, mu::Int, md::Int, mw::Int, ma::Int) where {DST <: DescriptorStateSpace}
    m = size(sys,2) 
    mu > m && error("number of control inputs $mu exceeds the number of system inputs $m") 
    mu < 0 && error("number of control inputs must be nonnegative")
    md < 0 && error("number of disturbance inputs must be nonnegative")
    mw < 0 && error("number of noise inputs must be nonnegative")
    ma < 0 && error("number of auxiliary inputs must be nonnegative")
    mu+md+mw+ma > m && error("the specified total number of inputs exceeds the number of system inputs")
end

function mdmodel_validation(sys::Vector{<:DescriptorStateSpace}, mu::Int, md::Vector{Int}, 
                            mw::Vector{Int}, ma::Vector{Int}) 
    N = length(sys) 
    N == length(md) || error("improper length of disturbance input vector dimensions")      
    N == length(mw) || error("improper length of noise input vector dimensions")  
    N == length(ma) || error("improper length of auxliary input vector dimensions")  
    m = minimum(size.(sys,2)) 
    mu > m && error("number of control inputs $mu exceeds the minimum number of system inputs $m") 
    mu < 0 && error("number of control inputs must be nonnegative")
    any(md .< 0) && error("number of disturbance inputs must be nonnegative")
    any(mw .< 0) && error("number of noise inputs must be nonnegative")
    any(ma .< 0) && error("number of auxiliary inputs must be nonnegative")
    mi = size.(sys,2)
    any(md+mw+ma-mi .> -mu) && error("specified number of inputs exceeds the number of system inputs")
end
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, sysm::MDMModel)
    summary(io, sysm); println(io)
    mu = sysm.mu
    for i = 1:length(sysm.sys)
        println(io, "\nModel #$i:")
        display(sysm.sys[i])
        mdi = sysm.md[i]
        mwi = sysm.mw[i]
        mauxi = sysm.ma[i]
        if mu+mdi+mwi+mauxi > 0
            println(io, "Input groups:")
            println(io, "Name         Channels")
            mu > 0 &&    println(io,"controls     $(1:mu)")
            mdi > 0 &&   println(io,"disturbances $(mu+1:mu+mdi)")
            mwi > 0 &&   println(io,"noise        $(mu+mdi+1:mu+mdi+mwi)")
            mauxi > 0 && println(io,"aux          $(mu+mdi+mwi+1:mu+mdi+mwi+mauxi)")
         end
           
    end
end
"""
    mdmodset(sys; controls, c, disturbances, d, noise, n, aux) -> sysm::MDMModel

Build for a vector of linear time-invariant system models `sys`, with the `i`-th component model 
`sys[i] = (Ai-λEi,Bi,Ci,Di)`, a multiple synthesis model object `sysm::MDMModel` 
to be used in conjunction with the analysis and synthesis functions of model detection filters. 

The information on the partition of the input components in control, disturbance, noise and auxiliary
inputs can be specified using the following keyword arguments:

`controls = inpu` or `c = inpu` specifies the indices `inpu` of the control inputs (default: void)

`disturbances = inpd` or `d = inpd` specifies the indices or a vector of indices `inpd` of the disturbance inputs (default: void)

`noise = inpn` or  `noise = inpn` specifies the indices or a vector of indices `inpn` of the noise inputs (default: void)

`aux = inpa` specifies the indices or a vector of indices `inpa` of the auxiliary inputs (default: void)

The indices of inputs can be specified as integer vectors, integer scalars or integer `UnitRange`s.
For disturbance, noise and auxiliary inputs, vectors of integer vectors or 
vectors of integer `UnitRange`s can be used to specify possibly different set of indices for each component model.

The resulting `sysm` contains the vector `sysm.sys` of partitioned systems, with the `i`-th model  
`sysm.sys[i] = (Ai-λEi,[Bui Bdi Bwi Bvi],Ci,[Dui Ddi Dwi Dvi])`, where 
`Bui`, `Bdi`, `Bwi` and `Bvi` are the input matrices from the control inputs `u`, disturbance inputs `di`, 
noise inputs `wi` and auxiliary inputs `vi`, respectively, and `Dui`, `Ddi`, `Dwi` and `Dvi` are the feedthrough matrices from those inputs.
The dimensions of control, disturbance, fault, noise and auxiliary input vectors are contained in 
`sysm.mu`, `sysm.md[i]`, `sysm.mw[i]` and `sysm.ma[i]`, respectively.  

_Method:_ If `Gi(λ)` is the `p x mi` transfer function matrix of `sys[i]`, then the resulting component system `sysm.sys[i]` has an 
equivalent input output form `[Gui(λ) Gdi(λ) Gwi(λ) Gvi(λ)]`, where the following relations define the component matrices:
`Gui(λ) = Gi(λ)*Su`, `Gdi(λ) = Gi(λ)*Sdi`,  `Gwi(λ) = G(λ)*Swi`,  `Gvi(λ) = Gi(λ)*Svi`, with
the selection matrices `Su`, `Sdi`, `Swi` and `Svi` formed from the columns of the `mi`-th order identity matrix 
(note that `Su` is the same for all component models). 
"""
function mdmodset(sys::Vector{<:DescriptorStateSpace}; 
                   controls::VRS = Int[], c::VRS = Int[], 
                   disturbances::Union{VRS,Vector{Vector{Int}},Vector{UnitRange{Int}}} = fill(Int[],length(sys)), 
                   d::Union{VRS,Vector{Vector{Int}},Vector{UnitRange{Int}}} = fill(Int[],length(sys)), 
                   noise::Union{VRS,Vector{Vector{Int}},Vector{UnitRange{Int}}} = fill(Int[],length(sys)), 
                   n::Union{VRS,Vector{Vector{Int}},Vector{UnitRange{Int}}} = fill(Int[],length(sys)), 
                   aux::Union{VRS,Vector{Vector{Int}},Vector{UnitRange{Int}}} = fill(Int[],length(sys))) 
    N = length(sys)  
    if N > 0
       p = size(sys[1],1) 
       Ts = sys[1].Ts
    end            
    m = maximum(size.(sys,2)) 
    sysn = similar(sys,N)
    inpu = sort(unique([controls;c]))
    mu = length(inpu)
    mu == 0 || inpu[mu] <= m || error("selected index/indices of control inputs larger than the minimum number of system inputs $m")   
    vecdist = isa(disturbances,Vector)
    vecd = isa(d,Vector)
    vecnoise = isa(noise,Vector)
    vecn = isa(n,Vector)
    vecaux = isa(aux,Vector)
    vecdist && length(disturbances) != N && error("the vector of disturbance input indices must have length $N") 
    vecd && length(disturbances) != N && error("the vector of disturbance input indices must have length $N") 
    vecnoise && length(noise) != N && error("the vector of noise input indices must have length $N") 
    vecn && length(n) != N && error("the vector of noise input indices must have length $N") 
    vecaux && length(aux) != N && error("the vector of auxiliary input indices must have length $N") 
    md = similar(Vector{Int},N)
    mw = similar(Vector{Int},N)
    ma = similar(Vector{Int},N)
    for i = 1:N
        i == 1 || (Ts = DescriptorSystems.promote_Ts(Ts,sys[i].Ts))
    end
    for i = 1:N
        i == 1 || p == size(sys[i],1) || error("all systems must have the same number of outputs")
        if vecdist && vecd
            inpdi = unique([disturbances[i];d[i]])
        elseif vecn
            inpdi = unique([disturbances;d[i]])
        elseif vecnooise
            inpdi = unique([disturbances[i];d])
        else
            inpdi = unique([disturbances;d])
        end
        isempty(intersect(inpu,inpdi)) || error("control and disturbance inputs must be distinct")
        if vecnoise && vecn
            inpwi = unique([noise[i];n[i]])
        elseif vecn
            inpwi = unique([noise;n[i]])
        elseif vecnooise
            inpwi = unique([noise[i];n])
        else
            inpwi = unique([noise;n])
        end
        inpai = vecaux ? unique([aux[i];Int[]]) : unique([aux;Int[]])
        mdi = length(inpdi)
        mwi = length(inpwi)
        mai = length(inpai)
        mdi == 0 || maximum(inpdi) <= m || error("selected index/indices of disturbance inputs larger than the minimum number of system inputs $m")   
        mwi == 0 || maximum(inpwi) <= m || error("selected index/indices of noise inputs larger than the minimum number of system inputs $m")   
        mai == 0 || maximum(inpauxi) <= m || error("selected index/indices of auxiliary inputs larger than the minimum number of system inputs $m")  
        sysn[i] = dss(dssdata(sys[i][:,[inpu;inpdi;inpwi;inpai]])...;Ts)
        md[i] = mdi
        mw[i] = mwi
        ma[i] = mai
    end
    return MDMModel(sysn; mu, md, mw, ma)          
end
"""
    mdmodset(sysc::Vector{<:MDModel}) -> sysm::MDMModel

Build for a vector of component synthesis models `sysc`, a multiple synthesis model object `sysm::MDMModel` 
to be used in conjunction with the analysis and synthesis functions of model detection filters. 
"""
function mdmodset(sysc::Vector{<:MDModel}) 
    N = length(sysc)  
    sysn = similar(Vector{DescriptorStateSpace},N)
    md = similar(Vector{Int},N)
    mw = similar(Vector{Int},N)
    ma = similar(Vector{Int},N)
    if N > 0
       p = size(sysc[1].sys,1) 
       mu = sysc[1].mu
       Ts = sysc[1].sys.Ts
       sysn[1] = sysc[1].sys
       md[1] = sysc[1].md
       mw[1] = sysc[1].mw
       ma[1] = sysc[1].ma
    end   
    for i = 2:N  
        p == size(sysc[i].sys,1) || error("all component models must have the same number of outputs")
        mu == sysc[i].mu || error("all component models must have the same number of control inputs")
        Ts = DescriptorSystems.promote_Ts(Ts,sysc[i].sys.Ts)
        sysn[i] = sysc[i].sys
        md[i] = sysc[i].md
        mw[i] = sysc[i].mw
        ma[i] = sysc[i].ma
    end
    return MDMModel(sysn; mu, md, mw, ma)          
end

"""
    MDFilter <: AbstractFDDObject

Type for model detection filters resulted as solutions of model detection problems.
    
If `filter::MDFilter` is the model detection filter object, 
the underlying `i`-th descriptor system model
can be obtained via `filter.sys[i]` and the dimensions of the partitioned filter input vectors as 
`outputs` and `controls`,  
can be accessed as the integers contained in `filter.ny` and `filter.mu`, respectively.
"""
struct MDFilter{T} <: AbstractFDDObject where T 
    sys::Vector{DescriptorStateSpace{T}}
    ny::Int
    mu::Int
    function MDFilter{T}(sys::Vector{DescriptorStateSpace{T}}, ny::Int, mu::Int) where T 
        N = length(sys)
        sysn = similar(sys,N)
        
        mf = ny+mu
        mf > size(sys[1],2) && error("total number of inputs exceeds the number of filter inputs $m") 
   
        sysn[1] = sys[1][:,1:mf]
        for i = 1:N
            mf > size(sys[i],2) && error("total number of inputs exceeds the number of filter inputs $m") 
            sysn[i] = sys[i][:,1:mf]
        end
        new{T}(sysn, ny, mu)
    end
end
# """
#     MDFilter(sys, p, mu) -> Q::MDFilter

# Build for a vector of linear time-invariant descriptor system models `sys[i] = (Ai-λEi,Bi,Ci,Di)`
# with the same number of inputs, a model detection filter object `Q`, 
# as determined with the synthesis functions of model detection filters. 
# `p` and `mu` are the number of measured outputs and the number of control inputs, respectively. 
# It is assumed that each `Bi = [Byi Bui Bvi]` and `Di = [Dyi Dui Dvi]` are partitioned matrices such that
# `Byi` and `Dyi` have `p` columns, and `Bui` and `Dui` have `mu` columns, 
# where `Byi` and `Bui` are the input matrices from the measured outputs `y` and 
# control inputs `u`, `Dyi` and `Dui` are the feedthrough matrices from the measured outputs `y` and 
# control inputs `u`. 

# The resulting `Q` contains the vector of partitioned systems
# `Q.sys[i] = (Ai-λEi,[Byi Bdi],Ci,[Dyi Dui])` and the indices of inputs corresponding 
# to the measured outputs and control inputs are contained in the associated 
# integer vectors `Q.outputs` and `Q.controls`. 
# """
# function MDFilter(sys::Vector{DescriptorStateSpace{T}}, p::Int, mu::Int) where T
#     m = size(sys[1],2) 
#     p < 0 && error("number of measured outputs must be non-negative")
#     mu < 0 && error("number of control inputs must be non-negative")
#     p+mu > m && error("number of measured outputs and control inputs exceeds the number of filter inputs $m")
#     return MDFilter{T}(sys, Vector(1:p), Vector(p+1:p+mu))          
# end
# function MDFilter(sys::Vector{DescriptorStateSpace}, p::Int, mu::Int)
#     m = size(sys[1],2) 
#     p < 0 && error("number of measured outputs must be non-negative")
#     mu < 0 && error("number of control inputs must be non-negative")
#     p+mu > m && error("number of measured outputs and control inputs exceeds the number of filter inputs $m")
#     return MDFilter{eltype(sys[1])}(sys, Vector(1:p), Vector(p+1:p+mu))          
# end

# function MDFilter(sys::Vector{DescriptorStateSpace{T}}; outputs::VRS = Int[], controls::VRS = Int[]) where T 
#     return MDFilter{T}(sys, vec([outputs; Int[]]), vec([controls; Int[]]))          
# end
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, filter::MDFilter)
    summary(io, filter); println(io)
    for i = 1:length(filter.sys)
        println(io, "Filter #$i:")
        display(filter.sys[i])
    end
    p = filter.ny
    mu = filter.mu
    if p+mu > 0
       println(io, "Input groups:")
       println(io,          "Name      Channels")
       p > 0 &&  println(io,"outputs   $(1:p)")
       mu > 0 && println(io,"controls  $(p+1:p+mu)")
   end
end

"""
    MDFilterIF <: AbstractFDDObject

Type for the internal form of model detection filters resulted as solutions of model detection problems.
    
If `filter::MDFilterIF` is the model detection filter internal form object, 
the underlying `(i,j)`-th descriptor system model
can be obtained via `filter.sys[i,j]` and the corresponding
dimensions of control, disturbance, fault, noise and auxiliary input vectors are contained in 
`sysm.mu`, `sysm.md[j]`, `sysm.mw[j]` and `sysm.ma[j]`, respectively.  
"""
struct MDFilterIF{T} <: AbstractFDDObject where T 
    sys::Matrix{DescriptorStateSpace{T}}
    mu::Int
    md::Vector{Int}
    mw::Vector{Int}
    ma::Vector{Int}
    function MDFilterIF{T}(sys::Matrix{DescriptorStateSpace{T}}, mu::Int, md::Vector{Int}, 
                           mw::Vector{Int}, ma::Vector{Int}) where T 
        M, N = size(sys)
        sysn = similar(sys,M,N)
        for j = 1:N
            mj = mu+md[j]+mw[j]+ma[j]
            mj > size(sys[1,j],2) && error("number of inputs exceeds the number of system inputs $m") 
            sysn[1,j] = sys[1,j][:,1:mj]
            for i = 2:M
                mj > size(sys[i,j],2) && error("number of inputs exceeds the number of system inputs $m") 
                sysn[i,j] = sys[i,j][:,1:mj]
            end
        end
        new{T}(sysn, mu, md, mw, ma)
    end
end
# function mdfilterIF_validation(sys::DescriptorStateSpace{T}, mu::Int, md::Int, mw::Int, ma::Int) where T
#     m = size(sys,2) 
    
#     inpu = unique(controls); 
#     inpd = unique(disturbances); 
#     inpw = unique(noise); 
#     inpaux = unique(aux); 
    
#     isempty(intersect(inpu,inpd)) || error("control and disturbance inputs must be distinct")
#     isempty(intersect(inpd,inpw)) || error("disturbance and noise inputs must be distinct")
#     isempty(intersect(inpw,inpaux)) || error("noise and aux inputs must be distinct")  
    
#     mu+md+mw+ma > m && error("number of inputs exceeds the number of system inputs $m") 
    
#     return inpu, inpd, inpw, inpaux
# end
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, filter::MDFilterIF)
    summary(io, filter); println(io)
    M, N = size(filter.sys)
    mu = filter.mu
    for j = 1:N
        md = filter.md[j]
        mw = filter.mw[j]
        maux = filter.ma[j]
        for i = 1:M
            println(io, "\n(Filter,Model) = (#$i,#$j):")
            display(filter.sys[i,j])
            if mu+md+mw+maux > 0
               println(io, "Input groups:")
               println(io, "Name         Channels")
               mu > 0 && println(io,"controls     $(1:mu)")
               md > 0 && println(io,"disturbances $(mu+1:mu+md)")
               mw > 0 && println(io,"noise        $(mu+md+1:mu+md+mw)")
               maux > 0 && println(io,"aux         $(mu+md+mw+1:mu+md+mw+maux)")
            end
        end
    end
end
"""
    mdIFeval(sysQ::MDFilter, sysm::MDMModel; minimal = false, atol, atol1 = atol, atol2 = atol, rtol, fast = true) -> sysR::MDFilterIF

Compute the internal form `sysR` of the model detection filter `sysQ` applied to the multiple synthesis model `sysm`. 
If the `j`-th component model `sysm.sys[j]` has the partitioned transfer function matrix `Gj(λ) = [ Guj(λ)  Gdj(λ) Gwj(λ) Gvj(λ) ]` in accordance with
the partitioned system inputs as `controls`, `disturbances`, `noise` and `auxiliary` inputs, respectively,
and `Qi(λ) = [ Qyi(λ) Qui(λ) ]` is the partitioned transfer function matrix of the `i`-th filter `sysQ.sys[i]` 
in accordance with the partitioned filter inputs as `outputs` and `controls`, then 
the transfer function matrix `Rij(λ)` corresponding to the (`i`-th filter,`j`-th model) pair in the resulting internal form `sysR.sys[i,j]` is given by
     
     Rij(λ) = | Qyi(λ)  Qui(λ) | * | Guj(λ)  Gdj(λ) Gwj(λ) Gvj(λ) |
                                   |  I      0      0      0      |

Minimal descriptor realizations are computed if `minimal = true` and a (possibly) non-minimal 
realization is determined if `minimal = false` (default). 

The minimal realization computation relies on pencil manipulation algorithms which 
employ rank determinations based on either the use of 
rank revealing QR-decomposition with column pivoting, if `fast = true`, or the SVD-decomposition.
The rank decision based on the SVD-decomposition is generally more reliable, but the involved computational effort is higher.

If `(Arij-λErij,Brij,Crij,Drij)` is the full order descriptor realization of `sysR.sys[i,j]`, then 
the keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Arij`, `Brij`, `Crij`, `Drij`, the absolute tolerance for the nonzero elements of `Erij`,  
and the relative tolerance for the nonzero elements of `Arij`, `Brij`, `Crij`, `Drij` and `Eirj`.
The default relative tolerance is `nij*ϵ`, where `ϵ` is the working _machine epsilon_ 
and `nij` is the order of the system `sysR.sys[i,j]`.  
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 

"""
function mdIFeval(Q::MDFilter, sysm::MDMModel; minimal::Bool = false, 
                  atol::Real = 0., atol1::Real = atol, atol2::Real = atol, rtol::Real = 0., fast::Bool = true)
   # compute the internal form as
   #   R[i,j] = Q[i] * [ Gu[j] Gd[j] Gw[j] Gaux[j]]
   #                   [ I     0     0     0      ]
   M = length(Q.sys)
   N = length(sysm.sys)
   T1 = eltype(Q.sys[1])
   sysR = similar(Matrix{DescriptorStateSpace{T1}},M,N)
   p = size(sysm.sys[1],1)
   mu = sysm.mu
   (Q.ny == p && Q.mu == mu) || error("filter Q is incompatible with the given multiple model")
   for j = 1:N
       m = size(sysm.sys[j],2)
       sysej = [sysm.sys[j];  I zeros(mu,m-mu)]
       for i = 1:M
           sysR[i,j] = minimal ? gminreal(Q.sys[i]*sysej; atol1, atol2, rtol, fast) : Q.sys[i]*sysej
       end
   end
   return MDFilterIF{T1}(sysR, sysm.mu, sysm.md, sysm.mw, sysm.ma) 
end

