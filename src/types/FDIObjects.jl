"""
    FDIModel <: AbstractFDDObject

Type for synthesis models employed to solve fault detection and isolation problems.
    
If `sysf::FDIModel` is the synthesis model object, the underlying descriptor system model
can be obtained via `sysf.sys` and the dimensions of the control, disturbance, fault, noise and auxiliary
vectors are contained in the integers `sysf.mu`, `sysf.md`, `sysf.mf`, `sysf.mw` and `sysf.ma`, respectively.
The ranges of indices of control, disturbance, fault, noise and auxiliary inputs can be accessed as
`sysf.controls`, `sysf.disturbances`, `sysf.faults`, `sysf.noise` and `sysf.aux`, respectively.
"""
struct FDIModel{T} <: AbstractFDDObject where T 
    sys::DescriptorStateSpace{T}
    mu::Int
    md::Int
    mf::Int
    mw::Int
    ma::Int
    function FDIModel{T}(sys::DescriptorStateSpace{T}, mu::Int, md::Int, mf::Int, mw::Int, ma::Int)  where T 
        fdimodel_validation(sys, mu, md, mf, mw, ma)
        new{T}(sys[:,1:mu+md+mf+mw+ma], mu, md, mf, mw, ma)
    end
end
"""
    FDIModel(sys; mu, md, mf, mw, ma) -> sysf::FDIModel

Build for a linear time-invariant descriptor system `sys = (A-λE,B,C,D)` 
a synthesis model object `sysf::FDIModel` 
to be used in conjunction with the analysis and synthesis functions of fault detection filters. 

The resulting synthesis model object `sysf` contains the component model with partitioned inputs, 
`sysc.sys = (A-λE,[Bu Bd Bf Bw Bv],C,[Du Dd Df Dw Dv])`, where 
`Bu`, `Bd`, `Bf`, `Bw` and `Bv` are formed from the successive columns of `B` and are
 the input matrices from the control inputs `u`, disturbance inputs `d`, fault inputs `f`,
noise inputs `w` and auxiliary inputs `v`, respectively, 
and `Du`, `Dd`, `Df`, `Dw` and `Dv` are formed from the successive columns of `D` and 
are the feedthrough matrices from those inputs.
The dimensions of control, disturbance, fault, noise and auxiliary input vectors are contained in 
`sysf.mu`, `sysf.md`, `sysf.mf`, `sysf.mw` and `sysf.ma`, respectively.  

The information on the partition of the input components 
in control, disturbance, fault, noise and auxiliary inputs can be specified using the following keyword arguments:

`mu = nu` specifies the dimension `nu` of the control input vector `u` (default: `nu = 0`)

`md = nd` specifies the dimension `nd` of the disturbance input vector `d` (default: `nd = 0`)

`mf = nf` specifies the dimension `nf` of the fault input vector `f` (default: `nf = 0`)

`mw = nw` specifies the dimension `nw` of the noise input vector `w` (default: `nw = 0`)

`ma = na` specifies the dimension `na` of the auxiliary input vector `v` (default: `na = 0`)
"""
function FDIModel{T}(sys::DescriptorStateSpace{T}; mu::Int = 0, md::Int = 0, mf::Int = 0, mw::Int = 0, ma::Int = 0)  where T
    FDIModel{T}(sys, mu, md, mf, mw, ma)
end
function fdimodel_validation(sys::DescriptorStateSpace{T}, mu::Int, md::Int, mf::Int, mw::Int, ma::Int) where T
    m = size(sys,2) 
    mu < 0 && error("number of control inputs must be nonnegative")
    md < 0 && error("number of disturbance inputs must be nonnegative")
    mf < 0 && error("number of disturbance inputs must be nonnegative")
    mw < 0 && error("number of noise inputs must be nonnegative")
    ma < 0 && error("number of auxiliary inputs must be nonnegative")
    mu+md+mf+mw+ma > m && error("the specified total number of inputs exceeds the number of system inputs")
end
function FDIModel(sys::DescriptorStateSpace, mu::Int, md::Int, mf::Int, mw::Int, ma::Int)
    FDIModel{eltype(sys)}(sys; mu, md, mf, mw, ma)
end
"""
    FDFilter <: AbstractFDDObject

Type for fault detection filters resulted as solutions of fault detection problems.
    
If `filter::FDFilter` is the fault detection filter object, the underlying descriptor system model
can be obtained via `filter.sys` and the dimensions of the partitioned filter input vectors as 
`measured outputs` and `control inputs`,  
can be accessed as the integers contained in `filter.ny` and `filter.mu`, respectively.
The ranges of the indices of output and control inputs 
can be accessed as `filter.outputs` and `filter.controls`, respectively.
"""
struct FDFilter{T} <: AbstractFDDObject where T 
    sys::DescriptorStateSpace{T}
    ny::Int
    mu::Int
    function FDFilter{T}(sys::DescriptorStateSpace{T}, ny::Int, mu::Int) where T 
        ny < 0 && error("number of measured outputs must be non-negative")
        mu < 0 && error("number of control inputs must be non-negative")
        m = ny+mu
        m > size(sys,2) && error("total number of inputs exceeds the number of filter inputs") 
        new{T}(sys[:,1:m], ny, mu)
    end
end
"""
    FDFilter(sys, ny, mu) -> Q::FDFilter

Build for a given linear time-invariant descriptor system model `sys = (A-λE,B,C,D)`, 
a fault detection filter object `Q`, as determined with the synthesis functions of FDI filters. 
`ny` and `mu` are the number of measured outputs and the number of control inputs, respectively. 
It is assumed that `B = [By Bu Bv]` and `D = [Dy Du Dv]` are partitioned matrices such that
`By` and `Dy` have `ny` columns, and `Bu` and `Du` have `mu` columns, 
where `By` and `Bu` are the input matrices from the measured outputs `y` and 
control inputs `u`, `Dy` and `Du` are the feedthrough matrices from the measured outputs `y` and 
control inputs `u`. 

The resulting `Q` contains the partitioned system 
`Q.sys = (A-λE,[By Bd],C,[Dy Du])` and the dimensions of the 
partitioned filter input vectors as 
`measured outputs` and `control inputs`,  
can be accessed as the integers contained in `Q.ny` and `Q.mu`, respectively. 
"""
function FDFilter(sys::DescriptorStateSpace, ny::Int, mu::Int)
    return FDFilter{eltype(sys)}(sys, ny, mu)
end    
function *(filter1::FDFilter{T1}, filter2::FDFilter{T2}) where {T1,T2}
    return FDFilter(filter1.sys*filter2.sys, filter2.ny, filter2.mu)
end
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, filter::FDFilter)
    summary(io, filter); println(io)
    display(filter.sys)
    ny = filter.ny
    mu = filter.mu
    if ny+mu > 0
       println(io, "Input groups:")
       println(io,          "Name      Channels")
       ny > 0 && println(io,"outputs   $(1:ny)")
       mu > 0 && println(io,"controls  $(ny+1:ny+mu)")
    end
end
"""
    FDIFilter <: AbstractFDDObject

Type for fault detection and isolation filters resulted as solutions of fault detection and isolation problems.
    
If `filter::FDIFilter` is the fault detection and isolation filter object, 
the underlying `i`-th descriptor system model
can be obtained via `filter.sys[i]` and the dimensions of the partitioned filter input vectors as 
`measured outputs` and `control inputs`,  
can be accessed as the integers contained in `filter.ny` and `filter.mu`, respectively.
The ranges of the indices of output and control inputs 
can be accessed as `filter.outputs` and `filter.controls`, respectively.
"""
struct FDIFilter{T} <: AbstractFDDObject where T 
    sys::Vector{DescriptorStateSpace{T}}
    ny::Int
    mu::Int
    function FDIFilter{T}(sys::Vector{DescriptorStateSpace{T}}, ny::Int, mu::Int) where T 
        N = length(sys)
        sysn = similar(sys,N)
        ny < 0 && error("number of measured outputs must be non-negative")
        mu < 0 && error("number of control inputs must be non-negative")  
        N == 0 &&  (return new{T}(sysn, ny, mu))      
        m = ny+mu
        m > size(sys[1],2) && error("total number of inputs exceeds the number of filter inputs") 
        sysn[1] = sys[1][:,1:m]
        for i = 2:N
            m > size(sys[i],2) && error("total number of inputs exceeds the number of filter inputs") 
            sysn[i] = sys[i][:,1:m]
        end
        new{T}(sysn, ny, mu)
    end
end
"""
    FDIFilter(sys, ny, mu) -> Q::FDIFilter

Build for a vector of linear time-invariant descriptor system models `sys[i] = (Ai-λEi,Bi,Ci,Di)`
with the same number of inputs, a fault detection and isolation filter object `Q`, 
as determined with the synthesis functions of FDI filters. 
`ny` and `mu` are the number of measured outputs and the number of control inputs, respectively. 
It is assumed that each `Bi = [Byi Bui Bvi]` and `Di = [Dyi Dui Dvi]` are partitioned matrices such that
`Byi` and `Dyi` have `ny` columns, and `Bui` and `Dui` have `mu` columns, 
where `Byi` and `Bui` are the input matrices from the measured outputs `y` and 
control inputs `u`, `Dyi` and `Dui` are the feedthrough matrices from the measured outputs `y` and 
control inputs `u`. 

The resulting `Q` contains the vector of partitioned systems
`Q.sys[i] = (Ai-λEi,[Byi Bdi],Ci,[Dyi Dui])` and the dimensions of the 
partitioned filter input vectors as 
`measured outputs` and `control inputs`,  
can be accessed as the integers contained in `Q.ny` and `Q.mu`, respectively. 
"""
function FDIFilter(sys::Vector{DescriptorStateSpace{T}}, ny::Int, mu::Int) where T
    return FDIFilter{T}(sys, ny, mu)
end
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, filter::FDIFilter)
    summary(io, filter); println(io)
    for i = 1:length(filter.sys)
        println(io, "Filter# $i:")
        display(filter.sys[i])
    end
    ny = filter.ny
    mu = filter.mu
    if ny+mu > 0
       println(io, "Input groups:")
       println(io,          "Name      Channels")
       ny > 0 && println(io,"outputs   $(1:ny)")
       mu > 0 && println(io,"controls  $(ny+1:ny+mu)")
   end
end

"""
    FDFilterIF <: AbstractFDDObject

Type for the internal form of fault detection filters resulted as solutions of fault detection problems.
    
If `filter::FDFilterIF` is the fault detection filter internal form object, 
the underlying descriptor system model
can be obtained via `filter.sys` and the dimensions of the control, disturbance, fault, noise and auxiliary
vectors are contained in the integers `filter.mu`, `filter.md`, `filter.mf`, `filter.mw` and `filter.ma`, respectively. 
The ranges of indices of control, disturbance, fault, noise and auxiliary inputs can be accessed as
`filter.controls`, `filter.disturbances`, `filter.faults`, `filter.noise` and `filter.aux`, respectively.
"""
struct FDFilterIF{T} <: AbstractFDDObject where T 
    sys::DescriptorStateSpace{T}
    mu::Int
    md::Int
    mf::Int
    mw::Int
    ma::Int
    function FDFilterIF{T}(sys::DescriptorStateSpace{T}, mu::Int, md::Int, mf::Int, mw::Int, ma::Int)  where T 
        fdimodel_validation(sys, mu, md, mf, mw, ma)
        new{T}(sys[:,1:mu+md+mf+mw+ma], mu, md, mf, mw, ma)
    end
end
FDFilterIF(sys::DescriptorStateSpace,  mu::Int, md::Int, mf::Int, mw::Int, ma::Int)  = 
           FDFilterIF{eltype(sys)}(sys, mu, md, mf, mw, ma)
"""
    FDFilterIF(sys; mu = 0, md = 0, mf = 0, mw = 0, ma = 0, moff = 0) -> R::FDFilterIF

Build for a given linear time-invariant descriptor system model `sys = (A-λE,B,C,D)`, 
a fault detection filter internal form object `R`, as determined with the synthesis functions of FDI filters. 
`mu`, `md`, `mf`, `mw` and `maux` are the dimensions of control, disturbance, fault, noise and auxiliary input vectors, respectively.
It is assumed that `B = [Boff Bu Bd Bf Bw Bv]` and `D = [Doff Du Dd Df Dw Dv]` are partitioned matrices such that
`Boff` and `Doff` have `moff` columns, `Bu` and `Du` have `mu` columns, `Bd` and `Dd` have `md` columns, 
`Bf` and `Df` have `mf` columns,  `Bw` and `Dw` have `mw` columns, and `Bv` and `Dv` have `maux` columns.   

The resulting `R` contains the partitioned system 
`R.sys = (A-λE,[Bu Bd Bf Bw Bv],C,[Du Dd Df Dw Dv])` and 
the dimensions of control, disturbance, fault, noise and 
auxiliary input vectors are contained in 
`R.mu`, `R.md`, `R.mf`, `R.mw` and `R.ma`, respectively.  
"""
function FDFilterIF(sys::DescriptorStateSpace; mu::Int = 0, md::Int = 0, mf::Int = 0, mw::Int = 0, ma::Int = 0, moff::Int = 0) 
    m = moff+mu+md+mf+mw+ma
    m > size(sys,2) && error("number of selected inputs exceeds the number of system inputs")
    return FDFilterIF{eltype(sys)}(sys[:,moff+1:m], mu, md, mf, mw, ma)
end
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, sysf::Union{FDIModel,FDFilterIF})
    summary(io, sysf); println(io)
    display(sysf.sys)
    mu = sysf.mu
    md = sysf.md
    mf = sysf.mf
    mw = sysf.mw
    ma = sysf.ma
    m = mu+md+mf+mw+ma
    if m > 0
       println(io, "Input groups:")
       println(io,          "Name         Channels")
       mu > 0 && println(io,"controls     $(1:mu)")
       md > 0 && println(io,"disturbances $(mu+1:mu+md)")
       mf > 0 && println(io,"faults       $(mu+md+1:mu+md+mf)")
       mw > 0 && println(io,"noise        $(mu+md+mf+1:mu+md+mf+mw)")
       ma > 0 && println(io,"aux          $(mu+md+mf+mw+1:m)")
    end
end
"""
    FDIFilterIF <: AbstractFDDObject

Type for the internal form of fault detection and isolation filters resulted as solutions of fault detection and isolation problems.
    
If `filter::FDIFilterIF` is the fault detection and isolation filter internal form object, 
the underlying `i`-th descriptor system model
can be obtained via `filter.sys[i]` and the dimensions of the control, disturbance, fault, noise and auxiliary
vectors are contained in the integers `filter.mu`, `filter.md`, `filter.mf`, `filter.mw` and `filter.ma`, respectively. 
The ranges of indices of control, disturbance, fault, noise and auxiliary inputs can be accessed as
`filter.controls`, `filter.disturbances`, `filter.faults`, `filter.noise` and `filter.aux`, respectively.
.
"""
struct FDIFilterIF{T} <: AbstractFDDObject where T 
    sys::Vector{DescriptorStateSpace{T}}
    mu::Int
    md::Int
    mf::Int
    mw::Int
    ma::Int
    function FDIFilterIF{T}(sys::Vector{DescriptorStateSpace{T}}, mu::Int, md::Int, mf::Int, mw::Int, ma::Int) where T 
        N = length(sys)
        sysn = similar(sys,N)
        inps = 1:mu+md+mf+mw+ma
        for i = 1:N
            fdimodel_validation(sys[i], mu, md, mf, mw, ma)
            sysn[i] = sys[i][:,inps]
        end
        new{T}(sysn, mu, md, mf, mw, ma)
    end
end
FDIFilterIF(sys::Vector{DescriptorStateSpace{T}}, mu::Int, md::Int, mf::Int, mw::Int, ma::Int) where T  = 
           FDIFilterIF{T}(sys, mu, md, mf, mw, ma)
"""
    FDIFilterIF(sys; mu = 0, md = 0, mf = 0, mw = 0, ma = 0, moff = 0 ) -> R::FDIFilterIF

Build for a vector of linear time-invariant descriptor system models `sys[i] = (Ai-λEi,Bi,Ci,Di)`, 
a fault detection and isolation filter internal form object `R`, as determined with the synthesis functions of FDI filters. 
`mu`, `md`, `mf`, `mw` and `ma` are the dimensions of control, disturbance, fault, noise and auxiliary input vectors, respectively.
It is assumed that each `Bi = [Boffi Bui Bdi Bfi Bwi Bvi]` and `Di = [Doffi Dui Ddi Dfi Dwi Dvi]` are partitioned matrices such that
`Boffi` and `Doffi` have `moff` columns, `Bui` and `Dui` have `mu` columns, `Bdi` and `Ddi` have `md` columns, 
`Bfi` and `Dfi` have `mf` columns,  `Bwi` and `Dwi` have `mw` columns, and `Bvi` and `Dvi` have `ma` columns.   

The resulting `R` contains the vector of partitioned systems 
`R.sys[i] = (A-λE,[Bui Bdi Bfi Bwi Bvi],C,[Dui Ddi Dfi Dwi Dvi])` and 
the dimensions of control, disturbance, fault, noise and auxiliary input vectors are contained in 
`R.mu`, `R.md`, `R.mf`, `R.mw` and `R.ma`, respectively.  
"""
function FDIFilterIF(sys::Vector{DescriptorStateSpace{T}}; mu::Int = 0, md::Int = 0, mf::Int = 0, mw::Int = 0, ma::Int = 0, moff::Int = 0) where T
    mu < 0 && error("number of control inputs must be nonnegative")
    moff < 0 && error("the offset must be nonnegative")
    md < 0 && error("number of disturbance inputs must be nonnegative")
    mf < 0 && error("number of disturbance inputs must be nonnegative")
    mw < 0 && error("number of noise inputs must be nonnegative")
    ma < 0 && error("number of auxiliary inputs must be nonnegative")
    m = moff+mu+md+mf+mw+ma
    N = length(sys)
    sysn = similar(sys,N)
    inps = moff+1:m
    for i = 1:N
        m > size(sys[i],2) && error("the specified total number of inputs exceeds the number of system inputs")
        sysn[i] = sys[i][:,inps]
    end
    return FDIFilterIF{T}(sysn, mu, md, mf, mw, ma)
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, filter::FDIFilterIF)
    summary(io, filter); println(io)
    for i = 1:length(filter.sys)
        println(io, "Filter# $i:")
        display(filter.sys[i])
    end
    mu = filter.mu
    md = filter.md
    mf = filter.mf
    mw = filter.mw
    ma = filter.ma
    m = mu+md+mf+mw+ma
    if m > 0
        println(io, "Input groups:")
        println(io,          "Name         Channels")
        mu > 0 && println(io,"controls     $(1:mu)")
        md > 0 && println(io,"disturbances $(mu+1:mu+md)")
        mf > 0 && println(io,"faults       $(mu+md+1:mu+md+mf)")
        mw > 0 && println(io,"noise        $(mu+md+mf+1:mu+md+mf+mw)")
        ma > 0 && println(io,"aux          $(mu+md+mf+mw+1:m)")
     end
end

"""
    fdimodset(sys; controls, c, disturbances, d, faults, f, fa, faults_sen, fs, noise, n, aux) -> sysf::FDIModel

Build for a given linear time-invariant system model `sys = (A-λE,B,C,D)`, a synthesis model object `sysf::FDIModel` 
to be used in conjunction with the analysis and synthesis functions of FDI filters. 
If `sys` is a vector of system models, then `sysf` results as a vector of synthesis model objects.

The information on the partition of the input components in control, disturbance, fault, noise and auxiliary
inputs can be specified using the following keyword arguments:

`controls = inpu` or `c = inpu` specifies the indices `inpu` of the control inputs (default: void)

`disturbances = inpd` or `d = inpd` specifies the indices `inpd` of the disturbance inputs (default: void)

`faults = inpf` or `f  = inpf` specifies the indices `inpf` of the fault inputs (default: void)

`fa = inpfa ` specifies the indices `inpfa` of control inputs subject to actuator fault (default: void)

`faults_sen = inpfs` or `fs = inpfs` specifies the indices `inpfs` of the system outputs subject to sensor fault inputs (default: void)

`noise = inpn` or  `noise = inpn` specifies the indices `inpn` of the noise inputs (default: void)

`aux = inpa` specifies the indices `inpa` of the auxiliary inputs (default: void)

The indices of inputs or outputs can be specified as integer vectors, integer scalars or integer `UnitRange`s.

The resulting `sysf` contains the partitioned system 
`sysf.sys = (A-λE,[Bu Bd Bf Bw Bv],C,[Du Dd Df Dw Dv])`, where 
`Bu`, `Bd`, `Bf`, `Bw` and `Bv` are the input matrices from the control inputs `u`, disturbance inputs `d`, fault inputs `f`, 
noise inputs `w` and auxiliary inputs `v`, respectively, and `Du`, `Dd`, `Df`, `Dw` and `Dv` are the feedthrough matrices from those inputs.
The dimensions of control, disturbance, fault, noise and auxiliary input vectors are contained in 
`sysm.mu`, `sysm.md`, `sysm.mf`, `sysm.mw` and `sysm.ma`, respectively.  

_Method:_ If `G(λ)` is the `p x m` transfer function matrix of `sys`, then the resulting system `sysf` has an 
equivalent input output form `[Gu(λ) Gd(λ) Gf(λ) Gw(λ) Gv(λ)]`, where the following relations define the component matrices:
`Gu(λ) = G(λ)*Su`, `Gd(λ) = G(λ)*Sd`,  `Gf(λ) = [G(λ)*Sf Ss]`, `Gw(λ) = G(λ)*Sw`,  `Gv(λ) = G(λ)*Sv`, with
the selection matrices `Su`, `Sd`, `Sf`, `Sw` and `Sv` formed from the columns of the `m`-th order identity matrix and 
`Ss` is formed  from the columns of the `p`-th order identity matrix. 
"""
function fdimodset(sys::DescriptorStateSpace{T}; 
                   controls::VRS = Int[], c::VRS = Int[], 
                   disturbances::VRS = Int[], d::VRS = Int[], 
                   faults::VRS = Int[], f::VRS = Int[], fa::VRS = Int[], 
                   faults_sen::VRS = Int[], fs::VRS = Int[],
                   noise::VRS = Int[], n::VRS = Int[], 
                   aux::VRS = Int[]) where T
    p, m = size(sys)
    nx = order(sys)
    inpu = unique([controls;c])
    inpd = unique([disturbances;d])
    isempty(intersect(inpu,inpd)) || error("control and disturbance inputs must be distinct")
    inpf1 = unique([faults;f;fa]) 
    isempty(intersect(inpd,inpf1)) || error("disturbance and fault inputs must be distinct")
    inpfs = unique([faults_sen;fs]) 
    inpw = unique([noise;n])
    isempty(intersect(inpw,inpf1)) || error("noise and fault inputs must be distinct")
    inpaux = unique([aux;Int[]])
    mu = length(inpu)
    md = length(inpd)
    mf1 = length(inpf1)
    mf2 = length(inpfs)
    mf = mf1+mf2
    mw = length(inpw)
    ma = length(inpaux)
    mu == 0 || maximum(inpu) <= m || error("selected index/indices of control inputs larger than the number of system inputs $m")   
    md == 0 || maximum(inpd) <= m || error("selected index/indices of disturbance inputs larger than the number of system inputs $m")   
    mf1 == 0 || maximum(inpf1) <= m || error("selected index/indices of fault inputs larger than the number of system inputs $m")   
    mf2 == 0 || maximum(inpfs) <= p || error("selected index/indices of sensor fault inputs larger than the number of system outputs $p")   
    mw == 0 || maximum(inpw) <= m || error("selected index/indices of noise inputs larger than the number of system inputs $m")   
    ma == 0 || maximum(inpaux) <= m || error("selected index/indices of auxiliary inputs larger than the number of system inputs $m")   
    Dsf = mf2 > 0 ? eye(T, p, p)[:,inpfs] : zeros(T, p, 0)
    
    Be = [sys.B[:,inpu] sys.B[:,inpd] sys.B[:,inpf1] zeros(T, nx, mf2) sys.B[:,inpw] sys.B[:,inpaux]];
    De = [sys.D[:,inpu] sys.D[:,inpd] sys.D[:,inpf1] Dsf sys.D[:,inpw] sys.D[:,inpaux]];
    return FDIModel{T}(dss(sys.A, sys.E, Be, sys.C, De, Ts = sys.Ts); mu, md, mf, mw, ma)          
end
function fdimodset(sys::Vector{DescriptorStateSpace}; kwargs...)
    N = length(sys)
    sysf = similar(Vector{FDIModel},N)
    for i = 1:N
        sysf[i] = fdimodset(sys[i]; kwargs...)
    end   
    return sysf
end
"""
    fdIFeval(sysQ::FDFilter, sysf::FDIModel; minimal = false, atol, atol1 = atol, atol2 = atol, rtol, fast = true) -> sysR::FDFilterIF

Compute the internal form `sysR` of the fault detection filter `sysQ` applied to the synthesis model `sysf`. 
If `sysf` has the partitioned transfer function matrix `G(λ) = [ Gu(λ)  Gd(λ) Gf(λ) Gw(λ) Gv(λ) ]` in accordance with
the partitioned system inputs as `controls`, `disturbances`, `faults`, `noise` and `auxiliary` inputs, respectively,
and `Q(λ) = [ Qy(λ) Qu(λ) ]` is the partitioned transfer function matrix of the fault detection filter `sysQ` 
in accordance with the partitioned filter inputs as `measurable outputs` and `control inputs`, then 
the transfer function matrix `R(λ)` of the resulting internal form `sysR` is given by
     
     R(λ) = | Qy(λ)  Qu(λ) | * | Gu(λ)  Gd(λ) Gf(λ) Gw(λ) Gv(λ) |
                               |  I       0     0     0     0   |

A minimal descriptor realization is computed if `minimal = true` and a possibly non-minimal 
realization is determined if `minimal = false` (default). 

The minimal realization computation relies on pencil manipulation algorithms which 
employ rank determinations based on either the use of 
rank revealing QR-decomposition with column pivoting, if `fast = true`, or the SVD-decomposition.
The rank decision based on the SVD-decomposition is generally more reliable, but the involved computational effort is higher.

If `(Ar-λEr,Br,Cr,Dr)` is the full order descriptor realization of `sysR.sys`, then 
the keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Ar`, `Br`, `Cr`, `Dr`, the absolute tolerance for the nonzero elements of `Er`,  
and the relative tolerance for the nonzero elements of `Ar`, `Br`, `Cr`, `Dr` and `Er`.
The default relative tolerance is `nϵ`, where `ϵ` is the working _machine epsilon_ 
and `n` is the order of the system `sysR`.  
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 

"""
function fdIFeval(Q::FDFilter, sysf::FDIModel; minimal::Bool = false, 
                  atol::Real = 0., atol1::Real = atol, atol2::Real = atol, rtol::Real = 0., fast::Bool = true)
   # compute the internal form as
   #   R = Q * [ Gu Gd Gf Gw Gaux]
   #           [ I  0  0  0  0   ]
   p, m = size(sysf.sys)
   mu = sysf.mu
   (Q.ny == p && Q.mu == mu) || error("filter Q is incompatible with the given system")
   sysR = Q.sys*[sysf.sys;  I zeros(eltype(sysf.sys),mu,m-mu)]
   return FDFilterIF(minimal ? gminreal(sysR; atol1, atol2, rtol, fast) : sysR, 
                 mu, sysf.md, sysf.mf, sysf.mw, sysf.ma) 
   
end
"""
    fdIFeval(sysQ::FDIFilter, sysf::FDIModel; minimal = false, atol, atol1 = atol, atol2 = atol, rtol, fast = true) -> sysR::FDIFilterIF

Compute the internal form `sysR` of the fault detection and isolation filter `sysQ` applied to the synthesis model `sysf`. 
If `sysf` has the partitioned transfer function matrix `G(λ) = [ Gu(λ)  Gd(λ) Gf(λ) Gw(λ) Gv(λ) ]` in accordance with
the partitioned system inputs as `controls`, `disturbances`, `faults`, `noise` and `auxiliary` inputs, respectively,
and `Qi(λ) = [ Qyi(λ) Qui(λ) ]` is the partitioned transfer function matrix of the `i`-th filter `sysQ.sys[i]` 
in accordance with the partitioned filter inputs as `measurable outputs` and `control inputs`, then 
the transfer function matrix `Ri(λ)` of the `i`-th filter in the resulting internal form `sysR.sys[i]` is given by
     
     Ri(λ) = | Qyi(λ)  Qui(λ) | * | Gu(λ)  Gd(λ) Gf(λ) Gw(λ) Gv(λ) |
                                  |  I       0     0     0     0   |

Minimal descriptor realizations are computed if `minimal = true` and a possibly non-minimal 
realization is determined if `minimal = false` (default). 

The minimal realization computation relies on pencil manipulation algorithms which 
employ rank determinations based on either the use of 
rank revealing QR-decomposition with column pivoting, if `fast = true`, or the SVD-decomposition.
The rank decision based on the SVD-decomposition is generally more reliable, but the involved computational effort is higher.

If `(Ari-λEri,Bri,Cri,Dri)` is the full order descriptor realization of `sysR.sys[i]`, then 
the keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Ari`, `Bri`, `Cri`, `Dri`, the absolute tolerance for the nonzero elements of `Eri`,  
and the relative tolerance for the nonzero elements of `Ari`, `Bri`, `Cri`, `Dri` and `Eir`.
The default relative tolerance is `ni*ϵ`, where `ϵ` is the working _machine epsilon_ 
and `ni` is the order of the system `sysR.sys[i]`.  
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 

"""
function fdIFeval(Q::FDIFilter, sysf::FDIModel; minimal::Bool = false, 
                  atol::Real = 0., atol1::Real = atol, atol2::Real = atol, rtol::Real = 0., fast::Bool = true)
   # compute the internal form as
   #   R[i] = Q[i] * [ Gu Gd Gf Gw Gaux]
   #                 [ I  0  0  0  0   ]
   N = length(Q.sys)
   sysR = similar(Q.sys,N)
   p, m = size(sysf.sys)
   mu = sysf.mu
   (Q.ny == p && Q.mu == mu) || error("filter Q is incompatible with the given system")
   for i = 1:N
       sysR[i] = minimal ? gminreal(Q.sys[i]*[sysf.sys;  I zeros(eltype(sysf.sys),mu,m-mu)]; atol1, atol2, rtol, fast) : 
                           Q.sys[i]*[sysf.sys;  I zeros(eltype(sysf.sys),mu,m-mu)]
   end
   return FDIFilterIF{eltype(Q.sys[1])}(sysR, mu, sysf.md, sysf.mf, sysf.mw, sysf.ma) 
end
gbalmr(Q::FDFilter{T}; kwargs...) where T = FDFilter(gbalmr(Q.sys; kwargs...)[1], Q.ny, Q.mu)
gminreal(Q::FDFilter{T}; kwargs...) where T = FDFilter(gminreal(Q.sys; kwargs...), Q.ny, Q.mu)
function gbalmr(Q::FDIFilter{T}; kwargs...) where T 
    for i = 1:length(Q.sys)
        Q.sys[i] = gbalmr(Q.sys[i]; kwargs...)[1]
    end
    return Q
end
function gminreal(Q::FDIFilter{T}; kwargs...) where T 
    for i = 1:length(Q.sys)
        Q.sys[i] = gminreal(Q.sys[i]; kwargs...)
    end
    return Q
end
gpole(Q::FDFilter{T}; kwargs...) where T = gpole(Q.sys;  kwargs...)
gpole(Q::FDIFilter{T}; kwargs...) where T = gpole.(Q.sys;  kwargs...)

# left multiplication of FD objects with descriptor system objects
function *(sys::DescriptorStateSpace{T1}, sysr::FDFilterIF{T2}) where {T1,T2}
    return FDFilterIF(sys*sysr.sys, sysr.mu, sysr.md, sysr.mf, sysr.mw, sysr.ma) 
end
function *(sys::DescriptorStateSpace{T1}, sysr::FDIModel{T2}) where {T1,T2}
    return FDIModel(sys*sysr.sys, sysr.mu, sysr.md, sysr.mf, sysr.mw, sysr.ma) 
end
function *(sys::Vector{DescriptorStateSpace{T1}}, sysr::FDIFilterIF{T2}) where {T1,T2}
    return FDIFilterIF(sys .* sysr.sys, sysr.mu, sysr.md, sysr.mf, sysr.mw, sysr.ma) 
end
function Base.getproperty(sys::Union{FDIModel,FDFilterIF,FDIFilterIF}, d::Symbol)  
    if d === :controls
        return 1:sys.mu
    elseif d === :disturbances
        return sys.mu+1:sys.mu+sys.md
    elseif d === :faults
        return sys.mu+sys.md+1:sys.mu+sys.md+sys.mf
    elseif d === :noise
        return sys.mu+sys.md+sys.mf+1:sys.mu+sys.md+sys.mf+sys.mw
    elseif d === :aux
        return sys.mu+sys.md+sys.mf+sys.mw+1:sys.mu+sys.md+sys.mf+sys.mw+sys.ma
    else
        getfield(sys, d)
    end
end
Base.propertynames(sys::Union{FDIModel,FDFilterIF,FDIFilterIF}) =
    (fieldnames(typeof(sys))..., :controls, :disturbances, :faults, :noise, :aux)
function Base.getproperty(sys::Union{FDFilter,FDIFilter}, d::Symbol)  
    if d === :outputs
        return 1:sys.ny
    elseif d === :controls
        return sys.ny+1:sys.ny+sys.mu
    else
        getfield(sys, d)
    end
end
Base.propertynames(sys::Union{FDFilter,FDIFilter}) =
    (fieldnames(typeof(sys))..., :outputs, :controls)
