"""
    FDIModel <: AbstractDescriptorStateSpace

Type for synthesis models employed to solve fault detection and isolation problems.
    
If `sysf::FDIModel` is the synthesis model object, the underlying descriptor system model
can be obtained via `sysf.sys` and the indices of control, disturbance, fault, noise 
and auxiliary inputs can be accessed as the integer vectors 
contained in `sysf.controls`, `sysf.disturbances`, `sysf.faults`, 
`sysf.noise` and `sysf.aux`, respectively.
"""
struct FDIModel{T} <: AbstractDescriptorStateSpace where T 
    sys::DescriptorStateSpace{T}
    controls::Vector{Int}
    disturbances::Vector{Int}
    faults::Vector{Int}
    noise::Vector{Int}
    aux::Vector{Int}
    function FDIModel{T}(sys::DescriptorStateSpace{T}; controls::Vector{Int} = Int[], disturbances::Vector{Int} = Int[], 
                         faults::Vector{Int} = Int[], noise::Vector{Int} = Int[], aux::Vector{Int} = Int[]) where T 
        fdimodel_validation(sys, controls, disturbances, faults, noise, aux)
        new{T}(sys, controls, disturbances, faults, noise, aux)
    end
end
function fdimodel_validation(sys::DescriptorStateSpace{T}, controls::Vector{Int}, disturbances::Vector{Int}, 
                             faults::Vector{Int}, noise::Vector{Int}, aux::Vector{Int}) where T
    m = size(sys,2) 
    isempty(intersect(controls,disturbances)) || error("control and disturbance inputs must be distinct")
    isempty(intersect(disturbances,faults)) || error("disturbance and fault inputs must be distinct")
    isempty(intersect(noise,faults)) || error("noise and fault inputs must be distinct")
     
    any(controls .> m) && error("control inputs indices larger than the number of system inputs $m") 
    any(disturbances .> m) && error("disturbance inputs indices larger than the number of system inputs $m") 
    any(faults .> m) && error("fault inputs indices larger than the number of system inputs $m") 
    any(noise .> m) && error("noise inputs indices larger than the number of system inputs $m") 
    any(aux .> m) && error("auxiliary inputs indices larger than the number of system inputs $m") 
end
"""
    FDFilter <: AbstractDescriptorStateSpace

Type for fault detection filters resulted as solutions of fault detection problems.
    
If `filter::FDFilter` is the fault detection filter object, the underlying descriptor system model
can be obtained via `filter.sys` and the indices of output and control inputs 
can be accessed as the integer vectors 
contained in `filter.outputs` and `filter.controls`, respectively.
"""
struct FDFilter{T} <: AbstractDescriptorStateSpace where T 
    sys::DescriptorStateSpace{T}
    outputs::Vector{Int}
    controls::Vector{Int}
    function FDFilter{T}(sys::DescriptorStateSpace{T}, outputs::Vector{Int}, controls::Vector{Int}) where T 
        inpout, inpu = fdfilter_validation(sys, outputs, controls)
        mout = length(inpout)
        mu = length(inpu)
        new{T}(sys[:,[inpout; inpu]], Vector(1:mout), mout .+ Vector(1:mu))
    end
end
function fdfilter_validation(sys::DescriptorStateSpace{T}, outputs::Vector{Int}, controls::Vector{Int}) where T
    m = size(sys,2) 
    isempty(intersect(controls,outputs)) || error("output and control inputs must be distinct")

    inpout = unique(outputs) 
    inpu = unique(controls)
   
    any(inpout .> m) && error("output inputs indices larger than the number of system inputs $m") 
    any(inpu .> m) && error("control inputs indices larger than the number of system inputs $m") 

    return inpout, inpu
end
FDFilter(sys::DescriptorStateSpace, outputs::Vector{Int}, controls::Vector{Int})  = FDFilter{eltype(sys)}(sys, outputs, controls)
"""
    FDFilter(sys, p, mu) -> Q::FDFilter

Build for a given linear time-invariant descriptor system model `sys = (A-λE,B,C,D)`, 
a fault detection filter object `Q`, as determined with the synthesis functions of FDI filters. 
`p` and `mu` are the number of measured outputs and the number of control inputs, respectively. 
It is assumed that `B = [By Bu Bv]` and `D = [Dy Du Dv]` are partitioned matrices such that
`By` and `Dy` have `p` columns, and `Bu` and `Du` have `mu` columns, 
where `By` and `Bu` are the input matrices from the measured outputs `y` and 
control inputs `u`, `Dy` and `Du` are the feedthrough matrices from the measured outputs `y` and 
control inputs `u`. 

The resulting `Q` contains the partitioned system 
`Q.sys = (A-λE,[By Bd],C,[Dy Du])` and the indices of inputs corresponding 
to the measured outputs and control inputs are contained in the associated 
integer vectors `Q.outputs` and `Q.controls`. 
"""
function FDFilter(sys::DescriptorStateSpace{T}, p::Int, mu::Int) where T 
    m = size(sys,2) 
    p < 0 && error("number of measured outputs must be non-negative")
    mu < 0 && error("number of control inputs must be non-negative")
    p+mu > m && error("number of measured outputs and control inputs exceeds the number of filter inputs $m")
    return FDFilter{T}(sys[:,1:p+mu], Vector(1:p), Vector(p+1:p+mu))          
end
function FDFilter(sys::DescriptorStateSpace{T}; outputs::VRS = Int[], controls::VRS = Int[]) where T 
    return FDFilter{T}(sys, vec([outputs; Int[]]), vec([controls; Int[]]))          
end
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, filter::FDFilter)
    summary(io, filter); println(io)
    display(filter.sys)
    p = length(filter.outputs)
    mu = length(filter.controls)
    if p+mu > 0
       println(io, "Input groups:")
       println(io,          "Name      Channels")
       p > 0 &&  println(io,"outputs   $(filter.outputs)")
       mu > 0 && println(io,"controls  $(filter.controls)")
   end
end
"""
    FDFilterIF <: AbstractDescriptorStateSpace

Type for the internal form of fault detection filters resulted as solutions of fault detection problems.
    
If `filter::FDFilterIF` is the fault detection filter internal form object, 
the underlying descriptor system model
can be obtained via `filter.sys` and the indices of control, disturbance, fault, noise 
and auxiliary inputs can be accessed as the integer vectors 
contained in `filter.controls`, `filter.disturbances`, `filter.faults`, 
`filter.noise` and `filter.aux`, respectively.
"""
struct FDFilterIF{T} <: AbstractDescriptorStateSpace where T 
    sys::DescriptorStateSpace{T}
    controls::Vector{Int}
    disturbances::Vector{Int}
    faults::Vector{Int}
    noise::Vector{Int}
    aux::Vector{Int}
    function FDFilterIF{T}(sys::DescriptorStateSpace{T}, controls::Vector{Int}, disturbances::Vector{Int}, 
                           faults::Vector{Int}, noise::Vector{Int}, aux::Vector{Int}) where T 
        inpu, inpd, inpf, inpw, inpaux = fdfilterIF_validation(sys, controls, disturbances, faults, noise, aux)
        mu = length(inpu)
        md = length(inpd)
        mf = length(inpf)
        mw = length(inpw)
        maux = length(inpaux)
        new{T}(sys[:,[inpu; inpd; inpf; inpw; inpaux]], Vector(1:mu), mu .+ Vector(1:md) , (mu+md) .+ Vector(1:mf), 
               (mu+md+mf) .+ Vector(1:mw), (mu+md+mf+mw) .+ Vector(1:maux))
    end
end
function fdfilterIF_validation(sys::DescriptorStateSpace{T}, controls::Vector{Int}, disturbances::Vector{Int}, 
                             faults::Vector{Int}, noise::Vector{Int}, aux::Vector{Int}) where T
    m = size(sys,2) 
    isempty(intersect(controls,disturbances)) || error("control and disturbance inputs must be distinct")
    isempty(intersect(disturbances,faults)) || error("disturbance and fault inputs must be distinct")
    isempty(intersect(faults,noise)) || error("fault and noise inputs must be distinct")
    isempty(intersect(noise,aux)) || error("noise and aux inputs must be distinct")

    inpu = unique(controls); 
    inpd = unique(disturbances); 
    inpf = unique(faults); 
    inpw = unique(noise); 
    inpaux = unique(aux); 

    isempty(intersect(inpu,inpd)) || error("control and disturbance inputs must be distinct")
    isempty(intersect(inpd,inpf)) || error("disturbance and fault inputs must be distinct")
    isempty(intersect(inpw,inpf)) || error("noise and fault inputs must be distinct")
    isempty(intersect(inpw,inpaux)) || error("noise and aux inputs must be distinct")  
   
    any(inpu .> m) && error("control inputs indices larger than the number of system inputs $m") 
    any(inpd .> m) && error("disturbance inputs indices larger than the number of system inputs $m") 
    any(inpf .> m) && error("fault inputs indices larger than the number of system inputs $m") 
    any(inpw .> m) && error("noise inputs indices larger than the number of system inputs $m") 
    any(inpaux .> m) && error("auxiliary inputs indices larger than the number of system inputs $m") 

    return inpu, inpd, inpf, inpw, inpaux
end
FDFilterIF(sys::DescriptorStateSpace, controls::Vector{Int}, disturbances::Vector{Int}, faults::Vector{Int}, noise::Vector{Int}, aux::Vector{Int})  = 
           FDFilterIF{eltype(sys)}(sys, controls, disturbances, faults, noise, aux)
"""
    FDFilterIF(sys, mu, md, mf, mw = 0, maux = 0; moff = 0 ) -> R::FDFilterIF

Build for a given linear time-invariant descriptor system model `sys = (A-λE,B,C,D)`, 
a fault detection filter internal form object `R`, as determined with the synthesis functions of FDI filters. 
`mu`, `md`, `mf`, `mw` and `maux` are the dimensions of control, disturbance, fault, noise and auxiliary input vectors, respectively.
It is assumed that `B = [Boff Bu Bd Bf Bw Bv]` and `D = [Doff Du Dd Df Dw Dv]` are partitioned matrices such that
`Boff` and `Doff` have `moff` columns, `Bu` and `Du` have `mu` columns, `Bd` and `Dd` have `md` columns, 
`Bf` and `Df` have `mf` columns,  `Bw` and `Dw` have `mw` columns, and `Bv` and `Dv` have `maux` columns.   

The resulting `R` contains the partitioned system 
`R.sys = (A-λE,[Bu Bd Bf Bw Bv],C,[Du Dd Df Dw Dv])` and the indices of inputs corresponding 
to the control, disturbance, fault, noise and auxiliary inputs are contained in the associated 
integer vectors `R.controls`, `R.disturbances`, `R.faults`, `R.noise` and `R.aux`.  
"""
function FDFilterIF(sys::DescriptorStateSpace{T}, mu::Int, md::Int, mf::Int, mw::Int = 0, maux::Int = 0; moff::Int = 0) where T
    m = moff+mu+md+mf+mw+maux
    m > size(sys,2) && error("number of selected inputs exceeds the number of system inputs")
    return FDFilterIF{T}(sys, moff .+ Vector(1:mu), (moff+mu) .+ Vector(1:md), (moff+mu+md) .+ Vector(1:mf), 
                         (moff+mu+md+mf) .+  Vector(1:mw), (moff+mu+md+mf+mw) .+ Vector(1:maux))
end
function FDFilterIF(sys::DescriptorStateSpace{T}; controls::VRS = Int[], disturbances::VRS = Int[], 
                faults::VRS = Int[], noise::VRS = Int[], aux::VRS = Int[]) where T 
    return FDFilterIF{T}(sys, vec([controls; Int[]]), vec([disturbances; Int[]]), 
                         vec([faults; Int[]]), vec([noise; Int[]]), vec([aux; Int[]]))
end
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, sysf::Union{FDIModel,FDFilterIF})
    summary(io, sysf); println(io)
    display(sysf.sys)
    mu = length(sysf.controls)
    md = length(sysf.disturbances)
    mf = length(sysf.faults)
    mw = length(sysf.noise)
    maux = length(sysf.aux)
    if mu+md+mf+mw+maux > 0
       println(io, "Input groups:")
       println(io, "Name         Channels")
       mu > 0 && println(io,"controls     $(sysf.controls)")
       md > 0 && println(io,"disturbances $(sysf.disturbances)")
       mf > 0 && println(io,"faults       $(sysf.faults)")
       mw > 0 && println(io,"noise        $(sysf.noise)")
       maux > 0 && println(io,"aux          $(sysf.aux)")
    end
end
 
"""
    fdimodset(sys; controls, c, disturbances, d, faults, f, fa, faults_sen, fs, noise, n, aux) -> sysf::FDIModel

Build for a given linear time-invariant system model `sys = (A-λE,B,C,D)`, a synthesis model object `sysf` 
to be used in conjunction with the analysis and synthesis functions of FDI filters. 
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
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysf.controls`, `sysf.disturbances`, `sysf.faults`, `sysf.noise` and `sysf.aux`, respectively.

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
    inpu = sort(unique([controls;c]))
    inpd = sort(unique([disturbances;d]))
    isempty(intersect(inpu,inpd)) || error("control and disturbance inputs must be distinct")
    inpf1 = sort(unique([faults;f;fa])) 
    isempty(intersect(inpd,inpf1)) || error("disturbance and fault inputs must be distinct")
    inpfs = sort(unique([faults_sen;fs])) 
    inpw = sort(unique([noise;n]))
    isempty(intersect(inpw,inpf1)) || error("noise and fault inputs must be distinct")
    inpaux = sort(unique([aux;Int[]]))
    mu = length(inpu)
    md = length(inpd)
    mf1 = length(inpf1)
    mf2 = length(inpfs)
    mf = mf1+mf2
    mw = length(inpw)
    maux = length(inpaux)
    mu == 0 || inpu[mu] <= m || error("selected index/indices of control inputs larger than the number of system inputs $m")   
    md == 0 || inpd[md] <= m || error("selected index/indices of disturbance inputs larger than the number of system inputs $m")   
    mf1 == 0 || inpf1[mf1] <= m || error("selected index/indices of fault inputs larger than the number of system inputs $m")   
    mf2 == 0 || inpfs[mf2] <= p || error("selected index/indices of sensor fault inputs larger than the number of system outputs $p")   
    mw == 0 || inpw[mw] <= m || error("selected index/indices of noise inputs larger than the number of system inputs $m")   
    maux == 0 || inpaux[maux] <= m || error("selected index/indices of auxiliary inputs larger than the number of system inputs $m")   
    Dsf = mf2 > 0 ? eye(T, p, p)[:,inpfs] : zeros(T, p, 0)
    
    Be = [sys.B[:,inpu] sys.B[:,inpd] sys.B[:,inpf1] zeros(T, nx, mf2) sys.B[:,inpw] sys.B[:,inpaux]];
    De = [sys.D[:,inpu] sys.D[:,inpd] sys.D[:,inpf1] Dsf sys.D[:,inpw] sys.D[:,inpaux]];
    return FDIModel{T}(dss(sys.A, sys.E, Be, sys.C, De, Ts = sys.Ts),  
            controls = Vector(1:mu), disturbances = mu .+ Vector(1:md) , faults = (mu+md) .+ Vector(1:mf), 
            noise = (mu+md+mf) .+ Vector(1:mw), aux = (mu+md+mf+mw) .+ Vector(1:maux))          
end
