"""
    FDIModel <: AbstractFDDObject

Type for synthesis models employed to solve fault detection and isolation problems.
    
If `sysf::FDIModel` is the synthesis model object, the underlying descriptor system model
can be obtained via `sysf.sys` and the indices of control, disturbance, fault, noise 
and auxiliary inputs can be accessed as the integer vectors 
contained in `sysf.controls`, `sysf.disturbances`, `sysf.faults`, 
`sysf.noise` and `sysf.aux`, respectively.
"""
struct FDIModel{T} <: AbstractFDDObject where T 
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
    FDFilter <: AbstractFDDObject

Type for fault detection filters resulted as solutions of fault detection problems.
    
If `filter::FDFilter` is the fault detection filter object, the underlying descriptor system model
can be obtained via `filter.sys` and the indices of output and control inputs 
can be accessed as the integer vectors 
contained in `filter.outputs` and `filter.controls`, respectively.
"""
struct FDFilter{T} <: AbstractFDDObject where T 
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
function FDFilter(sys::DescriptorStateSpace, p::Int, mu::Int)
    m = size(sys,2) 
    p < 0 && error("number of measured outputs must be non-negative")
    mu < 0 && error("number of control inputs must be non-negative")
    p+mu > m && error("number of measured outputs and control inputs exceeds the number of filter inputs $m")
    return FDFilter{eltype(sys)}(sys[:,1:p+mu], Vector(1:p), Vector(p+1:p+mu))          
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
    FDIFilter <: Vector{AbstractFDDObject}

Type for fault detection and isolation filters resulted as solutions of fault detection and isolation problems.
    
If `filter::FDIFilter` is the fault detection and isolation filter object, 
the underlying i-th descriptor system models
can be obtained via `filter.sys[i]` and the indices of output and control inputs 
can be accessed as the integer vectors 
contained in `filter.outputs` and `filter.controls`, respectively.
"""
struct FDIFilter{T} <: AbstractFDDObject where T 
    sys::Vector{DescriptorStateSpace{T}}
    outputs::Vector{Int}
    controls::Vector{Int}
    function FDIFilter{T}(sys::Vector{DescriptorStateSpace{T}}, outputs::Vector{Int}, controls::Vector{Int}) where T 
        N = length(sys)
        sysn = similar(sys,N)
        inpout, inpu = fdfilter_validation(sys[1], outputs, controls)
        sysn[1] = sys[1][:,[inpout; inpu]]
        #sys[1] = sys[1][:,[inpout; inpu]]
        for i = 2:N
            inpouti, inpui = fdfilter_validation(sys[i], outputs, controls)
            (inpout == inpouti && inpu == inpui) || error("all component filters must have the same inputs from output and control inputs ")
            sysn[i] = sys[i][:,[inpout; inpu]]
            #sys[i] = sys[i][:,[inpout; inpu]]
        end
        mout = length(inpout)
        mu = length(inpu)
        new{T}(sysn, Vector(1:mout), mout .+ Vector(1:mu))
        #new{T}(sys, Vector(1:mout), mout .+ Vector(1:mu))
    end
end
"""
    FDIFilter(sys, p, mu) -> Q::FDIFilter

Build for a vector of linear time-invariant descriptor system models `sys[i] = (Ai-λEi,Bi,Ci,Di)`
with the same number of inputs, a fault detection and isolation filter object `Q`, 
as determined with the synthesis functions of FDI filters. 
`p` and `mu` are the number of measured outputs and the number of control inputs, respectively. 
It is assumed that each `Bi = [Byi Bui Bvi]` and `Di = [Dyi Dui Dvi]` are partitioned matrices such that
`Byi` and `Dyi` have `p` columns, and `Bui` and `Dui` have `mu` columns, 
where `Byi` and `Bui` are the input matrices from the measured outputs `y` and 
control inputs `u`, `Dyi` and `Dui` are the feedthrough matrices from the measured outputs `y` and 
control inputs `u`. 

The resulting `Q` contains the vector of partitioned systems
`Q.sys[i] = (Ai-λEi,[Byi Bdi],Ci,[Dyi Dui])` and the indices of inputs corresponding 
to the measured outputs and control inputs are contained in the associated 
integer vectors `Q.outputs` and `Q.controls`. 
"""
function FDIFilter(sys::Vector{DescriptorStateSpace{T}}, p::Int, mu::Int) where T
    m = size(sys[1],2) 
    p < 0 && error("number of measured outputs must be non-negative")
    mu < 0 && error("number of control inputs must be non-negative")
    p+mu > m && error("number of measured outputs and control inputs exceeds the number of filter inputs $m")
    return FDIFilter{T}(sys, Vector(1:p), Vector(p+1:p+mu))          
end
function FDIFilter(sys::Vector{DescriptorStateSpace}, p::Int, mu::Int)
    m = size(sys[1],2) 
    p < 0 && error("number of measured outputs must be non-negative")
    mu < 0 && error("number of control inputs must be non-negative")
    p+mu > m && error("number of measured outputs and control inputs exceeds the number of filter inputs $m")
    return FDIFilter{eltype(sys[1])}(sys, Vector(1:p), Vector(p+1:p+mu))          
end

function FDIFilter(sys::Vector{DescriptorStateSpace{T}}; outputs::VRS = Int[], controls::VRS = Int[]) where T 
    return FDIFilter{T}(sys, vec([outputs; Int[]]), vec([controls; Int[]]))          
end
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, filter::FDIFilter)
    summary(io, filter); println(io)
    for i = 1:length(filter.sys)
        println(io, "Filter number $i:")
        display(filter.sys[i])
    end
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
    FDFilterIF <: AbstractFDDObject

Type for the internal form of fault detection filters resulted as solutions of fault detection problems.
    
If `filter::FDFilterIF` is the fault detection filter internal form object, 
the underlying descriptor system model
can be obtained via `filter.sys` and the indices of control, disturbance, fault, noise 
and auxiliary inputs can be accessed as the integer vectors 
contained in `filter.controls`, `filter.disturbances`, `filter.faults`, 
`filter.noise` and `filter.aux`, respectively.
"""
struct FDFilterIF{T} <: AbstractFDDObject where T 
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
    FDIFilterIF <: AbstractFDDObject

Type for the internal form of fault detection and isolation filters resulted as solutions of fault detection and isolation problems.
    
If `filter::FDIFilterIF` is the fault detection and isolation filter internal form object, 
the underlying `i`-th descriptor system models
can be obtained via `filter.sys[i]` and the indices of control, disturbance, fault, noise 
and auxiliary inputs can be accessed as the integer vectors 
contained in `filter.controls`, `filter.disturbances`, `filter.faults`, 
`filter.noise` and `filter.aux`, respectively.
"""
struct FDIFilterIF{T} <: AbstractFDDObject where T 
    sys::Vector{DescriptorStateSpace{T}}
    controls::Vector{Int}
    disturbances::Vector{Int}
    faults::Vector{Int}
    noise::Vector{Int}
    aux::Vector{Int}
    function FDIFilterIF{T}(sys::Vector{DescriptorStateSpace{T}}, controls::Vector{Int}, disturbances::Vector{Int}, 
                           faults::Vector{Int}, noise::Vector{Int}, aux::Vector{Int}) where T 
        N = length(sys)
        sysn = similar(sys,N)
        inpu, inpd, inpf, inpw, inpaux = fdfilterIF_validation(sys[1], controls, disturbances, faults, noise, aux)
        sysn[1] = sys[1][:,[inpu; inpd; inpf; inpw; inpaux]]
        #sys[1] = sys[1][:,[inpu; inpd; inpf; inpw; inpaux]]
        for i = 2:N
            inpui, inpdi, inpfi, inpwi, inpauxi = fdfilterIF_validation(sys[i], controls, disturbances, faults, noise, aux)
            (inpu == inpui && inpd == inpdi && inpf == inpfi && inpw == inpwi && inpaux == inpauxi) ||
                 error("all component filters must have the same inputs groups")
            sysn[i] = sys[i][:,[inpu; inpd; inpf; inpw; inpaux]]
            #sysn[i] = sys[i][:,[inpu; inpd; inpf; inpw; inpaux]]
            #sys[i] = sys[i][:,[inpu; inpd; inpf; inpw; inpaux]]
        end
        mu = length(inpu)
        md = length(inpd)
        mf = length(inpf)
        mw = length(inpw)
        maux = length(inpaux)
        new{T}(sysn, Vector(1:mu), mu .+ Vector(1:md) , (mu+md) .+ Vector(1:mf), 
               (mu+md+mf) .+ Vector(1:mw), (mu+md+mf+mw) .+ Vector(1:maux))
        # new{T}(sys, Vector(1:mu), mu .+ Vector(1:md) , (mu+md) .+ Vector(1:mf), 
        #        (mu+md+mf) .+ Vector(1:mw), (mu+md+mf+mw) .+ Vector(1:maux))
    end
end
function FDIFilterIF(sys::Vector{DescriptorStateSpace{T}}; controls::VRS = Int[], disturbances::VRS = Int[], 
                     faults::VRS = Int[], noise::VRS = Int[], aux::VRS = Int[]) where T 
    return FDIFilterIF{T}(sys, vec([controls; Int[]]), vec([disturbances; Int[]]), 
                          vec([faults; Int[]]), vec([noise; Int[]]), vec([aux; Int[]]))
end
"""
    FDIFilterIF(sys, mu, md, mf, mw = 0, maux = 0; moff = 0 ) -> R::FDIFilterIF

Build for a a vector of linear time-invariant descriptor system models `sys[i] = (Ai-λEi,Bi,Ci,Di)` 
with the same number of inputs, , 
a fault detection and isolation filter internal form object `R`, as determined with the synthesis functions of FDI filters. 
`mu`, `md`, `mf`, `mw` and `maux` are the dimensions of control, disturbance, fault, noise and auxiliary input vectors, respectively.
It is assumed that each `Bi = [Boffi Bui Bdi Bfi Bwi Bvi]` and `Di = [Doffi Dui Ddi Dfi Dwi Dvi]` are partitioned matrices such that
`Boffi` and `Doffi` have `moff` columns, `Bui` and `Dui` have `mu` columns, `Bdi` and `Ddi` have `md` columns, 
`Bfi` and `Dfi` have `mf` columns,  `Bwi` and `Dwi` have `mw` columns, and `Bvi` and `Dvi` have `maux` columns.   

The resulting `R` contains the vector of partitioned systems 
`R.sys[i] = (A-λE,[Bui Bdi Bfi Bwi Bvi],C,[Dui Ddi Dfi Dwi Dvi])` and the indices of inputs corresponding 
to the control, disturbance, fault, noise and auxiliary inputs are contained in the associated 
integer vectors `R.controls`, `R.disturbances`, `R.faults`, `R.noise` and `R.aux`.  
"""
function FDIFilterIF(sys::Vector{DescriptorStateSpace{T}}, mu::Int, md::Int, mf::Int, mw::Int = 0, maux::Int = 0; moff::Int = 0) where T
    m = moff+mu+md+mf+mw+maux
    m > size(sys[1],2) && error("number of selected inputs exceeds the number of system inputs")
    return FDIFilterIF{T}(sys, moff .+ Vector(1:mu), (moff+mu) .+ Vector(1:md), (moff+mu+md) .+ Vector(1:mf), 
                         (moff+mu+md+mf) .+  Vector(1:mw), (moff+mu+md+mf+mw) .+ Vector(1:maux))
end

function FDIFilterIF(sys::Vector{DescriptorStateSpace}, mu::Int, md::Int, mf::Int, mw::Int = 0, maux::Int = 0; moff::Int = 0) 
    m = moff+mu+md+mf+mw+maux
    m > size(sys[1],2) && error("number of selected inputs exceeds the number of system inputs")
    return FDIFilterIF{eltype(sys[1])}(sys, moff .+ Vector(1:mu), (moff+mu) .+ Vector(1:md), (moff+mu+md) .+ Vector(1:mf), 
                         (moff+mu+md+mf) .+  Vector(1:mw), (moff+mu+md+mf+mw) .+ Vector(1:maux))
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, filter::FDIFilterIF)
    summary(io, filter); println(io)
    for i = 1:length(filter.sys)
        println(io, "Filter number $i:")
        display(filter.sys[i])
    end
    mu = length(filter.controls)
    md = length(filter.disturbances)
    mf = length(filter.faults)
    mw = length(filter.noise)
    maux = length(filter.aux)
    if mu+md+mf+mw+maux > 0
       println(io, "Input groups:")
       println(io, "Name         Channels")
       mu > 0 && println(io,"controls     $(filter.controls)")
       md > 0 && println(io,"disturbances $(filter.disturbances)")
       mf > 0 && println(io,"faults       $(filter.faults)")
       mw > 0 && println(io,"noise        $(filter.noise)")
       maux > 0 && println(io,"aux          $(filter.aux)")
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
"""
    fdIFeval(sysQ::FDFilter, sysf::FDIModel; minimal = false, atol, atol1 = atol, atol2 = atol, rtol, fast = true) -> sysR::FDFilterIF

Compute the internal form `sysR` of the fault detection filter `sysQ` applied to the synthesis model `sysf`. 
If `sysf` has the partitioned transfer function matrix `G(λ) = [ Gu(λ)  Gd(λ) Gf(λ) Gw(λ) Gv(λ) ]` in accordance with
the partitioned system inputs as `controls`, `disturbances`, `faults`, `noise` and `auxiliary` inputs, respectively,
and `Q(λ) = [ Qy(λ) Qu(λ) ]` is the partitioned transfer function matrix of the fault detection filter `sysQ` 
in accordance with the partitioned filter inputs as `outputs` and `controls`, then 
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
   mu = length(sysf.controls)
   (length(Q.outputs) == p && length(Q.controls) == mu) || error("filter Q is incompatible with the given system")
   # return fdRset(gminreal(Q.sys*[sysf.sys;  I zeros(eltype(sysf.sys),mu,m-mu)]; atol1, atol2, rtol); 
   #               sysf.controls, sysf.disturbances, sysf.faults, sysf.noise, sysf.aux) 
   sysR = Q.sys*[sysf.sys;  I zeros(eltype(sysf.sys),mu,m-mu)]
   return FDFilterIF(minimal ? gminreal(sysR; atol1, atol2, rtol, fast) : sysR, 
                 sysf.controls, sysf.disturbances, sysf.faults, sysf.noise, sysf.aux) 
   
end
"""
    fdIFeval(sysQ::FDIFilter, sysf::FDIModel; minimal = false, atol, atol1 = atol, atol2 = atol, rtol, fast = true) -> sysR::FDIFilterIF

Compute the internal form `sysR` of the fault detection and isolation filter `sysQ` applied to the synthesis model `sysf`. 
If `sysf` has the partitioned transfer function matrix `G(λ) = [ Gu(λ)  Gd(λ) Gf(λ) Gw(λ) Gv(λ) ]` in accordance with
the partitioned system inputs as `controls`, `disturbances`, `faults`, `noise` and `auxiliary` inputs, respectively,
and `Qi(λ) = [ Qyi(λ) Qui(λ) ]` is the partitioned transfer function matrix of the i-th filter `sysQ.sys[i]` 
in accordance with the partitioned filter inputs as `outputs` and `controls`, then 
the transfer function matrix `Ri(λ)` of the i-th filter in the resulting internal form `sysR.sys[i]` is given by
     
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
   mu = length(sysf.controls)
   (length(Q.outputs) == p && length(Q.controls) == mu) || error("filter Q is incompatible with the given system")
   for i = 1:N
       sysR[i] = minimal ? gminreal(Q.sys[i]*[sysf.sys;  I zeros(eltype(sysf.sys),mu,m-mu)]; atol1, atol2, rtol, fast) : 
                           Q.sys[i]*[sysf.sys;  I zeros(eltype(sysf.sys),mu,m-mu)]
   end
   return FDIFilterIF{eltype(Q.sys[1])}(sysR, sysf.controls, sysf.disturbances, sysf.faults, sysf.noise, sysf.aux) 
end
gbalmr(Q::FDFilter; kwargs...) = FDFilter(gbalmr(Q.sys; kwargs...)[1], Q.outputs, Q.controls)
function gbalmr(Q::FDIFilter; kwargs...)
    for i = 1:length(Q.sys)
        Q.sys[i] = gbalmr(Q.sys[i]; kwargs...)[1]
    end
    return Q
end
