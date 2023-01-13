module CS2_2
using FaultDetectionTools, DescriptorSystems, LinearAlgebra, 
      JLD2, Polynomials, Test

## CS2_2  - Case-study example: Monitoring air data sensor faults
#           Robust least order LPV synthesis
println("Case study CS2_2 with Fig8.8")

## Part 1 - Model setup
# load matrices of the aircraft multiple-model SYSACM 
# and percentages of mass variations massi

cd(joinpath(pkgdir(FaultDetectionTools), "Examples"))
SYSACSM, massi = load("cs2data.jld2", "SYSACSM","massi")

N = size(SYSACSM,1)        # number of models
p, m = size(SYSACSM[1].D)  # output and input dimensions
n = order(SYSACSM[1])      # state dimension
nom = 7;                   # index of nominal system 
atol = 1.e-6               # absolute tolerance
maxdeg = 4                 # maximum interpolation order
T = Float64

# set sensor indices and set system dimensions of inputs
#     [ AoA VCAS ]
sen = [  2,    4  ]; 
md = 2; mu = m-md; mf = length(sen); 

# build minimal realizations of AC-models with sensor faults
idm = eye(p);
syssen = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    syssen[i] = [ SYSACSM[i] idm[:,sen]]
end

# build synthesis models with sensor faults
syssenf = fdimodset(syssen; c = Vector(1:mu), 
                    d = mu .+ Vector(1:md), f = (mu+md) .+ (1:mf));

## Part 2 - Setup of the synthesis specifications

# compute achievable strong specifications for constant faults
FDtol = 1.e-5      # weak fault detection threshold
FDGainTol = 0.0001  # strong fault detection threshold
FDfreq = 0         # frequency for strong detectability checks 
S_strong = fdigenspec(syssenf[nom]; atol, FDtol, FDGainTol, 
                      FDfreq, sdeg = -0.05);

# check number of maximum possible specifications
@test size(S_strong,1) == 2^mf-1

## Part 3 - LPV synthesis using up to 4th order polynomial gain interpolation

# form extended system [ Gu Gd Gf; I 0 0] with outputs extended with u
# and determine constant solutions based on Remark 8.1 of the FDI book
syse = similar(Vector{DescriptorStateSpace},N)
D = zeros(T,mf,p+mu,N)
Dref = [zeros(T,mf,n+mu+md) eye(T,mf)]
for i = 1:N
    syse[i] = [syssenf[i].sys;eye(mu,mu+md+mf)]
    cde = [syse[i].C syse[i].D] 
    D[:,:,i] = Dref/cde
end     

Nintp = similar(Matrix{T},N,maxdeg+1)
Qpol = Polynomial.(zeros(T,mf,p+mu))
y = []
for deg = 0:maxdeg
    for i = 1:mf
        for j = 1:p+mu
            Qpol[i,j] = fit(view(massi,:), view(D,i,j,:), deg)
        end
    end
    for i = 1:N
        Ri = (massi[i] .|> Qpol)*syse[i]
        rezi = Ri - [zeros(mf,mu+md) eye(2)]
        Nintp[i,deg+1]  = glinfnorm(rezi)[1]
        deg == maxdeg && push!(y,stepresp(Ri,10)[1])
    end     
end
@test all(maximum(Nintp,dims = 1) .> [545.  138.  16.  0.88  0.12 0.022 0.002][:,1:maxdeg+1])

println(" Table 8.5  Robustness analysis results of interpolation based approximations")
println(prod([" N$k        " for k = 0:maxdeg]))
display(Nintp)

tout = Vector(0:0.1:10.)
include("Fig8_8.jl")
Fig8_8 = fig
#display(CS2_2.Fig8_8) 

end  # module
using Main.CS2_2

 
