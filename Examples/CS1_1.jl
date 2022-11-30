module CS1_1
using FaultDetectionTools, DescriptorSystems, LinearAlgebra, Test
using JLD2

## CS1_1  - Case-study example: Monitoring flight actuator faults
#           No measurements of control surface angles are used
println("Case study CS1_1 with Fig8.2 and Fig8.3")

## Part 1 - Model setup
# load matrices of the aircraft multiple-model SYSACM, 
# actuator model ACT and output-feedback gain K
cd(joinpath(pkgdir(FaultDetectionTools), "Examples"))
SYSACM, K, ACT = load("cs1data.jld2", "SYSACM", "K", "ACT")

# set dimensions
N = size(SYSACM,1)        # number of models
p, m = size(SYSACM[1].D)  # output and input dimensions
atol = 1.e-6              # absolute tolerance
nom = 7;                  # index of nominal system 

# set primary actuator indices
#           [ ailerons,  elevators, stabilizer, ruder ]
act_prim = [ [1,2,15,16]; [17,19];     18;       20   ]; 
mu = size(ACT,1); md = m-mu; mf = length(act_prim);
indu = 1:mu      # system control input indicess
indy = 1:p       # system measured output indices

# build minimal realizations of AC-models with actuator faults
sysact = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    sysact[i] = gir(SYSACM[i] * append(ACT,dss(eye(md))); atol)
end

# build synthesis models with actuator faults 
sysactf = fdimodset(sysact; c = Vector(1:mu), 
          d = mu .+ Vector(1:md), f = act_prim)

# determine closed-loop stability margin
sdeg_cloop = -5
for i = 1:N
    syst = feedback(sysact[i],K,indu,indy; negative = false)
    global sdeg_cloop = max(sdeg_cloop, 
                            maximum(real(eigvals(syst.A))))
end
@test sdeg_cloop < 0
println("sdeg_cloop = $sdeg_cloop")

## Part 2 - Setup of the synthesis specifications

# compute the achievable weak specifications
FDtol = 1.e-5      # weak fault detection threshold
S = fdigenspec(sysactf[nom]; atol, FDtol);

# compute achievable strong specifications for constant faults
FDGainTol = 0.001  # strong fault detection threshold
FDfreq = 0         # frequency strong detectability checks 
S_strong = fdigenspec(sysactf[nom]; atol, FDtol, FDGainTol, 
                      FDfreq, sdeg = -0.05)

# check inclusion of strong specification
@test issubset(binmat2dec(S_strong), binmat2dec(S)) 

# define SFDI, the signatures for isolation of single faults
SFDI = [0 1 1 1 1 0 1 0
        1 0 1 1 0 1 1 0
        1 1 0 1 1 0 1 0
        1 1 1 0 1 1 0 0
        1 1 1 1 0 0 0 0
        0 0 0 0 0 0 0 1] .> 0;

## Part 3 - Synthesis using Procedure EFDI 
# set options for least order synthesis with EFDISYN
Q, Rf, info = efdisyn(sysactf[nom], SFDI; atol, FDGainTol, FDfreq,  
                      rdim = 1, sdeg = -5, smarg = -0.2)  
println("orders = $(order.(Q.sys))")
 
##  Part 4 - Assessment of synthesis results
# form extended open-loop systems [ Gu Gd Gf; I 0 0] 
syse = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    syse[i] = [sysactf[i].sys;eye(mu,mu+md+mf)]
end     

# build overall detector filter and its internal form
Qtot = vcat(Q.sys...)
Rftilde = vcat(Rf.sys...)

# check open-loop synthesis conditions:
# (1) Q*[Gu Gd;I 0] = 0
# (2) Q*[Gf; 0] = Rftilde
# (3) achieved structure matrix = SFDI
@test iszero(Qtot*syse[nom][:,1:mu+md]; atol) && 
      iszero(Qtot*syse[nom][:,mu+md+1:end]-Rftilde; atol) &&
      isequal(fditspec_(Rftilde; FDfreq, atol)[:,:,1],SFDI)

# form extended closed-loop system with output feedback u = K*y+v
# and evaluate [Ru Rd Rf] for the closed-loop setup

sysefb = similar(Vector{DescriptorStateSpace},N)
Rtot = similar(Vector{DescriptorStateSpace},N)

for i = 1:N
    sysefb[i] = feedback(syse[i], K, indu, indy; negative = false)
    Rtot[i] = Qtot*sysefb[i]
end     

# check robustness by computing the H-inf norms 
NormRu = similar(Vector{Float64},N)
NormRd = similar(Vector{Float64},N)
NormRfmRfnom = similar(Vector{Float64},N)
indc = sysactf[1].controls       # indices of control inputs
indd = sysactf[1].disturbances   # indices of disturbances inputs
indf = sysactf[1].faults         # indices of fault inputs
for i = 1:N
    NormRu[i] = glinfnorm(Rtot[i][:,indc])[1]
    NormRd[i] = glinfnorm(Rtot[i][:,indd])[1]
    NormRfmRfnom[i] = glinfnorm(Rtot[i][:,indf]-Rftilde)[1]
end  
# comment out for a more accurate evaluation of norm  
# rtolinf = 1.e-6
# for i = 1:N
#     NormRu[i] = glinfnorm(Rtot[i][:,indc];rtolinf)[1]
#     NormRd[i] = glinfnorm(Rtot[i][:,indd];rtolinf)[1]
#     NormRfmRfnom[i] = glinfnorm(Rtot[i][:,indf]-Rftilde;rtolinf)[1]
# end

println(" Table 8.2  Robustness analysis results for the nominal synthesis")
println(" NormRu        NormRd        NormRfmRfnom")
display([NormRu NormRd NormRfmRfnom])
@test NormRu[nom] < 1.e-6 &&  NormRd[nom] < 1.e-6 &&   NormRfmRfnom[nom] < 1.e-6 

# generate figures
y, tout, x = stepresp(Rftilde,10);
include("Fig8_2.jl")
Fig8_2 = f

y = [stepresp(Rtot[i][:,indf],10)[1] for i in 1:N];
include("Fig8_3.jl")
Fig8_3 = f

end  # module

using Main.CS1_1
CS1_1.Fig8_2
CS1_1.Fig8_3
 
