module CS1_2
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Test
using JLD2

## CS1_2  - Case-study example: Monitoring flight actuator faults
#           Local measurements of control surface angles are used
println("Case study CS1_2")

## Part 1 - Model setup
# load matrices of the aircraft multiple-model SYSACM, actuator model ACT,
# output-feedback gain K, percentages of mass variations massi

cd(joinpath(pkgdir(FaultDetectionTools), "Examples"))
@load "cs1data.jld2"

N = size(SYSACM,1)
p, m = size(SYSACM[1].D)


# build minimal realizations of AC-models with actuator faults
# set dimensions
nom = 7;             # nominal system 
# set primary actuator indices
#           [ ailerons,  elevators, stabilizer, ruder ]
act_prim = [ [1,2,15,16]; [17,19];     18;       20   ]; 
mu = size(ACT,1); md = m-mu; mf = length(act_prim);

# form augmented aircraft model with extended measurement set
ee = eye(m); p1 = p+mf;


sysact = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    sysact[i] = gir([ee[act_prim,:]; SYSACM[i]] * append(ACT,dss(eye(md))); atol = 1.e-6)
end

sysactf = fdimodset(sysact; c = Vector(1:mu), d = mu .+ Vector(1:md), f = act_prim)

# determine closed-loop stability margin
sdegcl = -5
for i = 1:N
    global sdegcl = max(sdegcl, maximum(real(eigvals(feedback(sysact[i],K,1:mu, mf.+(1:p); negative = false).A))))
end
@test sdegcl < 0
println("sdegcl = $sdegcl")

## Part 2 - Setup of the synthesis specifications

# compute the achievable strong specifications for constant faults
S_strong = fdigenspec(sysactf[nom], atol = 1.e-6, atol3 = 1.e-6, FDtol = 1.e-5, FDGainTol = 0.001, FDfreq = 0, sdeg = -0.05);

@test size(S_strong,1) == 2^mf-1

# define SFDI, the signatures for isolation of single faults
SFDI = eye(mf) .> 0; 

## Part 3 - Synthesis using Procedure EFDI 
# set options for least order synthesis with EFDISYN
Q, Rf, info = efdisyn( sysactf[nom], SFDI;  atol = 1.e-6, sdeg = -5, smarg = -1, simple = false, FDfreq = 0, FDGainTol = 0.001, rdim = 1);  
 
##  Part 4 - Assesment of synthesis results
# form extended system [ Gu Gd Gf; I 0 0] with outputs extended with u
syse = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    syse[i] = [sysactf[i].sys;eye(mu,mu+md+mf)]
end  

# build overall detector and fault detection system for open-loop system
Qtot = vcat(Q.sys...)
Rftilde = vcat(Rf.sys...)
# check synthesis conditions
# open-loop checks: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rftilde
# check of achieved structure matrix 
@test iszero(Qtot*syse[nom][:,1:mu+md]; atol=1.e-5) && 
      iszero(Qtot*syse[nom][:,mu+md+1:end]-Rftilde; atol=1.e-5) &&
      isequal(fditspec_(Rftilde; FDfreq = 0, atol = 1.e-5)[:,:,1],SFDI)


# form extended closed-loop system with output feedback u = K*y+v
# and evaluate [Ru Rd Rf] for the closed-loop setup

sysefb = similar(Vector{DescriptorStateSpace},N)
Rtot = similar(Vector{DescriptorStateSpace},N)

for i = 1:N
    sysefb[i] = feedback(syse[i], K, 1:mu, mf.+(1:p); negative = false)
    Rtot[i] = Qtot*sysefb[i]
end     

# check robustness by computing the H-inf norms 
NormRu = similar(Vector{Float64},N)
NormRd = similar(Vector{Float64},N)
NormRfmRfnom = similar(Vector{Float64},N)
for i = 1:N
    NormRu[i] = glinfnorm(Rtot[i][:,sysactf[i].controls])[1]
    NormRd[i] = glinfnorm(Rtot[i][:,sysactf[i].disturbances])[1]
    NormRfmRfnom[i] = glinfnorm(Rtot[i][:,sysactf[i].faults]-Rftilde)[1]
end     

norms = [NormRu NormRd NormRfmRfnom]
println(" Table 8.3  Robustness analysis results for the nominal synthesis with position measurements")
println(" NormRu      NormRd       NormRfmRfnom")
display(norms)
@test norm(NormRu,Inf) < 1.e-5 &&  norm(NormRd,Inf) < 1.e-5 &&   norm(NormRfmRfnom,Inf) < 1.e-5 
println("orders = $(order.(Q.sys))")

end  # module

 
