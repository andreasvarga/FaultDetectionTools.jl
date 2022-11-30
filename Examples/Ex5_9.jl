module Ex5_9

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.9 - Solution of an AFDP using AFDSYN with rdim = 2
println("Example 5.9")

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gw = [1/(s+2); 0];                   # enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       # enter dimensions

# setup the synthesis model 
sysf = fdimodset(dss([Gu Gf Gw],minimal = true),c =1:mu, f = mu .+ (1:mf), n = (mu+mf) .+ (1:mw));

# call of AFDSYN with the options for stability degree -3 and the synthesis 
# of a filter with 2 outputs
Q, Rfw, info = afdsyn(sysf, atol = 1.e-7, minimal = true, nullspace = false, smarg = -2, sdeg = -2, poles = [-2,-3],
                      rdim = 2, HDesign = [1 0]); info

# permute rows of Q and R to match example
ip = [2,1]
println("Q = $(dss2rm(Q.sys[ip,:], atol = 1.e-7))")
println("Rf = $(dss2rm(Rfw.sys[ip,Rfw.faults], atol = 1.e-7))")
println("Rw = $(dss2rm(Rfw.sys[ip,Rfw.noise], atol = 1.e-7))")
gap = fdif2ngap(Rfw)[1]    
println("gap = $gap")
β = fdiscond(Rfw)[1]
println("scond = $β")

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7)
@test iszero(R.sys[:,[R.controls; R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,[R.faults;R.noise]]-Rfw.sys[:,[Rfw.faults;Rfw.noise]]) && 
      fdisspec(Rfw, block = false) == [info.S; info.S2] && 
      gap ≈ info.gap


end # module
using Main.Ex5_9
