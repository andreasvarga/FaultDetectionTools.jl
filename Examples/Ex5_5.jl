module Ex5_5
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.5 - Solution of an AFDP using EFDSYN
println("Example 5.5")

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gw = [1/(s+2); 0];                   # enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       # enter dimensions

# setup the synthesis model 
sysf = fdimodset(dss([Gu Gf Gw],minimal = true),c =1:mu, f = mu .+ (1:mf), n = (mu+mf) .+ (1:mw));

# call of EFDSYN with the options for stability degree -3 and the synthesis 
# of a scalar output filter
Q, Rfw, info = efdsyn(sysf, atol = 1.e-7, minimal = false, nullspace = false, smarg = -2, sdeg = -3); info

# normalize Q and Rf to match example
println("Q = $(dss2rm(Q.sys, atol = 1.e-7))")
println("Rf = $(dss2rm(Rfw.sys[:,Rfw.faults], atol = 1.e-7))")
println("Rw = $(dss2rm(Rfw.sys[:,Rfw.noise], atol = 1.e-7))")
gap = fdif2ngap(Rfw)[1]
println("gap = $gap")

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7)
@test iszero(R.sys[:,[R.controls; R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,[R.faults;R.noise]]-Rfw.sys[:,[Rfw.faults;Rfw.noise]]) && 
      fdisspec(Rfw, block = true) == Bool[1 1] && 
      maximum(real(gpole(Q))) >= -2 && minimum(real(gpole(Q))) <= -3


end # module
using Main.Ex5_5
