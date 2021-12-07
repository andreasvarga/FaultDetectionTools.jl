module Ex5_11c
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.11c - Solution of an AFDIP 
println("Example 5.11c")

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gw = [1/(s+2); 0];                   # enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       # enter dimensions
SFDI = eye(mf) .> 0

# build the synthesis model with additive faults 
sysf = fdimodset(dss([Gu Gf Gw]), c = 1:mu,f = mu .+ (1:mf), n = (mu+mf) .+ (1:mw));


Q, R, info = afdisyn(sysf, SFDI; smarg =-3, sdeg = -3)

# scale with -1 to match example
Q.sys[1] = -Q.sys[1]; R.sys[1] = -R.sys[1];

for i = 1:size(SFDI,1)
    println("Q[$i] = $(dss2rm(Q.sys[i], atol = 1.e-7))")
    println("Rf[$i] = $(dss2rm(R.sys[i][:,R.faults], atol = 1.e-7))")
    println("Rw[$i] = $(dss2rm(R.sys[i][:,R.noise], atol = 1.e-7))")
end

# check synthesis conditions: Qt*[Gu Gd;I 0] = 0 and Qt*[Gf; 0] = Rf
Rt = fdIFeval(Q,sysf) # form Qt*[Gu Gd Gf;I 0 0];
@test iszero(vcat(Rt.sys...)[:,[Rt.controls;Rt.disturbances]],atol=1.e-7) &&
      iszero(vcat(Rt.sys...)[:,[Rt.faults;Rt.noise]]-vcat(R.sys...),atol=1.e-7) &&
      fdif2ngap(R,SFDI)[1] ≈ info.gap


println("Q = ")
display(Q)
println("R = ")
display(R)
println("gap = $(fdif2ngap(R,SFDI)[1])")
   
end # module
