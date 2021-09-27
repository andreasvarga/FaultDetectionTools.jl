module Ex5_4c
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.4c - Solution of an EFDP using EFDSYN
println("Example 5.4c")

# define s as an improper transfer function
s = rtf('s');
# define Gu(s) and Gd(s)
Gu = [(s+1)/(s-2); (s+2)/(s-3)];     # enter Gu(s)
Gd = [(s-1)/(s+2); 0];               # enter Gd(s)
Gf = [(s+1)/(s-2) 0; (s+2)/(s-3) 1]; # enter Gf(s)
p = 2; mu = 1; md = 1; mf = 2;       # set dimensions

# setup the synthesis model with faults with Gf(s) = [ Gu(s) [0;1]]
sysf = fdimodset(dss([Gu Gd]),c =1,d = 2,f = 1,fs = 2);

# call of EFDSYN with the options for stability degree -3 and the synthesis 
# of a scalar output filter
Qt, Rft, info = efdsyn(sysf, sdeg = -3, rdim = 1); 

# normalize Q and Rf to match example
scale = evalfr(Rft.sys[1,1],Inf)[1,1];
Q = FDFilter(Qt.sys/scale,Qt.outputs,Qt.controls);
Rf = FDFilterIF(Rft.sys/scale,faults = Rft.faults);
println("Q = $(dss2rm(Q.sys, atol = 1.e-7))")
println("Rf = $(dss2rm(Rf.sys, atol = 1.e-7))")


# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7)
@test iszero(R.sys[:,[R.controls; R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys) && 
      fdisspec(Rf) == Bool[1 1] && 
      gpole(Q.sys) ≈ [-3] && gpole(Rf.sys) ≈ [-3] 


end # module
