module Ex7_3c
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 7.3c - Nullspace-based synthesis 

println("Example 7.3c")

A = diagm([ -2, -3, -1, 2, 3]);
Bu = [1 1 0 0 0]'; Bd = [0 0 2 0 0]';
Bf = [0 0 0 2 2;0 0 0 0 0]';
C = [-1 0 -1 1.5 0; 0 -1 0 0 2.5];
Du = [ 1 1]'; Dd = [1 0]'; Df = [1 0; 1 1];
p = 2; mu = 1; md = 1; mf = 2;       # enter dimensions

# setup the synthesis model 
sysf = fdimodset(dss(A,[Bu Bd Bf],C,[Du Dd Df]),c =1:mu, d = mu .+ (1:md), f = (mu+md) .+ (1:mf));


# call of EFDSYN with the options for stability degree -3 and the synthesis 
# of a scalar output filter
Qt, Rft, info = efdsyn(sysf, sdeg = -3, rdim = 1); 

# compute Q and Rf; scale to meet example
scale = -1
Q = FDFilter(Qt.sys*scale,Qt.outputs,Qt.controls);
Rf = FDFilterIF(Rft.sys*scale,faults = Rft.faults);

# display results
println("Q = $(dss2rm(Q.sys))")
println("Rf = $(dss2rm(Rf.sys))")

end # module
