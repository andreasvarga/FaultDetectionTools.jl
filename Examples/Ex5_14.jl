module Ex5_14
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.14 - Solution of an H∞ AMMP 
println("Example 5.14")

s = rtf('s'); # define the Laplace variable s
Gf = 1/(s+1);        # enter $G_f(s)$
Gw = 1/(s+2);        # enter $G_w(s)$
mu = 0; mf = 1; mw = 1; p = 1; # set dimensions

# build the synthesis model with additive faults 
sysf = fdimodset(dss([Gf Gw]), f = 1:mf, n = mf .+ (1:mw));

# define Mr(s) = 1/(s+3)
Mr = FDFilterIF(dss(1/(s+3)); mf)

Q, R, info = ammsyn(sysf, Mr; atol = 1.e-7, reltol = 1.e-5, sdeg = -1, normalize = "gain");

println("Q = $(zpk(dss2rm(Q.sys)[1,1]))")

# check suboptimal solution
Ge = [sysf.sys[:,[sysf.faults;sysf.noise]]; zeros(mu,mf+mw)]; Me = [info.M*Mr.sys zeros(mf,mw)];
@test iszero(R.sys-Q.sys*Ge, atol = 1.e-7) && 
      fdimmperf(R,info.M*Mr) ≈ info.gammasub && glinfnorm(Me-Q.sys*Ge)[1] ≈ info.gammasub

println("gamma_opt0 = $(info.gammaopt0)")  # optimal performance for the initial problem 
println("gamma_opt  = $(info.gammaopt)")   # optimal performance 
println("gamma_sub  = $(info.gammasub)")   # suboptimal performance 


end # module
