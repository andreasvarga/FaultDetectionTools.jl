module Ex5_16c
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.16c - Solution of an H∞ AMMP 
println("Example 5.16c")

# define system with control, noise and actuator fault inputs 
A = [-.8 0 0; 0 -.5 .6; 0 -.6 -.5];
Bu = [1 1;1 0;0 1]; Bw = 0.25*[0 0;0 1;1 0]; Bf = Bu;
C = [0 1 1; 1 1 0]; Du = zeros(2,2); Dw = zeros(2,2);  
p, mu = size(Du); mf = mu; mw = size(Dw,2); 

# setup synthesis model with additive actuator faults with Bf = Bu, Df = Du
sysf = fdimodset(dss(A,[Bu Bw],C,[Du Dw]), c = 1:mu, n = mu .+ (1:mw),f = 1:mu); 

# define Mr(s) = I
Mr = FDFilterIF(dss(eye(mf)),0,0,mf)

Q, R, info = ammsyn(sysf,Mr; atol = 1.e-7, nullspace = false, reltol = 5.e-4, sdeg = -10, normalize = "dcgain");

gamma_opt0 = info.gammaopt0  # optimal performance for the initial problem 
gamma_opt  = info.gammaopt   # optimal performance 
gamma_sub  = info.gammasub   # suboptimal performance 

# check suboptimal solution
Ge = [sysf.sys[:,[sysf.faults;sysf.noise]]; zeros(mu,mf+mw)]; Me = [info.M*Mr.sys zeros(mf,mw)];
@test iszero(R.sys-Q.sys*Ge, atol = 1.e-7) && 
      fdimmperf(R,info.M*Mr) ≈ gamma_sub && glinfnorm(Me-Q.sys*Ge)[1] ≈ gamma_sub


end # module
