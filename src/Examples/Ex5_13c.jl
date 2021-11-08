module Ex5_13c
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.13c - Solution of an EMMP using EMMSYN
println("Example 5.13c")

# define s as an improper transfer function
s = rtf('s');
# enter Gu(s), Gf(s) and Mr(s)
Gu = [s/(s^2+3*s+2) 1/(s+2);
     s/(s+1) 0;
      0 1/(s+2)];
Mr = dss(eye(2));                  # enter Mr(s)
p, mu = size(Gu); mf = mu

# build the synthesis model with additive faults 
sysf = fdimodset(dss(Gu), c = 1:mu, f = 1:mu);

# enter reference model
Mr = fdimodset(dss(eye(mf)), f = 1:mf)

atol = 1.e-7                # tolerance for rank tests
sdeg = -1                   # set stability degree

# solve an exact model-matching problem using EMMSYN
Q, R, info = emmsyn(sysf, Mr; atol, sdeg, minimal = false); info

# display results
println("Q = $(dss2rm(Q.sys, atol = 1.e-7))")
println("M = $(dss2rm(info.M, atol = 1.e-7))")

# check solution
G = [sysf.sys; eye(mu,mu+mf)]; F = [zeros(mf,mu) info.M*Mr.sys];
@test iszero(Q.sys*G-F; atol)

end # module
