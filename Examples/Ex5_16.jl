module Ex5_16
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.16 - Solution of an H∞ AMMP 
println("Example 5.16")

# define system with control, noise and actuator fault inputs 
A = [-.8 0 0; 0 -.5 .6; 0 -.6 -.5];
Bu = [1 1;1 0;0 1]; Bw = 0.25*[0 0;0 1;1 0]; Bf = Bu;
C = [0 1 1; 1 1 0]; Du = zeros(2,2);
# define Gu(s), Gw(s), Gf(s) and Mr(s)
Gu = dss(A,Bu,C,0); Gw = dss(A,Bw,C,0); Gf = Gu;
mf = size(Gf,2); Mr = dss(eye(mf)); 
p, mu = size(Gu); mw = size(Gw,2); n = order(Gu); 

# compute left nullspace basis as Q1(s) = [ I -Gu(s) ]
Q1 = dss(A,[zeros(n,p) -Bu],C,[eye(p) -Du]); Rf = Gf; Rw = Gw;

# check solvability condition
rank(evalfr(Rf,rand())) == mf || error("No solution exist")

# check for unstable or infinite zeros
Rf_Rw = dss(A,[Bf Bw],C,0);
gzero(Rf_Rw)   # two infinite zeros 

atol = 1.e-7                # tolerance for rank tests
sdeg = -10                  # set stability degree

# compute the quasi-co-outer-co-inner factorization [Rf Rw] = [Go 0]*Gi
Gi, Go = goifac(Rf_Rw; atol); 

# compute Q = inv(Go)*Q1 using explicit formulas 
Qbar = dss([Go.A Go.B; Go.C Go.D], [eye(n,n+mf); zeros(mf,n+mf)],
           [Q1.B; Q1.D], [ zeros(mf,n) -eye(mf)], zeros(mf,p+mu));
    

# compute [F1 F2 ] = [Mr 0]*Gi'
F1_F2 = [Mr zeros(mf,mw)]*Gi'; 

# solve the L∞-LDP
Q4, mininf = glinfldp(F1_F2, mw; atol, reltol = 5.e-4);

# update Q to make it stable and proper
Qtilde = gminreal(Q4*Qbar; atol); 

# compute stable and proper Q = Q4*Qtilde with suitable diagonal Q4 = M
Q = dss(zeros(0,p+mu)); M = dss(zeros(0,0));
for i = 1:mf
    Qi, Mi = glcf(Qtilde[i,:]; atol, sdeg);
    # normalize Mi to unit H∞ norm to match example
    sc = ghinfnorm(Mi)[1]*sign(dcgain(Mi)[1,1]);  
    global Q = [Q; Qi/sc]; 
    global M = append(M, Mi/sc);
end

# convert to standard state space representation
Q = gss2ss(Q)[1]; M = gss2ss(M)[1] 

# display results
println("Q = $(dss2rm(Q, atol = 1.e-7))")
println("M = $(dss2rm(M, atol = 1.e-7))")

# check solution
G = [Gf Gw Gu; zeros(mu,mf+mw) eye(mu)]; F = M*Mr*eye(mf,mu+mw+mf);
err_sub = glinfnorm(Q*G-F)[1]

# compare with the (improper) optimal solution
Yopt = glinfldp(M*F1_F2, mw; atol, reltol = 5.e-4)[1];
Qopt = Yopt*Qbar;
err_opt = glinfnorm(Qopt*G-F)[1]
@test err_sub-err_opt < 0.002

end # module
