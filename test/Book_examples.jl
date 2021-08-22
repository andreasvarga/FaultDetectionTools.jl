module Book_examples

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

## Synthesis examples

# Example 5.3 - Solution of an EFDP using EFDSYN
# Uses DescriptorSystems.jl v1.0 (or later) FaultDetectionTools.jl V0.1 (or later)

# define s as an improper transfer function
s = rtf('s');
# define Gu(s) and Gd(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gd = [(s-1)/(s+2); 0];               # enter Gd(s)
p = 2; mu = 1; md = 1; mf = 2;       # set dimensions

# setup the synthesis model with additive faults with Gf(s) = [ Gu(s) [0;1]]
sysf = fdimodset(dss([Gu Gd]),c =1,d = 2,f = 1,fs = 2);

# call of EFDSYN with the options for stability degree -3 and the synthesis 
# of a scalar output filter
Q, Rf, info = efdsyn(sysf, smarg =-3, sdeg = -3, rdim = 1); 

# normalize Q and Rf to match example
scale = evalfr(Rf.sys[1,1],Inf)[1,1];
dss2rm(Q.sys/scale, atol = 1.e-7)
dss2rm(Rf.sys/scale, atol = 1.e-7) 

# check synthesis results
R = fdRcomp(Q, sysf; minimal = true, atol = 1.e-7)
@test ghinfnorm(R.sys[:,R.controls])[1] < 1.e-7
@test ghinfnorm(R.sys[:,R.disturbances])[1] < 1.e-7
@test iszero(R.sys[:,R.faults]-Rf.sys) 
t, gains = fdisspec(Rf);
@test t[:,:,1] == Bool[1 1] && gains ≈ abs.(dcgain(Rf.sys))

# Example 5.4m - Solution of an EFDP using EFDSYN
# Uses DescriptorSystems.jl v1.0 (or later) FaultDetectionTools.jl V0.1 (or later)

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
Q, Rf, info = efdsyn(sysf, sdeg = -3, rdim = 1); 

# normalize Q and Rf to match example
scale = evalfr(Rf.sys[1,1],Inf)[1,1];
dss2rm(Q.sys/scale, atol = 1.e-7)
dss2rm(Rf.sys/scale, atol = 1.e-7) 

# check synthesis results
R = fdRcomp(Q, sysf; minimal = true, atol = 1.e-7)
@test ghinfnorm(R.sys[:,R.controls])[1] < 1.e-7
@test ghinfnorm(R.sys[:,R.disturbances])[1] < 1.e-7
@test iszero(R.sys[:,R.faults]-Rf.sys) 
t, gains = fdisspec(Rf);
@test t[:,:,1] == Bool[1 1] && gains ≈ abs.(dcgain(Rf.sys))



# Example 5.4 - Solution of an EFDP
# Uses only DescriptorSystems.jl v1.0 (or later) 

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gd(s), Gf(s)
Gu = [(s+1)/(s-2); (s+2)/(s-3)];     # enter Gu(s)
Gd = [(s-1)/(s+2); 0];               # enter Gd(s)
Gf = [(s+1)/(s-2) 0; (s+2)/(s-3) 1]; # enter Gf(s)
p = 2; mu = 1; md = 1; mf = 2;       # set dimensions

# compute a left nullspace basis Q of [Gu Gd; I 0]
Q1 = glnull(dss([Gu Gd;eye(mu,mu+md)]))[1];

# compute Rf1 = Q1[Gf;0]
Rf1 = gir(Q1*dss([Gf;zeros(mu,mf)]));

# check solvability using a random frequency
if minimum(abs.(evalfr(Rf1,rand()))) > 0.01
   # compute a stable left coprime factorization [Q1 Rf1]=inv(Q3)*[Q,Rf]
   # enforce stability degree -3
   Q_Rf, Q3 = glcf([Q1 Rf1];sdeg = -3);
   # extract Q and Rf
   Q = Q_Rf[:,1:p+mu]; Rf = Q_Rf[:,p+mu+1:end]; 
   # normalize Q and Rf to match example
   scale = evalfr(Rf[1,1],Inf)[1,1]
   Q = dss2rm(Q/scale,atol=1.e-7)
   Rf = dss2rm(Rf/scale,atol=1.e-7)
else
   @info "No solution exists"
end

end
