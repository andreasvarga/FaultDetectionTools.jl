module Ex5_12c
using FaultDetectionTools, DescriptorSystems, Test

# Example 5.12c - Solution of an EMMP 
println("Example 5.12c")

# enter output and fault vector dimensions
p = 3; mf = 3;
# generate random dimensions for system order and input vectors
nu = Int(floor(1+4*rand())); mu = Int(floor(1+4*rand()));
nd = Int(floor(1+4*rand())); md = Int(floor(1+4*rand()));

# define random Gu(s) and Gd(s) with triplex sensor redundancy
# and Gf(s) for triplex sensor faults
Gu = ones(3,1)*rss(nu,1,mu); # enter Gu(s) in state-space form
Gd = ones(3,1)*rss(nd,1,md); # enter Gd(s) in state-space form
Gf = eye(3);                 # enter Gf(s) for sensor faults
atol = 1.e-7;                # tolerance for rank tests

# build model with faults sysf = [Gu Gd Gf]
sysf = fdimodset([Gu Gd Gf], c = 1:mu, d = mu.+(1:md), 
                 f = (mu+md).+(1:mf))   

# enter reference model
Mr = FDFilterIF(dss([ 0 1 -1; -1 0 1; 1 -1 0]); mf)

# solve an exact model-matching problem using EMMSYN
Q, R, info = emmsyn(sysf,Mr; atol); 

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Mr
Rt = fdIFeval(Q,sysf; atol) # form Q*[Gu Gd Gf;I 0 0];
@test iszero(Rt.sys[:,[Rt.controls;Rt.disturbances]]; atol) &&
      iszero(Rt.sys[:,Rt.faults]-Mr.sys; atol)

# one step solution
# solve Qbar*Ge = Me, where Ge = [Gu Gd Gf; I 0 0] and Me = [0 0 Mr ].
Ge = [sysf.sys; eye(mu,mu+md+mf)]; Me = [zeros(p,mu+md) Mr.sys];
Qbar = glsol(Ge, Me; atol)[1]

# compare solutions 
@test iszero(Q.sys-Qbar; atol)

   
end # module
using Main.Ex5_12c
