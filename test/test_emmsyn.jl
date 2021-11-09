module Test_emmsyn

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Random
using Test

# Testing examples for AFDSYN 
println("Test_emmsyn")
@testset "emmsyn" begin
#rand()

##
p = 1; mf = 0; mw = 0;
sysf = fdimodset(rss(1,p,mf+mw),faults = 1:mf,noise = mf.+Vector(1:mw))
sysr = FDFilterIF(sysf.sys,0,0,mf,mw)

# solve an EMMP for a single reference model
Q, R, info = emmsyn(sysf,sysr; minimal = false); info
@test iszero(R.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,sysr) == 0

# solve an EMMP for a bank of two reference models
sysr = FDIFilterIF([sysf.sys,sysf.sys],0,0,mf,mw)
Q, R, info = emmsyn(sysf,sysr; minimal = false); info
@test iszero(vcat(R.sys...)-vcat(Q.sys...)*sysf.sys, atol = 1.e-7) && fdimmperf(R,sysr) == [0,0]

##
p = 2; mf = 2;
sysf = fdimodset(rss(1,p,mf,stable = true),faults = 1:mf)
sysr = FDFilterIF(sysf.sys,0,0,mf)

# solve an EMMP with sysr = sysf: solution Q = I
Q, Rf, info = emmsyn(sysf,sysr,atol = 1.e-7); info
@test iszero(Rf.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(Rf,info.M*sysr) == 0

# solve an EMMP for a bank of two reference models
sysr = FDIFilterIF([sysf.sys,sysf.sys],0,0,mf)
Q, Rf, info = emmsyn(sysf,sysr,atol = 1.e-7); info
@test iszero(vcat(Rf.sys...)-vcat(Q.sys...)*sysf.sys, atol = 1.e-7) && 
      isapprox(fdimmperf(Rf,info.M*sysr), [0,0], atol = 1.e-7)

##
p = 3; mf = 2;
sysf = fdimodset(rss(5,p,mf,stable=true),faults = 1:mf)

sysr = FDFilterIF(sysf.sys,0,0,mf)

# solve an EMMP with sysr = sysf: solution Q = I
Q, Rf, info = emmsyn(sysf, sysr, atol = 1.e-7); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rf.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,info.M*sysr) < 1.e-7 &&
      iszero(R.sys-info.M*sysr.sys, atol = 1.e-7)

# solve an EMMP with sysr = sysf: solution Q = I
Q, Rf, info = emmsyn(sysf, sysr, atol = 1.e-7, minimal = true); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rf.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,info.M*sysr) < 1.e-7

# solve an EMMP for a bank of two reference models
sysr = FDIFilterIF([sysf.sys,sysf.sys],0,0,mf)
Q, Rf, info = emmsyn(sysf,sysr,atol = 1.e-7); info
@test iszero(vcat(Rf.sys...)-vcat(Q.sys...)*sysf.sys, atol = 1.e-7) && 
      isapprox(fdimmperf(Rf,info.M*sysr), [0,0], atol = 1.e-7)

##
p = 3; mf = 2;
sysf = fdimodset(rss(1,p,mf,stable=true),faults = 1:mf)
sysr = FDFilterIF(sysf.sys[1:mf,:],0,0,mf)

# solve an EMMP with sysr = sysf: solution Q = [I 0]
Q, Rf, info = emmsyn(sysf, sysr, atol = 1.e-7, minimal = false); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rf.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,info.M*sysr) < 1.e-7 &&
      iszero(R.sys-info.M*sysr.sys, atol = 1.e-7)

# solve an EMMP using the computed design matrix
Q1, Rf1, info1 = emmsyn(sysf,sysr, atol = 1.e-7, minimal = false, HDesign = info.HDesign); 
@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(Rf.sys-Rf1.sys,atol=1.e-7)


# solve an EMMP with sysr = sysf: solution Q = [I 0]
Q, Rf, info = emmsyn(sysf, sysr, atol = 1.e-7, minimal = true); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rf.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,info.M*sysr) < 1.e-7 &&
      iszero(R.sys-info.M*sysr.sys, atol = 1.e-7) && ismissing(info.HDesign)

Q1, Rf1, info1 = emmsyn(sysf,sysr, atol = 1.e-7, minimal = true, HDesign = eye(2,3)); 
@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(Rf.sys-Rf1.sys,atol=1.e-7)


##
p = 3; mf = 2;
sysf = fdimodset(rss(1,p,mf,stable=true),faults = 1:mf)
sysr = FDFilterIF(sysf.sys[1:mf,:],0,0,mf)

# solve an EMMP with sysr = sysf: solution Q = [I 0]
Q, Rf, info = emmsyn(sysf, sysr, atol = 1.e-7, minimal = false, simple = true); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rf.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,info.M*sysr) < 1.e-5 &&
      iszero(R.sys-info.M*sysr.sys, atol = 1.e-7)


## 
p = 3; mf = 2; mw = 2;
sysf = fdimodset(rss(1,p,mf+mw,stable=true),faults = 1:mf,noise = mf.+Vector(1:mw))
sysr = FDFilterIF(sysf.sys[1:mf,sysf.faults],0,0,mf)

# solve an EMMP with reference model without noise input
Q, Rfw, info = emmsyn(sysf,sysr; atol = 1.e-7, minimal = false); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rfw.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,info.M*sysr) ≈ fdimmperf(R)

sysr = FDFilterIF(sysf.sys[1:mf,:],0,0,mf,mw)

# solve an EMMP with reference model with noise input
Q, Rfw, info = emmsyn(sysf,sysr; atol = 1.e-7, minimal = false); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rfw.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,info.M*sysr) < 1.e-7 

##
p = 3; mu = 2; md = 1; mf = 2; mw = 0;
m = mu+md+mf+mw;
sysf = fdimodset(rss(3,p,m), c = 1:mu,d = mu .+ (1:md),f = (mu+md) .+ (1:mf), n = (mu+md+mf) .+ (1:mw))
sysr = fdimodset(dss(zeros(3,mu+md)),c = 1:mu,d = mu .+ (1:md)) 

Q, Rfw, info = emmsyn(sysf,sysr; atol = 1.e-7); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rfw.sys-Q.sys*[sysf.sys[:,[sysf.controls;sysf.disturbances]]; eye(mu,mu+md)], atol = 1.e-7) && 
       fdimmperf(R,info.M*sysr) < 1.e-7 


# Example 5.12c - Solution of an EMMP 

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

# build model with faults
#sysf = [Gu Gd Gf]
sysf = fdimodset([Gu Gd Gf], c = 1:mu, d = mu.+(1:md), f = (mu+md).+(1:mf))   

# enter reference model
Mr = FDFilterIF(dss([ 0 1 -1; -1 0 1; 1 -1 0]), 0, 0, mf)

# solve an exact model-matching problem using EMMSYN
Q, R, info = emmsyn(sysf,Mr; atol); info

# one step solution
# solve Qbar*Ge = Me, where Ge = [Gu Gd Gf; I 0 0] and Me = [0 0 Mr ].
Ge = [sysf.sys; eye(mu,mu+md+mf)]; Me = [zeros(p,mu+md) Mr.sys];
Qbar = glsol(Ge, Me; atol)[1]

# compare solutions 
@test iszero(Q.sys-Qbar; atol)

# solve an EMMP without explicit nullspace computation
Mre = FDFilterIF([zeros(mf,mu) Mr.sys],mu,0,mf); 

Q1, R1, info1 = emmsyn(sysf,Mre; atol); info1
Rt1 = fdIFeval(Q1, sysf; atol = 1.e-7, minimal = true); 

@test iszero(Q.sys-Q1.sys; atol) && iszero(R1.sys[:,R1.faults]-R.sys[:,R.faults]; atol) 


# Example 5.13c - Solution of an EMMP using EMMSYN

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

# check solution
G = [sysf.sys; eye(mu,mu+mf)]; F = [zeros(mf,mu) info.M*Mr.sys];
@test iszero(Q.sys*G-F; atol)

# Example 5.4c - Solution of an EFDP using EFDSYN
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

#  solve an EMMP
QM, RM, info = emmsyn(sysf,Rf)
@test iszero(Q.sys-QM.sys,atol=1.e-7) && iszero(Rf.sys-RM.sys,atol=1.e-7) 

## Model Niemann 1998, Optim. Appl. Meth. 
n = 5; mu = 2; md = 1; mf = 2; m = mu+md+mf; p = 5;
a = [
-0.0782 0.2939 0.0220 0.0208 -0.0291
0.0077 -0.0278 0.0015 -0.0015 0.0026
-1.3251 5.32608 -0.5263 0.2214 -0.4777
1.0809 -4.4452 0.3770 -0.4631 0.4032
2.1532 -8.6386 0.7811 -0.5745 0.7816]*1000;
bd = [
    0.0443
   -0.0042
    0.7910
   -0.6598
   -1.2881]*1000;
bu = [
    -0.007197 0.003005
0.003473 0.000284
1.218502 -0.032866
1.322502 0.020147
-0.082265 0.024388 ];
bf = bu;
c = eye(n); dd = zeros(n,md); du = zeros(n,mu); df = zeros(n,mf);

sys = dss(a,[bu bd bf],c,[du dd df]);
sysf = fdimodset(sys, c = 1:mu, d = mu .+ (1:md) , f = (mu+md) .+ (1:mf));

SFDI = fdigenspec(sysf) # determine achievable structure matrix
nb = size(SFDI,1)
@time Q, Rf = efdisyn(sysf, SFDI; atol = 1.e-7, rdim = 1);
@time Q, Rf = efdisyn(sysf, SFDI; atol = 1.e-7, rdim = 1, sdeg = -5, smarg = -0.05,
                                  FDfreq = 0, FDGainTol = 0.0001);

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
R = fdIFeval(Q,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      order.(Q.sys) == ones(Int,nb)

#  solve an EMMP
QM, RM, info = emmsyn(sysf,Rf; atol = 1.e-7, minimal = true);
# @test iszero(vcat(Q.sys...)-vcat(QM.sys...),atol=1.e-7) && iszero(Rf.sys-RM.sys,atol=1.e-7) 

R1 = fdIFeval(QM,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R1.sys...)[:,[R1.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R1.sys...)[:,R1.faults]-vcat(RM.sys...),atol=1.e-7) &&
      order.(QM.sys) == ones(Int,nb) && 
      iszero(vcat((info.M*Rf).sys...)[:,Rf.faults]-vcat(RM.sys...),atol=1.e-7) 


## Example 2, Xiong & Saif 2000, Int. J. Robust Nonlinear Control
n = 5; mu = 3; md = 2; mf = 2; m = mu+md+mf; p = 3;
a = [
0 0 -0.0034 0 0
0 -0.041 0.0013 0 0
0 0 -1.1471 0 0
0 0 -0.0036 0 0
0 0.094 0.0057 0 -0.051];
bu = [
    -1 0 0
     0 0 0
     0 0 0.948
     0.916 -1 0
    -0.598 0 0];
bd = [
     0 1
     0.062 -0.132
     0 -7.189
     0 0
     0 0 ];
bf = bu[:,1:2];
c = [
1 0 0 0 0
0 0 1 0 0
0 0 0 1 0]; 
dd = zeros(p,md); du = zeros(p,mu); df = zeros(p,mu-1);

sys = dss(a,[bu bd bf],c,[du dd df]);

sysf = fdimodset(sys, c = 1:mu, d = mu .+ (1:md) , f = (mu+md) .+ (1:mf));

SFDI = fdigenspec(sysf) # determine achievable structure matrix
nb = size(SFDI,1)
@time Q, Rf = efdisyn(sysf, SFDI; atol = 1.e-7, rdim = 1, sdeg =-2, smarg = -1,
                                  FDfreq = 0, FDGainTol = 0.0001);

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
R = fdIFeval(Q,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      order.(Q.sys) == [1, 1, 1]

# check weak and strong fault detectability
@test fditspec(Rf) == fdisspec(Rf) == SFDI 


#  solve an EMMP
QM, RM, info = emmsyn(sysf,Rf; atol = 1.e-7, minimal = true);
# @test iszero(vcat(Q.sys...)-vcat(QM.sys...),atol=1.e-7) && iszero(Rf.sys-RM.sys,atol=1.e-7) 

R1 = fdIFeval(QM,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R1.sys...)[:,[R1.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R1.sys...)[:,R1.faults]-vcat(RM.sys...),atol=1.e-7) &&
      order.(QM.sys) == ones(Int,nb) && 
      iszero(vcat((info.M*Rf).sys...)[:,Rf.faults]-vcat(RM.sys...),atol=1.e-7) 


end # test


end   # module
