module Test_ammsyn

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Random
using Test

# Testing examples for AMMSYN 
println("Test_ammsyn")
@testset "ammsyn" begin

p = 1; mf = 0; mw = 0;
sysf = fdimodset(rss(1,p,mf+mw),faults = 1:mf,noise = mf.+Vector(1:mw))
sysr = FDFilterIF(sysf.sys,0,0,mf,mw)

# solve an EMMP for a single reference model
@time Q, R, info = ammsyn(sysf, sysr); info
@test iszero(R.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,sysr) == 0

# solve an EMMP for a bank of two reference models
sysr = FDIFilterIF([sysf.sys,sysf.sys],0,0,mf,mw)
@time Q, R, info = ammsyn(sysf,sysr); info
@test iszero(vcat(R.sys...)-vcat(Q.sys...)*sysf.sys, atol = 1.e-7) && fdimmperf(R,sysr) == [0,0]

##
p = 2; mf = 2;
sysf = fdimodset(rss(1,p,mf,stable = true),faults = 1:mf)
sysr = FDFilterIF(sysf.sys,0,0,mf)

# solve an EMMP with sysr = sysf: solution Q = I
@time Q, Rf, info = ammsyn(sysf,sysr,atol = 1.e-7); info
@test iszero(Rf.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(Rf,info.M*sysr) < 1.e-7

# solve an EMMP for a bank of two reference models
sysr = FDIFilterIF([sysf.sys,sysf.sys],0,0,mf)
Q, Rf, info = ammsyn(sysf,sysr,atol = 1.e-7); info
@test iszero(vcat(Rf.sys...)-vcat(Q.sys...)*sysf.sys, atol = 1.e-7) && 
      isapprox(fdimmperf(Rf,info.M*sysr), [0,0], atol = 1.e-7)


p = 3; mf = 2;
sysf = fdimodset(rss(5,p,mf,stable=true),faults = 1:mf)

sysr = FDFilterIF(sysf.sys,0,0,mf)

# solve an EMMP with sysr = sysf: solution Q = I
@time Q, Rf, info = ammsyn(sysf, sysr, atol = 1.e-7); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rf.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,info.M*sysr) < 1.e-7 &&
      iszero(R.sys-info.M*sysr.sys, atol = 1.e-7)

# solve an EMMP with sysr = sysf: solution Q = I
Q, Rf, info = ammsyn(sysf, sysr, atol = 1.e-7, mindeg = true); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rf.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,info.M*sysr) < 1.e-7

# solve an EMMP for a bank of two reference models
sysr = FDIFilterIF([sysf.sys,sysf.sys],0,0,mf)
@time Q, Rf, info = ammsyn(sysf,sysr,atol = 1.e-7); info
@test iszero(vcat(Rf.sys...)-vcat(Q.sys...)*sysf.sys, atol = 1.e-7) && 
      isapprox(fdimmperf(Rf,info.M*sysr), [0,0], atol = 1.e-7)

##
p = 3; mf = 2;
sysf = fdimodset(rss(1,p,mf,stable=true),faults = 1:mf)
sysr = FDFilterIF(sysf.sys[1:mf,:],0,0,mf)

# solve an EMMP with sysr = sysf: solution Q = [I 0]
@time Q, Rf, info = ammsyn(sysf, sysr, atol = 1.e-7, mindeg = false); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rf.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,info.M*sysr) < 1.e-7 &&
      iszero(R.sys-info.M*sysr.sys, atol = 1.e-7)

# solve an EMMP using the computed design matrix
Q1, Rf1, info1 = ammsyn(sysf,sysr, atol = 1.e-7, mindeg = false, HDesign = info.HDesign); 
@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(Rf.sys-Rf1.sys,atol=1.e-7)


# solve an EMMP with sysr = sysf: solution Q = [I 0]
@time Q, Rf, info = ammsyn(sysf, sysr, atol = 1.e-7, mindeg = true); info
R = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(Rf.sys-Q.sys*sysf.sys, atol = 1.e-7) && fdimmperf(R,info.M*sysr) < 1.e-7 &&
      iszero(R.sys-info.M*sysr.sys, atol = 1.e-7) 

Q1, Rf1, info1 = ammsyn(sysf,sysr, atol = 1.e-7, mindeg = true, HDesign = info.HDesign); 
@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(Rf.sys-Rf1.sys,atol=1.e-7)

# Example 5.6c - Solution of an AFDP using AMMSYN

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gw = [1/(s+2); 0];                   # enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       # enter dimensions

# setup the synthesis model 
sysf = fdimodset(dss([Gu Gf Gw],minimal = true),c =1:mu, f = mu .+ (1:mf), n = (mu+mf) .+ (1:mw));

Mr = FDFilterIF(dss(ones(1,2)),0,0,2);                 

# call of AMMSYN with the options for stability degree -3 and the synthesis 
# of a scalar output filter
@time Q, R, info = ammsyn(sysf, Mr, atol = 1.e-7, mindeg = true, smarg = -3, sdeg = -3, nullspace = false, HDesign = [1 1]); info

R1 = fdIFeval(Q, sysf; atol = 1.e-7, minimal = true); 
@test iszero(R1.sys[:,[R1.controls;R1.disturbances]],atol = 1.e-7) &&
      iszero(R.sys[:,[R.faults;R.noise]]-R1.sys[:,[R1.faults;R1.noise]], atol = 1.e-7) && 
      fdimmperf(R,info.M*Mr) ≈ info.gammaopt 

# check with known performance measures
@test fdiscond(R)[1] ≈ 0.7817359599705505 && fdimmperf(R,Mr) ≈ 0.7385533646476293 &&
      fdif2ngap(R)[1] ≈ 2.345207879911668

# Example 5.16c - Solution of an H∞ AMMP 

# define system with control, noise and actuator fault inputs 
A = [-.8 0 0; 0 -.5 .6; 0 -.6 -.5];
Bu = [1 1;1 0;0 1]; Bw = 0.25*[0 0;0 1;1 0]; Bf = Bu;
C = [0 1 1; 1 1 0]; Du = zeros(2,2); Dw = zeros(2,2);  
p, mu = size(Du); mf = mu; mw = size(Dw,2); 

# setup synthesis model with additive actuator faults with Bf = Bu, Df = Du
sysf = fdimodset(dss(A,[Bu Bw],C,[Du Dw]), c = 1:mu, n = mu .+ (1:mw),f = 1:mu); 

# define Mr(s) = I
Mr = FDFilterIF(dss(eye(mf)),0,0,mf)

@time Q, R, info = ammsyn(sysf,Mr; atol = 1.e-7, nullspace = false, reltol = 5.e-4, sdeg = -10, normalize = "dcgain");

gamma_opt0 = info.gammaopt0  # optimal performance for the initial problem 
gamma_opt  = info.gammaopt   # optimal performance 
gamma_sub  = info.gammasub   # suboptimal performance 

# check suboptimal solution
Ge = [sysf.sys[:,[sysf.faults;sysf.noise]]; zeros(mu,mf+mw)]; Me = [info.M*Mr.sys zeros(mf,mw)];
@test iszero(R.sys-Q.sys*Ge, atol = 1.e-7) && 
      fdimmperf(R,info.M*Mr) ≈ gamma_sub && glinfnorm(Me-Q.sys*Ge)[1] ≈ gamma_sub

@time Q, R, info = ammsyn(sysf,Mr; H2syn = true, atol = 1.e-7, nullspace = false, reltol = 5.e-4, sdeg = -10, normalize = "dcgain");

# check suboptimal solution
Ge = [sysf.sys[:,[sysf.faults;sysf.noise]]; zeros(mu,mf+mw)]; Me = [info.M*Mr.sys zeros(mf,mw)];
@test iszero(R.sys-Q.sys*Ge, atol = 1.e-7) && 
      fdimmperf(R,info.M*Mr,2; atolinf = 1.e-7) ≈ info.gammasub && gl2norm(Me-Q.sys*Ge) ≈ info.gammasub


# Example 5.11c - Solution of an AFDIP 

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gw = [1/(s+2); 0];                   # enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       # enter dimensions

# build the synthesis model with additive faults 
sysf = fdimodset(dss([Gu Gf Gw]), c = 1:mu,f = mu .+ (1:mf), n = (mu+mf) .+ (1:mw));

# define Mr(s) = I
Mr = FDFilterIF(dss(eye(mf)),0,0,mf)

@time Q, R, info = ammsyn(sysf,Mr; atol = 1.e-7);

# check suboptimal solution
Ge = [sysf.sys[:,[sysf.faults;sysf.noise]]; zeros(mu,mf+mw)]; Me = [info.M*Mr.sys zeros(mf,mw)];
@test iszero(R.sys-Q.sys*Ge, atol = 1.e-7) && 
      fdimmperf(R,info.M*Mr) ≈ info.gammasub && glinfnorm(Me-Q.sys*Ge)[1] ≈ info.gammasub


# Example 5.14 of (V,2017)
s = rtf('s'); # define the Laplace variable s
Gf = 1/(s+1);        # enter $G_f(s)$
Gw = 1/(s+2);        # enter $G_w(s)$
mu = 0; mf = 1; mw = 1; p = 1; # set dimensions

# build the synthesis model with additive faults 
sysf = fdimodset(dss([Gf Gw]), f = 1:mf, n = mf .+ (1:mw));

# define Mr(s) = 1/(s+3)
Mr = FDFilterIF(dss(1/(s+3)),0,0,mf)

@time Q, R, info = ammsyn(sysf, Mr; atol = 1.e-7, reltol = 1.e-5, sdeg = -1, normalize = "gain");
# check suboptimal solution
Ge = [sysf.sys[:,[sysf.faults;sysf.noise]]; zeros(mu,mf+mw)]; Me = [info.M*Mr.sys zeros(mf,mw)];
@test iszero(R.sys-Q.sys*Ge, atol = 1.e-7) && 
      fdimmperf(R,info.M*Mr) ≈ info.gammasub && glinfnorm(Me-Q.sys*Ge)[1] ≈ info.gammasub


Me = FDFilterIF([zeros(mf,mu) Mr.sys zeros(mf,mw)], 0, 0, mf, mw )

@time Q1, R1, info1 = ammsyn(sysf, Me; regmin = false, atol = 1.e-7, reltol = 1.e-5, sdeg = -1, normalize = "gain");
@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(R.sys-R1.sys,atol=1.e-7)


# Example 5.15 - Solution of an H2 AMMP 

s = rtf('s'); # define the Laplace variable s
Gf = 1/(s+1);        # enter $G_f(s)$
Gw = 1/(s+2);        # enter $G_w(s)$
mu = 0; mf = 1; mw = 1; p = 1; # set dimensions

# build the synthesis model with additive faults 
sysf = fdimodset(dss([Gf Gw]), f = 1:mf, n = mf .+ (1:mw));

# define Mr(s) = 1/(s+3)
Mr = FDFilterIF(dss(1/(s+3)),0,0,mf)

Q, R, info = ammsyn(sysf, Mr; H2syn = true, atol = 1.e-7, reltol = 1.e-5, sdeg = -1, normalize = "gain");

# check suboptimal solution
Ge = [sysf.sys[:,[sysf.faults;sysf.noise]]; zeros(mu,mf+mw)]; Me = [info.M*Mr.sys zeros(mf,mw)];
@test iszero(R.sys-Q.sys*Ge, atol = 1.e-7) && 
      fdimmperf(R, info.M*Mr,2; atolinf = 1.e-7) ≈ info.gammasub && gl2norm(Me-Q.sys*Ge)[1] ≈ info.gammasub

# Example 5.17 - Solution of an H2 AMMP 

# define system with control, noise and actuator fault inputs 
A = [-.8 0 0; 0 -.5 .6; 0 -.6 -.5];
Bu = [1 1;1 0;0 1]; Bw = 0.25*[0 0;0 1;1 0]; Bf = Bu;
C = [0 1 1; 1 1 0]; Du = zeros(2,2); Dw = zeros(2,2);  
p, mu = size(Du); mf = mu; mw = size(Dw,2); 

# setup synthesis model with additive actuator faults with Bf = Bu, Df = Du
sysf = fdimodset(dss(A,[Bu Bw],C,[Du Dw]), c = 1:mu, n = mu .+ (1:mw),f = 1:mu); 

# define Mr(s) = 10/(s+10)
Mr = FDFilterIF(dss([10/(s+10) 0; 0 10/(s+10)]),0,0,mf)

@time Q, R, info = ammsyn(sysf, Mr; H2syn = true, atol = 1.e-7, nullspace = false, reltol = 5.e-4, sdeg = -10, normalize = "dcgain");

# check suboptimal solution
Ge = [sysf.sys[:,[sysf.faults;sysf.noise]]; zeros(mu,mf+mw)]; Me = [info.M*Mr.sys zeros(mf,mw)];
@test iszero(R.sys-Q.sys*Ge, atol = 1.e-7) && 
      fdimmperf(R, info.M*Mr,2; atolinf = 1.e-7) ≈ info.gammasub && gl2norm(Me-Q.sys*Ge)[1] ≈ info.gammasub


# Example 5.11c - Solution of an AFDIP 

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gw = [1/(s+2); 0];                   # enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       # enter dimensions
SFDI = eye(mf) .> 0
nb = mf

# build the synthesis model with additive faults 
sysf = fdimodset(dss([Gu Gf Gw]), c = 1:mu,f = mu .+ (1:mf), n = (mu+mf) .+ (1:mw));


Q, R, info = afdisyn(sysf, SFDI; smarg =-3, sdeg = -3)

Mr = FDFilterIF([R.sys[1][:,R.faults]; R.sys[2][:,R.faults]],0,0,mf);

@time QM, RM, info1 = ammsyn(sysf, Mr; atol = 1.e-7, reltol = 1.e-5, sdeg = -3);

R1 = fdIFeval(QM,sysf); # form Q*[Gu Gd Gf;I 0 0];

# check suboptimal solution
Ge = [sysf.sys[:,[sysf.faults;sysf.noise]]; zeros(mu,mf+mw)]; Me = [info1.M*Mr.sys zeros(mf,mw)];
@test iszero(RM.sys-QM.sys*Ge, atol = 1.e-7) && 
      fdimmperf(R1, info1.M*Mr) ≈ info1.gammasub && glinfnorm(Me-QM.sys*Ge)[1] ≈ info1.gammasub &&
      fdimmperf(R1,fditspec(Mr)) <= info1.gammasub
      isapprox(fdif2ngap(R1,fditspec(Mr))[1], info.gap; atol=1.e-7)

Mr = FDIFilterIF([R.sys[1][:,R.faults], R.sys[2][:,R.faults]],0,0,mf);  # choose targeted reference model

@time QM, RM, info1 = ammsyn(sysf, Mr; atol = 1.e-7, reltol = 1.e-4, sdeg = -3);

R1 = fdIFeval(QM,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R1.sys...)[:,[R1.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R1.sys...)[:,[R1.faults;R1.noise]]-vcat(RM.sys...),atol=1.e-7) &&
      fdimmperf(R1,fditspec(Mr))[1] .<= info1.gammasub[1] &&
      isapprox(fdif2ngap(R1,fditspec(Mr))[1], info.gap; atol=1.e-5)


# Example 5.3 - Solution of an AFDP with strong synthesis 

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
mu = 1; mw = 1; p = 2; mf = mu+p; # set dimensions
Gu = [(s+1)/(s+2); (s+2)/(s+3)];  # enter $G_u(s)$
Gf = [Gu eye(p)];                 # enter $G_f(s)$
Gw = [(s-1)/(s+2); 0];            # enter $G_w(s)$

# build the synthesis model with additive faults 
sysf = fdimodset(dss([Gu Gf Gw]; minimal = true), c = 1:mu, f = mu .+ (1:mf), n = (mu+mf) .+ (1:mw));

Q1, R1, info1 = afdsyn(sysf; FDfreq = [0, 0.01], sdeg = -1)

Mr = FDFilterIF(R1.sys[:,R1.faults],0,0,mf);

@time Q, R, info = ammsyn(sysf, Mr; atol = 1.e-7, reltol = 1.e-5, sdeg = -1);

R2 = fdIFeval(Q,sysf); # form Q*[Gu Gd Gf;I 0 0];

# check suboptimal solution
Ge = [sysf.sys; eye(mu,mu+mf+mw)]; Me = [zeros(1,mu) info.M*Mr.sys zeros(1,mw)];
@test iszero(R2.sys-Q.sys*Ge, atol = 1.e-7) && 
      fdimmperf(R, info.M*Mr) ≈ info.gammasub && glinfnorm(Me-Q.sys*Ge)[1] ≈ info.gammasub 

Mre = FDFilterIF(Me,mu,0,mf);  # choose targeted reference model

@time QM, RM, info2 = ammsyn(sysf, Mre; atol = 1.e-7, reltol = 1.e-4, sdeg = -1);

R3 = fdIFeval(QM,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(R3.sys-RM.sys, atol=1.e-7) &&
      fdimmperf(RM,Mre) ≈ info2.gammasub 

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

Q, Rf, info = efdsyn(sysf, atol = 1.e-7, sdeg =-5, smarg = -0.05, minimal = true,
                     FDfreq = 0, FDGainTol = 0.00001, rdim = 1, simple = false);  info

gammasub = fdimmperf(Rf)
fscond = fdiscond(Rf)[1]
@test gammasub == 0 

sysr = FDFilterIF(dss([1 1]),0,0,mf);

Q1, Rf1, info1 = ammsyn(sysf,sysr, atol = 1.e-7, sdeg =-5, reltol = 1.e-5)

gammasub1 = fdimmperf(Rf1,info1.M*sysr)
fscond1 = fdiscond(Rf1)[1]
@test gammasub <= gammasub1 && fscond <= fscond1

sysr2 = FDFilterIF(dss([1 0]),0,0,mf);

Q2, Rf2, info2 = ammsyn(sysf, sysr2, atol = 1.e-7, sdeg =-5, reltol = 1.e-5)

@test fdimmperf(Rf2,info2.M*sysr2) < 1.e-7 && all(fdiscond(Rf2,fditspec(sysr2))[1] .≈ 1)

sysr2 = FDFilterIF(dss([0 1]),0,0,mf);

Q2, Rf2, info2 = ammsyn(sysf, sysr2, atol = 1.e-7, sdeg =-5, reltol = 1.e-5)

@test fdimmperf(Rf2,info2.M*sysr2) < 1.e-7 && all(fdiscond(Rf2,fditspec(sysr2))[1] .≈ 1)



end # test
end # module
