module Test_amdsyn

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Random
using Test

# Testing examples for EMDSYN 
@testset "amdsyn" begin
#rand()

## Actuator faults
s = rtf('s');
k = 14; N = 4; mu = 1;
sys1 = dss(k/(s+k));           # nominal model
sys2 = dss(0.5*k/(s+0.5*k));   # 50# loss of efficiency
sys3 = dss(10*k/(s+10*k));     # disconnected actuator
sys4 = dss(0.01*k/(s+0.01*k)); # stall load fault model

# setup synthesis model
sysm = mdmodset([sys1, sys2, sys3, sys4], controls = 1:mu);

@time Q, R, info = amdsyn(sysm; sdeg = -15, poles = [-20]); info.MDperf
@time distgap, fpeak = mddist(sysm)
@test sortperm(distgap[1,:]) == sortperm(info.MDperf[1,:]) && 
      sortperm(distgap[2,:]) == sortperm(info.MDperf[2,:]) && 
      sortperm(distgap[3,:]) == sortperm(info.MDperf[3,:]) && 
      sortperm(distgap[4,:]) == sortperm(info.MDperf[4,:])

@time Q1, R1, info1 = amdsyn(sysm; sdeg = -15, poles = [-20], MDfreq = 1); info.MDperf
@time distgap1, fpeak1 = mddist(sysm, MDfreq = 1)
@test sortperm(distgap1[1,:]) == sortperm(info1.MDperf[1,:]) && 
      sortperm(distgap1[2,:]) == sortperm(info1.MDperf[2,:]) && 
      sortperm(distgap1[3,:]) == sortperm(info1.MDperf[3,:]) && 
      sortperm(distgap1[4,:]) == sortperm(info1.MDperf[4,:])

@test all(mdsspec(R) .== 0) && (mdsspec(R,1) .== 1) == mdsspec(R,1)
@time mdgain,fpeak,mind = mdmatch(Q,MDModel(sys2;mu))
@test mind == 2 && argmin(mdgain) == mind


# Example 6.1c - Solution of a EMDP

# Lateral aircraft model without faults
A = [-.4492 0.046 .0053 -.9926;
     0 0 1 0.0067;
     -50.8436 0 -5.2184 .722;
     16.4148 0 .0026 -.6627];
Bu = [0.0004 0.0011; 0 0; -1.4161 .2621; -0.0633 -0.1205]; 
C = eye(4); p = size(C,1); mu = size(Bu,2); 
# define the LOE faults Gamma_i
Gamma = 1 .- [ 0  0 0 .5 .5 .5 1  1 1;
            0 .5 1  0 .5  1 0 .5 1 ]';
N = size(Gamma,1);
# define multiple physical fault model Gui = Gu*Gamma_i
sysu = similar(Vector{DescriptorStateSpace},N)

for i = 1:N
    sysu[i] = dss(A,Bu*diagm(Gamma[i,:]),C,0);
end

# setup synthesis model
sysm = mdmodset(sysu, controls = 1:mu);


# call of AMDSYN with the options for stability degree -1 and pole -1 for
# the filters, tolerance and a design matrix H to form a linear combination
# of the left nullspace basis vectorsH = [ 0.7645 0.8848 0.5778 0.9026 ];
H = [ 0.7645 0.8848 0.5778 0.9026 ];
@time Q, R, info = amdsyn(sysm, sdeg = -1, poles = [-1], HDesign = H); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7)) &&
      all(info.MDgap .== Inf) && info.MDperf ≈ mdperf(R)[1]

@time Q1, R1, info1 = amdsyn(sysm, sdeg = -1, poles = [-1], HDesign = info.HDesign); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))


@time Q, R, info = amdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], HDesign = H); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q, R, info = amdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], rdim = 1); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = amdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], rdim = 1, HDesign = info.HDesign); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))

@time Q, R, info = amdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1]); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = amdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], HDesign = info.HDesign); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))


@time Q, R, info = amdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], minimal = false); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = amdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], minimal = false, HDesign = info.HDesign); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))

@time Q, R, info = amdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], minimal = true, MDfreq = [0,1], HDesign = H); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q, R, info = amdsyn(sysm, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = true, MDfreq = [0,1],simple = true); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && maximum(mdperf(R,MDfreq = [0,1])[1]-info.MDperf) < 1.e-7 &&
      all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q, R, info = amdsyn(sysm, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = false, MDfreq = [0,1]); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && maximum(mdperf(R,MDfreq = [0,1])[1]-info.MDperf) < 1.e-7 &&
      all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q, R, info = amdsyn(sysm, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = true, MDfreq = [0,1], MDGainTol=0.0001); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && maximum(mdperf(R,MDfreq = [0,1])[1]-info.MDperf) < 1.e-7 &&
      all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = amdsyn(sysm, atol = 1.e-7, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = true, MDfreq = [0,1], 
                       HDesign = info.HDesign, MDGainTol=0.0001); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))

@time Q, R, info = amdsyn(sysm, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = true); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && maximum(mdperf(R,MDfreq = [0,1])[1]-info.MDperf) < 1.e-7 &&
      all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = amdsyn(sysm, atol = 1.e-7, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = true, HDesign = info.HDesign); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))


@time Q, R, info = amdsyn(sysm, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = false, MDfreq = [0,1], MDGainTol=0.0001); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && maximum(mdperf(R,MDfreq = [0,1])[1]-info.MDperf) < 1.e-7 &&
      all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = amdsyn(sysm, atol = 1.e-7, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = false, MDfreq = [0,1], 
                       HDesign = info.HDesign, MDGainTol=0.0001); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))


# Example 6.2c - Solution of an AMDP

# Lateral aircraft model without faults
A = [-.4492 0.046 .0053 -.9926;
     0 0 1 0.0067;
     -50.8436 0 -5.2184 .722;
     16.4148 0 .0026 -.6627];
Bu = [0.0004 0.0011; 0 0; -1.4161 .2621; -0.0633 -0.1205]; 
n, mu = size(Bu); p = 2; mw = n+p; m = mu+mw; 
Bw = eye(n,mw);
C = 180/pi*eye(p,n); Du = zeros(p,mu); Dw = [zeros(p,n) eye(p)];
# define the LOE faults Gamma_i
Gamma = 1 .- [ 0  0 0 .5 .5 .5 1  1 1;
            0 .5 1  0 .5  1 0 .5 1 ]';
N = size(Gamma,1);
# define multiple physical fault model Gui = Gu*Gamma_i
sysuw = similar(Vector{DescriptorStateSpace},N)

for i = 1:N
  sysuw[i] = dss(A,[Bu*diagm(Gamma[i,:]) Bw],C,[Du Dw]);
  #sysu[i] = gir(dss(A,Bu*diagm(Gamma[i,:]),C,0),atol = 1.e-7);
end

# setup synthesis model
sysm = mdmodset(sysuw, controls = 1:mu, noise = mu+1:mu+mw);

# use nonminimal design with  AMDSYN
@time Q, R, info = amdsyn(sysm,sdeg = -1, poles = [-1], minimal = false); 

R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true);
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7)) &&
      info.MDgap ≈ mdgap(R)[1] && info.MDperf ≈ mdperf(R)[1] && isstable(Q)

@time Q2 = gbalmr(Q)
R2 = mdIFeval(Q2, sysm, atol=1.e-7, minimal = true);
@test  all(iszero.(Q.sys .- Q2.sys,atol=1.e-7)) && all(iszero.(R.sys .- R2.sys,atol=1.e-7)) 

@time Q2 = gminreal(Q,atol=1.e-7)
R2 = mdIFeval(Q2, sysm, atol=1.e-7, minimal = true);
@test  all(iszero.(Q.sys .- Q2.sys,atol=1.e-7)) && all(iszero.(R.sys .- R2.sys,atol=1.e-7)) 

# use strong nonminimal design with  AMDSYN 
@time Q, R, info = amdsyn(sysm, sdeg = -1, poles = [-1], minimal = false, MDfreq = [0,1]); 

R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true);
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7)) &&
      info.MDgap ≈ mdgap(R, MDfreq = [0,1])[1] && info.MDperf ≈ mdperf(R, MDfreq = [0,1])[1] && isstable(Q)
 

# tests for noise input with varying dimensions
mw = Int.(floor.((n+p+1)*rand(N)))
for i = 1:N
  sysuw[i] = dss(A,[Bu*diagm(Gamma[i,:]) Bw[:,1:mw[i]]],C,[Du Dw[:,1:mw[i]]]);
  #sysu[i] = gir(dss(A,Bu*diagm(Gamma[i,:]),C,0),atol = 1.e-7);
end

sysm = MDMModel(sysuw; mu, mw); 

# use nonminimal design with  AMDSYN
@time Q, R, info = amdsyn(sysm, atol = 1.e-7, sdeg = -1, poles = [-1], minimal = false); 

R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true);
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7)) &&
      isapprox(info.MDgap,mdgap(R)[1],atol=1.e-5) && all(abs.(info.MDperf - mdperf(R,atol=1.e-7)[1]) .< 0.001) 

@time Q, R, info = amdsyn(sysm, atol = 1.e-7, sdeg = -1, poles = [-1], minimal = true); 

R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true);
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7)) &&
isapprox(info.MDgap,mdgap(R)[1],atol=1.e-5) && all(abs.(info.MDperf - mdperf(R,atol=1.e-7)[1]) .< 0.001) 


for nonstd = 2:5
   @time Q, R, info = amdsyn(sysm, atol = 1.e-7, sdeg = -1, poles = [-1], minimal = true, nonstd = nonstd); 

   R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true);
   @test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7)) &&
      isapprox(info.MDgap,mdgap(R)[1],atol=1.e-5) && 
      all(abs.(info.MDperf - mdperf(R,atol=1.e-7)[1]) .< 0.001) 
end

end # test amdsyn



end # module
