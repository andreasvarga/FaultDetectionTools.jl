module Test_emdsyn

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Random
using Test

# Testing examples for EMDSYN 
@testset "emdsyn" begin
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
display(sysm)

@time distgap, fpeak = mddist(sysm)
@time distgapc, fpeakc = mddist2c(sysm,sysm)
@test norm(distgap-distgapc,Inf) < 1.e-7 

@time distgap1, fpeak1 = mddist(sysm, MDfreq = 1)
@time distgap1c, fpeak1c = mddist2c(sysm, sysm, MDfreq = 1)
@test norm(distgap1-distgap1c,Inf) < 1.e-7 

@time distgap2, fpeak2 = mddist(sysm, distance = "2")
@time distgap2c, fpeak2c = mddist2c(sysm, sysm, distance = "2")
@test norm(distgap2-distgap2c,Inf) < 1.e-7 

@time distgap3, fpeak3 = mddist(sysm, distance = "Inf")
@time distgap3c, fpeak3c = mddist2c(sysm, sysm, distance = "Inf")
@test norm(distgap3-distgap3c,Inf) < 1.e-7 

@time Q, R, info = emdsyn(sysm; sdeg = -15, poles = [-20]); info.MDperf
display(Q)
display(R)
gpole(Q)
@test sortperm(distgap[1,:]) == sortperm(info.MDperf[1,:]) && 
      sortperm(distgap[2,:]) == sortperm(info.MDperf[2,:]) && 
      sortperm(distgap[3,:]) == sortperm(info.MDperf[3,:]) && 
      sortperm(distgap[4,:]) == sortperm(info.MDperf[4,:])

@time Q1, R1, info1 = emdsyn(sysm; sdeg = -15, poles = [-20], MDfreq = 1); info.MDperf
@test sortperm(distgap1[1,:]) == sortperm(info1.MDperf[1,:]) && 
      sortperm(distgap1[2,:]) == sortperm(info1.MDperf[2,:]) && 
      sortperm(distgap1[3,:]) == sortperm(info1.MDperf[3,:]) && 
      sortperm(distgap1[4,:]) == sortperm(info1.MDperf[4,:])

@test (mdspec(R,atol = 1.e-7) .== 1) == mdsspec(R,1)
@test all(mdsspec(R) .== 0) && (mdsspec(R,1) .== 1) == mdsspec(R,1)
@time mdgain,fpeak,mind = mdmatch(Q,MDModel(sys2;mu))
@test mind == 2 && argmin(mdgain) == mind

@time Q, R, info = emdsyn(sysm; sdeg = -15, poles = [-20], MDSelect = [1]); info.MDperf
@time Q, R, info = emdsyn(sysm; sdeg = -15, poles = [-20], MDSelect = [2]); info.MDperf 
@time Q, R, info = emdsyn(sysm; sdeg = -15, poles = [-20], normalize = true); info.MDperf 
@test info.MDperf[1,:] == info.MDperf[:,1]

sysc1=MDModel(rss(3,2,6,stable = true); mu = 2, md = 1,mw = 2, ma = 1)
display(sysc1)
sysc2=MDModel(rss(3,2,6,stable = true); mu = 2, md = 2, mw = 1)
sysm = mdmodset([sysc1,sysc2])
@test mddist2c([sysc1,sysc2],[sysc1,sysc2])[1] == mddist2c([sysc1,sysc2],sysm)[1] == mddist2c(sysm,[sysc1,sysc2])[1]
@test mddist([sysc1,sysc2],distance = "Inf")[1] == mddist(sysm,distance = "Inf")[1] 

# check models with variable number of disturbances
p = 4; mu = 2; md = [0,3,1,2];
sys1 = rss(3,p,mu+md[1],stable = true);          
sys2 = rss(1,p,mu+md[2],stable = true);
sys3 = rss(2,p,mu+md[3],stable = true);
sys4 = rss(0,p,mu+md[4],stable = true);

# setup synthesis model
sysm = mdmodset([sys1, sys2, sys3, sys4], controls = 1:mu, 
        disturbances = [mu .+ Vector(1:md[1]), mu .+ Vector(1:md[2]),mu .+ Vector(1:md[3]),mu .+ Vector(1:md[4])])
sysm1 = mdmodset([sys1, sys2, sys3, sys4], controls = 1:mu, 
        disturbances = [mu .+ (1:md[1]), mu .+ (1:md[2]),mu .+ (1:md[3]),mu .+ (1:md[4])]);
sysm2 = MDMModel([sys1, sys2, sys3, sys4]; mu, md);

@test all(iszero.(sysm.sys .- sysm1.sys,atol=1.e-7)) && sysm.mu == sysm1.mu &&
      sysm.md == sysm1.md && sysm.mw == sysm1.mw && sysm.ma == sysm1.ma

@test all(iszero.(sysm.sys .- sysm2.sys,atol=1.e-7)) && sysm.mu == sysm2.mu &&
      sysm.md == sysm2.md && sysm.mw == sysm2.mw && sysm.ma == sysm2.ma


@time distgap, fpeak = mddist(sysm)
@time distgap1, fpeak1 = mddist(sysm, MDfreq = fpeak[:])
@test distgap ≈ distgap1 && fpeak ≈ fpeak1

@time distgapd, fpeakd = mddist(sysm, cdinp = true)
@time distgapd1, fpeakd1 = mddist(sysm, cdinp = true, MDfreq = fpeakd[:])
@test distgapd ≈ distgapd1 && fpeakd ≈ fpeakd1 && all(abs.(distgapd-distgap+1.e-7*I) .>= 0)

@time distgapc, fpeakc = mddist2c(sysm,sysm)
@time distgapc1, fpeakc1 = mddist2c(sysm, sysm, MDfreq = fpeakc[:])
@test norm(distgapc-distgapc1,Inf) < 1.e-7 

@time distgapcd, fpeakcd = mddist2c(sysm,sysm, cdinp = true)
@time distgapcd1, fpeakcd1 = mddist2c(sysm, sysm, cdinp = true, MDfreq = fpeakcd[:])
@test norm(distgapcd-distgapcd1,Inf) < 1.e-7 && all(abs.(distgapcd-distgapc+1.e-7*I) .>= 0)


@time distgap, fpeak = mddist(sysm, distance= "2")
@time distgapc, fpeakc = mddist2c(sysm, sysm, distance= "2")
@test isapprox(distgapc,distgap,atol = 1.e-7) && fpeakc ≈ fpeak && all(diag(distgapc) .< 1.e-7)

@time distgapd, fpeakd = mddist(sysm, cdinp = true, distance= "2")
@time distgapcd, fpeakcd = mddist2c(sysm, sysm, cdinp = true, distance= "2")
@test isapprox(distgapcd,distgapd,atol = 1.e-7) && fpeakcd ≈ fpeakd && all(diag(distgapcd) .< 1.e-7) 


@time distgap, fpeak = mddist(sysm, distance = "Inf")
@time distgap1, fpeak1 = mddist(sysm, distance = "Inf", MDfreq = fpeak[:])
@test distgap ≈ distgap1 && fpeak ≈ fpeak1

@time distgapd, fpeakd = mddist(sysm, distance = "Inf", cdinp = true)
@time distgapd1, fpeakd1 = mddist(sysm, distance = "Inf", cdinp = true, MDfreq = fpeakd[:])
@test distgapd ≈ distgapd1 && fpeakd ≈ fpeakd1 && all(abs.(distgapd-distgap+1.e-7*I) .>= 0)

@time distgapc, fpeakc = mddist2c(sysm,sysm, distance = "Inf")
@time distgapc1, fpeakc1 = mddist2c(sysm, sysm, distance = "Inf", MDfreq = fpeakc[:])
@test norm(distgapc-distgapc1,Inf) < 1.e-7 

@time distgapcd, fpeakcd = mddist2c(sysm,sysm, cdinp = true, distance = "Inf")
@time distgapcd1, fpeakcd1 = mddist2c(sysm, sysm, cdinp = true, distance = "Inf", MDfreq = fpeakcd[:])
@test norm(distgapcd-distgapcd1,Inf) < 1.e-7 && all(abs.(distgapcd-distgapc+1.e-7*I) .>= 0)


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


# call of EMDSYN with the options for stability degree -1 and pole -1 for
# the filters, tolerance and a design matrix H to form a linear combination
# of the left nullspace basis vectorsH = [ 0.7645 0.8848 0.5778 0.9026 ];
H = [ 0.7645 0.8848 0.5778 0.9026 ];
@time Q, R, info = emdsyn(sysm, sdeg = -1, poles = [-1], HDesign = H); 
display(Q)
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = emdsyn(sysm, sdeg = -1, poles = [-1], HDesign = info.HDesign); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))


@time Q, R, info = emdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], HDesign = H); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q, R, info = emdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], rdim = 1); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = emdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], rdim = 1, HDesign = info.HDesign); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))

@time Q, R, info = emdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1]); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = emdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], HDesign = info.HDesign); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))


@time Q, R, info = emdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], minimal = false); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = emdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], minimal = false, HDesign = info.HDesign); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))

@time Q, R, info = emdsyn(sysm, nullspace = true, sdeg = -1, poles = [-1], minimal = true, MDfreq = [0,1], HDesign = H); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q, R, info = emdsyn(sysm, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = true, MDfreq = [0,1],simple = true); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && maximum(mdperf(R,MDfreq = [0,1])[1]-info.MDperf) < 1.e-7 &&
      all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q, R, info = emdsyn(sysm, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = false, MDfreq = [0,1]); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && maximum(mdperf(R,MDfreq = [0,1])[1]-info.MDperf) < 1.e-7 &&
      all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q, R, info = emdsyn(sysm, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = true, MDfreq = [0,1], MDGainTol=0.0001); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && maximum(mdperf(R,MDfreq = [0,1])[1]-info.MDperf) < 1.e-7 &&
      all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = emdsyn(sysm, atol = 1.e-7, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = true, MDfreq = [0,1], 
                       HDesign = info.HDesign, MDGainTol=0.0001); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))

@time Q, R, info = emdsyn(sysm, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = true); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && maximum(mdperf(R,MDfreq = [0,1])[1]-info.MDperf) < 1.e-7 &&
      all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = emdsyn(sysm, atol = 1.e-7, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = true, HDesign = info.HDesign); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))


@time Q, R, info = emdsyn(sysm, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = false, MDfreq = [0,1], MDGainTol=0.0001); 
R1 = mdIFeval(Q, sysm, atol=1.e-7, minimal = true)
@test norm(diag(info.MDperf)) < 1.e-7 && maximum(mdperf(R,MDfreq = [0,1])[1]-info.MDperf) < 1.e-7 &&
      all(iszero.(R.sys .- R1.sys,atol=1.e-7))

@time Q1, R1, info1 = emdsyn(sysm, atol = 1.e-7, nullspace = true, sdeg = -1, smarg = -1, poles = [-1], minimal = false, MDfreq = [0,1], 
                       HDesign = info.HDesign, MDGainTol=0.0001); 
R11 = mdIFeval(Q1, sysm, atol=1.e-7, minimal = true)
@test norm(info.MDperf-info1.MDperf) < 1.e-7 && all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) &&
all(iszero.(R.sys .- R11.sys,atol=1.e-7))




end # test emdsyn



end # module
