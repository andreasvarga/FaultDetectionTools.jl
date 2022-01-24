module Test_fditools

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Testing examples for EFDSYN 
println("Test_fdiutils")
@testset "fdiutils" begin

# setup the system with additive faults
# [Gu(s) Gd(s) Gf(s)], where Gf(s) = [ Gu(s) I]
s = rtf('s'); # define the Laplace variable s
Gu = [(s+1)/(s+2); (s+2)/(s+3)] # enter Gu(s)
Gd = [(s-1)/(s+2); 0]; # enter Gd(s)
# build state space model of [Gu(s) Gd(s) Gf(s)] and set input groups
sysf = fdimodset(dss([Gu Gd]),c =1,d = 2,f = 1,fs = 1:2)
display(sysf)

# tests FDFilter
filter = rss(3,1,3)
Q = FDFilter(filter,1,2)
display(Q)

# tests FDIFilter
filter = [rss(3,1,3), rss(2,2,3)]
Q = FDIFilter(filter,1,2)
display(Q)

# tests FDFilterIF
filteri = rss(3,1,6)
R = FDFilterIF(filteri,2,1,1,1,1);
display(R)

R1 = FDFilterIF(filteri; mu = 2, md = 1, mf = 1, mw = 1, ma = 1);
@test iszero(R1.sys-R.sys,atol=1.e-7)
           
R1 = FDFilterIF(filteri; mf = 1, mw = 1, ma = 1, moff = 3);
R2 = FDFilterIF(filteri[:,4:end],0,0,1,1,1);
@test iszero(R1.sys-R2.sys,atol=1.e-7)


# tests FDIFilterIF
filteri = [rss(3,1,6), rss(2,2,6)]
R = FDIFilterIF(filteri,2,1,1,1,1);
display(R)

R1 = FDIFilterIF(filteri; mu = 2, md = 1, mf = 1, mw = 1, ma = 1);
@test all(iszero.(R1.sys .- R.sys, atol=1.e-7))


R2 = FDIFilterIF(filteri; mf = 1, mw = 1, ma = 1, moff = 3);
display(R2)

sysf1 = fdimodset(rss(3,1,6,stable=true),c = 1:2, d = 3, f = 4:5, n = 6, aux = 2:6)
Q = FDFilter(rss(3,1,3,stable=true),1,2)
R = fdIFeval(Q, sysf1)
syse = [sysf1.sys;eye(2,11)]
@test iszero(Q.sys*syse[:,1:2]-R.sys[:,R.controls],atol=1.e-7) && 
      iszero(Q.sys*syse[:,3]-R.sys[:,R.disturbances],atol=1.e-7) &&
      iszero(Q.sys*syse[:,4:5]-R.sys[:,R.faults],atol=1.e-7) &&
      iszero(Q.sys*syse[:,6]-R.sys[:,R.noise],atol=1.e-7) &&
      iszero(Q.sys*syse[:,7:end]-R.sys[:,R.aux],atol=1.e-7)
@test fditspec(R) == fditspec_(R.sys[:,R.faults])
@time scond, β, γ = fdiscond(R)
@test β == fdhinfminus(R.sys[:,R.faults])[1] && γ == fdhinfmax(R.sys[:,R.faults])[1] 
freq = 0
@time scond, β, γ = fdiscond(R,freq)
@test β == fdhinfminus(R.sys[:,R.faults],freq)[1] && γ == fdhinfmax(R.sys[:,R.faults],freq)[1] 
freq = [0, 3, 5]; 
@time scond, β, γ = fdiscond(R,freq)
@test β == fdhinfminus(R.sys[:,R.faults],freq)[1] && γ == fdhinfmax(R.sys[:,R.faults],freq)[1] 

sysf1 = fdimodset(rss(3,1,6,stable=true),c = 1:2, d = 3, f = 4:5, n = 6, aux = 2:6)
Q = FDIFilter([rss(3,1,3,stable=true),rss(2,2,3,stable=true)],1,2)
R = fdIFeval(Q, sysf1)
syse = [sysf1.sys;eye(2,11)]
@test all(iszero.(Q.sys .* [syse]-R.sys,atol=1.e-7))

@time S = fdisspec(R,[0,1]; FDGainTol = 1.e-5)
@test S ==  ([1  1 ; 1  1] .> 0) 
@time S1 = fditspec(R; FDfreq = [0,1], FDtol = 1.e-5)
@test S == S1
@time S2 = fditspec(R)
@test S2 == ([1  1 ; 1  1] .> 0) 



@test fditspec_(sysf.sys[:,sysf.faults],block=true,FDtol = 1.e-3) == ([1  1  1] .> 0)
@test fditspec_(sysf.sys[:,sysf.faults],block=true,FDfreq=0)[:,:,1] == ([1  1  1] .> 0)
@test fditspec_(sysf.sys[:,sysf.faults],block=true,FDfreq=0)[:,:,1] == ([1  1  1] .> 0)
@test fditspec_(sysf.sys[:,sysf.faults],block=false) == ([1  1  0; 1  0  1] .> 0)
@test fditspec_(sysf.sys[:,sysf.faults],block=false,FDfreq=0)[:,:,1] == ([1  1  0; 1  0  1] .> 0)
@test fditspec_(sysf.sys[:,sysf.faults],block=false,FDfreq=[0,1,Inf])[:,:,1] == ([1  1  0; 1  0  1] .> 0)

sysf1 = fdimodset(rss(3,1,6,stable=true),c = 1:2, d = 3, f = 4:5, n = 6, aux = 2:6)
sysf2 = fdimodset(rss(2,2,6,stable=true),c = 1:2, d = 3, f = 4:5, n = 6, aux = 2:6)


@time S, gains = fdisspec_(sysf.sys[:,sysf.faults],[0,1,Inf]; block=false)
@test S[:,:,1] == ([1  1  0; 1  0  1] .> 0) && 
     gains[:,:,1] ≈  [0.5  1.0  0.0; 2/3  0.0  1.0] && gains[:,:,3] ≈  [1.0 1.0 0.0; 1.0 0.0 1.0] && gains[:,:,2] ≈  [0.6324555320336759 1.0 0.0; 0.7071067811865475 0.0 1.0]

@time S1, gains1 = fdisspec_(sysf.sys[:,sysf.faults],[0,1,Inf]; block=true)
@test S1[:,:,1] == Bool[1 1 1] && gains1[:,:,1] ≈  [0.8333333333333334 1.0 1.0]

# some non-trivial examples for frequency values equal to poles
# S1 = [1 0]
s = rtf('s')
g = dss([1/(s^2+1) 1/s]);
S1, g1 = fdisspec_(g, [0,1]; block=true,stabilize = true)
S2 = fditspec_(g,FDfreq = [0,1],block=true,poleshift = true)
@test S1 == S2

S1, g1 = fdisspec_(g, [0,1]; block=false,stabilize = true)
S2 = fditspec_(g,FDfreq =[0,1],block=false,poleshift = true)
@test S1 == S2

# S1 = [1 0]
s = rtf('s')
g = dss([1/s 1/(s+2)]);
S1, g1 = fdisspec_(g,block=true,stabilize = true)
S2 = fditspec_(g,FDfreq = 0,block=true,poleshift = true)
@test S1 == S2

S1, g1 = fdisspec_(g,block=false,stabilize = true)
S2 = fditspec_(g,FDfreq = 0,block=false,poleshift = true)
@test S1 == S2

# S1 = [1 1]
s = rtf('s')
g = dss([1/s 1/s/(s+2)]);
S1, g1 = fdisspec_(g,block=true,stabilize = true)
S2 = fditspec_(g,FDfreq = 0,block=true,poleshift = true)
@test S1 == S2

S1, g1 = fdisspec_(g,block=false,stabilize = true)
S2 = fditspec_(g,FDfreq = 0,block=false,poleshift = true)
@test S1 == S2

s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [s/(s+1) 0; (s^2+4)/(s+3)^2 1]; # enter Gf(s)
# new information in smat and zero local H-minus index values

smat, gains = fdisspec_(dss(Gf), [0, 2]; FDGainTol = 1.e-3) 
@test smat[:,:,1] == Bool[0 0; 1 1] && smat[:,:,2] == Bool[1 0; 0 1] && 
      gains[:,:,1] ≈ [0.0 0.0; 0.4444444444444444 1.0] && gains[:,:,2] ≈ [0.894427190999916 0.0; 0.0 1.0]
ssmat = fditspec_(dss(Gf), atol = 1.e-7,FDtol = 1.e-3, FDfreq = [0, 2]) 
@test smat == ssmat
swmat = fditspec_(dss(Gf), atol=1.e-7,FDtol = 1.e-3) 
@test swmat == Bool[1 0; 1 1] 

s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [1/s; 1/(s^2+4)]; # enter Gf(s)

smat, gains = fdisspec_(dss(Gf), [0, 2]; FDGainTol = 1.e-3, stabilize = true) 
@test smat[:,:,1] == trues(2,1) && smat[:,:,2] == trues(2,1) && all(maximum(gains,dims=1) .> 0.001)
ssmat = fditspec_(dss(Gf), atol = 1.e-7,FDtol = 1.e-3, FDfreq = [0, 2],poleshift = true) 
@test smat == ssmat
swmat = fditspec_(dss(Gf), atol=1.e-7,FDtol = 1.e-3) 
@test swmat == trues(2,1) 

s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [s; s^2+4]; # enter Gf(s)

smat, gains = fdisspec_(dss(Gf), [0, 2]; FDGainTol = 1.e-3,stabilize = true) 
@test smat[:,:,1] == eye(2)[:,2:2] && smat[:,:,2] == eye(2)[:,1:1] && 
      gains[1,1,1] < 1.e-7 && gains[2,1,2] < 1.e-7
ssmat = fditspec_(dss(Gf), atol = 1.e-7,FDtol = 1.e-3, FDfreq = [0, 2],poleshift = true) 
@test smat == ssmat
swmat = fditspec_(dss(Gf), atol=1.e-7,FDtol = 1.e-3) 
@test swmat ==trues(2,1) 



end # fditools


end