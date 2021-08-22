module Test_fditools

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Testing examples for EFDSYN 
println("Test_fditools")
@testset "fditools" begin

# setup the system with additive faults
# [Gu(s) Gd(s) Gf(s)], where Gf(s) = [ Gu(s) I]
s = rtf('s'); # define the Laplace variable s
Gu = [(s+1)/(s+2); (s+2)/(s+3)] # enter Gu(s)
Gd = [(s-1)/(s+2); 0]; # enter Gd(s)
# build state space model of [Gu(s) Gd(s) Gf(s)] and set input groups
sysf = fdimodset(dss([Gu Gd]),c =1,d = 2,f = 1,fs = 1:2)

# tests FDFilter
filter = rss(3,1,3)
Q = FDFilter(filter,1,2)
Q2 = FDFilter(filter; outputs = [1],controls = [2,3])
@test iszero(Q.sys-Q2.sys,atol=1.e-7) && 
            Q.controls == Q2.controls && Q.outputs == Q2.outputs
Q3 = FDFilter(filter,outputs = 1,controls = 2:3)
@test iszero(Q.sys-Q3.sys,atol=1.e-7) && 
            Q.controls == Q3.controls && Q.outputs == Q3.outputs


# tests FDFilterIF
filteri = rss(3,1,6)
R = FDFilterIF(filteri,controls=1:2,disturbances=[3],aux=6,noise=5,faults=[4])
R1 = FDFilterIF(filteri,2,1,1,1,1)
@test iszero(R.sys-R1.sys,atol=1.e-7) && 
       R.controls == R1.controls && R.disturbances == R1.disturbances &&
       R.faults == R1.faults && R.noise == R1.noise &&
       R.aux == R1.aux
            
R = FDFilterIF(filteri,aux=[6],noise=[5],faults=[4])
R1 = FDFilterIF(filteri,0,0,1,1,1;moff = 3)
@test iszero(R.sys-R1.sys,atol=1.e-7) && 
       R.controls == R1.controls && R.disturbances == R1.disturbances &&
       R.faults == R1.faults && R.noise == R1.noise &&
       R.aux == R1.aux
R2 = FDFilterIF(filteri[:,4:end],0,0,1,1,1)
@test iszero(R.sys-R2.sys,atol=1.e-7) && 
       R.controls == R2.controls && R.disturbances == R2.disturbances &&
       R.faults == R2.faults && R.noise == R2.noise &&
       R.aux == R2.aux

sysf1 = fdimodset(rss(3,1,6,stable=true),c = 1:2, d = 3, f = 4:5, n = 6, aux = 2:6)
Q = FDFilter(rss(3,1,3,stable=true),1,2)
R = fdIFeval(Q, sysf1)
syse = [sysf1.sys;eye(2,11)]
@test iszero(Q.sys*syse[:,1:2]-R.sys[:,R.controls],atol=1.e-7) && 
      iszero(Q.sys*syse[:,3]-R.sys[:,R.disturbances],atol=1.e-7) &&
      iszero(Q.sys*syse[:,4:5]-R.sys[:,R.faults],atol=1.e-7) &&
      iszero(Q.sys*syse[:,6]-R.sys[:,R.noise],atol=1.e-7) &&
      iszero(Q.sys*syse[:,7:end]-R.sys[:,R.aux],atol=1.e-7)
@test fditspec(R) == fditspec(R.sys[:,R.faults])
@time scond, β, γ = fdscond(R)
@test β == fdhinfminus(R.sys[:,R.faults])[1] && γ == fdhinfmax(R.sys[:,R.faults])[1] 
freq = 0
@time scond, β, γ = fdscond(R,freq)
@test β == fdhinfminus(R.sys[:,R.faults],freq)[1] && γ == fdhinfmax(R.sys[:,R.faults],freq)[1] 
freq = [0, 3, 5]; 
@time scond, β, γ = fdscond(R,freq)
@test β == fdhinfminus(R.sys[:,R.faults],freq)[1] && γ == fdhinfmax(R.sys[:,R.faults],freq)[1] 



@test fditspec(sysf.sys[:,sysf.faults],block=true,FDtol = 1.e-3) == ([1  1  1] .> 0)
@test fditspec(sysf.sys[:,sysf.faults],block=true,FDfreq=0)[:,:,1] == ([1  1  1] .> 0)
@test fditspec(sysf.sys[:,sysf.faults],block=false) == ([1  1  0; 1  0  1] .> 0)
@test fditspec(sysf.sys[:,sysf.faults],block=false,FDfreq=0)[:,:,1] == ([1  1  0; 1  0  1] .> 0)
@test fditspec(sysf.sys[:,sysf.faults],block=false,FDfreq=[0,1,Inf])[:,:,1] == ([1  1  0; 1  0  1] .> 0)

@time S, gains = fdisspec(sysf.sys[:,sysf.faults],[0,1,Inf]; block=false)
@test S[:,:,1] == ([1  1  0; 1  0  1] .> 0) && 
     gains[:,:,1] ≈  [0.5  1.0  0.0; 2/3  0.0  1.0] && gains[:,:,3] ≈  [1.0 1.0 0.0; 1.0 0.0 1.0] && gains[:,:,2] ≈  [0.6324555320336759 1.0 0.0; 0.7071067811865475 0.0 1.0]

@time S1, gains1 = fdisspec(sysf.sys[:,sysf.faults],[0,1,Inf]; block=true)
@test S1[:,:,1] == Bool[1 1 1] && gains1[:,:,1] ≈  [0.8333333333333334 1.0 1.0]

end # fditools


end