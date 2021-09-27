module Test_fdigenspec

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Random
using Test

# Testing examples for fdigenspec
println("Test_fdigenspec")
@testset "fdigenspec & fdichkspec" begin

p = 4; mu = 0; md = 0; mf = 4; n = 0;
sysf = fdimodset(rss(n,p,mf),faults = 1:mf);
@time S_weak = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5)
@test size(S_weak,1) == 15

@time S_strong = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5, FDGainTol = 1.e-3, FDfreq=[0,1])  
@test S_weak == S_strong

@time rdims, orders, leastorders = fdichkspec(sysf, S_weak; atol = 1.e-7, FDtol = 1.e-5)
@test all(rdims .> 0) && all(orders .== 0)
@time rdims, orders, leastorders = fdichkspec(sysf, S_weak; atol = 1.e-7, FDtol = 1.e-5, minimal = true)
@test all(rdims .> 0) && all(orders .== 0) && all(leastorders .== 0)

## Example without control and disturbance inputs
p = 4; mu = 0; md = 0; mf = 4; n = 2;
sysf = fdimodset(rss(n,p,mf),faults = 1:mf);
@time S_weak = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5)
@test size(S_weak,1) == 15

@time S_strong = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5, FDGainTol = 1.e-3, FDfreq=[0,1])  
@test S_weak == S_strong

@time rdims, orders, leastorders = fdichkspec(sysf, S_weak; atol = 1.e-7, FDtol = 1.e-5)
@test all(rdims .> 0) && all(orders .== n)

@time rdims, orders, leastorders = fdichkspec(sysf, S_weak; atol = 1.e-7, FDtol = 1.e-5, minimal = true)
@test all(rdims .> 0) && all(orders .== n) && all(leastorders .<= n)

## Yuan et al. IJC (1997)  weak synthesis with pole assignment and
#  enforcing least global order of 6!
p = 3; mu = 1; md = 0; mf = 8; n = 4;
A = [
-1 1 0 0
1 -2 1 0
0 1 -2 1
0 0 1 -2  ];
Bu = [1 0 0 0]';
Bf = [
1 0 0 0 1 0 0 0
0 1 0 0 -1 1 0 0
0 0 1 0 0 -1 1 0
0 0 0 1 0 0 -1 1];    
C = [
 1 0 0 0
0 0 1 0
0 0 0 1];
Du = zeros(p,mu); Df = zeros(p,mf);

sys = dss(A,[Bu Bf],C,[Du Df]);         
sysf = fdimodset(sys, c = 1:mu, d = mu .+ (1:md) , f = (mu+md) .+ (1:mf));

@time S_weak = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5)
@test size(S_weak,1) == 18 && 
sort(binmat2dec(S_weak)) == sort(binmat2dec([
0   0   0   1   0   0   1   1
0   1   1   0   1   1   1   0
0   1   1   1   1   1   0   1
0   1   1   1   1   1   1   1
1   0   1   0   1   1   1   0
1   0   1   1   1   1   0   1
1   0   1   1   1   1   1   1
1   1   0   0   1   1   0   0
1   1   0   1   1   1   1   1
1   1   1   0   0   1   1   0
1   1   1   0   1   0   1   0
1   1   1   0   1   1   1   0
1   1   1   1   0   1   0   1
1   1   1   1   0   1   1   1
1   1   1   1   1   0   0   1
1   1   1   1   1   0   1   1
1   1   1   1   1   1   0   1
1   1   1   1   1   1   1   1
]))

@time S_strong = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5, FDGainTol = 1.e-3, FDfreq=[0,1])  
@test size(S_strong,1) == 12 && 
sort(binmat2dec(S_strong)) == sort(binmat2dec([
      0   0   0   1   0   0   1   1
      0   1   1   0   1   1   1   0
      0   1   1   1   1   1   0   1
      0   1   1   1   1   1   1   1
      1   0   1   0   1   1   1   0
      1   0   1   1   1   1   0   1
      1   0   1   1   1   1   1   1
      1   1   0   0   1   1   0   0
      1   1   0   1   1   1   1   1
      1   1   1   0   1   1   1   0
      1   1   1   1   1   1   0   1
      1   1   1   1   1   1   1   1
   ]))

@time rdims, orders, leastorders = fdichkspec(sysf, S_weak; atol = 1.e-7, FDtol = 1.e-5, minimal = true)
@test all(rdims .> 0) && all(orders .<= n)  

@time rdims, orders, leastorders = fdichkspec(sysf, S_weak; FDfreq=[0,1], atol = 1.e-7, FDtol = 1.e-5, minimal = true)
@test sort(binmat2dec(S_strong)) == sort(binmat2dec(S_weak[rdims .> 0,:])) 

# some timing results
@time S_weak = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5);
Stest = dec2binmat(Vector(1:2^mf-1));
@time rdims, orders, leastorders = fdichkspec(sysf, Stest; atol = 1.e-7, FDtol = 1.e-5, minimal = true);
@test size(S_weak,1) == 18 && sort(binmat2dec(S_weak)) == sort(binmat2dec(Stest[rdims .> 0,:]))

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

@time S_weak = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5)
@test S_weak == Bool[1 1]

@time S_strong = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5, FDGainTol = 1.e-3, FDfreq=[0,1])  
@test S_strong == Bool[1 1]

Stest = dec2binmat(Vector(1:2^mf-1));

@time rdims, orders, leastorders = fdichkspec(sysf, Stest; atol = 1.e-7, FDtol = 1.e-5, minimal = true)
@test rdims == [0, 0, 1] && orders == [-1, -1, 1] && leastorders == [-1, -1, 1]


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

@time S_weak = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5)
@test S_weak == Bool[1 1; 0 1; 1 0]

@time S_strong = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5, FDGainTol = 1.e-3, FDfreq=[0,1])  
@test S_strong == Bool[1 1; 0 1; 1 0]


Stest = dec2binmat(Vector(1:2^mf-1));

@time rdims, orders, leastorders = fdichkspec(sysf, Stest; atol = 1.e-7, FDtol = 1.e-5, minimal = true)
@test rdims == [3, 3, 4] && orders == [3, 3, 4] && leastorders == [1, 1, 1]

## Example 1, Xiong & Saif 2000, Int. J. Robust Nonlinear Control
n = 4; mu = 2; md = 2; mf = 2; m = mu+md+mf; p = 4;
a = [ -9.9477 -0.7476 0.2632 -5.0337
      52.1659 2.7452 5.5532 -24.4221
      26.0922 2.6361 -4.1975 -19.2774
      0.0 0.0 1.0 0.0];
bu = [ 0.4422 0.1761
       3.5446 -7.5922
        -5.52 4.49
        0.0 0.0];
bd = [ 0 0; 0 1 ; 1 0; 0 0 ];
bf = [bu[:,1] zeros(n,1)];
c = [ 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 1 1 1]; 
dd = zeros(n,md); du = zeros(n,mu); df = [zeros(n,mf-1) [0 0 1 0]'];

sys = dss(a,[bu bd bf],c,[du dd df]);
sysf = fdimodset(sys, c = 1:mu, d = mu .+ (1:md) , f = (mu+md) .+ (1:mf));


@time S_weak = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5)
@test S_weak == Bool[1 1; 0 1; 1 0]

@time S_strong = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5, FDGainTol = 1.e-3, FDfreq=[0,1])  
@test S_strong == Bool[1 1; 0 1; 1 0]


Stest = dec2binmat(Vector(1:2^mf-1));

@time rdims, orders, leastorders = fdichkspec(sysf, Stest; atol = 1.e-7, FDtol = 1.e-5, minimal = true)
@test rdims == [1, 1, 2] && orders == [1, 2, 2] && leastorders == [1, 2, 1]


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

@time S_weak = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5)
@test S_weak == Bool[1 1; 0 1; 1 0]

@time S_strong = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5, FDGainTol = 1.e-3, FDfreq=[0,1])  
@test S_strong == Bool[1 1; 0 1; 1 0]


Stest = dec2binmat(Vector(1:2^mf-1));

@time rdims, orders, leastorders = fdichkspec(sysf, Stest; atol = 1.e-7, FDtol = 1.e-5, minimal = true)
@test rdims == [1, 1, 2] && orders == [1, 1, 2] && leastorders == [1, 1, 1]


## Linearized Boeing-747 (Varga,2007)

A = [
-0.4861 0.000317 -0.5588 0 -2.04e-6
0 -0.0199 3.0796 -9.8048 8.98e-5
1.0053 -0.0021 -0.5211 0 9.30e-6
1 0 0 0 0
0 0 -92.6 92.6 0];

Bu = [
-0.1455 -0.1455 -0.1494 -0.1494 -1.2860 0.0013 0.0035 0.0035 0.0013
0 0 0 0 -0.3122 0.1999 0.1999 0.1999 0.1999
-0.0071 -0.0071 -0.0074 -0.0074 -0.0676 -0.0004 -0.0004 -0.0004 -0.0004
0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0];

C = [
0 0 1 0 0
0 -0.0199 3.0796 -9.8048 8.98e-5
0 0 0 1 0
1 0 0 0 0
0 0 -92.6 92.6 0
0 0 0 0 1 ];   

Du = [
 0 0 0 0 0 0 0 0 0
0 0 0 0 -0.3122 0.1999 0.1999 0.1999 0.1999
0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0  ];

Bf = Bu[:,1:5]; Df = Du[:,1:5];

p = 6; n = 5; mu = 9; mf = 5; md = 0;

sys = dss(A,[Bu Bf],C,[Du Df]);    

sysf = fdimodset(sys, c = 1:mu, d = mu .+ (1:md) , f = (mu+md) .+ (1:mf));

@time S_weak = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5)
@test sort(binmat2dec(S_weak)) == sort(binmat2dec([
      0   0   0   0   1
      0   0   1   1   0
      0   0   1   1   1
      1   1   0   0   0
      1   1   0   0   1
      1   1   1   1   0
      1   1   1   1   1
      ]))

@time S_strong = fdigenspec(sysf, atol = 1.e-8, FDtol = 1.e-5, FDGainTol = 0.00001, FDfreq=[0,1])  
@test sort(binmat2dec(S_strong)) == sort(binmat2dec([
    1  1  1  1  1
    0  0  1  1  1
    1  1  0  0  1
    1  1  1  1  0
         ]))


Stest = dec2binmat(Vector(1:2^mf-1));

@time rdims, orders, leastorders = fdichkspec(sysf, Stest; atol = 1.e-8, FDtol = 0.00001, minimal = true)
@test rdims[rdims .> 0] == [4, 4, 5, 4, 5, 5, 6] && orders[rdims .> 0] == [3, 4, 4, 4, 4, 5, 5] && 
leastorders[rdims .> 0] == [1, 2, 1, 2, 1, 1, 1]

@time rdims, orders, leastorders = fdichkspec(sysf, Stest; atol = 1.e-8, FDtol = 0.00001, FDGainTol = 0.00001, FDfreq=[0,1], minimal = true)
@test rdims[rdims .> 0] == [5, 5, 5, 6] && orders[rdims .> 0] == [4, 4, 5, 5] && 
leastorders[rdims .> 0] == [1, 1, 1, 1]

@time rdims, orders, leastorders = fdichkspec(sysf, S_weak; atol = 1.e-8, FDtol = 0.000001, FDGainTol = 0.00001, FDfreq=0)
@test rdims[rdims .> 0] == [6, 5, 5, 5] && orders[rdims .> 0] == [5, 4, 4, 5] 

@time rdims, orders, leastorders = fdichkspec(sysf, S_weak; atol = 1.e-8, FDtol = 0.000001, FDGainTol = 0.00001, FDfreq=0, minimal = true)
@test rdims[rdims .> 0] == [6, 5, 5, 5] && orders[rdims .> 0] == [5, 4, 4, 5] && leastorders[rdims .> 0] == [1, 1, 1, 1]


S_weak[sortperm(binmat2dec(S_weak)),:]   # reordering according to increasing decimal values of rows
S_weak[sortperm(binmat2dec(S_weak[sum(S_weak,dims=2)[:],:])),:]  # reording according to increasing numbers of zeros in rows


# Lecture 6 example ()
s = rtf('s'); # define the Laplace variable s
Gu = [1/s; 1/s] # enter Gu(s)
Gd = [0; s/(s+3)]; # enter Gd(s)
Gf = [(s+1)/(s+2); 1/(s+2)]; # enter Gd(s)
# build state space model of [Gu(s) Gd(s) Gf(s)] and set input groups
sys = dss([Gu Gd Gf],minimal=true)
sysf = fdimodset(dss([Gu Gd Gf],minimal=true),c =1,d = 2,f = 3)

@time S_weak = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5)
@test S_weak ==  dec2binmat(1)

@time S_strong = fdigenspec(sysf, atol = 1.e-7, FDtol = 1.e-5, FDGainTol = 0.00001, FDfreq=[0,1])  
@test isempty(S_strong)

@time rdims, ords, lords = fdichkspec(sysf,trues(1,1), FDfreq = 0)
@test rdims[1] == 0


s = rtf('s')
sys = dss([1/(s^2+1) 1/s]);
@time S_weak = fdigenspec(sys, atol = 1.e-7, FDtol = 1.e-5)
@test S_weak ==  Bool[1 1]

@time S_strong = fdigenspec(sys, atol = 1.e-7, FDtol = 1.e-5, FDGainTol = 0.00001, FDfreq=[0,1])  
@test isempty(S_strong)

@time rdims, ords, lords = fdichkspec(sys, S_weak, FDfreq=[1])
@test rdims[1] == 0

s = rtf('s')
sys = dss([1/s/(s^2+1) 1/s/(s^2+1)/(s+2)]);
@time S_weak = fdigenspec(sys, atol = 1.e-7, FDtol = 1.e-5)
@test S_weak ==  Bool[1 1]

@time S_strong = fdigenspec(sys, atol = 1.e-7, FDtol = 1.e-5, FDGainTol = 0.00001, FDfreq=[0,1])  
@test S_weak == S_strong

@time rdims, ords, lords = fdichkspec(sys, S_weak, FDfreq=[0,1])
@test rdims[1] == 1

## Example with solvable strong synthesis (frequencies are also zeros) 
#
s = rtf('s');
# define 2x2 Gf(s) without zeros 
Gf = [s/(s+1) 0; (s^2+4)/(s+3)^2 1]; # enter Gf(s)

@time sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1:2);
S_weak = fdigenspec(sysf)
S_strong = fdigenspec(sysf, FDfreq = [0], FDtol = 1.e-7, FDGainTol = 1.e-3)
r, ord, lord = fdichkspec(sysf, S_weak; FDfreq = [0], FDGainTol = 1.e-3)
@test S_strong == S_weak[r.> 0,:]

S_strong = fdigenspec(sysf, FDfreq = [2], FDtol = 1.e-7, FDGainTol = 1.e-3)
r, ord, lord = fdichkspec(sysf, S_weak; FDfreq = [2], FDGainTol = 1.e-3)
@test S_strong == S_weak[r.> 0,:]

S_strong = fdigenspec(sysf, FDfreq = [0,2], FDtol = 1.e-7, FDGainTol = 1.e-3)
r, ord, lord = fdichkspec(sysf, S_weak; FDfreq = [0,2], FDGainTol = 1.e-3)
@test S_strong == S_weak[r.> 0,:]


end # test fdigenspec



end # module
