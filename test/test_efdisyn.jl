module Test_efdisyn

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Random
using Test

# Testing examples for EFDSYN 
println("Test_efdisyn")
@testset "efdisyn" begin
rand()

## Example without control and disturbance inputs
p = 3; mf = 2;
sysf = fdimodset(rss(1,p,mf),faults = 1:mf);


# solve an EFDIP using the nullspace based approach
@time Q, Rf, info = efdisyn(sysf; atol3 = 1.e-7, minimal = false, rdim = [2]); info
R = fdIFeval(Q, sysf)
@test iszero(Rf.sys[1][:,Rf.faults]-R.sys[1][:,R.faults],atol=1.e-7) && info.HDesign[1] == [0.0 1.0 0.0; 0.0 0.0 1.0]

@time Q, Rf, info = efdisyn(sysf; atol3 = 1.e-7, minimal = false, rdim = [2], separate = true); info
R = fdIFeval(Q, sysf)
@test iszero(Rf.sys[1][:,Rf.faults]-R.sys[1][:,R.faults],atol=1.e-7) && info.HDesign[1] == [0.0 1.0 0.0; 0.0 0.0 1.0]


## Example without disturbance inputs
p = 3; mf = 2; mu = 2
sysf = fdimodset(rss(1,p,mf+mu),controls = 1:mu, faults = mu+1:mu+mf);

# solve an EFDP using the nullspace based approach
@time Q, Rf, info = efdisyn(sysf; FDtol = 0.00001, atol = 1.e-7, minimal = false, rdim = [2]); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys[1][:,Rf.faults]-R.sys[1][:,R.faults],atol=1.e-7) && iszero(R.sys[1][:,R.controls],atol=1.e-7) &&
      info.HDesign[1] == [0.0 1.0 0.0; 0.0 0.0 1.0]

# solve using observer based nullspace
@time Q, Rf, info = efdisyn(sysf; nullspace = false, atol3 = 1.e-7, minimal = false, rdim = [2]); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys[1][:,Rf.faults]-R.sys[1][:,R.faults],atol=1.e-7) && iszero(R.sys[1][:,R.controls],atol=1.e-7) &&
      info.HDesign[1] == [0.0 1.0 0.0; 0.0 0.0 1.0]

@time Q, Rf, info = efdisyn(sysf; nullspace = false, atol3 = 1.e-7, minimal = true, rdim = [1]); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys[1][:,Rf.faults]-R.sys[1][:,R.faults],atol=1.e-7) && iszero(R.sys[1][:,R.controls],atol=1.e-7) &&
      info.HDesign[1] == [0.0 1.0 0.0]

SFDI = fdigenspec(sysf, atol = 1.e-7)
@time Q, Rf, info = efdisyn(sysf, SFDI; nullspace = false, atol3 = 1.e-7, minimal = true, rdim = ones(Int,3)); info
R = fdIFeval(Q, sysf);
@test iszero(vcat(Rf.sys...)-vcat(R.sys...)[:,R.faults],atol=1.e-7) && iszero(vcat(R.sys...)[:,R.controls],atol=1.e-7) &&
      info.HDesign == [[0.0 1.0 0.0], [0.0 1.0], [0.0 1.0]] &&
      fditspec(Rf, FDtol = 1.e-6) == SFDI  && fdisspec(Rf, FDGainTol = 1.e-3) == SFDI

SFDI = fdigenspec(sysf, atol = 1.e-7)
@time Q, Rf, info = efdisyn(sysf, SFDI[1:1,:]; nullspace = false, atol3 = 1.e-7, minimal = true, rdim = ones(Int,1),HDesign = [[0.03964565583901385 0.8670124660000067 -1.7144387091257782]]); info
#@time Q, Rf, info = efdisyn(sysf, SFDI[1:1,:]; nullspace = false, atol3 = 1.e-7, minimal = true, rdim = ones(Int,1),HDesign = [randn(1,3)]); info
fdscond(Rf,SFDI[1:1,:])


# Example 5.10c - Solution of an EFDIP

# enter output and fault vector dimensions
p = 3; mf = 4;
# generate random dimensions for system order and input vectors
nu = Int(floor(1+4*rand())); mu = Int(floor(1+4*rand()));
nd = Int(floor(1+4*rand())); md = Int(floor(1+4*rand()));

# define random Gu(s) and Gd(s) with mf-times sensor redundancy
# and Gf(s) for triplex sensor faults
Gu = ones(mf,1)*rss(nu,1,mu); # enter Gu(s) in state-space form
Gd = ones(mf,1)*rss(nd,1,md); # enter Gd(s) in state-space form

# build model with faults
sysf = fdimodset([Gu Gd], c = 1:mu, d = mu.+(1:md), fs = 1:mf);     

SFDI = fdigenspec(sysf) # determine achievable structure matrix
nb = size(SFDI,1)
@time Q, Rf = efdisyn(sysf, SFDI; atol = 1.e-7, rdim = 1);

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
R = fdIFeval(Q,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      order.(Q.sys) == zeros(Int,nb) 

# check weak and strong fault detectability
@test fditspec(Rf) == fdisspec(Rf) == SFDI 


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

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
R = fdIFeval(Q,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      order.(Q.sys) == ones(Int,nb)

# check weak and strong fault detectability
@test fditspec(Rf) == fdisspec(Rf) == SFDI 


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

SFDI = fdigenspec(sysf) # determine achievable structure matrix
nb = size(SFDI,1)
@time Q, Rf = efdisyn(sysf, SFDI; atol = 1.e-7, rdim = 1);

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
R = fdIFeval(Q,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      order.(Q.sys) == [1, 1, 2]

# check weak and strong fault detectability
@test fditspec(Rf) == fdisspec(Rf) == SFDI 


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
@time Q, Rf = efdisyn(sysf, SFDI; atol = 1.e-7, rdim = 1);

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
R = fdIFeval(Q,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      order.(Q.sys) == [1, 1, 1]

# check weak and strong fault detectability
@test fditspec(Rf) == fdisspec(Rf) == SFDI 




## Yuan et al. IJC (1997)  strong synthesis
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

SFDI = fdigenspec(sysf; FDfreq = 0, atol = 1.e-7) # determine achievable strong structure matrix
nb = size(SFDI,1)

@time Q, Rf = efdisyn(sysf, SFDI; FDfreq = 0, atol = 1.e-7, rdim = 1);

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
R = fdIFeval(Q,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      order.(Q.sys) == [2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 2, 2]

# check strong fault detectability
@test fdisspec(Rf) == SFDI 


## Yuan et al. IJC (1997)  weak synthesis with stabillization
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

SFDI = fdigenspec(sysf; atol = 1.e-7) # determine achievable weak structure matrix
nb = size(SFDI,1)

@time Q, Rf = efdisyn(sysf, SFDI; atol = 1.e-7, rdim = 1);

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
R = fdIFeval(Q,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      order.(Q.sys) == [2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2] &&
      order(vcat(Q.sys...)) == 32

# check weak  fault detectability
@test fditspec(Rf) == SFDI 



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

SFDI = fdigenspec(sysf; atol = 1.e-7) # determine achievable weak structure matrix
nb = size(SFDI,1)

@time Q, Rf = efdisyn(sysf, SFDI; smarg = -5, poles = [-5, -5], atol = 1.e-7, rdim = 1);

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
R = fdIFeval(Q,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      order.(Q.sys) == [2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2] &&
      order(vcat(Q.sys...)) == 32 && order(gbalmr(vcat(Q.sys...),atol=1.e-7)[1]) == 6 &&
      count(ghanorm(vcat(Q.sys...))[2] .> 1.e-7) == 6 &&
      fditspec(Rf) == SFDI

@test all(iszero.(gbalmr(Q; balance = true).sys - Q.sys, atol=1.e-7) )

## Example with solvable strong synthesis  
#
s = rtf('s');
# define 2x2 Gf(s) without zeros 
Gf = [s/(s+1) 0; (s^2+4)/(s+3)^2 1]; # enter Gf(s)
sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1:2);

# weak achievable specifications
SFDI = fdigenspec(sysf)   

@time Q, Rf, info = efdisyn(sysf, SFDI; sdeg = -3, rdim = 1); info

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
R = fdIFeval(Q,sysf); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      fditspec(Rf) == SFDI

# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdisyn(sysf, SFDI; sdeg = -3, rdim = 1, HDesign = info.HDesign); info1 

@test all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) && all(iszero.(Rf.sys .- Rf1.sys,atol=1.e-7))



# build minimal realization of the model with faults
sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1:2);

# strong achievable specifications
SFDI = fdigenspec(sysf, FDfreq = [0, 2], FDGainTol = 1.e-3)

@time Q, Rf, info = efdisyn(sysf, SFDI; sdeg = -3, rdim = 1, FDfreq = [0, 2]); info

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      fdisspec(Rf,[0,2]) == SFDI


# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdisyn(sysf, SFDI; sdeg = -3, rdim = 1, FDfreq = [0, 2], FDGainTol = 0.001, HDesign = info.HDesign); info1 

@test all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) && all(iszero.(Rf.sys .- Rf1.sys,atol=1.e-7))

## Example with solvable strong synthesis (frequencies are also zeros) 
#  Using simple basis based synthesis
#
s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [s/(s+1) 0; (s^2+4)/(s+3)^2 1]; # enter Gf(s)

# build minimal realzation of the model with faults
sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1:2);

# strong achievable specifications
SFDI = fdigenspec(sysf, FDfreq = [0, 2], FDGainTol = 1.e-3)


@time Q, Rf, info = efdisyn(sysf, SFDI; sdeg = -3, rdim = 1, FDfreq = [0, 2], simple = true); info

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      fdisspec(Rf,[0,2]) == SFDI

Q1, Rf1, info1 = efdisyn(sysf, SFDI; sdeg = -3, rdim = 1, FDfreq = [0, 2], simple = true, HDesign = info.HDesign); info1 

@test all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) && all(iszero.(Rf.sys .- Rf1.sys,atol=1.e-7))


## Example with solvable strong synthesis (frequencies are also poles)  
#
s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [1/s; 1/(s^2+4)]; # enter Gf(s)

# build minimal realzation of the model with faults
sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1);

# strong achievable specifications
SFDI = fdigenspec(sysf, FDfreq = [0, 2], FDGainTol = 1.e-3);
SFDI = [SFDI;SFDI];

@time Q, Rf, info = efdisyn(sysf, SFDI; sdeg = -3, rdim = [1,2], FDfreq = [0, 2], simple = true); info

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      fdisspec(Rf,[0,2]) == SFDI

Q1, Rf1, info1 = efdisyn(sysf, SFDI; sdeg = -3, rdim = [1,2], FDfreq = [0, 2], simple = true, HDesign = info.HDesign); info1 

@test all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) && all(iszero.(Rf.sys .- Rf1.sys,atol=1.e-7))




## Example with solvable strong synthesis (frequencies are also poles)  
#  Using simple basis based synthesis
s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [1/s; 1/(s^2+4)]; # enter Gf(s)

# build minimal realzation of the model with faults
sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1);

# strong achievable specifications
SFDI = fdigenspec(sysf, FDfreq = [0, 2], FDGainTol = 1.e-3);
SFDI = [SFDI;SFDI];

@time Q, Rf, info = efdisyn(sysf, SFDI; sdeg = -3, rdim = [1,2], FDfreq = [0, 2], simple = true); info

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      fdisspec(Rf,[0,2]) == SFDI

Q1, Rf1, info1 = efdisyn(sysf, SFDI; sdeg = -3, rdim = [1,2], FDfreq = [0, 2], simple = true, HDesign = info.HDesign); info1 

@test all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) && all(iszero.(Rf.sys .- Rf1.sys,atol=1.e-7))


## Example with solvable strong synthesis (frequencies are also zeros)  

s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [s; s^2+4]; # enter Gf(s)

# build minimal realzation of the model with faults
sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1);

# strong achievable specifications
SFDI = fdigenspec(sysf, FDfreq = [0, 2], FDGainTol = 1.e-3);
SFDI = [SFDI;SFDI];

separate = true
@time Q, Rf, info = efdisyn(sysf, SFDI; atol = 1.e-7, separate, smarg = -3, sdeg = -3, rdim = [1,2], FDGainTol=1.e-5, FDfreq = [0, 2], simple = true); info

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      fdisspec(Rf,[0,2]; FDGainTol=1.e-5) == SFDI

Q1, Rf1, info1 = efdisyn(sysf, SFDI; atol = 1.e-7, separate, smarg = -3, sdeg = -3, rdim = [1,2], FDGainTol=1.e-5, FDfreq = [0, 2], simple = true, HDesign = info.HDesign); info1 
#Q1, Rf1, info1 = efdisyn(sysf, SFDI; atol = 1.e-7, smarg = -3, sdeg = -3, rdim = [1,2], FDGainTol=1.e-5, FDfreq = [0, 2], simple = false, HDesign = info.HDesign); info1 
R1 = fdIFeval(Q1, sysf; minimal = true, atol = 1.e-7);
@test iszero(vcat(R1.sys...)[:,[R1.controls;R1.disturbances]],atol=1.e-7) &&
      iszero(vcat(R1.sys...)[:,R1.faults]-vcat(Rf1.sys...),atol=1.e-7) &&
      fdisspec(Rf1,[0,2]; FDGainTol=1.e-5) == SFDI

@test all(iszero.(Q.sys .- Q1.sys,atol=1.e-7)) && all(iszero.(Rf.sys .- Rf1.sys,atol=1.e-7))



end # test efdisyn



end # module
