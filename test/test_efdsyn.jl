module Test_efdsyn

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Random
using Test

# Testing examples for EFDSYN 
println("Test_efdsyn")
@testset "efdsyn" begin
rand()

## Example without control and disturbance inputs
p = 3; mf = 2;
sysf = fdimodset(rss(1,p,mf),faults = 1:mf);


# solve an EFDP using the nullspace based approach
@time Q, Rf, info = efdsyn(sysf; atol3 = 1.e-7, minimal = false, rdim = 2); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys[:,Rf.faults]-R.sys[:,R.faults],atol=1.e-7) && info.S == Bool[1 1; 1 1; 1 1] && info.HDesign == [0.0 1.0 0.0; 0.0 0.0 1.0]

## Example without disturbance inputs
p = 3; mf = 2; mu = 2
sysf = fdimodset(rss(1,p,mf+mu),controls = 1:mu, faults = mu+1:mu+mf);

# solve an EFDP using the nullspace based approach
@time Q, Rf, info = efdsyn(sysf; FDtol = 0.00001, atol = 1.e-7, minimal = false, rdim = 2); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && iszero(R.sys[:,R.controls],atol=1.e-7) &&
      info.S == Bool[1 1; 1 1; 1 1] && info.HDesign == [0.0 1.0 0.0; 0.0 0.0 1.0]

# solve using observer based nullspace
@time Q, Rf, info = efdsyn(sysf; nullspace = false, atol3 = 1.e-7, minimal = false, rdim = 2); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && iszero(R.sys[:,R.controls],atol=1.e-7) &&
      info.S == Bool[1 1; 1 1; 1 1] && info.HDesign == [1.0 0.0 0.0; 0.0 1.0 0.0]

@time Q, Rf, info = efdsyn(sysf; nullspace = false, atol3 = 1.e-7, minimal = true, rdim = 1); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && iszero(R.sys[:,R.controls],atol=1.e-7) &&
      info.S == Bool[1 1; 1 1; 1 1] && info.HDesign == [1.0 0.0 0.0]


## Example with solvable strong synthesis (frequencies are also zeros) 
#
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

# build minimal realzation of the model with faults
sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1);

@time Q, Rf, info = efdsyn(sysf; sdeg = -3, rdim = 1, FDfreq = [0, 2]); info
# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(R.sys[:,[R.controls;R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys) && 
      fdisspec(Rf, [0, 2])[:] == Bool[1] 

# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdsyn(sysf; sdeg = -3, rdim = 1, FDfreq = [0, 2], FDGainTol = 0.001, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(Rf.sys-Rf1.sys,atol=1.e-7)

sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1:2);
@time Q, Rf, info = efdsyn(sysf; sdeg = -3, rdim = 2, FDfreq = [0, 2]); info
# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(R.sys[:,[R.controls;R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys) && 
      fdisspec(Rf, [0, 2]; block = false) == Bool[0 1; 0 0] &&
      fdisspec(Rf, [0]; block = false) == Bool[1 1; 0 0] &&
      fdisspec(Rf, [2]; block = false) == Bool[0 1; 1 0] &&
      fdisspec(Rf, [0, 2]; block = true) == Bool[1 1] 


## Example with solvable strong synthesis (frequencies are also zeros) 
#  Using simple basis based synthesis
#
s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [s/(s+1) 0; (s^2+4)/(s+3)^2 1]; # enter Gf(s)

# build minimal realzation of the model with faults
sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1);

@time Q, Rf, info = efdsyn(sysf; sdeg = -3, rdim = 1, FDfreq = [0, 2], simple = true); info

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(R.sys[:,[R.controls;R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys) && 
      fdisspec(Rf, [0, 2])[:] == Bool[1] 

Q1, Rf1, info1 = efdsyn(sysf; sdeg = -3, rdim = 1, FDfreq = [0, 2], FDGainTol = 0.001, simple = true, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(Rf.sys-Rf1.sys,atol=1.e-7)


## Example with solvable strong synthesis (frequencies are also poles)  
#
s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [1/s; 1/(s^2+4)]; # enter Gf(s)

# build minimal realzation of the model with faults
sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1);

@time Q, Rf, info = efdsyn(sysf; sdeg = -3, rdim = 1, FDfreq = [0, 2]); info
# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(R.sys[:,[R.controls;R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys) && 
      fdisspec(Rf, [0, 2])[:] == Bool[1] 

Q1, Rf1, info1 = efdsyn(sysf; sdeg = -3, rdim = 1, FDfreq = [0, 2], FDGainTol = 0.001, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(Rf.sys-Rf1.sys,atol=1.e-7)


## Example with solvable strong synthesis (frequencies are also poles)  
#  Using simple basis based synthesis
s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [1/s; 1/(s^2+4)]; # enter Gf(s)

@time Q, Rf, info = efdsyn(sysf; sdeg = -3, rdim = 1, FDfreq = [0, 2], simple = true); info
# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(R.sys[:,[R.controls;R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys) && 
      fdisspec(Rf, [0, 2])[:] == Bool[1] 

Q1, Rf1, info1 = efdsyn(sysf; sdeg = -3, rdim = 1, FDfreq = [0, 2], FDGainTol = 0.001, simple = true, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(Rf.sys-Rf1.sys,atol=1.e-7)


## Example with solvable strong synthesis (frequencies are also zeros)  
#
s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [s; s^2+4]; # enter Gf(s)

# build minimal realzation of the model with faults
sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1);

@time Q, Rf, info = efdsyn(sysf; sdeg = -3, smarg = -3, rdim = 1, FDfreq = [0, 2]); info

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(R.sys[:,[R.controls;R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys) && 
      fdisspec(Rf, [0, 2])[:] == Bool[1] 

Q1, Rf1, info1 = efdsyn(sysf; sdeg = -3, smarg = -3, rdim = 1, FDfreq = [0, 2], HDesign = info.HDesign); info1 
@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(Rf.sys-Rf1.sys,atol=1.e-7)

@time Q, Rf, info = efdsyn(sysf; sdeg = -3, smarg = -3, rdim = 1, FDfreq = [0, 2], simple = true); info

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(R.sys[:,[R.controls;R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys) && 
      fdisspec(Rf, [0, 2])[:] == Bool[1] 

Q1, Rf1, info1 = efdsyn(sysf; sdeg = -3, smarg = -3, rdim = 1, FDfreq = [0, 2], HDesign = info.HDesign); info1 
@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(Rf.sys-Rf1.sys,atol=1.e-7)


## Example with solvable strong synthesis (frequencies are also zeros)  
#  Using simple basis based synthesis
#
s = rtf('s');
# define 2x1 Gf(s) without zeros 
Gf = [s; s^2+4]; # enter Gf(s)

# build minimal realzation of the model with faults
sysf = fdimodset(dss(Gf; minimal = true, atol = 1.e-7),faults = 1);

@time Q, Rf, info = efdsyn(sysf; sdeg = -3, smarg = -3, rdim = 1, FDfreq = [0, 2], simple = true); info

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(R.sys[:,[R.controls;R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys) && 
      fdisspec(Rf, [0, 2])[:] == Bool[1] 

Q1, Rf1, info1 = efdsyn(sysf; sdeg = -3, smarg = -3, rdim = 1, FDfreq = [0, 2], simple = true, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys,atol=1.e-7) && iszero(Rf.sys-Rf1.sys,atol=1.e-7)

# Example 5.2 - EFDP has no solution 

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gd(s), Gf(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gd = [1/(s+2); 0];                   # enter Gd(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; md = 1; mf = 2;       # set dimensions

# build model with additive faults 
#sysf = fdimodset(ss([Gu Gd Gf]),struct('c',1:mu,'d',mu+(1:md),'f',mu+md+(1:mf)));
sysf = fdimodset(dss([Gu Gd Gf]),c = 1:mu, d = mu.+(1:md), f = (mu+md).+(1:mf))

# call of EFDSYN with default options 
try
  Q, Rf, info = efdsyn(sysf)
  @test false
catch err
  display(err)
  # compute achievable structure matrix
  S_achievable = fdigenspec(sysf)  
  @test true
end

# Example 5.3 - Solution of an EFDP using EFDSYN

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
@time Q, Rf, info = efdsyn(sysf, smarg =-3, sdeg = -3, rdim = 1); 

# normalize Q and Rf to match example
scale = evalfr(Rf.sys[1,1],Inf)[1,1];
dss2rm(Q.sys/scale, atol = 1.e-7)
dss2rm(Rf.sys/scale, atol = 1.e-7) 

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(R.sys[:,[R.controls;R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys) && 
      fdisspec(Rf, [0]) == Bool[1 1] 


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
@time Q, Rf, info = efdsyn(sysf, sdeg = -3, rdim = 1); 

# normalize Q and Rf to match example
scale = evalfr(Rf.sys[1,1],Inf)[1,1];
dss2rm(Q.sys/scale, atol = 1.e-7)
dss2rm(Rf.sys/scale, atol = 1.e-7) 

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(R.sys[:,[R.controls;R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys) && 
      fdisspec(Rf, [0]) == Bool[1 1] 

@test iszero(gbalmr(Q; balance = true).sys-Q.sys, atol=1.e-7) 


## Test covering design matrices
##
p = 3; mf = 2;
sysf = fdimodset(rss(1,p,mf),faults = 1:mf);


# solve an EFDP using the nullspace based approach
@time Q, Rf, info = efdsyn(sysf; atol3 = 1.e-7, minimal = false, rdim = 2); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && info.S == Bool[1 1; 1 1; 1 1] && info.HDesign == [0.0 1.0 0.0; 0.0 0.0 1.0]

# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdsyn(sysf; atol3 = 1.e-7, minimal = false, simple = true, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys, atol = 1.e-7) && iszero(Rf.sys-Rf1.sys, atol = 1.e-7)


## solve an EFP using the nullspace based approach
@time Q, Rf, info = efdsyn(sysf; atol3 = 1.e-7, minimal = false, rdim = 2, simple = true); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && info.S == Bool[1 1; 1 1; 1 1] && info.HDesign == [0.0 1.0 0.0; 0.0 0.0 1.0]

# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdsyn(sysf; atol3 = 1.e-7, minimal = false, simple = true, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys, atol = 1.e-7) && iszero(Rf.sys-Rf1.sys, atol = 1.e-7)


## solve an EFP using the nullspace based approach
@time Q, Rf, info = efdsyn(sysf; atol3 = 1.e-7, minimal = false, rdim = 1, simple = true); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && info.S == Bool[1 1; 1 1; 1 1] && info.HDesign == [0.0  1.0  0.0]

# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdsyn(sysf; atol3 = 1.e-7, minimal = false, simple = true, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys, atol = 1.e-7) && iszero(Rf.sys-Rf1.sys, atol = 1.e-7)

## solve an EFD using the nullspace based approach
@time Q, Rf, info = efdsyn(sysf; atol3 = 1.e-7, minimal = true, simple = true); info
R = fdIFeval(Q, sysf; atol = 1.e-7);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && iszero(R.sys[:,R.controls],atol=1.e-7) && 
      info.S == Bool[1 1; 1 1; 1 1] && info.HDesign ==[0.0 1.0 0.0]

# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdsyn(sysf; atol3 = 1.e-7, minimal = false, simple = true, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys, atol = 1.e-7) && iszero(Rf.sys-Rf1.sys, atol = 1.e-7)


## solve an EFD using the nullspace based approach
@time Q, Rf, info = efdsyn(sysf; atol3 = 1.e-7, minimal = true, simple = false); info
R = fdIFeval(Q, sysf; atol = 1.e-7);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && iszero(R.sys[:,R.controls],atol=1.e-7) && 
      info.S == Bool[1 1; 1 1; 1 1] && info.HDesign ==[0.0 1.0 0.0]

# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdsyn(sysf; atol3 = 1.e-7, minimal = false, simple = false, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys, atol = 1.e-7) && iszero(Rf.sys-Rf1.sys, atol = 1.e-7)


##
## Example with fault and noise inputs
p = 3; mf = 2; mw = 2; mu = 0
sysf = fdimodset(rss(1,p,mf+mw), noise = (mu+mf).+ (1:mw), faults = mu+1:mu+mf);


## solve an AMMP using the nullspace based approach
Q, Rfw, info = efdsyn(sysf; atol = 1.e-7, atol3 = 1.e-7); info
R = fdIFeval(Q, sysf; atol = 1.e-7);
@test iszero(Rfw.sys[:,Rfw.faults]-R.sys[:,R.faults],atol=1.e-7) && iszero(R.sys[:,R.controls],atol=1.e-7) && 
      iszero(Rfw.sys[:,Rfw.noise]-R.sys[:,R.noise],atol=1.e-7) &&
      info.S == Bool[1 1; 1 1; 1 1] && info.HDesign ==[0.0 1.0 0.0]

# the same results  must result using the resulting design matrix H
Q1, Rfw1, info1 = efdsyn(sysf; atol3 = 1.e-7, minimal = false, simple = false, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys, atol = 1.e-7) && iszero(Rfw.sys-Rfw1.sys, atol = 1.e-7)


## special cases which require forming rdim combinations of nout vecotors, with nout > rdim  

# Gf = [ I;0]
p = 4; mf = 2;
sys = dss([eye(2);zeros(2,2)]);
sysf = fdimodset(sys,faults = 1:mf);


# solve an EFDP using the nullspace based approach
@time Q, Rf, info = efdsyn(sysf; atol3 = 1.e-7, minimal = true, rdim = 1); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && info.S == Bool[0 0; 0 0; 0 1; 1 0] 

# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdsyn(sysf; atol3 = 1.e-7, minimal = true, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys, atol = 1.e-7) && iszero(Rf.sys-Rf1.sys, atol = 1.e-7)



# solve an EFDP using the nullspace based approach
@time Q, Rf, info = efdsyn(sysf; atol3 = 1.e-7, minimal = false, rdim = 1); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && info.S == Bool[0 0; 0 0; 0 1; 1 0] 

# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdsyn(sysf; atol3 = 1.e-7, minimal = false, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys, atol = 1.e-7) && iszero(Rf.sys-Rf1.sys, atol = 1.e-7)

# solve an EFDP using the nullspace based approach
@time Q, Rf, info = efdsyn(sysf; atol3 = 1.e-7, minimal = true, rdim = 1, simple = true); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && info.S == Bool[0 0; 0 0; 0 1; 1 0] 

# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdsyn(sysf; atol3 = 1.e-7, minimal = true, simple = true, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys, atol = 1.e-7) && iszero(Rf.sys-Rf1.sys, atol = 1.e-7)

# solve an EFDP using the nullspace based approach
@time Q, Rf, info = efdsyn(sysf; atol3 = 1.e-7, minimal = false, rdim = 1, simple = true); info
R = fdIFeval(Q, sysf);
@test iszero(Rf.sys-R.sys[:,R.faults],atol=1.e-7) && info.S == Bool[0 0; 0 0; 0 1; 1 0] 

# the same results  must result using the resulting design matrix H
Q1, Rf1, info1 = efdsyn(sysf; atol3 = 1.e-7, minimal = false, simple = true, HDesign = info.HDesign); info1 

@test iszero(Q.sys-Q1.sys, atol = 1.e-7) && iszero(Rf.sys-Rf1.sys, atol = 1.e-7)


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

# opt = struct('tol',1.e-7,'FDGainTol',.01,...
#     'FDFreq',0,'sdeg',-0.05,'m1',mu+md);
# apply fdigenspec to [Gu Gd Gf; I 0 0]
#S_strong = fdigenspec([sys; eye(mu,mu+md+mf)],opt)
S_strong = Bool[1 0;0 1;1 1]
SFDI = S_strong

nb = size(SFDI,1)
orders = zeros(Int,nb)

Qtot = dss(zeros(0,p+mu)); Rftilde = dss(zeros(0,mf));
info = Any[]; HDesign = Any[]
for i = 1:nb
   # set input groups for i-th design
   indd = Vector(1:mf)[SFDI[i,:] .== false] 
   indf = Vector(1:mf)[SFDI[i,:] .== true] 
   sysc = fdimodset(sysf.sys, c = 1:mu, d = [Vector(mu .+ (1:md)); (mu+md) .+ indd], 
                    f = (mu+md) .+ indf, aux = (mu+md) .+ (1:mf));
   #
   Qi, Rfi, infoi = efdsyn(sysc, atol = 1.e-7, sdeg =-5, smarg = -0.05, minimal = true,
                           FDfreq = 0, FDGainTol = 0.0001, rdim = 1, simple = false);  infoi
   info = [info; infoi]       
   HDesign = [HDesign; [infoi.HDesign]]                 
   Qtot = [Qtot; Qi.sys]; 
   Rftilde = [Rftilde; Rfi.sys[:,Rfi.aux]];
   orders[i] = order(Qi.sys)
end

# checks: Q[Gu Gd;I 0] = 0 and Q[Gf; 0] = Rftilde
syse = [sysf.sys; eye(mu,mu+md+mf)];  
@test iszero(Qtot*syse[:,[sysf.controls; sysf.disturbances]], atol = 1.e-7)
@test iszero(Qtot*syse[:,sysf.faults]-Rftilde,atol=1.e-7)

@test all(orders .== 1)
#step(Rftilde,10)

S, gains = fdisspec_(Rftilde, 0;  FDGainTol = 0.0001)
@test isequal(SFDI,S[:,:,1])
S1 = fditspec_(Rftilde; atol = 1.e-7,FDStol = 0.0001, FDfreq = 0)
@test isequal(SFDI,S1[:,:,1])

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

# opt = struct('tol',1.e-7,'FDGainTol',.01,...
#     'FDFreq',0,'sdeg',-0.05,'m1',mu+md);
# apply fdigenspec to [Gu Gd Gf; I 0 0]
#S_strong = fdigenspec([sys; eye(mu,mu+md+mf)],opt)
S_strong = Bool[1 0;0 1;1 1]
SFDI = S_strong

nb = size(SFDI,1)
orders = zeros(Int,nb)

Qtot = dss(zeros(0,p+mu)); Rftilde = dss(zeros(0,mf));
info = Any[]; HDesign = Any[]
for i = 1:nb
   # set input groups for i-th design
   indd = Vector(1:mf)[SFDI[i,:] .== false] 
   indf = Vector(1:mf)[SFDI[i,:] .== true] 
   sysc = fdimodset(sysf.sys, c = 1:mu, d = [Vector(mu .+ (1:md)); (mu+md) .+ indd], 
                    f = (mu+md) .+ indf, aux = (mu+md) .+ (1:mf));
   #
   Qi, Rfi, infoi = efdsyn(sysc, atol = 1.e-7, sdeg =-5, smarg = -0.05, minimal = true,
                           FDfreq = 0, FDGainTol = 0.0001, rdim = 1, simple = false);  infoi
   info = [info; infoi]       
   HDesign = [HDesign; [infoi.HDesign]]                 
   Qtot = [Qtot; Qi.sys]; 
   Rftilde = [Rftilde; Rfi.sys[:,Rfi.aux]];
   orders[i] = order(Qi.sys)
end

# checks: Q[Gu Gd;I 0] = 0 and Q[Gf; 0] = Rftilde
syse = [sysf.sys; eye(mu,mu+md+mf)];  
@test iszero(Qtot*syse[:,[sysf.controls; sysf.disturbances]], atol = 1.e-7)
@test iszero(Qtot*syse[:,sysf.faults]-Rftilde,atol=1.e-7)

@test orders == [2, 1, 1]
#step(Rftilde,10)

S, gains = fdisspec_(Rftilde, 0;  FDGainTol = 0.0001)
@test isequal(SFDI,S[:,:,1])
S1 = fditspec_(Rftilde; atol = 1.e-7,FDStol = 0.0001, FDfreq = 0)
@test isequal(SFDI,S1[:,:,1])


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

# opt = struct('tol',1.e-7,'FDGainTol',.01,...
#     'FDFreq',0,'sdeg',-0.05,'m1',mu+md);
# apply fdigenspec to [Gu Gd Gf; I 0 0]
#S_strong = fdigenspec([sys; eye(mu,mu+md+mf)],opt)
S_strong = Bool[1 0;0 1;1 1]
SFDI = S_strong

nb = size(SFDI,1)
orders = zeros(Int,nb)

Qtot = dss(zeros(0,p+mu)); Rftilde = dss(zeros(0,mf));
info = Any[]; HDesign = Any[]
for i = 1:nb
   # set input groups for i-th design
   indd = Vector(1:mf)[SFDI[i,:] .== false] 
   indf = Vector(1:mf)[SFDI[i,:] .== true] 
   sysc = fdimodset(sysf.sys, c = 1:mu, d = [Vector(mu .+ (1:md)); (mu+md) .+ indd], 
                    f = (mu+md) .+ indf, aux = (mu+md) .+ (1:mf));
   #
   Qi, Rfi, infoi = efdsyn(sysc, atol = 1.e-7, sdeg =-5, smarg = -0.05, minimal = true,
                           FDfreq = 0, FDGainTol = 0.0001, rdim = 1, simple = false);  infoi
   info = [info; infoi]       
   HDesign = [HDesign; [infoi.HDesign]]                 
   Qtot = [Qtot; Qi.sys]; 
   Rftilde = [Rftilde; Rfi.sys[:,Rfi.aux]];
   orders[i] = order(Qi.sys)
end

# checks: Q[Gu Gd;I 0] = 0 and Q[Gf; 0] = Rftilde
syse = [sysf.sys; eye(mu,mu+md+mf)];  
@test iszero(Qtot*syse[:,[sysf.controls; sysf.disturbances]], atol = 1.e-7)
@test iszero(Qtot*syse[:,sysf.faults]-Rftilde,atol=1.e-7)

@test orders == [1, 1, 1]
#step(Rftilde,10)

S, gains = fdisspec_(Rftilde, 0;  FDGainTol = 0.0001)
@test isequal(SFDI,S[:,:,1])
S1 = fditspec_(Rftilde; atol = 1.e-7,FDStol = 0.0001, FDfreq = 0)
@test isequal(SFDI,S1[:,:,1])


# use individual synthesis with specified structure vectors 
nb = size(SFDI,1)
orders = zeros(Int,nb)

Qtot = dss(zeros(0,p+mu)); Rftilde = dss(zeros(0,mf));
info = Any[]; HDesign = Any[]
for i = 1:nb
   Qi, Rfi, infoi = efdsyn(sysf, SFDI[i,:]; atol = 1.e-7, sdeg =-5, smarg = -0.05, minimal = true,
                           FDfreq = 0, FDGainTol = 0.0001, rdim = 1, simple = false);  infoi
   info = [info; infoi]       
   HDesign = [HDesign; [infoi.HDesign]]                 
   Qtot = [Qtot; Qi.sys]; 
   Rftilde = [Rftilde; Rfi.sys];
   orders[i] = order(Qi.sys)
end

# checks: Q[Gu Gd;I 0] = 0 and Q[Gf; 0] = Rftilde
syse = [sysf.sys; eye(mu,mu+md+mf)];  
@test iszero(Qtot*syse[:,[sysf.controls; sysf.disturbances]], atol = 1.e-7)
@test iszero(Qtot*syse[:,sysf.faults]-Rftilde,atol=1.e-7)

@test orders == [1, 1, 1]
#step(Rftilde,10)

S, gains = fdisspec_(Rftilde, 0;  FDGainTol = 0.0001)
@test isequal(SFDI,S[:,:,1])
S1 = fditspec_(Rftilde; atol = 1.e-7,FDStol = 0.0001, FDfreq = 0)
@test isequal(SFDI,S1[:,:,1])


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

# opt = struct('tol',1.e-7,'FDGainTol',.01,...
#     'FDFreq',0,'sdeg',-0.05,'m1',mu+md);
# # apply fdigenspec to [Gu Gd Gf; I 0 0]
# S_strong = fdigenspec([sys; eye(mu,mu+md+mf)],opt)
# 

S_strong = Bool[   0   0   0   1   0   0   1   1
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
1   1   1   1   1   1   1   1]


SFDI = S_strong

nb = size(SFDI,1)
orders = zeros(Int,nb)

Qtot = dss(zeros(0,p+mu)); Rftilde = dss(zeros(0,mf));
info = Any[]; HDesign = Any[]
for i = 1:nb
   # set input groups for i-th design
   indd = Vector(1:mf)[SFDI[i,:] .== false] 
   indf = Vector(1:mf)[SFDI[i,:] .== true] 
   sysc = fdimodset(sysf.sys, c = 1:mu, d = [Vector(mu .+ (1:md)); (mu+md) .+ indd], 
                    f = (mu+md) .+ indf, aux = (mu+md) .+ (1:mf));
   #
   Qi, Rfi, infoi = efdsyn(sysc, atol = 1.e-7, sdeg =-5, smarg = -0.05, minimal = true,
                           FDfreq = 0, FDGainTol = 0.0001, rdim = 1, simple = false);  infoi
   info = [info; infoi]       
   HDesign = [HDesign; [infoi.HDesign]]                 
   Qtot = [Qtot; Qi.sys]; 
     Rftilde = [Rftilde; Rfi.sys[:,Rfi.aux]];
   orders[i] = order(Qi.sys)
end

# checks: Q[Gu Gd;I 0] = 0 and Q[Gf; 0] = Rftilde
syse = [sysf.sys; eye(mu,mu+md+mf)];  
@test iszero(Qtot*syse[:,[sysf.controls; sysf.disturbances]], atol = 1.e-7)
@test iszero(Qtot*syse[:,sysf.faults]-Rftilde,atol=1.e-7)

@test orders == [1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2]
#step(Rftilde,10)

S, gains = fdisspec_(Rftilde, 0;  FDGainTol = 0.0001)
@test isequal(SFDI,S[:,:,1])
S1 = fditspec_(Rftilde; atol = 1.e-7,FDStol = 0.0001, FDfreq = 0)
@test isequal(SFDI,S1[:,:,1])


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

# opt = struct('tol',1.e-7,'FDGainTol',.01,...
#     'FDFreq',0,'sdeg',-0.05,'m1',mu+md);
# # apply fdigenspec to [Gu Gd Gf; I 0 0]
# S_strong = fdigenspec([sys; eye(mu,mu+md+mf)],opt)
# 

S_weak = Bool[   0   0   0   1   0   0   1   1
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
1   1   1   1   1   1   1   1]

SFDI = S_weak;
nb = size(SFDI,1)
orders = zeros(Int,nb)

Qtot = dss(zeros(0,p+mu)); Rftilde = dss(zeros(0,mf));
info = Any[]; HDesign = Any[]
for i = 1:nb
   # set input groups for i-th design
   indd = Vector(1:mf)[SFDI[i,:] .== false] 
   indf = Vector(1:mf)[SFDI[i,:] .== true] 
   sysc = fdimodset(sysf.sys, c = 1:mu, d = [Vector(mu .+ (1:md)); (mu+md) .+ indd], 
                    f = (mu+md) .+ indf, aux = (mu+md) .+ (1:mf));
   #
   Qi, Rfi, infoi = efdsyn(sysc, atol = 1.e-7, sdeg =-5, smarg = -0.05, minimal = true,
                           FDtol = 0.0001, rdim = 1, simple = false);  infoi
   info = [info; infoi]       
   HDesign = [HDesign; [infoi.HDesign]]                 
   Qtot = [Qtot; Qi.sys]; 
   Rftilde = [Rftilde; Rfi.sys[:,Rfi.aux]];
   orders[i] = order(Qi.sys)
end

# checks: Q[Gu Gd;I 0] = 0 and Q[Gf; 0] = Rftilde
syse = [sysf.sys; eye(mu,mu+md+mf)];  
@test iszero(Qtot*syse[:,[sysf.controls; sysf.disturbances]], atol = 1.e-7)
@test iszero(Qtot*syse[:,sysf.faults]-Rftilde,atol=1.e-7)

@test orders == [1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
#step(Rftilde,10)

S = fditspec_(Rftilde;  atol = 1.e-5, FDtol = 0.0001)
@test isequal(SFDI,S)  

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

# opt = struct('tol',1.e-7,'FDGainTol',.01,...
#     'FDFreq',0,'sdeg',-0.05,'m1',mu+md);
# # apply fdigenspec to [Gu Gd Gf; I 0 0]
# S_strong = fdigenspec([sys; eye(mu,mu+md+mf)],opt)
# 

S_weak = Bool[   0   0   0   1   0   0   1   1
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
1   1   1   1   1   1   1   1]

SFDI = S_weak;
nb = size(SFDI,1)
orders = zeros(Int,nb)

Qtot = dss(zeros(0,p+mu)); Rftilde = dss(zeros(0,mf));
info = Any[]; HDesign = Any[]
for i = 1:nb
   # set input groups for i-th design
   indd = Vector(1:mf)[SFDI[i,:] .== false] 
   indf = Vector(1:mf)[SFDI[i,:] .== true] 
   sysc = fdimodset(sysf.sys, c = 1:mu, d = [Vector(mu .+ (1:md)); (mu+md) .+ indd], 
                    f = (mu+md) .+ indf, aux = (mu+md) .+ (1:mf));
   #
   Qi, Rfi, infoi = efdsyn(sysc, atol = 1.e-7, sdeg =-5, smarg = -5, poles = [-5, -6], minimal = true,
                            FDtol = 0.0001, rdim = 1, simple = false);  infoi
   info = [info; infoi]       
   HDesign = [HDesign; [infoi.HDesign]]                 
   Qtot = [Qtot; Qi.sys]; 
   Rftilde = [Rftilde; Rfi.sys[:,Rfi.aux]];
   orders[i] = order(Qi.sys)
end

# checks: Q[Gu Gd;I 0] = 0 and Q[Gf; 0] = Rftilde
syse = [sysf.sys; eye(mu,mu+md+mf)];  
@test iszero(Qtot*syse[:,[sysf.controls; sysf.disturbances]], atol = 1.e-7)
@test iszero(Qtot*syse[:,sysf.faults]-Rftilde,atol=1.e-7)

@test orders == [1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
#step(Rftilde,10)


S = fditspec_(Rftilde;  atol = 1.e-7, FDtol = 0.0001);
@test isequal(SFDI,S)  

@test  count(ghanorm(Qtot)[2] .> 1.e-7) == 6


# ## DC Motor (Ding, 2013)
# s = rtf('s');
# Gu = 8.75/(1+1.225*s)/(1+0.03*s)/(1+0.005*s); 
# Gd = -31.07/s/(1+0.005*s);
# Gc = 1.6*(1 + 1.225*s)/s;

# G = [Gu*Gc Gu 1 Gd ]./(1+Gu*Gc);
# # set input groups
# sysf = fdimodset(dss(G,minimal=true,atol=1.e-7),c = 1, f = 2:4);
# S_weak = fdigenspec(sysf,struct('tol',1.e-7))
# S_strong = fdigenspec(sysf,struct('tol',1.e-7,'FDFreq',0))


# # # check with expressions from Ding (2013)
# # gu = 14/(s*(1+0.03*s)*(1+0.005*s)+14);
# # gd = -31.07*(1+0.03*s)/(s*(1+0.03*s)*(1+0.005*s)+14);
# # gfa = 8.75*s/(1+1.225*s)/(s*(1+0.03*s)*(1+0.005*s)+14);
# # gfs = s*(1+0.03*s)*(1+0.005*s)/(s*(1+0.03*s)*(1+0.005*s)+14);
# # norm(gu-G(:,1),inf), norm(gd-G(:,4),inf), norm(gfa-G(:,2),inf), norm(gfs-G(:,3),inf),
# # 

# ## case 0: 
# f = [0 0 1]; T = [0:0.1:5]'; N = length(T); ucontrol = ones(N,1); 
# Uf = [ucontrol, f(1)*ones(N,1) f(2)*ones(N,1) f(3)*ones(N,1)];
# lsim(G,Uf,T); 




# ## case 1 (Echeveria, DC Motor Benchmark): 
# f = [0.87 -0.12 0.53]; T = [0:0.1:5]'; N = length(T); ucontrol = rand(N,1); 
# Uf = [ucontrol, f(1)*ones(N,1) f(2)*ones(N,1) f(3)*ones(N,1)];
# yf = lsim(G,Uf,T); 
# yfp = yf.*(ones(N,1)+0.02*rand(N,1));
# yfp2 = yf.*(ones(N,1)+0.08*rand(N,1));

# x = rand(3,1); 
# err = dcmotor( x, G, yf, T, ucontrol);
# lb = -ones(3,1); ub = ones(3,1); 


# [x,err] = fmincon(@(x) dcmotor(x,G,yf,T,ucontrol),zeros(3,1),[],[],[],[],lb,ub)

# [xp,err1] = fmincon(@(x) dcmotor(x,G,yfp,T,ucontrol),zeros(3,1),[],[],[],[],lb,ub)
# [xp2,err2] = fmincon(@(x) dcmotor(x,G,yfp2,T,ucontrol),zeros(3,1),[],[],[],[],lb,ub)

end # test efdsyn



end # module
