module Test_evaluation

using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Random
using Test
# using control inputss
using Measurements

# Testing evaluation examples 
println("Test_evaluation")
@testset "FD & FDI" begin

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
scale = evalfr(Rf.sys[1,Rf.mu+Rf.md+1],Inf)[1,1];
dss2rm(Q.sys/scale, atol = 1.e-7)
dss2rm(Rf.sys/scale, atol = 1.e-7) 

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7);
@test iszero(R.sys[:,[R.controls;R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,R.faults]-Rf.sys[:,Rf.faults]) && 
      fdisspec(Rf, [0]) == Bool[1 1] 

Ts = 0.01;  
Qd = c2d(Q,Ts); 


t = 0:Ts:1; ns = length(t); tfault = 0.05; kfault = round(Int,tfault/Ts); fnorm_min = 0.02; 
commands = sin.(2pi*t).*ones(ns,1); disturbances = randn(ns); faults = 0.5*[zeros(kfault,2); ones(ns-kfault,2)]

# no faults
inputs = [commands disturbances 0*faults] # 
y, = timeresp(sysf.sys, inputs, t); 

# Fault detection
# theta = zeros(ns)
# reval = zeros(ns)
fdsys = FDSystem(Qd, α = [0], β = [10], γ = [.9])
isig = zeros(Int,ns)
for i = 1:ns
    tstep!(fdsys, y[i,:], commands[i,1:1])
    # theta[i] = fdsys.θ[1]
    isig[i] = fdsys.isig[1]
    # reval[i] = fdsys.re[1]
end
@test isempty(fdsys.indfault) && iszero(isig) # no false alarms

# with faults
inputs = [commands disturbances faults] # 
y, = timeresp(sysf.sys, inputs, t); 

# Fault detection
# theta = zeros(ns)
# reval = zeros(ns)
fdsys = FDSystem(Qd, α = [0], β = [10], γ = [.9])
isig = zeros(Int,ns)
for i = 1:ns
    tstep!(fdsys, y[i,:], commands[i,1:1])
    # theta[i] = fdsys.θ[1]
    isig[i] = fdsys.isig[1]
    # reval[i] = fdsys.re[1]
end
@test fdsys.indfault == [1]
# plot(t,y) # y1, y2
# plot!(t,reval) # y3
# plot!(t,theta) # y4
# plot!(t,isig) # y5
# plot!(t,commands_withd_jamming)


# Fault detection with measurement noise
t = 0:Ts:1; ns = length(t); tfault = 0.05; kfault = round(Int,tfault/Ts); fnorm_min = 0.02; 
commands = sin.(2pi*t).*ones(ns,1); disturbances = randn(ns); faults = 0.5*[zeros(kfault,2); ones(ns-kfault,2)]


# Fault detection with measurement noise

# no faults
inputs = [commands disturbances 0*faults] 
y, = timeresp(sysf.sys, inputs, t); 
fdsys = FDSystem(Qd, x0 = [measurement(0.,0.)], α = [0], β = [10], γ = [.9]) 
isig = zeros(Int,ns)
for i = 1:ns
    tstep!(fdsys, y[i,:]*measurement(1,0.2), commands[i,1:1])
    isig[i] = fdsys.isig[1]
end
@test isempty(fdsys.indfault) && iszero(isig)


# with faults
inputs = [commands disturbances faults] 
y, = timeresp(sysf.sys, inputs, t); 
fdsys = FDSystem(Qd, x0 = [measurement(0.,0.)], α = [0], β = [10], γ = [.9]) 
isig = zeros(Int,ns)
for i = 1:ns
    tstep!(fdsys, y[i,:]*measurement(1,0.2), commands[i,1:1])
    isig[i] = fdsys.isig[1]
end
@test fdsys.indfault == [1]



# Fault detection with parametric uncertainties

s = rtf('s');
# define Gup(s) and Gdp(s) with parametric uncertainties
Gup = [(s+measurement(1,0.01))/(s+measurement(2,0.02)); (s+2)/(s+measurement(3,0.03))];     # enter Gup(s)
Gdp = [(s-measurement(1,0.01))/(s+measurement(2,0.02)); 0];               # enter Gdp(s)
p = 2; mu = 1; md = 1; mf = 2;       # set dimensions
sysfp = fdimodset(dss([Gup Gdp]),c =1,d = 2,f = 1,fs = 2);

# no faults
inputs = [commands disturbances 0*faults] 
yp, = timeresp(sysfp.sys, inputs, t; interpolation = "Tustin");

fdsys = FDSystem(Qd, x0 = [measurement(0.,0.)], α = [0], β = [10], γ = [.9]) 
isig = zeros(Int,ns)
for i = 1:ns
    tstep!(fdsys, yp[i,:], commands[i,1:1])
    isig[i] = fdsys.isig[1]
end
@test isempty(fdsys.indfault) && iszero(isig)


# with faults
inputs = [commands disturbances faults] 
yp, = timeresp(sysfp.sys, inputs, t; interpolation = "Tustin");

fdsys = FDSystem(Qd, x0 = [measurement(0.,0.)], α = [0], β = [10], γ = [.9]) 
isig = zeros(Int,ns)
for i = 1:ns
    tstep!(fdsys, yp[i,:], commands[i,1:1])
    isig[i] = fdsys.isig[1]
end
@test fdsys.indfault == [1]


# Example 5.9 - Solution of an AFDP using AFDSYN with rdm = 2

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gw = [1/(s+2); 0];                   # enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       # enter dimensions

# setup the synthesis model 
sysf = fdimodset(dss([Gu Gf Gw],minimal = true),c =1:mu, f = mu .+ (1:mf), n = (mu+mf) .+ (1:mw));

# call of AFDSYN with the options for stability degree -3 and the synthesis 
# of a scalar output filter
Q, Rfw, info = afdsyn(sysf, atol = 1.e-7, minimal = true, nullspace = false, smarg = -2, sdeg = -2, poles = [-2,-3],
                      rdim = 2, HDesign = [1 0]); info

# check synthesis results
R = fdIFeval(Q, sysf; minimal = true, atol = 1.e-7)
@test iszero(R.sys[:,[R.controls; R.disturbances]], atol = 1.e-7) &&
      iszero(R.sys[:,[R.faults;R.noise]]-Rfw.sys[:,[Rfw.faults;Rfw.noise]]) && 
      fdisspec(Rfw, block = false) == [info.S; info.S2] && 
      fdif2ngap(Rfw)[1] ≈ info.gap

Ts = 0.01;  
Qd = c2d(Q,Ts); 

t = 0:Ts:1; ns = length(t); tfault = 0.05; kfault = round(Int,tfault/Ts); fnorm_min = 0.02; 
commands = sin.(2pi*t).*ones(ns,1); noise = 0.1*randn(ns); faults = 0.5*[zeros(kfault,2); ones(ns-kfault,2)]


# no faults
inputs = [commands 0*faults noise] 
y, = timeresp(sysf.sys, inputs, t); 

# Fault detection
theta = zeros(ns)
# reval = zeros(ns)
fdsys = FDSystem(Qd, α = [0], β = [5], γ = [.8])
isig = zeros(Int,ns)
for i = 1:ns
    tstep!(fdsys, y[i,:], commands[i,1:1])
    theta[i] = fdsys.θ[1]
    isig[i] = fdsys.isig[1]
    #reval[i] = fdsys.re[1]
end
@test isempty(fdsys.indfault) && iszero(isig) # no false alarms

# with faults
inputs = [commands faults noise] 
y, = timeresp(sysf.sys, inputs, t); 

# Fault detection
# theta = zeros(ns)
# reval = zeros(ns)
fdsys = FDSystem(Qd, α = [0], β = [5], γ = [.8])
isig = zeros(Int,ns)
for i = 1:ns
    tstep!(fdsys, y[i,:], commands[i,1:1])
    # theta[i] = fdsys.θ[1]
    isig[i] = fdsys.isig[1]
    # reval[i] = fdsys.re[1]
end
@test fdsys.indfault == [1]


# Fault detection and isolation

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
Ts = 0.01;
sysd = c2d(dss(a,[bu bd bf],c,[du dd df]),Ts)[1];
sysfd = fdimodset(sysd, c = 1:mu, d = mu .+ (1:md) , f = (mu+md) .+ (1:mf));

# SFDI = fdigenspec(sysf,FDtol=1.e-5) # determine achievable structure matrix
# SFDI = fdigenspec(sysfd,FDtol=1.e-5) # determine achievable structure matrix
SFDI = Matrix(I(2) .> 0)
nb = size(SFDI,1)
@time Qd, Rf = efdisyn(sysfd, SFDI; FDtol=1.e-5, atol = 1.e-7, rdim = 1);

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
R = fdIFeval(Qd,sysfd); # form Q*[Gu Gd Gf;I 0 0];
@test iszero(vcat(R.sys...)[:,[R.controls;R.disturbances]],atol=1.e-7) &&
      iszero(vcat(R.sys...)[:,R.faults]-vcat(Rf.sys...),atol=1.e-7) &&
      order.(Qd.sys) == ones(Int,nb)

# check weak and strong fault detectability
@test fditspec(Rf,FDtol=1.e-5) == fdisspec(Rf) == SFDI 

Ts = 0.01
t = 0:Ts:1; ns = length(t); tfault = 0.05; kfault = round(Int,tfault/Ts); fnorm_min = 0.02; f1 = 1; f2 = 1
commands = sin.(2pi*t).*ones(ns,mu); disturbances = 0.1*randn(ns,md); faults = 10*[zeros(5,2); ones(ns-5,2)]*[f1 0; 0 f2]

# no faults
inputs = [commands disturbances 0*faults]
y,  = timeresp(sysfd.sys, inputs, t); 

Ne = length(Qd.sys)
#SFDI = fditspec(Rf)
fdisys = FDISystem(Qd,SFDI; α = .0*ones(Ne), β = 1.1*ones(Ne), γ = .7*ones(Ne))
theta = zeros(ns,Ne);
isig = zeros(Int,ns,Ne);
for i = 1:ns
    tstep!(fdisys, y[i,:], commands[i,:])
    theta[i,:] = fdisys.θ
    isig[i,:] = fdisys.isig
end
@test isempty(fdisys.indfault) && iszero(isig)

# with fault 1
f1 = 1; f2 = 0
faults = 0.5*[zeros(5,2); ones(ns-5,2)]*[f1 0; 0 f2]

inputs = [commands disturbances faults]
y,  = timeresp(sysfd.sys, inputs, t); 

#Qd = c2d(Q,Ts)
Ne = length(Qd.sys)
#SFDI = fditspec(Rf)
fdisys = FDISystem(Qd,SFDI; α = .0*ones(Ne), β = 1.1*ones(Ne), γ = .7*ones(Ne))
theta = zeros(ns,Ne);
isig = zeros(Int,ns,Ne);
for i = 1:ns
    tstep!(fdisys, y[i,:], commands[i,:])
    theta[i,:] = fdisys.θ
    isig[i,:] = fdisys.isig
end
@test fdisys.indfault == [1]

# with fault 2
f1 = 0; f2 = 1
faults = 0.5*[zeros(5,2); ones(ns-5,2)]*[f1 0; 0 f2]

inputs = [commands disturbances faults]
y,  = timeresp(sysfd.sys, inputs, t); 

Ne = length(Qd.sys)
fdisys = FDISystem(Qd,SFDI; α = .0*ones(Ne), β = 1.1*ones(Ne), γ = .7*ones(Ne))
theta = zeros(ns,Ne);
isig = zeros(Int,ns,Ne);
for i = 1:ns
    tstep!(fdisys, y[i,:], commands[i,:])
    theta[i,:] = fdisys.θ
    isig[i,:] = fdisys.isig
end
@test fdisys.indfault == [2]

# with fault 1 and fault 2
f1 = 1; f2 = 1
faults = 0.5*[zeros(5,2); ones(ns-5,2)]*[f1 0; 0 f2]

inputs = [commands disturbances faults]
y,  = timeresp(sysfd.sys, inputs, t); 

Ne = length(Qd.sys)
fdisys = FDISystem(Qd,SFDI; α = .0*ones(Ne), β = 1.1*ones(Ne), γ = .7*ones(Ne))
theta = zeros(ns,Ne);
isig = zeros(Int,ns,Ne);
for i = 1:ns
    tstep!(fdisys, y[i,:], commands[i,:])
    theta[i,:] = fdisys.θ
    isig[i,:] = fdisys.isig
end
@test fdisys.indfault == [1,2]


## FDI example with multiple fault detection 
p = 3; mf = 2; mu = 2
sysf = fdimodset(rss(1,p,mf+mu;stable = true),controls = 1:mu, faults = mu+1:mu+mf);

SFDI = fdigenspec(sysf, atol = 1.e-7)

#SFDI = Matrix(I(2) .> 0)

@time Q, R, info = efdisyn(sysf, SFDI; nullspace = false, atol = 1.e-7, minimal = true); info

@test iszero(vcat(R.sys...)-vcat(R.sys...)[:,Rf.faults],atol=1.e-7) && iszero(vcat(R.sys...)[:,R.controls],atol=1.e-7) &&
      info.HDesign == [[0.0 1.0 0.0], [0.0 1.0], [0.0 1.0]] &&
      fditspec(R, FDtol = 1.e-6) == SFDI  && fdisspec(R, FDGainTol = 1.e-3) == SFDI


t = 0:0.01:1; ns = length(t); f1 = 1; f2 = 1
commands = sin.(2pi*t).*ones(ns,2); disturbances = randn(ns,0); faults = 0.5*[zeros(5,2); ones(ns-5,2)]*[f1 0; 0 f2]
inputs = [commands disturbances faults]
y,  = timeresp(sysf.sys, inputs, t); 

Qd = c2d(Q,0.01)
Ne = size(SFDI,1)
fdisys = FDISystem(Qd,SFDI; α = 15*ones(Ne), β = .5*ones(Ne), γ = .8*ones(Ne))
nsys = length(Qd.sys);
#theta = zeros(ns,nsys);
isig = zeros(Int,ns,nsys);
# define an extension of the structure matrix and the corresponding vector with the indices of multiple faults  
multifaults = [[1,2]]
SFDIext = sum(SFDI[:,multifaults[1]],dims=2) .> 0
for i = 1:ns
    tstep!(fdisys, y[i,:], commands[i,:])
    #theta[i,:] = fdisys.θ
    isig[i,:] = fdisys.isig
    if isempty(fdisys.indfault) 
       ind = findfirst(!iszero,([SFDIext[:,i] == fdisys.isig for i in 1:size(SFDIext,2)] .== 1)) 
       fdisys.indfault = isnothing(ind) ? Int[] : multifaults[ind]
    end
end
@test fdisys.indfault == [1,2]
# plot(t,y)
# plot!(t,theta)
# plot!(t,isig)


# Example 5.11c - Solution of an AFDIP 

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gw = [1/(s+2); 0];                   # enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       # enter dimensions
SFDI = eye(mf) .> 0

# build the synthesis model with additive faults 
sysf = fdimodset(dss([Gu Gf Gw]), c = 1:mu,f = mu .+ (1:mf), n = (mu+mf) .+ (1:mw));

Q, R, info = afdisyn(sysf, SFDI; smarg =-3, sdeg = -3)

# check synthesis conditions: Qt*[Gu Gd;I 0] = 0 and Qt*[Gf; 0] = Rf
Rt = fdIFeval(Q,sysf) # form Qt*[Gu Gd Gf;I 0 0];
@test iszero(vcat(Rt.sys...)[:,[Rt.controls;Rt.disturbances]],atol=1.e-7) &&
      iszero(vcat(Rt.sys...)[:,[Rt.faults;Rt.noise]]-vcat(R.sys...),atol=1.e-7) &&
      fdif2ngap(R,SFDI)[1] ≈ info.gap

# check weak and strong fault detectability
@test fditspec(R) == fdisspec(R) == SFDI 


Ts = 0.01;  
Qd = c2d(Q,Ts); 

t = 0:Ts:1; ns = length(t); tfault = 0.05; kfault = round(Int,tfault/Ts); fnorm_min = 0.02; 
commands = sin.(2pi*t).*ones(ns,1); noise = 0.1*randn(ns); faults = 0.5*[zeros(kfault,2); ones(ns-kfault,2)]


# no faults
inputs = [commands 0*faults noise] 
y, = timeresp(sysf.sys, inputs, t); 

# Fault detection
#theta = zeros(ns)
Ne = size(SFDI,1)
fdisys = FDISystem(Qd,SFDI; α = 15*ones(Ne), β = .5*ones(Ne), γ = .8*ones(Ne))
nsys = length(Qd.sys);
#theta = zeros(ns,nsys);
isig = zeros(Int,ns,Ne)
for i = 1:ns
    tstep!(fdisys, y[i,:], commands[i,1:1])
    #theta[i] = fdisys.θ[1]
    isig[i] = fdisys.isig[1]
end
@test isempty(fdisys.indfault) && iszero(isig) # no false alarms

# two faults
inputs = [commands faults noise] 
y, = timeresp(sysf.sys, inputs, t); 

# Fault detection
#theta = zeros(ns)
Ne = size(SFDI,1)
fdisys = FDISystem(Qd,SFDI; α = 15*ones(Ne), β = .5*ones(Ne), γ = .8*ones(Ne))
nsys = length(Qd.sys);
#theta = zeros(ns,nsys);
isig = zeros(Int,ns,Ne)
for i = 1:ns
    tstep!(fdisys, y[i,:], commands[i,1:1])
    #theta[i] = fdisys.θ[1]
    isig[i] = fdisys.isig[1]
end
@test fdisys.indfault == [1,2]





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
@time Q, R, info = emmsyn(sysf, Mr; atol, sdeg, minimal = false); info

# check solution
G = [sysf.sys; eye(mu,mu+mf)]; F = [zeros(mf,mu) info.M*Mr.sys];
@test iszero(Q.sys*G-F; atol)

# setup with additive faults
commands = sin.(2pi*t).*ones(ns,2); disturbances = randn(ns,0); faults = 0.5*[zeros(kfault,2); ones(ns-kfault,2)]*[1 0;0 1]
inputs = [commands disturbances faults] # 
y, = timeresp(sysf.sys, inputs, t); 

Ts = 0.01;  
Qd = c2d(Q,Ts); 
fdsys = FDSystem(Qd,fditspec(Mr), x0 = zeros(Measurement{Float64},order(Qd.sys)), α = [0;0], β = [10;10], γ = [.9,.9])
isig = zeros(Int,ns,2)
for i = 1:ns
    tstep!(fdsys, y[i,:]*measurement(1,0.2), commands[i,:])
    isig[i,:] = fdsys.isig
end
@test fdsys.indfault == [1,2]


# Example 5.16c - Solution of an H∞ AMMP 

# define system with control, noise and actuator fault inputs 
A = [-.8 0 0; 0 -.5 .6; 0 -.6 -.5];
Bu = [1 1;1 0;0 1]; Bw = 0.25*[0 0;0 1;1 0]; Bf = Bu;
C = [0 1 1; 1 1 0]; Du = zeros(2,2); Dw = zeros(2,2);  
p, mu = size(Du); mf = mu; mw = size(Dw,2); 

# setup synthesis model with additive actuator faults with Bf = Bu, Df = Du
sysf = fdimodset(dss(A,[Bu Bw],C,[Du Dw]), c = 1:mu, n = mu .+ (1:mw),f = 1:mu); 

# define Mr(s) = I
Mr = FDFilterIF(dss(eye(mf)); mf)

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


SFDI = fdisspec(Mr)

Ts = 0.01;  
Qd = c2d(Q,Ts); 

t = 0:Ts:1; ns = length(t); tfault = 0.05; kfault = round(Int,tfault/Ts); fnorm_min = 0.02; 
commands = sin.(2pi*t).*ones(ns,2); noise = 0.1*randn(ns,2); faults = 0.5*[zeros(kfault,2); ones(ns-kfault,2)]


# no faults
inputs = [commands 0*faults noise] 
y, = timeresp(sysf.sys, inputs, t); 

# Fault detection
#theta = zeros(ns)
Ne = size(SFDI,1)
fdsys = FDSystem(Qd,SFDI; α = 1*ones(Ne), β = .5*ones(Ne), γ = .8*ones(Ne))
nsys = length(Qd.sys);
#theta = zeros(ns,nsys);
isig = zeros(Int,ns,Ne)
for i = 1:ns
    tstep!(fdsys, y[i,:], commands[i,:])
    #theta[i] = fdisys.θ[1]
    isig[i] = fdsys.isig[1]
end
@test isempty(fdsys.indfault) && iszero(isig) # no false alarms

# two faults
inputs = [commands faults noise] 
y, = timeresp(sysf.sys, inputs, t); 

# Fault detection
#theta = zeros(ns)
Ne = size(SFDI,1)
fdsys = FDSystem(Qd,SFDI; α = 1*ones(Ne), β = .5*ones(Ne), γ = .8*ones(Ne))
nsys = length(Qd.sys);
#theta = zeros(ns,nsys);
isig = zeros(Int,ns,Ne)
for i = 1:ns
    tstep!(fdsys, y[i,:], commands[i,:])
    #theta[i] = fdisys.θ[1]
    isig[i] = fdsys.isig[1]
end
@test fdsys.indfault == [1,2]






end    




@testset "Model detection" begin

## Actuator faults - exact models
s = rtf('s');
k = 14; N = 4; mu = 1;
sys1 = dss(k/(s+k));           # nominal model
sys2 = dss(0.5*k/(s+0.5*k));   # 50% loss of efficiency
sys3 = dss(10*k/(s+10*k));     # disconnected actuator
sys4 = dss(0.01*k/(s+0.01*k)); # stall load fault model

# setup synthesis model
sysm = mdmodset([sys1, sys2, sys3, sys4], controls = 1:mu);
@time distgap, fpeak = mddist(sysm)

@time Q, R, info = emdsyn(sysm; sdeg = -15, poles = [-20]); 
@test sortperm(distgap[1,:]) == sortperm(info.MDperf[1,:]) && 
      sortperm(distgap[2,:]) == sortperm(info.MDperf[2,:]) && 
      sortperm(distgap[3,:]) == sortperm(info.MDperf[3,:]) && 
      sortperm(distgap[4,:]) == sortperm(info.MDperf[4,:])


Qd = c2d(Q,0.001)
mdsys = MDSystem(Qd; α = 15*ones(4), β = .5*ones(4), γ = .9*ones(4))
t = 0:0.001:1; ns = length(t); 
commands = sin.(10pi*t).*ones(ns,1); 
inputs = commands
nsys = length(Qd.sys);
for k = 1:nsys
    y, = timeresp(sysm.sys[k], inputs, t, interpolation = "Tustin"); 
    for i = 1:ns
        tstep!(mdsys, y[i,:], commands[i,:])
    end
    @test mdsys.indminim == k   
    @test mdsys.indmodel == k   
end


## Actuator faults - models with gain uncertainties
kp = 14*measurement(1,0.1); 
sys1p = dss(kp/(s+kp));           # nominal model
sys2p = dss(0.5*kp/(s+0.5*kp));   # 50% loss of efficiency
sys3p = dss(10*kp/(s+10*kp));     # disconnected actuator
sys4p = dss(0.01*kp/(s+0.01*kp)); # stall load fault model

# setup synthesis model
sysmp = mdmodset([sys1p, sys2p, sys3p, sys4p], controls = 1:mu);
nsys = length(Qd.sys);

mdsysp = MDSystem(Qd; x0 = [zeros(Measurement{Float64},order(Qd.sys[i])) for i in 1:nsys],α = 15*ones(4), β = .5*ones(4), γ = .9*ones(4))
t = 0:0.001:1; ns = length(t); 
commands = sin.(10pi*t).*ones(ns,1); 
inputs = commands
for k = 1:nsys
    yp, = timeresp(sysmp.sys[k], inputs, t, interpolation = "Tustin"); 
    for i = 1:ns
        tstep!(mdsysp, yp[i,:], commands[i,:])
    end
    @test mdsysp.indminim == k   
    @test mdsysp.indmodel == k   
end

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


Qd = c2d(Q,0.001)
mdsys = MDSystem(Qd; α = 15*ones(N), β = .5*ones(N), γ = .9*ones(N))
t = 0:0.001:1; ns = length(t); 
commands = sin.(10pi*t).*ones(ns,2); 
inputs = commands
nsys = length(Qd.sys);
for k = 1:nsys
    y, = timeresp(sysm.sys[k], inputs, t); 
    for i = 1:ns
        tstep!(mdsys, y[i,:], commands[i,:])
    end
    @test mdsys.indminim == k   
    @test mdsys.indmodel == k   
end

end
end