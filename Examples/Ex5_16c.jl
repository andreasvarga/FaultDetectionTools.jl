module Ex5_16c
using FaultDetectionTools, DescriptorSystems, Measurements, Test

# Example 5.16c - Solution of an H∞ AMMP 
println("Example 5.16c with Fig5.2")

# define system with control, noise and actuator fault inputs 
A = [-.8 0 0; 0 -.5 .6; 0 -.6 -.5]
Bu = [1 1;1 0;0 1]; Bw = 0.25*[0 0;0 1;1 0]; Bf = Bu
C = [0 1 1; 1 1 0]; Du = zeros(2,2); Dw = zeros(2,2)  
p, mu = size(Du); mf = mu; mw = size(Dw,2) 

atol = 1.e-7;                # tolerance for rank tests

# set up synthesis model with additive actuator faults 
# with Bf = Bu, Df = Du
sysf = fdimodset(dss(A,[Bu Bw],C,[Du Dw]), c = 1:mu, 
                 n = mu .+ (1:mw), f = 1:mu) 

# define Mr(s) = I
Mr = FDFilterIF(dss(eye(mf)); mf)

Q, R, info = ammsyn(sysf, Mr; atol, nullspace = false, sdeg = -10,
                            reltol = 5.e-4, normalize = "dcgain")

gamma_opt0 = info.gammaopt0 # optimal initial problem performance 
gamma_opt  = info.gammaopt  # optimal performance 
gamma_sub  = info.gammasub  # suboptimal performance 

# check suboptimal solution
Ge = [sysf.sys[:,[sysf.faults;sysf.noise]]; zeros(mu,mf+mw)] 
Me = [info.M*Mr.sys zeros(mf,mw)]
@test iszero(R.sys-Q.sys*Ge; atol) && 
      fdimmperf(R,info.M*Mr) ≈ gamma_sub && 
      glinfnorm(Me-Q.sys*Ge)[1] ≈ gamma_sub

# compute step responses with uncertainty bounds
# define uncertain parameters p = (ρ1,ρ2) 
ρ1 = measurement(0, 0.25); ρ2 = measurement(0, 0.25);
# build uncertain state matrix A(p)
Ap = [-.8 0 0;0 -.5*(1+ρ1) .6*(1+ρ2); 0 -0.6*(1+ρ2) -0.5*(1+ρ1)];
# build the uncertain system with controls and faults
uncsysf = fdimodset(dss(Ap,Bu,C,Du), c = 1:mu, f = 1:mu); 

# determine the internal form with uncertainties 
Rt = fdIFeval(Q, uncsysf) 

# determine the step responses with uncertainty bounds
Rtd = c2d(Rt.sys[:,[Rt.faults;Rt.controls]], 0.01,"Tustin")[1]
uncy, tout, _ = stepresp(Rtd, 10) 
y = getfield.(uncy,:val)        # nominal step responses
yerr = getfield.(uncy,:err)     # error bounds
yl = y .- yerr; yu = y .+ yerr  # lower and upper bounds      

# plot step responses
include("Fig5_2.jl")
Fig5_2 = f

end
using Main.Ex5_16c
Ex5_16c.Fig5_2
