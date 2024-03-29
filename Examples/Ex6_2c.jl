module Ex6_2c

using FaultDetectionTools, DescriptorSystems, LinearAlgebra, Test

# Example 6.2c - Solution of an AMDP
println("Example 6.2c with Fig6.3")

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

# define multiple physical fault model Gui = Gu*Gamma_i with noise inputs
sysuw = similar(Vector{DescriptorStateSpace},N)

for i = 1:N
  sysuw[i] = dss(A,[Bu*diagm(Gamma[i,:]) Bw],C,[Du Dw]);
end

# setup synthesis model
sysm = mdmodset(sysuw, controls = 1:mu, noise = mu+1:mu+mw);


# use nonminimal design with  AMDSYN
Q, R, info = amdsyn(sysm,sdeg = -1, poles = [-1], minimal = false); 

distinf = info.MDperf
@test all(diag(distinf) .< 1.e-7) && all(distinf.+eps(2.)+I .>= 1)

println("Norms of final residual models")
display(distinf )

println("Resulting gaps")
display(info.MDgap)

# plot residual responses
include("Fig6_3.jl")
Fig6_3 = fig
display(Ex6_2c.Fig6_3)

end
using Main.Ex6_2c
