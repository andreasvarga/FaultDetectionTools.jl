module Ex2_4
using FaultDetectionTools, DescriptorSystems, Measurements, LinearAlgebra, Test

# Example 2.4 - Generation of a synthesis model with noise inputs from a LPV uncertain model 
println("Example 2.4")
# define ρ1 and ρ2 as uncertain parameters
ρ1 = measurement(0, rand()); ρ2 = measurement(0, rand());

# define A(ρ1,ρ2), Bu, C, Du
n = 3; mu = 2; p = 2;       # enter dimensions
A = [ -.8 0 0;
      0 -0.5*(1+ρ1) 0.6*(1+ρ2);
      0 -0.6*(1+ρ2) -0.5*(1+ρ1) ];
Bu = [ 1 1; 1 0; 0 1]; C = [0 1 1;1 1 0]; Du = zeros(p,mu);

# build ΔA and ΔC
ΔA = getfield.(A,:err)
ΔC = zeros(p,n)

# determine number of noise inputs 
mw = rank(ΔA)

# compute orthogonal range of [ΔA;0] 
U = svd([ΔA;ΔC]).U[:,1:mw]


# determin Bw and Dw 
Bw = U[1:n,:]
Dw = U[n+1:end,:]

# setup synthesis model with additive actuator faults with Bf = Bu, Df = Du
sysf = fdimodset(dss(getfield.(A,:val),[Bu Bw],C,[Du Dw]), c = 1:mu, n = mu .+ (1:mw), f = 1:mu); 
display(sysf)

end # module
using Main.Ex2_4
