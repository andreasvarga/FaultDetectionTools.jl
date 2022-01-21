module Ex6_2
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
#using Plots
using Test

# Example 6.2 - Solution of an AMDP
println("Example 6.2")

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
# define multiple physical fault model Gui = Gu*Gamma_i
sysuw = similar(Vector{DescriptorStateSpace},N)

for i = 1:N
  sysuw[i] = dss(A,[Bu*diagm(Gamma[i,:]) Bw],C,[Du Dw]);
  #sysu[i] = gir(dss(A,Bu*diagm(Gamma[i,:]),C,0),atol = 1.e-7);
end

# setup initial full order model detector Q1i = [I -Gui]
Q1 = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    Q1[i] = [eye(p) -sysuw[i][:,1:mu]];
end

# solve a minimum dynamic cover problem
# select row 1 of h*Q1i to combine it with rows 1:4
# of Q1i; the result is a least order Qi = Q2i*Q1i
atol = 1.e-7;                     # set tolerance
Q = similar(Vector{DescriptorStateSpace},N)
R = similar(Matrix{DescriptorStateSpace},N,N)
distinf = similar(Matrix{Float64},N,N)
for i = 1:N
    rwi = gir(Q1[i][:,1:p]*sysuw[i][:,mu+1:m]; atol);
    gi, go = goifac(rwi; atol);
    Q[i] = gminreal(go\Q1[i], atol1 = atol, atol2 = atol);
    for j = 1:N
        R[i,j] = gir(Q[i]*[sysuw[j]; eye(mu,m)];atol);
        distinf[i,j] = glinfnorm(R[i,j][:,1:mu])[1]
    end
end

# scale Qi and Rij; determine gap 
beta = similar(Vector{Float64},N)
for i = 1:N
  scale = minimum(distinf[i,[1:i-1; i+1:N]]);
  distinf[i,:] = distinf[i,:]/scale;
  Q[i] = Q[i]/scale;
  for j = 1:N
     R[i,j] = R[i,j]/scale;
  end
  beta[i] = scale;
end
gap = beta

@test norm(diag(distinf)) < 1.e-7


println("Norms of final residual models")
display(distinf)

println("Resulting gaps")
display(gap)

end # module
