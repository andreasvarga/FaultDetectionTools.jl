module Ex6_1
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
#using Plots
using Test

# Example 6.1 - Solution of an EMDP
println("Example 6.1")

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
  #sysu[i] = gir(dss(A,Bu*diagm(Gamma[i,:]),C,0),atol = 1.e-7);
end

# setup initial full order model detector Q1i = [I -Gui]
Q1 = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    Q1[i] = [eye(p) -sysu[i]];
end

# solve a minimum dynamic cover problem
# select row 1 of h*Q1i to combine it with rows 1:4
# of Q1i; the result is a least order Qi = Q2i*Q1i
h = [ 0.7645   0.8848   0.5778   0.9026];
atol = 1.e-7;                     # set tolerance
Q = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
  Q[i] = glmcover1([h;eye(p)]*Q1[i],1; atol)[1];  
  #Q[i] = glcf(Q[i]; atol, sdeg = -1); 
end

# compute Rij and their norms
R = similar(Matrix{DescriptorStateSpace},N,N)
distinf = similar(Matrix{Float64},N,N)
for i = 1:N
  for j = 1:N
    temp = Q[i]*[sysu[j]; eye(mu)];
    R[i,j] = gir(temp; atol);
    distinf[i,j] = glinfnorm(R[i,j])[1]
  end
end

# scale Qi and Rij
distinf1 = similar(Matrix{Float64},N,N)
for i=1:N
  scale = 1/minimum(distinf[i,[1:i-1; i+1:N]])
  Q[i] = scale*Q[i];
  for j = 1:N
     R[i,j] = scale*R[i,j];
     distinf1[i,j] = glinfnorm(R[i,j])[1]
  end
end
@test iszero(diag(distinf1)) 

distinf0 = similar(Matrix{Float64},N,N)
for i = 1:N
  for j = 1:N
    distinf0[i,j] = glinfnorm(Q1[i]*[sysu[j]; eye(mu)])[1]
  end
end

println("Norms of initial residual models")
display(distinf0)


println("Norms of final residual models")
display(distinf1)

# display(surface(distinf0, xlim = (0,10), ylim = (0,10), 
#                 title = "Norms of initial residual models",
#                 xlabel = "Model numbers", 
#                 ylabel = "Residual numbers"))

#display(contour(distinf0, fill=true))

# display(surface(distinf1, xlim = (0,10), ylim = (0,10), 
#                 title = "Norms of final residual models",
#                 xlabel = "Model numbers", 
#                 ylabel = "Residual numbers"))

end # module
