module Ex6_1c
using FaultDetectionTools, DescriptorSystems, LinearAlgebra, Test

# Example 6.1c - Solution of an EMDP
println("Example 6.1c with Fig6.1(fig1) and Fig6.2(fig2)")

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


# call of EMDSYN with the options for stability degree -1 and pole -1 for
# the filters, tolerance and a design matrix H to form a linear combination
# of the left nullspace basis vectorsH = [ 0.7645 0.8848 0.5778 0.9026 ];
H = [ 0.7645 0.8848 0.5778 0.9026 ];
Q, R, info = emdsyn(sysm, sdeg = -1, poles = [-1], HDesign = H); 
distinf = info.MDperf

# check model detectability 
@test all(diag(distinf) .< 1.e-7) && all(distinf.+eps(2.)+I .>= 1)

println("Norms of final residual models")
display(distinf)

# plot distance mappings 
include("Fig6_1.jl")
Fig6_1 = fig1

# plot residual responses
include("Fig6_2.jl")
Fig6_2 = fig2
#display(Ex6_1c.Fig6_1)
#display(Ex6_1c.Fig6_2)

end
using Main.Ex6_1c
