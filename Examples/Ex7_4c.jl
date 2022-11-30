module Ex7_4c
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 7.4c - Nullspace-based synthesis 

println("Example 7.4c")

# define the state-space realizations of Gu, Gd and Gf
n = 5; p = 3; mu = 3; md = 1; mf = 5;       # enter dimensions
# define matrices of the state space realization
A = [ 0 0 1.132 0 -1;
0 -0.0538 -0.1712 0 0.0705;
0 0 0 1 0;
0 0.0485 0 -0.8556 -1.013;
0 -0.2909 0 1.0532 -0.6859];
Bu = [0 0 0;-0.12 1 0;0 0 0;4.419 0 -1.665;1.575 0 -0.0732];
Bd = Bu[:,mu]; Bf = [zeros(n,p) Bu[:,1:mu-1]];
C=eye(p,n); Du=zeros(p,mu); Dd=zeros(p,md); Df=eye(p,mf);
sys = dss(A,[Bu Bd Bf],C,[Du Dd Df]);        # define system

# setup the synthesis model 
sysf = fdimodset(dss(A,[Bu Bd Bf],C,[Du Dd Df]),c =1:mu, d = mu .+ (1:md), f = (mu+md) .+ (1:mf));

# call of EFDSYN with the options for stability degree -1 and the synthesis 
# of a scalar output filter using two design matrices
Q1, Rf1, info1 = efdsyn(sysf, HDesign = [0 1]; sdeg = -1, smarg = -1, rdim = 1); 
Q2, Rf2, info2 = efdsyn(sysf, HDesign = [-1 1]; poles = [-1.0000 + 0.3992im, -1.0000 - 0.3992im], sdeg = -1, smarg = -1, rdim = 1); 

# check admissibility, i.e., finite reciprocal sensitivity conditions
rscond1 = 1/fdiscond(Rf1,1)[1]
rscond2 = 1/fdiscond(Rf2,1)[1]
@test !isinf(rscond1) && !isinf(rscond1)
println("rscond1 = $rscond1 ")
println("rscond2 = $rscond2")

end # module
using Main.Ex7_4c
