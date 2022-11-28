module Ex7_4
using DescriptorSystems, Test

# Example 7.4 - Nullspace-based synthesis 
# Uses only DescriptorSystems.jl 
println("Example 7.4")

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

# compute [Q Rf], where Q is a left nullspace basis of
# [Gu Gd;I 0] and Rf = Q*[Gf;0];
Q_Rf, info = glnull([sys;eye(mu,mu+md+mf)],mf);
info.degs

# determine 1-st and 2-nd order scalar output designs
Q1_Rf1 = glmcover1([0 1;eye(2)]*Q_Rf,1)[1];    #  [Q1 Rf1]
Q2_Rf2 = glmcover1([-1 1;eye(2)]*Q_Rf,1)[1];    #  [Q2 Rf2]

# compute stable left coprime factorizations
Q1_Rf1 = glcf(Q1_Rf1, sdeg=-1, smarg = -1)[1];
# choose poles to meet example
Q2_Rf2 = glcf(Q2_Rf2, evals = [-1. + 0.3992im, -1. - 0.3992im], 
              sdeg = -1, smarg = -1)[1];

# check admissibility
# compute Q1 and Rf1; 
Q1 = Q1_Rf1[:,1:p+mu]; Rf1 = Q1_Rf1[:,p+mu+1:end];
# compute reciprocal sensitivity condition
g1 = abs.(evalfr(Rf1,im)); rscond1 = maximum(g1)/minimum(g1)
# compute Q2 and Rf2;  
Q2 = Q2_Rf2[:,1:p+mu]; Rf2 = Q2_Rf2[:,p+mu+1:end];
# compute reciprocal sensitivity condition
g2 = abs.(evalfr(Rf2,im)); rscond2 = maximum(g2)/minimum(g2)
println("|Rf1(im)| = $g1")
println("rscond1 = $rscond1 ")
println("|Rf2(im)| = $g2")
println("rscond2 = $rscond2")
# check finite reciprocal sensitivity conditions
@test !isinf(rscond1) && !isinf(rscond2)

# display results
println("Q1 = $(dss2rm(Q1))")
println("Q2 = $(dss2rm(Q2))")

end # module
using Main.Ex7_4
