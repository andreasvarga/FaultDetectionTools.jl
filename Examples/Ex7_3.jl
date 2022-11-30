module Ex7_3
using DescriptorSystems, LinearAlgebra

# Example 7.3 - Illustrating nullspace-based synthesis paradigm
# Uses only DescriptorSystems.jl 
println("Example 7.3")

A = diagm([ -2, -3, -1, 2, 3]);
Bu = [1 1 0 0 0]'; Bd = [0 0 2 0 0]';
Bf = [0 0 0 2 2;0 0 0 0 0]';
C = [-1 0 -1 1.5 0; 0 -1 0 0 2.5];
Du = [ 1 1]'; Dd = [1 0]'; Df = [1 0; 1 1];
p = 2; mu = 1; md = 1; mf = 2;       # enter dimensions
sys = dss(A,[Bu Bd Bf],C,[Du Dd Df]); # define system

# compute [Q Rf], where Q is a left nullspace basis of
# [Gu Gd;I 0] and Rf = Q*[Gf;0];
Q_Rf = glnull([sys;eye(mu,mu+md+mf)],mf)[1];

# stabilize using left coprime factorization
Q_Rf, M = glcf(Q_Rf, sdeg = -3, smarg = -2);

# compute Q and Rf; scale to meet example
Q = -sqrt(2)*Q_Rf[:,1:p+mu]; Rf = -sqrt(2)*Q_Rf[:,p+mu+1:end];

# display results
println("Q = $(dss2rm(Q))")
println("Rf = $(dss2rm(Rf))")
println("M = $(dss2rm(M))")

end # module
using Main.Ex7_3
