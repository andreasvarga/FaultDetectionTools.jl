module Ex5_13
using FaultDetectionTools, DescriptorSystems, LinearAlgebra, Test

# Example 5.13 - Solution of an EMMP 
println("Example 5.13")

# define s as an improper transfer function
s = rtf('s');
# enter Gu(s), Gf(s) and Mr(s)
Gu = [s/(s^2+3*s+2) 1/(s+2);
     s/(s+1) 0;
      0 1/(s+2)];
Gf = [s/(s^2+3*s+2) 1/(s+2);
     s/(s+1) 0;
      0 1/(s+2)];
Mr = dss(eye(2));                  # enter Mr(s)
p, mf = size(Gf); mu = size(Gu,2);

# compute left nullspace basis as Q1(s) = [ I -Gu(s) ]
Q1 = dss([I -Gu]); Rf = dss(Gf)

# check solvability condition
rank(evalfr(Rf,rand())) == mf || error("No solution exist")

# check for unstable or infinite zeros
gzero(Rf)   # zeros at infinity and in the origin exist 

atol = 1.e-7                # tolerance for rank tests
sdeg = -1                   # set stability degree

# using minimal dynamic covers, compute Q2 such that Q2*Q1 is a lower order
# left annihilator and Q2*Rf full row rank; 
# select basis vectors [2 3] of Q1 to combine them with basis vector 1
# to obtain Rf_Q = Q2*[Rf Q1], with Q2*Rf in the leading mf columns and 
# Q2*Q1 in the trailing p+mu columns
cinp = [ 2, 3, 1]
Rf_Q = glmcover1([Rf[cinp,:] Q1[cinp,:]], mf; atol)[1]

# compute the irreducible realization of Qtilde = Mr*(inv(Rf)*Q) by 
# first solving the linear rational matrix equation Rf*X = Q
X = grsol(Rf_Q, p+mu; atol)[1]
Qtilde = gir(Mr*X; atol)


# compute stable and proper Q = Q4*Qtilde with suitable diagonal Q4 = M
Q = dss(zeros(0,p+mu)); M = dss(zeros(0,0));
for i = 1:mf
    Qi, Mi = glcf(Qtilde[i,:]; atol, sdeg);
    # scale with gain to fit example
    sc = zpk(dss2rm(Mi; atol1 = atol, atol2 = atol)[1,1])[3]
    global Q = [Q; Qi/sc]; 
    global M = append(M, Mi/sc);
end

# convert to standard state space representation
Q = gss2ss(Q)[1]; M = gss2ss(M)[1] 

# display results
println("Q = $(dss2rm(Q, atol = 1.e-7))")
println("M = $(dss2rm(M, atol = 1.e-7))")

# check solution
G = dss([Gu Gf;eye(mu,mu+mf)]); F = [zeros(mf,mu) M*Mr];
@test iszero(Q*G-F; atol)

end # module
using Main.Ex5_13
