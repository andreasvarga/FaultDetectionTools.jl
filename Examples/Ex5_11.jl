module Ex5_11
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.11 - Solution of an AFDIP 
println("Example 5.11")

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gw = [1/(s+2); 0];                   # enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       # enter dimensions
S = eye(mf) .> 0

# Procedure AFDI

# Step 1): choose the left nullspace as Q1 = [I -Gu] and 
# form Rf1 = Q1*[Gf;0] and Rw1 = Q1*[Gw;0]
Q1 = dss([eye(p) -Gu]); Rf1 = dss(Gf); Rw1 = dss(Gw);
atol = 1.e-7;


# Step 2): determine Q[i] and corresponding Rf[i] and Rw[i]
 
# initialize overall filter Q and fault detection system Rf
nb = size(S,1);      # number of necessary filters 
Qt = similar(Vector{typeof(Q1)},nb)
Rfwt = similar(Vector{typeof(Q1)},nb)
scale = [1, -1]  # set appropriate scaling to match example

for i = 1:nb
    # perform Procedure AFD or EFD to compute Q[i]
    indd = Vector(1:mf)[S[i,:] .== false] 
    Qi1 = glnull(Rf1[:,indd])[1]
    # initialize Q[i], Rf[i] and Rw[i]
    Qi = Qi1*Q1; Rfi = Qi1*Rf1; Rwi = Qi1*Rw1; 

    if norm(evalfr(Rwi,rand())) > 0.0001
       # compute the quasi-co-outer-co-inner factorization
       Rwi, Rwo, = goifac(Rwi;atol) 
       # update 
       Qi = Rwo\Qi; Rfi = Rwo\Rfi 
    end   
    # update the solution if [Q[i], Rf[i], Rw[i]] is improper or unstable
    Qi_Rfi_Rwi, M = glcf([Qi Rfi Rwi], sdeg = -3, smarg = -3, mindeg = true)
   
    # adjust denominator M to unit infinity norm to match example
    Mnorm = ghinfnorm(M)[1]
    Qt[i] = (scale[i]/Mnorm)*Qi_Rfi_Rwi[:,1:p+mu]
    Rfwt[i] = (scale[i]/Mnorm)*Qi_Rfi_Rwi[:,p+mu+1:end]

    println("Q[$i] = $(dss2rm(Qt[i], atol = 1.e-7))")
    println("Rf[$i] = $(dss2rm(Rfwt[i][:,1:mf], atol = 1.e-7))")
    println("Rw[$i] = $(dss2rm(Rfwt[i][:,mf+1:end], atol = 1.e-7))")
end
Q = FDIFilter(gminreal.(Qt; atol), p, mu)
R = FDIFilterIF(gminreal.(Rfwt; atol); mf, mw)
println("Q = ")
display(Q)
println("R = ")
display(R)
println("gap = $(fdif2ngap(R,S)[1])")
   
end # module
