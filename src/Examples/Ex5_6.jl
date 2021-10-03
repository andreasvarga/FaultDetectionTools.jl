module Ex5_6
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.6 - Solution of an AFDP 
println("Example 5.6")

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gw = [1/(s+2); 0];                   # enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       # enter dimensions

Q = dss([eye(p) -Gu]); Rf = dss(Gf); Rw = dss(Gw);
atol = 1.e-7;

# check solvability using a random frequency
if minimum(maximum(abs.(evalfr(Rf,rand())),dims=1)) > 0.01

    # compress Rw to a full row rank matrix
    rw = rank(evalfr(Rw,rand())); nb = size(Q,1);
    if rw < nb
       h = ones(rw,nb);  # usable only for rw = 1
       # use alternatively h = rand(rw,nb);
       Q = h*Q; Rf = h*Rf; Rw = h*Rw;
    end
 
    # compute the quasi-co-outer-co-inner factorization 
    Rwi, Rwo, info = goifac(Rw; atol) 
 
    # compute optimal filter (for standard case)
    Q = gir(Rwo\Q; atol )              # update Q
    Rf = gir(Rwo\Rf; atol); Rw = Rwi;  # update Rf and Rw
  
    # check for poles on the extended imaginary axis
    poles = gpole([Q Rf])
    if any(isinf.(poles)) || minimum(abs.(real(poles))) < 0.0001
        # compute a stable and proper left coprime factorization
        # of [Q Rf Rw] with desired stability degree -3
        # update if solution is improper or unstable
        Q_Rf_Rw, M = glcf(gir([Q Rf Rw]; atol); atol, sdeg = -3, smarg = -3);  
  
        # adjust denominator to unit infinity norm to match example
        Mnorm = ghinfnorm(M)[1]
        Q_Rf_Rw = Q_Rf_Rw/Mnorm
        Q = Q_Rf_Rw[:,1:p+mu]
        Rf = Q_Rf_Rw[:,p+mu+1:p+mu+mf]
        Rw = Q_Rf_Rw[:,p+mu+mf+1:end]
        println("Q = $(dss2rm(Q, atol = 1.e-7))")
        println("Rf = $(dss2rm(Rf, atol = 1.e-7))")
        println("Rw = $(dss2rm(Rw, atol = 1.e-7))")
        gap = fdhinfminus(Rf)[1]
        println("gap = $gap")
    end
else
    @info "No solution exists"
end
   
end # module
