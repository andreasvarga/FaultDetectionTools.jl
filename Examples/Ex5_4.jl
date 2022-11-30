module Ex5_4
using DescriptorSystems, Test

# Example 5.4 - Solution of an EFDP
println("Example 5.4")

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gd(s), Gf(s)
Gu = [(s+1)/(s-2); (s+2)/(s-3)];     # enter Gu(s)
Gd = [(s-1)/(s+2); 0];               # enter Gd(s)
Gf = [(s+1)/(s-2) 0; (s+2)/(s-3) 1]; # enter Gf(s)
p = 2; mu = 1; md = 1; mf = 2;       # set dimensions

# compute a left nullspace basis Q of [Gu Gd; I 0]
Q1 = glnull(dss([Gu Gd;eye(mu,mu+md)]))[1];

# compute Rf1 = Q1[Gf;0]
Rf1 = gir(Q1*dss([Gf;zeros(mu,mf)]));

# check solvability using a random frequency
if minimum(abs.(evalfr(Rf1,rand()))) > 0.01
   # compute a stable left coprime factorization 
   # [Q1 Rf1]=inv(Q3)*[Q,Rf]
   # enforce stability degree -3
   Q_Rf, Q3 = glcf([Q1 Rf1];sdeg = -3);
   # extract Q and Rf
   Q = Q_Rf[:,1:p+mu]; Rf = Q_Rf[:,p+mu+1:end]; 
   scale = evalfr(Rf[1,1],Inf)[1,1]
   Q = Q/scale; Rf = Rf/scale;
   @test gpole(Q) ≈ [-3] && gpole(Rf) ≈ [-3] && 
         all(abs.(evalfr(Rf,rand())) .> 0) && 
         iszero(Rf-Q*dss([Gf;zeros(mu,mf)]),atol=1.e-7) && 
         iszero(Q*dss([Gu Gd;eye(mu,mu+md)]),atol=1.e-7)
   # normalize Q and Rf to match example
   println(" Q = $(dss2rm(Q,atol=1.e-7))")
   println(" Rf = $(dss2rm(Rf,atol=1.e-7))")
else
   @info "No solution exists"
end

end # module
using Main.Ex5_4
