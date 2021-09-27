module Ex5_2
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.2 - EFDP has no solution 
println("Example 5.2")
# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gd(s), Gf(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gd = [1/(s+2); 0];                   # enter Gd(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; md = 1; mf = 2;       # set dimensions

# build model with additive faults 
#sysf = fdimodset(ss([Gu Gd Gf]),struct('c',1:mu,'d',mu+(1:md),'f',mu+md+(1:mf)));
sysf = fdimodset(dss([Gu Gd Gf]),c = 1:mu, d = mu.+(1:md), f = (mu+md).+(1:mf))

# call of EFDSYN with default options 
try
  Q, Rf, info = efdsyn(sysf)
  @test false
catch err
  @test true
  display(err)
  # compute achievable structure matrix
  println("S_achievable = $(fdigenspec(sysf))")
end

end # module
