module Ex5_13c
using FaultDetectionTools, DescriptorSystems, Test

# Example 5.13c - Solution of an EMMP 
println("Example 5.13c")

# define s as an improper transfer function
s = rtf('s');
# enter Gu(s), Gf(s) and Mr(s)
Gu = [s/(s^2+3*s+2) 1/(s+2);
     s/(s+1) 0;
      0 1/(s+2)];
Mr = dss(eye(2));                  # enter Mr(s)
p, mu = size(Gu); mf = mu

# build the synthesis model with additive faults 
sysf = fdimodset(dss(Gu), c = 1:mu, f = 1:mu);

# enter reference model
Mr = fdimodset(dss(eye(mf)), f = 1:mf)

atol = 1.e-7                # tolerance for rank tests
sdeg = -1                   # set stability degree

# solve an exact model-matching problem using EMMSYN
Q, R, info = emmsyn(sysf, Mr; atol, sdeg, minimal = false); 

# display results
println("Q = $(dss2rm(Q.sys; atol))")
println("M = $(dss2rm(info.M; atol))")

# check synthesis conditions: Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = M*Mr
Rt = fdIFeval(Q, sysf; atol) # form Q*[Gu Gd Gf;I 0 0];
@test iszero(Rt.sys[:,[Rt.controls;Rt.disturbances]]; atol) &&
      iszero(Rt.sys[:,Rt.faults]-info.M*Mr.sys; atol)

end # module
using Main.Ex5_13c
