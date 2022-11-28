module Ex5_11c
using FaultDetectionTools, DescriptorSystems, Test

# Example 5.11c - Solution of an AFDIP 
println("Example 5.11c")

# define s as an improper transfer function
s = rtf('s');
# define Gu(s), Gf(s) and Gw(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     # enter Gu(s)
Gw = [1/(s+2); 0];                   # enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           # enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       # enter dimensions
atol = 1.e-7;                        # tolerance for rank tests

SFDI = eye(mf) .> 0                  # set structure matrix to I

# build the synthesis model with additive faults 
sysf = fdimodset(dss([Gu Gf Gw]), c = 1:mu, 
                 f = mu .+ (1:mf), n = (mu+mf) .+ (1:mw));

# call of AFDISYN with the option for stability degree -3 
Q, R, info = afdisyn(sysf, SFDI; smarg =-3, sdeg = -3)

# scale with -1 to match example
Q.sys[1] = -Q.sys[1]; R.sys[1] = -R.sys[1];

for i = 1:size(SFDI,1)
    println("Q[$i] = $(dss2rm(Q.sys[i]; atol))")
    println("Rf[$i] = $(dss2rm(R.sys[i][:,R.faults]; atol))")
    println("Rw[$i] = $(dss2rm(R.sys[i][:,R.noise]; atol))")
end

println("Q = "); display(Q)
println("R = "); display(R)
println("gap = $(fdif2ngap(R,SFDI)[1])")
   
# check synthesis conditions: 
# Q*[Gu Gd;I 0] = 0 and Q*[Gf; 0] = Rf
Rt = fdIFeval(Q,sysf) # form Q*[Gu Gd Gf;I 0 0];
indcd = [Rt.controls;Rt.disturbances]; 
indfn = [Rt.faults;Rt.noise]
@test iszero(vcat(Rt.sys...)[:,indcd]; atol) &&
      iszero(vcat(Rt.sys...)[:,indfn]-vcat(R.sys...); atol) &&
      fdif2ngap(R,SFDI)[1] ≈ info.gap ≈ [3,Inf]
end # module
