module Ex5_10c
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test

# Example 5.10c - Solution of an EFDIP
println("Example 5.10c")

# enter output and fault vector dimensions
p = 3; mf = 3;
# generate random dimensions for system order and input vectors
nu = Int(floor(1+4*rand())); mu = Int(floor(1+4*rand()));
nd = Int(floor(1+4*rand())); md = Int(floor(1+4*rand()));

# define random Gu(s) and Gd(s) with triplex sensor redundancy
# and Gf(s) for triplex sensor faults
Gu = ones(3,1)*rss(nu,1,mu); # enter Gu(s) in state-space form
Gd = ones(3,1)*rss(nd,1,md); # enter Gd(s) in state-space form
Gf = eye(3);                 # enter Gf(s) for sensor faults
atol = 1.e-7;                # tolerance for rank tests

# build model with faults
sysf = fdimodset([Gu Gd Gf], c = 1:mu, d = mu.+(1:md), f = (mu+md).+(1:mf))       

SFDI = [ 0 1 1; 1 0 1; 1 1 0] .> 0; # enter structure matrix

# 
Qt, Rft = efdisyn(sysf, SFDI; atol, rdim = 1);

# normalize Q and Rf to match example
scale = sign.([ Rft.sys[1].D[1,2], Rft.sys[2].D[1,3], Rft.sys[3].D[1,1]])
Q = FDIFilter(scale .* Qt.sys, p, mu)
R = FDIFilterIF(scale .* Rft.sys; mf)

# check synthesis conditions: Qt*[Gu Gd;I 0] = 0 and Qt*[Gf; 0] = Rf
Rt = fdIFeval(Qt,sysf) # form Qt*[Gu Gd Gf;I 0 0];
@test iszero(vcat(Rt.sys...)[:,[Rt.controls;Rt.disturbances]],atol=1.e-7) &&
      iszero(vcat(Rt.sys...)[:,Rt.faults]-vcat(Rft.sys...),atol=1.e-7)

# check weak and strong fault detectability
@test fditspec(Rft) == fdisspec(Rft) 

println("Q = ")
display(Q)
println("R = ")
display(R)


end # module
