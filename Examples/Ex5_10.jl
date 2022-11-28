module Ex5_10
using FaultDetectionTools, DescriptorSystems, LinearAlgebra, Test

# Example 5.10 - Solution of an EFDIP
println("Example 5.10")

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

# build model with faults: sysf = [Gu Gd Gf]
sysf = fdimodset([Gu Gd Gf], c = 1:mu, d = mu.+(1:md), f = (mu+md).+(1:mf))       

S = [ 0 1 1; 1 0 1; 1 1 0] .> 0;  # enter structure matrix

## Step 1) of Procedure EFDI
# compute Q1, the left nullspace of [Gu Gd;I 0], and  Rf1 = Q1*[Gf;0]
# QRf contains Q1 in QRf[:,1:p+mu] and Rf1 in QRf[:,p+mu+1:p+mu+mf]
QRf, info = glnull([sysf.sys; eye(mu,mu+md+mf)],mf; atol) 

## Step 2): of Procedure EFDI 
# initialize overall filter Q and fault detection system Rf
nb = size(S,1);      # number of filters 
Qt = similar(Vector{typeof(QRf)},nb)
Rft = similar(Vector{typeof(QRf)},nb)
for i = 1:nb
    # Step 2.1): define  Rf11 = Rf1[:,indd], Rf12 = Rf1[:,indf]
    indd = Vector(1:mf)[S[i,:] .== false] 
    indf = Vector(1:mf)[S[i,:] .== true] 
    # pack [ Rf11 Rf12 [Q1 Rf1] ]
    sysc = fdimodset(QRf, d = (p+mu) .+ indd, f = (p+mu) .+ indf, aux = Vector(1:p+mu+mf))
    # Step 2.2): apply  Procedure EFD to {0,Rf11,Rf12,[Q1 Rf1]}
    # to determine a least order Q1i such that Q1i*Rf11 = 0
    # QRfauxi contains: [Q1i*Rf12 Q1i*Q1 Q1i*Rf1]
    _, QRfauxi, infoi = efdsyn(sysc; rdim = 1, atol) 
    # extract [Q1i*Q1 Q1i*Rf1]
    QRfi = QRfauxi.sys[:,QRfauxi.aux]
    # Step 2.3): extract Q[i] = Q1i*Q1 and Rf[i] = Q1i*Rf1
    Qt[i] = QRfi[:,1:p+mu]
    Rft[i] = QRfi[:,(p+mu).+(1:mf)]
end

# normalize Q and Rf to match example
scale = sign.([ Rft[1].D[1,2], Rft[2].D[1,3], Rft[3].D[1,1]])
Q = FDIFilter(scale .* Qt, p, mu)
R = FDIFilterIF(scale .* Rft; mf)

# check all gaps are infinite
@test all(fdif2ngap(R,S)[1] .== Inf)
end # module
