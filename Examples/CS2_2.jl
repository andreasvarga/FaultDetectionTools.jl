module CS2_2
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test
using JLD2

## CS2_2  - Case-study example: Monitoring air data sensor faults
#           Robust least order LPV synthesis
println("Case study CS2_2")

## Part 1 - Model setup
# load matrices of the aircraft multiple-model SYSACM, actuator model ACT,
# output-feedback gain K, percentages of mass variations massi

cd(joinpath(pkgdir(FaultDetectionTools), "Examples"))
@load "cs2data.jld2"

N = size(SYSACSM,1)
p, m = size(SYSACSM[1].D)
T = Float64

# build minimal realizations of AC-models with actuator faults
# set dimensions
nom = 7;             # nominal system 
# set sensor indices and set dimensions of inputs
#     [ AoA VCAS ]
sen = [  2,    4  ]; 
md = 2; mu = m-md; mf = length(sen); n = maximum(order.(SYSACSM));

idm = eye(T,p);

syssen = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    syssen[i] = [ SYSACSM[i] idm[:,sen]]
end

syssenf = fdimodset(syssen; c = Vector(1:mu), d = mu .+ Vector(1:md), f = (mu+md) .+ (1:mf));


## Part 2 - Setup of the synthesis specifications

# compute the achievable strong specifications for constant faults
S_strong = fdigenspec(syssenf[nom], atol = 1.e-6, atol3 = 1.e-6, FDtol = 1.e-5, FDGainTol = 0.001, FDfreq = 0, sdeg = -0.05);

@test size(S_strong,1) == 2^mf-1

## Part 3 - LPV synthesis using up to 4th order polynomial gain interpolation

# form extended system [ Gu Gd Gf; I 0 0] with outputs extended with u
# and determine constant solutions based on Remark 8.1 of the FDI book
syse = similar(Vector{DescriptorStateSpace},N)
D = zeros(T,mf,p+mu,N)
Dref = [zeros(T,mf,n+mu+md) eye(T,mf)]
for i = 1:N
    syse[i] = [syssenf[i].sys;eye(mu,mu+md+mf)]
    cde = [syse[i].C syse[i].D] 
    D[:,:,i] = Dref/cde
end     

N4 = similar(Matrix{Float64},N,5)
Qpol = Polynomial.(zeros(T,mf,p+mu))
for deg = 0:4
    for i = 1:mf
        for j = 1:p+mu
            Qpol[i,j] = fit(view(massi,:), view(D,i,j,:), deg)
        end
    end
    for i = 1:N
        rezi = (massi[i] .|> Qpol)*syse[i] - [zeros(mf,mu+md) eye(2)]
        N4[i,deg+1]  = glinfnorm(rezi)[1]
    end     
end
@test all(maximum(N4,dims = 1) .> [545.  138.  16.  0.88  0.12])

println(" Table 8.5  Robustness analysis results of interpolation based approximations")
println(" N0        N1         N2         N3         N4")
display(N4)

## Part 4 - LPV synthesis using 3rd order polynomial gain interpolation
#  setup of the LPV synthesis similar to using MATLAB's SYSTUNE (to be implemented)


end  # module

 
