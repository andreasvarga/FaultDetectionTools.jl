module CS2_1
using FaultDetectionTools
using DescriptorSystems
using LinearAlgebra
using Polynomials
using Test
using Optim
#using Evolutionary
using JLD2

## CS2_1  - Case-study example: Monitoring air data sensor faults
#           Robust least order LTI synthesis
println("Case study CS2_1")

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

idm = eye(p);

syssen = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    syssen[i] = [ SYSACSM[i] idm[:,sen]]
end

syssenf = fdimodset(syssen; c = Vector(1:mu), d = mu .+ Vector(1:md), f = (mu+md) .+ (1:mf));

## Part 2 - Setup of the synthesis specifications

# compute the achievable strong specifications for constant faults
S_strong = fdigenspec(syssenf[nom], atol = 1.e-6, atol3 = 1.e-6, FDtol = 1.e-5, FDGainTol = 0.001, FDfreq = 0, sdeg = -0.05);

@test size(S_strong,1) == 2^mf-1

# define SFDI, the signatures for isolation of single faults
SFDI = eye(mf) .> 0;

## Part 3 - Multiple filter synthesis using Procedure EFDI 

Q = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    Qi = efdisyn( syssenf[i], SFDI;  atol = 1.e-6, sdeg = -5, smarg = -0.05, simple = false, FDfreq = 0, FDGainTol = 0.0001, rdim = 1)[1];  
    Q[i] = vcat(Qi.sys...)
end
@test order.(Q) == zeros(Int,N)

##  Part 4 - Assesment of synthesis results
# form extended system [ Gu Gd Gf; I 0 0] with outputs extended with u
syse = similar(Vector{DescriptorStateSpace},N)
Nref = similar(Vector{Float64},N)
for i = 1:N
    syse[i] = [syssenf[i].sys;eye(mu,mu+md+mf)]
    rezi = Q[i]*syse[i] - [zeros(mf,mu+md) eye(2)]
    Nref[i]  = glinfnorm(rezi)[1]
end     


##  Part 5 - Assesment of synthesis results for nominal synthesis
# form extended system [ Gu Gd Gf; I 0 0] with outputs extended with u
Nnom = similar(Vector{Float64},N)
for i = 1:N
    rezi = Q[nom]*syse[i] - [zeros(mf,mu+md) eye(2)]
    Nnom[i]  = glinfnorm(rezi)[1]
end     


##  Part 6 - Assesment of synthesis results for constant interpolated gain
D0 = zeros(mf,p+mu)
D = zeros(T,mf,p+mu,N)
for i = 1:N
    D[:,:,i] = Q[i].D
end     

for i = 1:mf
    for j = 1:p+mu
        D0[i,j] = fit(view(massi,:), view(D,i,j,:), 0)(0)
    end
end

N0 = similar(Vector{Float64},N)
for i = 1:N
    rezi = D0*syse[i] - [zeros(mf,mu+md) eye(2)]
    N0[i]  = glinfnorm(rezi)[1]
end     

## Part 7 - Multiple model synthesis of a constant gain 
# to be implemented using optimization-based tools

# optimal solution computed with MATLAB's SYSTUNE
Dopt0 = [    
    -0.7998    1.0000   -0.0315   -0.0000    0.0069   -0.2488    4.1934    0.3858   -0.1468   -0.0773    0.0347
    0.0567         0   -0.0207    1.0000   -1.4180   -0.0067    4.5920    0.2673   -0.0255   -0.1018   -0.0228
]


# define parameterized constant filter gain: those set to false are fixed
mask = trues(mf,mu+p)
mask[1,2] = false
mask[1,4] = false
mask[2,2] = false
mask[2,4] = false
#Dinit = Dopt0  # this is the optimal solution 
#Dinit = D0  # this  the normal initalization, unfortunately not converging
#Dres = Dinit[.!mask]
# use a 10% perturbed initialization
#Dinit = Dopt0 .* (1 .+ .1*randn(mf,mu+p))

function maxerror(x)
    err = 0.0
    D = similar(Matrix{Float64},mf,mu+p)
    D[mask] = x
    D[.!mask] = Dres
    for i in 1:N
        rezi = D*syse[i] - [zeros(mf,mu+md) eye(2)]  
        err = max(err, glinfnorm(rezi)[1])
        #err += glinfnorm(rezi)[1]^2
        #err += norm(dcgain(rezi))^2
    end
    return err
end

# initialize with the optimal solution computed with MATLAB's SYSTUNE
Dinit = Dopt0  # this is the optimal solution 
Dres = Dinit[.!mask]

println("Ninit = $(maxerror(Dinit[mask]))")
println("Please wait a few seconds - this optimization will finish soon")

res = optimize(maxerror, Dinit[mask], NelderMead(), Optim.Options(iterations = 20000))

println("Nfinal = $(maxerror(Optim.minimizer(res)))")

# comment out the next text for a more realistic try
# Dinit = Dopt0 .* (1 .+ .1*randn(mf,mu+p))
# Dres = Dopt0[.!mask]
# println("Ninit = $(maxerror(Dinit[mask]))")
# println("Please wait - this optimization take some time")

# res = optimize(maxerror, Dinit[mask], NelderMead(), Optim.Options(iterations = 20000))

# println("Nfinal = $(maxerror(Optim.minimizer(res)))")


Dopt = similar(Matrix{Float64},mf,mu+p)
Dopt[mask] = Optim.minimizer(res)
Dopt[.!mask] = Dres

Nopt = similar(Vector{Float64},N)
for i = 1:N
    rezi = Dopt*syse[i] - [zeros(mf,mu+md) eye(2)]
    Nopt[i]  = glinfnorm(rezi)[1]
end     

norms = [Nref        Nnom  N0    Nopt]
println(" Table 8.4  Robustness analysis results for constant approximations")
println(" Nref         Nnom          N0        Nopt")
display(norms)
@test norm(Nref,Inf) < 1.e-7 &&  norm(Nnom,Inf) > 643 &&   norm(N0,Inf) > 504 &&   norm(Nopt,Inf) > 0.76


end  # module

 
