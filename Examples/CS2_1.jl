module CS2_1
using FaultDetectionTools, DescriptorSystems, LinearAlgebra, 
      JLD2, Optim, Test
#using Evolutionary

## CS2_1  - Case-study example: Monitoring air data sensor faults
#           Robust least order LTI synthesis
println("Case study CS2_1")

## Part 1 - Model setup
# load matrices of the aircraft multiple-model SYSACSM and 
# percentages of mass variations massi

cd(joinpath(pkgdir(FaultDetectionTools), "Examples"))
#@load "cs2data.jld2"
SYSACSM, massi = load("cs2data.jld2", "SYSACSM","massi")

# set dimensions
N = size(SYSACSM,1)        # number of models
p, m = size(SYSACSM[1].D)  # output and input dimensions
n = order(SYSACSM[1])      # state dimension
atol = 1.e-6               # absolute tolerance
nom = 7;                   # index of nominal system 
T = Float64

# set sensor indices and set system dimensions of inputs
#     [ AoA VCAS ]
sen = [  2,    4  ]; 
md = 2; mu = m-md; mf = length(sen); 

# build minimal realizations of AC-models with sensor faults
idm = eye(p);
syssen = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    syssen[i] = [ SYSACSM[i] idm[:,sen]]
end

# build synthesis models with sensor faults
syssenf = fdimodset(syssen; c = Vector(1:mu), 
                    d = mu .+ Vector(1:md), f = (mu+md) .+ (1:mf));

## Part 2 - Setup of the synthesis specifications

# compute achievable strong specifications for constant faults
FDtol = 1.e-5      # weak fault detection threshold
FDGainTol = 0.0001  # strong fault detection threshold
FDfreq = 0         # frequency for strong detectability checks 
S_strong = fdigenspec(syssenf[nom]; atol, FDtol, FDGainTol, 
                      FDfreq, sdeg = -0.05);

# check number of maximum possible specifications
@test size(S_strong,1) == 2^mf-1

# define SFDI, the signatures for isolation of single faults
SFDI = eye(mf) .> 0;

## Part 3 - Multiple filter synthesis using Procedure EFDI 

Q = similar(Vector{DescriptorStateSpace},N)
for i = 1:N
    Qi = efdisyn(syssenf[i], SFDI;   rdim = 1, FDfreq, FDGainTol,
                 atol, sdeg = -5, smarg = -0.05)[1] 
    Q[i] = vcat(Qi.sys...)
end
@test all(order.(Q) .== 0)

##  Part 4 - Assesment of synthesis results

# evaluate R[i] := Q[i]*[ Gu[i] Gd[i] Gf[i]; I 0 0] 
# and determine ||R[i]-[0 0 I]||∞ for i = 1, ..., N
syse = similar(Vector{DescriptorStateSpace},N)
Nref = similar(Vector{T},N)
y = []
for i = 1:N
    syse[i] = [syssenf[i].sys;eye(mu,mu+md+mf)]
    Ri = Q[i]*syse[i]
    rezi = Ri - [zeros(mf,mu+md) eye(2)]
    Nref[i]  = glinfnorm(rezi)[1]
    push!(y,stepresp(Ri,10)[1])
end  

tout = Vector(0:0.1:10.)
include("Fig8_5.jl")
Fig8_5 = f

##  Part 5 - Assesment of synthesis results for nominal synthesis

# evaluate R[i] := Q[nom]*[ Gu[i] Gd[i] Gf[i]; I 0 0] 
# and determine ||R[i]-[0 0 I]||∞ for i = 1, ..., N

Nnom = similar(Vector{T},N)
y = []
for i = 1:N
    Ri = Q[nom]*syse[i]
    rezi = Ri - [zeros(mf,mu+md) eye(2)]
    Nnom[i]  = glinfnorm(rezi)[1]
    push!(y,stepresp(Ri,10)[1])
end     

tout = Vector(0:0.1:10.)
include("Fig8_6.jl")
Fig8_6 = f


##  Part 6 - Assesment of synthesis results for mean value gains

# compute mean values of gains
D0 = similar(Matrix{T},mf,p+mu)
for i = 1:mf
    for j = 1:p+mu
        D0[i,j] = sum([Q[k].D[i,j] for k in 1:N])/N
    end
end

# evaluate R[i] := D0*[ Gu[i] Gd[i] Gf[i]; I 0 0] 
# and determine N0[i] = ||R[i]-[0 0 I]||∞ for i = 1, ..., N
N0 = similar(Vector{T},N)
for i = 1:N
    rezi = D0*syse[i] - [zeros(mf,mu+md) eye(2)]
    N0[i]  = glinfnorm(rezi)[1]
end     

## Part 7 - Multiple model synthesis of a constant gain 

# optimal solution computed with MATLAB's SYSTUNE
# Dopt0 = [    
# -0.7998    1.0000   -0.0315   -0.0000    0.0069   -0.2488    4.1934    0.3858   -0.1468   -0.0773    0.0347
#  0.0567         0   -0.0207    1.0000   -1.4180   -0.0067    4.5920    0.2673   -0.0255   -0.1018   -0.0228
# ]
# optimal solution computed with MATLAB's SYSTUNE 
# rounded to two decimal figures
Dopt0 = [    
-0.79 1 -0.03 0  0.01 -0.25 4.2 0.39 -0.15 -0.08  0.03
 0.06 0 -0.02 1 -1.42 -0.01 4.6 0.27 -0.03 -0.10 -0.02
]

# define parameterized constant filter gain: those set to false are fixed
mask = trues(mf,mu+p)
mask[1,2] = false
mask[1,4] = false
mask[2,2] = false
mask[2,4] = false

function maxerror(x)
    err = 0.0
    D = similar(Matrix{T},mf,mu+p)
    D[mask] = x
    D[.!mask] = Dres
    for i in 1:N
        rezi = D*syse[i] - [zeros(mf,mu+md) eye(2)]  
        err = max(err, glinfnorm(rezi)[1])
    end
    return err
end

# initialize with the optimal solution computed with the MATLAB's SYSTUNE
# solution rounded to two exact decimal digits
Dinit = Dopt0  
Dres = Dinit[.!mask]
res = optimize(maxerror, Dinit[mask], NelderMead(), Optim.Options(iterations = 20000))

println("Ninit = $(maxerror(Dinit[mask]))")
println("Nfinal = $(maxerror(Optim.minimizer(res)))")

# comment out the next text for a more realistic try
# Dinit = Dopt0 .* (1 .+ .1*randn(mf,mu+p))
# Dres = Dopt0[.!mask]

# Dinit = D0 
# Dres = D0[.!mask]

# res = optimize(maxerror, Dinit[mask], NelderMead(), Optim.Options(iterations = 20000))

# println("Ninit = $(maxerror(Dinit[mask]))")
# println("Nfinal = $(maxerror(Optim.minimizer(res)))")


Dopt = similar(Matrix{T},mf,mu+p)
Dopt[mask] = Optim.minimizer(res)
Dopt[.!mask] = Dres

Nopt = similar(Vector{T},N)
y = []
for i = 1:N
    Ri = Dopt*syse[i]
    rezi = Ri - [zeros(mf,mu+md) eye(2)]
    Nopt[i]  = glinfnorm(rezi)[1]
    push!(y,stepresp(Ri,10)[1])
end     

norms = [Nref        Nnom  N0    Nopt]
println(" Table 8.4  Robustness analysis results for constant approximations")
println(" Nref         Nnom          N0        Nopt")
display(norms)
@test norm(Nref,Inf) < 1.e-7 &&  norm(Nnom,Inf) > 643 &&   norm(N0,Inf) > 504 &&   norm(Nopt,Inf) > 0.76

tout = Vector(0:0.1:10.)
include("Fig8_7.jl")
Fig8_7 = f

export Fig8_5, Fig8_6, Fig8_7

end  # module
using Main.CS2_1

 
