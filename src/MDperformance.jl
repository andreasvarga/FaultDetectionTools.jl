"""
     S = mdspec(sysR::MDFilterIF; cdinp = false, atol, atol1 = atol, atol2 = atol, rtol = 0)

Compute the weak binary structure matrix `S` 
of a collection of model detection filters
using the model detection filter internal form object `sysR::MDFilterIF`.

For an `M × N` array of filters `sysR`, `S` is an `M × N` binary array determined as follows. 
Let the `(i,j)`-th component filter `sysR.sys[i,j]` have the input-output form
 
     rij = Ruij(λ)*u + Rdij(λ)*dj + Rwij(λ)*wj + Rvij(λ)*vj ,

with the Laplace- or Z-transformed residual output `rij`, control inputs `u`, 
disturbance inputs `dj`, noise inputs `wj`, and auxiliary inputs `vj`,  
and with `Ruij(λ)`, `Rdij(λ)`, `Rwij(λ)` and `Rvij(λ)`, the corresponding transfer function matrices. 
Then, `S[i,j] = 1` if `Ruij(λ)` is nonzero for `cdinp = false` (default), or
if `[Ruij(λs) Rdij(λs)]` is nonzero for `cdinp = true`. Otherwise, `S[i,j] = 0`.

If `(Arij-λErij,Brij,Crij,Drij)` is the descriptor realization of `sysR.sys[i,j]`, then 
the keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Arij`, `Brij`, `Crij`, `Drij`, the absolute tolerance for the nonzero elements of `Erij`,  
and the relative tolerance for the nonzero elements of `Arij`, `Brij`, `Crij`, `Drij` and `Eirj`.
The default relative tolerance is `nij*ϵ`, where `ϵ` is the working _machine epsilon_ 
and `nij` is the order of the system matrix of `sysR.sys[i,j]`.  
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 
"""    
function mdspec(sysR::MDFilterIF{T};  cdinp::Bool = false,
                  atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, 
                  rtol::Real = 0) where T
    M, N = size(sysR.sys)
    S = trues(M,N)
    mu = sysR.mu
    cdinp && (md = sysR.md)
    for j = 1:N
        inpj = cdinp ? (1:mu+md[j]) : (1:mu)
        for i = 1:M
            S[i,j] = !iszero(sysR.sys[i,j][:,inpj]; atol1, atol2, rtol)
        end         
    end
    return S
end
"""
     S = mdsspec(sysR::MDFilterIF, freq; cdinp = false, MDGainTol = 0.01, 
                 atol, atol1 = atol, atol2 = atol, rtol = 0, fast = true)

Compute, for a given set of real frequencies `freq`, the strong binary structure matrix `S` 
of a collection of model detection filters
using the model detection filter internal form object `sysR::MDFilterIF`.

`freq` must contain a real frequency value or a vector of `nf` real frequencies 
which characterize the classes of persistent control and disturbance signals 
(default: `freq = 0`, i.e., characterizing constant signals) and defines 
the set `Ω` of complex frequencies which characterize the classes of persistent signals
as follows: 
if `f` is a real frequency in `freq`, then the corresponding complex frequency `λ` in `Ω` 
is `λ := im*f`, for a continuous-time system,
or `λ := exp(im*f*abs(Ts))`, for a discrete-time system with sampling-time `Ts`. 

For an `M × N` array of filters `sysR`, `S` is an `M × N` binary array determined as follows. 
Let the `(i,j)`-th component filter `sysR.sys[i,j]` have the input-output form
 
     rij = Ruij(λ)*u + Rdij(λ)*dj + Rwij(λ)*wj + Rvij(λ)*vj ,

with the Laplace- or Z-transformed residual output `rij`, control inputs `u`, 
disturbance inputs `dj`, noise inputs `wj`, and auxiliary inputs `vj`,  
and with `Ruij(λ)`, `Rdij(λ)`, `Rwij(λ)` and `Rvij(λ)`, the corresponding transfer function matrices. 
Then, `S[i,j] = 1` if `Ruij(λs)` is nonzero for any `λs ∈ Ω` and `cdinp = false` (default), or
if `[Ruij(λs) Rdij(λs)]` is nonzero for any `λs ∈ Ω` and `cdinp = true`. Otherwise, `S[i,j] = 0`.

`MDGainTol = tol` specifies an absolute  threshold `tol` for the nonzero magnitudes of 
the frequency response gains (default: `tol = 0.01`). 

The computation of minimal realizations of individual input-output channels relies on pencil manipulation algorithms,
which employ rank determinations based on either the use of 
rank revealing QR-decomposition with column pivoting, if `fast = true` (default) or 
the more reliable SVD-decompositions if `fast = false`.

If `(Arij-λErij,Brij,Crij,Drij)` is the descriptor realization of `sysR.sys[i,j]`, then 
the keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Arij`, `Brij`, `Crij`, `Drij`, the absolute tolerance for the nonzero elements of `Erij`,  
and the relative tolerance for the nonzero elements of `Arij`, `Brij`, `Crij`, `Drij` and `Eirj`.
The default relative tolerance is `nij*ϵ`, where `ϵ` is the working _machine epsilon_ 
and `nij` is the order of the system `sysR.sys[i,j]`.  
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 
"""    
function mdsspec(sysR::MDFilterIF{T}, freq::Union{AbstractVector{<:Real},Real} = 0;  cdinp::Bool = false,
                  MDGainTol::Real = 0.01, atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, 
                  rtol::Real = 0, fast::Bool = true) where T
    isa(freq,Vector) ? tfreq = freq : tfreq = [freq] 
    M, N = size(sysR.sys)
    S = trues(M,N)
    mu = sysR.mu
    cdinp && (md = sysR.md)
    for j = 1:N
        inpj = cdinp ? (1:mu+md[j]) : (1:mu)
        for i = 1:M
            S[i,j] = any(fdisspec_(sysR.sys[i,j][:,inpj], tfreq; FDGainTol = MDGainTol, block = true, atol1, atol2, rtol, fast)[1])
        end         
    end
    return S
end
"""
     mdperf(sysR::MDFilterIF; MDfreq, cdinp = false, rtolinf = 0.00001, 
            offset, atol, atol1 = atol, atol2 = atol, rtol = 0, fast = true) -> (mdgain, fpeak)

Compute the distance-mapping performance of a collection of model detection filters
using the model detection filter internal form object `sysR::MDFilterIF`.  
For an `M × N` array of filters `sysR`, the `M × N` arrays of model detection performance gains 
`mdgain` and the corresponding peak frequencies `fpeak` are determined as follows. 
Let the `(i,j)`-th component filter `sysR.sys[i,j]` have the input-output form
 
     rij = Ruij(λ)*u + Rdij(λ)*dj + Rwij(λ)*wj + Rvij(λ)*vj ,

with the Laplace- or Z-transformed residual output `rij`, control inputs `u`, 
disturbance inputs `dj`, noise inputs `wj`, and auxiliary inputs `vj`,  
and with `Ruij(λ)`, `Rdij(λ)`, `Rwij(λ)` and `Rvij(λ)`, the corresponding transfer function matrices. 
Then, the `(i,j)`-th performance gain is evaluated as 
`mdgain[i,j] = ||Ruij(λ)||∞` if `cdinp = false` (default) or `mdgain[i,j] = ||[Ruij(λ) Rdij(λ)]||∞` if `cdinp = true` 
and `fpeak[i,j]` contains the corresponding peak frequency. 

If `MDfreq = ω`, where `ω` is a given vector of real frequency values, then each gain `mdgain[i,j]` represents
the maximum of 2-norm pointwise gains evaluated for all frequencies in `ω` and 
`fpeak[i,j]` is the corresponding peak frequency.

The stability boundary offset, `β`, to be used to assess the finite poles which belong to the
boundary of the stability domain can be specified via the keyword parameter `offset = β`.
Accordingly, for a continuous-time system, these are the finite poles having 
real parts within the interval `[-β,β]`, while for a discrete-time system, 
these are the finite pole having moduli within the interval `[1-β,1+β]`. 
The default value used for `β` is `sqrt(ϵ)`, where `ϵ` is the working machine precision. 

Pencil reduction algorithms are employed to compute the H∞-norms. 
These algorithms perform rank decisions based on rank 
revealing QR-decompositions with column pivoting 
if `fast = true` (default) or the more reliable SVD-decompositions if `fast = false`.

If `(Arij-λErij,Brij,Crij,Drij)` is the descriptor realization of `sysR.sys[i,j]`, then 
the keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Arij`, `Brij`, `Crij`, `Drij`, the absolute tolerance for the nonzero elements of `Erij`,  
and the relative tolerance for the nonzero elements of `Arij`, `Brij`, `Crij`, `Drij` and `Eirj`.
The default relative tolerance is `nij*ϵ`, where `ϵ` is the working _machine epsilon_ 
and `nij` is the order of the system `sysR.sys[i,j]`.  
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 

The keyword argument `rtolinf = tol` specifies the relative accuracy `tol` to be used 
to compute the infinity norms. The default value used is `tol = 0.00001`.
"""    
function mdperf(sysR::MDFilterIF{T}; MDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                cdinp::Bool = false, offset::Real = sqrt(eps(T)), rtolinf::Real = 0.00001, atol::Real = zero(float(real(T))), 
                atol1::Real = atol, atol2::Real = atol, rtol::Real = 0, fast::Bool = true) where T
    if ismissing(MDfreq)   
       lfreq = 0
    else         
       isa(MDfreq,Vector) || (MDfreq = [MDfreq]) 
       lfreq = length(MDfreq);
    end
    M, N = size(sysR.sys)
    mdgain = similar(Matrix{T}, M, N)
    fpeak = similar(Matrix{T}, M, N)
    inpud = 1:sysR.mu
    for j = 1:N
        cdinp && (inpud = 1:sysR.mu+sysR.md[j]) 
        for i = 1:M
            if lfreq == 0
                mdgain[i,j], fpeak[i,j] = ghinfnorm(sysR.sys[i,j][:,inpud]; rtolinf, atol1, atol2, rtol, fast, offset)
            else
                gmax = zero(T)
                fmax = MDfreq[1]
                for k = 1:lfreq
                    tval = opnorm(evalfr(sysR.sys[i,j][:,inpud]; fval = MDfreq[k], atol1, atol2, rtol, fast))
                    gmax > tval || (gmax = tval; fmax = MDfreq[k])
                end
                mdgain[i,j] = gmax
                fpeak[i,j] = fmax
            end
        end
    end
    return mdgain, fpeak
end
"""
    mdmatch(sysQ::MDFilter, sysc::MDModel; MDfreq, minimal = false, rtolinf, offset, atol, atol1 = atol, atol2 = atol, rtol, fast = true) -> (mdgain,fpeak,mind)

Compute the distance-mapping performance vector `mdgain` achieved using the model detection filter object `sysQ::MDFilter`
applied to a component model `sysc::MDModel`, the corresponding vector of peak frequencies `fpeak`,
and the index `mind` of the component of `mdgain` for which the minimum gain value is achieved. 

If the `i`-th filter `sysQ.sys[i]` has the transfer function matrix `Qi(λ)` and 
the component model `sysc::MDModel` has the partitioned transfer function matrix 
`G(λ) = [Gu(λ)  Gd(λ) Gw(λ) Gv(λ)]` in accordance with
the partitioned system inputs as `controls`, `disturbances`, `noise` and `auxiliary` inputs, respectively,
then the distance-mapping performance of the `i`-th filter applied to the given component model
is computed as  `mdgain[i] = || Ri(λ) ||∞`, where `Ri(λ)` 
is the corresponding internal form

     Ri(λ) = Qi(λ) * | Gu(λ)  Gd(λ) Gw(λ) Gv(λ) | .
                     |  I     0     0     0     |

Minimal descriptor realizations are computed for `Ri(λ)` if `minimal = true` and a (possibly) non-minimal 
realization is determined if `minimal = false` (default). 

The computation of minimal realizations of individual input-output channels relies on pencil manipulation algorithms,
which employ rank determinations based on either the use of 
rank revealing QR-decomposition with column pivoting, if `fast = true` (default) or 
the more reliable SVD-decompositions if `fast = false`.

If `(Ari-λEri,Bri,Cri,Dri)` is the full order descriptor realization of `Ri(λ)`, then 
the keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Ari`, `Bri`, `Cri`, `Dri`, the absolute tolerance for the nonzero elements of `Eri`,  
and the relative tolerance for the nonzero elements of `Ari`, `Bri`, `Cri`, `Dri` and `Eir`.
The default relative tolerance is `ni*ϵ`, where `ϵ` is the working _machine epsilon_ 
and `ni` is the order of the realitation of `Ri(λ)`.  
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 

"""
function mdmatch(Q::MDFilter{T}, sysc::MDModel; minimal::Bool = false, MDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                 offset::Real = sqrt(eps(T)), rtolinf::Real = 0.00001, 
                 atol::Real = 0., atol1::Real = atol, atol2::Real = atol, rtol::Real = 0., fast::Bool = true) where T
   # compute the internal form as
   #   R[i,j] = Q[i] * [ Gu[j] Gd[j] Gw[j] Gaux[j]]
   #                   [ I     0     0     0      ]
   if ismissing(MDfreq)   
      lfreq = 0
   else         
      isa(MDfreq,Vector) || (MDfreq = [MDfreq]) 
      lfreq = length(MDfreq);
   end
   N = length(Q.sys)
   p, m = size(sysc.sys)
   mu = sysc.mu
   (Q.ny == p && Q.mu == mu) || error("filter Q is incompatible with the given component model")
   syse = [sysc.sys;  I zeros(mu,m-mu)]
   mdgain = similar(Vector{T}, N)
   fpeak = similar(Vector{T}, N)
   N == 0 && (return mdgain, fpeak, nothing)
   tmin = Inf
   mind = 0
   for i = 1:N
       sysRi = minimal ? gminreal(Q.sys[i]*syse; atol1, atol2, rtol, fast) : Q.sys[i]*syse
       if lfreq == 0
          mdgain[i], fpeak[i] = ghinfnorm(sysRi; rtolinf, atol1, atol2, rtol, fast, offset)
       else
          gmax = zero(T)
          fmax = MDfreq[1]
          for k = 1:lfreq
              tval = opnorm(evalfr(sysRi; fval = MDfreq[k], atol1, atol2, rtol, fast))
              gmax > tval || (gmax = tval; fmax = MDfreq[k])
          end
          mdgain[i] = gmax
          fpeak[i] = fmax
       end
       tmin <= mdgain[i] || (tmin = mdgain[i]; mind = i)
   end
   return mdgain, fpeak, mind
end
"""
     mdgap(sysR::MDFilterIF; MDfreq, cdinp = false, rtolinf = 0.00001, 
            offset, atol, atol1 = atol, atol2 = atol, rtol = 0, fast = true) -> (gap, β, γ)

Compute the noise gaps performance of a collection of model detection filters
using the model detection filter internal form object `sysR::MDFilterIF`.  
For an `M × N` array of filters `sysR`, let the `(i,j)`-th component filter 
`sysR.sys[i,j]` have the input-output form
 
     rij = Ruij(λ)*u + Rdij(λ)*dj + Rwij(λ)*wj + Rvij(λ)*vj ,

with the Laplace- or Z-transformed residual output `rij`, control inputs `u`, 
disturbance inputs `dj`, noise inputs `wj`, and auxiliary inputs `vj`,  
and with `Ruij(λ)`, `Rdij(λ)`, `Rwij(λ)` and `Rvij(λ)`, the corresponding transfer function matrices. 
Then, `gap`, `β` and `γ` are `M`-dimensional vectors, such that the `i`-th noise gap is evaluated as 
`gap[i] = β[i]/γ[i]`, where `β[i] = min(||Rij(λ)||∞` for `i` ``\\neq`` `j`) 
and `γ[i] = ||Rwii(λ)||∞`.  `Rij(λ)` is defined as `Rij(λ) = Ruij(λ)`  if `cdinp = false` (default) 
or `Rij(λ) = [Ruij(λ) Rdij(λ)]` if `cdinp = true`.  

If `MDfreq = ω`, where `ω` is a given vector of real frequency values, 
then each gain `β[i]` represents the minimum of 
the maximum of 2-norm pointwise gains evaluated for all frequencies in `ω`.

The stability boundary offset, `β`, to be used to assess the finite poles which belong to the
boundary of the stability domain can be specified via the keyword parameter `offset = β`.
Accordingly, for a continuous-time system, these are the finite poles having 
real parts within the interval `[-β,β]`, while for a discrete-time system, 
these are the finite pole having moduli within the interval `[1-β,1+β]`. 
The default value used for `β` is `sqrt(ϵ)`, where `ϵ` is the working machine precision. 

Pencil reduction algorithms are employed to compute the H∞-norms. 
These algorithms perform rank decisions based on rank 
revealing QR-decompositions with column pivoting 
if `fast = true` (default) or the more reliable SVD-decompositions if `fast = false`.

If `(Arij-λErij,Brij,Crij,Drij)` is the descriptor realization of `sysR.sys[i,j]`, then 
the keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Arij`, `Brij`, `Crij`, `Drij`, the absolute tolerance for the nonzero elements of `Erij`,  
and the relative tolerance for the nonzero elements of `Arij`, `Brij`, `Crij`, `Drij` and `Eirj`.
The default relative tolerance is `nij*ϵ`, where `ϵ` is the working _machine epsilon_ 
and `nij` is the order of the system `sysR.sys[i,j]`.  
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 

The keyword argument `rtolinf = tol` specifies the relative accuracy `tol` to be used 
to compute the infinity norms. The default value used is `tol = 0.00001`.
"""    
function mdgap(sysR::MDFilterIF{T}; MDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                cdinp::Bool = false, rtolinf::Real = 0.00001, offset::Real = sqrt(eps(T)), atol::Real = zero(float(real(T))), 
                atol1::Real = atol, atol2::Real = atol, rtol::Real = 0, fast::Bool = true) where T
    if ismissing(MDfreq)   
       lfreq = 0
    else         
       isa(MDfreq,Vector) || (MDfreq = [MDfreq]) 
       lfreq = length(MDfreq);
    end
    M, N = size(sysR.sys)
    β = similar(Vector{T}, M)
    γ = similar(Vector{T}, M)
    inpud = 1:sysR.mu
    for i = 1:M
        mi = sysR.mu + sysR.md[i] 
        inpw = mi : mi + sysR.mw[i] 
        γ[i] = ghinfnorm(sysR.sys[i,i][:,inpw]; rtolinf, atol1, atol2, rtol, fast, offset)[1]
        gmin = Inf
        for j = 1:N
            j == i && continue
            cdinp && (inpud = 1:sysR.mu+sysR.md[j]) 
            if lfreq == 0
                gmin = min(gmin,ghinfnorm(sysR.sys[i,j][:,inpud]; rtolinf, atol1, atol2, rtol, fast, offset)[1])
            else
                gmax = zero(T)
                for k = 1:lfreq
                   gmax = max(gmax,opnorm(evalfr(sysR.sys[i,j][:,inpud]; fval = MDfreq[k], atol1, atol2, rtol, fast)))
                end
                gmin = min(gmin,gmax)
            end
        end
        β[i] = gmin
    end

    return β./γ, β, γ
end
