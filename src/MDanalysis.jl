"""
    S = mdgenspec(sysm::MDMModel; emtest = false, fast = true, MDtol, MDGainTol, MDfreq, 
                  sdeg, atol, atol1, atol2, atol3, rtol) 
                           
Generate for a given multiple synthesis model `sysm::MDMModel` the achievable binary structure matrix `S` 
of the internal form corresponding to a set of stable model detection filters, chosen such that the `i`-th residual is 
insensitive to the `i`-th model for all control and disturbance inputs. 
If `N` is the number of component models, then `S` is an `N×N` binary matrix, whose `(i,i)`-th 
element is zero for `i = 1, ..., N`, 
and its `(i,j)`-th element is nonzero if there exists a model detection filter such that the corresponding 
`i`-th residual is sensitive to the `j`-th model for certain control inputs (if `emtest = false`) or 
for certain control and disturbance inputs (if `emtest = true`).   
   
`MDFreq = freq` specifies a vector of real frequency values or a scalar real frquency value
for strong model detectability checks (default: `MDFreq = missing`).

`MDtol = tol1` specifies the threshold `tol1` for model detectability checks
   (default: `tol1 = 0.0001`).

`MDGainTol = tol2` specifies the threshold `tol2` for strong model detectability checks
   (default: `tol2 = 0.01`). 

The keyword argument `sdeg = β` specifies a prescribed stability degree `β` for the poles of the internally 
generated candidate filters, such that the real parts of
filters poles must be less than or equal to `β`, in the continuous-time case, and 
the magnitudes of filter poles must be less than or
equal to `β`, in the discrete-time case. If `sdeg = missing` then no stabilization is performed if `FDFreq = missing`.
If `sdeg = missing` and `FDFreq = freq`, then the following default values are employed : `β = -0.05`, in continuous-time case, and  `β = 0.95`, 
in discrete-time case. 

The rank determinations in the performed reductions
are based on rank revealing QR-decompositions with column pivoting 
if `fast = true` (default) or the more reliable SVD-decompositions if `fast = false`.

The keyword arguments `atol1`, `atol2` and `rtol`, specify, respectively, 
the absolute tolerance for the nonzero elements of the state-space matrices (see below)
`Ai`, `Bi`, `Ci` and `Di`, the absolute tolerance for the nonzero elements of `Ei`,   
and the relative tolerance for the nonzero elements of all above matrices.  
The default relative tolerance is `ni*ϵ`, where `ϵ` is the working machine epsilon 
and `ni` is the orders of the system `sysm.sys[i]`. 
The keyword argument `atol` can be used 
to simultaneously set `atol1 = atol` and `atol2 = atol`. 

   
_Method:_ For a multiple synthesis model `sysm` containing `N` stable component models, 
the `i`-th component system `sysm.sys[i]` has a descriptor realization of the form 
`sysm.sys[i] = (Ai-λEi,Bi,Ci,Di)` and the corresponding input-output form is

      yi = Gui(λ)*u + Gdi(λ)*di + Gwi(λ)*wi + Gvi(λ)*vi 

with the Laplace- or Z-transformed plant outputs `yi`, control inputs `u`, 
disturbance inputs `di`, noise inputs `wi`, and auxiliary inputs `vi`,  
and with `Gui(λ)`, `Gdi(λ)`, `Gwi(λ)` and 
`Gvi(λ)`, the corresponding transfer function matrices. 

To generate the achievable structure matrix `S`, it is assumed that `N` model detection filters are used, 
`Q_i(λ)`, for `i = 1, ..., N`, where the `i`-th model detection filter `Q_i(λ)` generates the i-th residual

            r_i = Q_i(λ)*| y |
                         | u |               

and is determined such that it fulfills the nullspace synthesis condition
 
            Q_i(λ)*|Gui(λ) Gdi(λ)| = 0.
                   |  I     0    |

To generate the elements of the `i`-th row of `S`, the corresponding internal forms of the `i`-th filter 
are determined for `j = 1, ..., N` as                    

       | Ruij(λ) Rdij(λ) | := Q_i(λ)*| Guj(λ) Gdj(λ) | 
                                     |   I      0    |
                   
and `S[i,j]` is set as follows: 
   - if `MDfreq = missing` and `emtest = false`, 
     then `S[i,j] = 1` if `Ruij(λ) ̸= 0` or `S[i,j] = 0` if `Ruij(λ) = 0`;
   - if `MDfreq = missing` and `emtest = true`, 
     then `S[i,j] = 1` if `[Ruij(λ) Rdij(λ)] ̸= 0` or `S[i,j] = 0` if `[Ruij(λ) Rdij(λ)] = 0`;
   - if `MDfreq = freq` and `emtest = false`, 
     then `S[i,j] = 1` if `Ruij(λs) ̸= 0` for any `λs ∈ freq` or 
          `S[i,j] = 0` if `Ruij(λs) = 0` for all `λs ∈ freq` ;
   - if `MDfreq = freq` and `emtest = true`, 
     then `S[i,j] = 1` if `[Ruij(λs) Rdij(λs)] ̸= 0` for any `λs ∈ freq` or 
          `S[i,j] = 0` if `[Ruij(λs) Rdij(λs)] = 0` for all `λs ∈ freq` .


_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017; sec. 5.2.

[2] A. Varga, Least order fault and model detection using multi-models. 
       IEEE CDC'09, Shanghai, China, 2009.
"""
function mdgenspec(sysm::MDMModel; sdeg::Union{Real,Missing} = missing, emdtest::Bool = false,
                      MDtol::Real = 0.0001, MDGainTol::Real = 0.01, MDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                      atol::Real = zero(Float64), atol1::Real = atol, atol2::Real = atol, 
                      rtol::Real = (size(sysm.sys[1].A,1)+1)*eps(Float64)*iszero(max(atol1,atol2)), 
                      fast::Bool = true) 
   N = length(sysm.sys)   
   T1 = Float64
   # check stability of input model
   isstable(sysm) || error("the input multiple model must be stable")

   disc = (sysm.sys[1].Ts != 0);  # system type (continuous- or discrete-time)
  
   # decode options
    
   strongMD = !ismissing(MDfreq)
   strongMD && (freq = isa(MDfreq,Vector) ? MDfreq : [MDfreq]) 
  
   # set default stability degree
   sdegdefault = disc ? 0.95 : -0.05
   
   ismissing(sdeg) && (sdeg = sdegdefault)  # set desired stability degree to default value
  
 
   # decode input information
   mu = sysm.mu  
   md = sysm.md 
   m = mu .+ md      # total number of inputs
   sdegNS = strongMD ? sdegdefault : missing   
   S = falses(N,N)  
   for i = 1:N
       # solve the exact MD problem for the i-th system (with noise inputs)
       # 
       # compute a minimal left nullspace basis Q1i for the i-th system
      # form [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
      syse = [sysm.sys[i][:,1:m[i]]; eye(T1,mu,m[i])];
      #
      # compute a left nullspace basis Q = Q1 of G1 = [Gu Gd; I 0] = 0 and
      # obtain QR = [ Q R ], where R = [ Rw Ra] = Q*[Gw Ga;0 0 0]
      qtemp = glnull(syse; atol1, atol2, rtol, fast, sdeg = sdegNS)[1] 
      # check solvability conditions
      size(qtemp,1) == 0 && (@warn "empty nullspace basis for the $i-th model")


      #strongMD && (qtemp = glcf(qtemp; atol1, atol2, rtol, fast, sdeg = sdegNS)[1])
      for j = 1:N
         # check j-th model detectability
         if i != j
            mj = emdtest ? m[j] : mu
            temp = gir(qtemp*[sysm.sys[j][:,1:mj]; eye(T1,mu,mj)]; atol1, atol2, rtol)
            if strongMD
               S[i,j] = any(fdisspec_(temp, freq; FDGainTol = MDGainTol, atol1, atol2, rtol = 0, fast, stabilize = false)[1])
            else
               S[i,j] = any(fditspec_(temp; atol1, atol2, rtol, FDtol = MDtol))
            end
         end
      end
   end

   return S
   
   # end MDGENSPEC
end

"""
     mddist2c(sysm1, sysm2; MDfreq, cdinp = false, distance = "nugap", 
              rtolinf = 0.00001, offset, atol, atol1 = atol, atol2 = atol, rtol = 0, fast = true) -> (dist, fpeak)

Compute the pairwise distances between two sets of multiple component synthesis models 
`sysm1::MDMModel` and `sysm2::MDMModel`.  
If `sysm1` contains `N` component models and `sysm2` contains `M` component models, 
then the resulting `dist` and `fpeak` are `N × M` real matrices, whose `[i,j]`-th entries contain the 
distance beetwen the systems `sysm1.sys[i]` and `sysm2.sys[j]`, and, respectively, 
the corresponding peak frequency
for which the value of the computed distance is achieved.  
Both `sysm1` and  `sysm2` can be alternatively specified as vectors of component synthesis models. 

The `i`-th component system `sysm1.sys[i]` has a descriptor realization of the form 
`sysm1.sys[i] = (A1i-λE1i,B1i,C1i,D1i)` and the corresponding input-output form is

      y1i = Gu1i(λ)*u + Gd1i(λ)*di + Gw1i(λ)*wi + Gv1i(λ)*vi 

with the Laplace- or Z-transformed plant outputs `y1i`, control inputs `u`, 
disturbance inputs `di`, noise inputs `wi`, and auxiliary inputs `vi`,  
and with `Gu1i(λ)`, `Gd1i(λ)`, `Gw1i(λ)` and 
`Gv1i(λ)`, the corresponding transfer function matrices. Similarly, 
the `j`-th component system `sysm2.sys[j]` has a descriptor realization of the form 
`sysm2.sys[j]= (A2j-λE2j,B2j,C2j,D2j)` and the corresponding input-output form is

      y2j = Gu2j(λ)*u + Gd2j(λ)*dj + Gw2j(λ)*wj + Gv2j(λ)*vj 

with the Laplace- or Z-transformed plant outputs `y2j`, control inputs `u`, 
disturbance inputs `dj`, noise inputs `wj`, and auxiliary inputs `vj`,  
and with `Gu2j(λ)`, `Gd2j(λ)`, `Gw2j(λ)` and 
`Gv2j(λ)`, the corresponding transfer function matrices. 

The distance beetwen the systems `sysm1.sys[i]` and `sysm2.sys[j]` is evaluated as

      dist[i,j] = distf(G1i(λ),G2j(λ)),

where `distf(⋅,⋅)` is a distance function specified via the option parameter `distance` (see below)
and `G1i(λ)` and `G2j(λ)` are suitably defined transfer function matrices 
via the option parameter `cdinp` (see below). The corresponding peak frequency `fpeak[i,j]`
is the frequency value for which the distance is achieved.

`distance = job` specifies the distance function to be used as follows:

    job = "nugap"  - use the ν-gap distance ν(G1i(λ),G2j(λ)) defined in [1] (default)
    job = "inf"    - use the H∞-norm of the difference (i.e., ||G1i(λ)-G2j(λ)||_∞)
    job = "2"      - use the H2-norm of the difference (i.e., ||G1i(λ)-G2j(λ)||_2)

If `cdinp = false` (default), then `G1i(λ) = Gu1i(λ)` and `G2j(λ) = Gu2j(λ)`, while if `cdinp = true`
then `G1i(λ) = [Gu1i(λ) Gd1i(λ)]` and `G2j(λ) = [Gu2j(λ) Gd2j(λ)]`. 
If `Gd1i(λ)` has less columns than `Gd2j(λ)`, then `[Gd1i(λ) 0]` (with suitably padded zero columns) 
is used instead `Gd1i(λ)`, while 
if `Gd1i(λ)` has more columns than `Gd2j(λ)`, then `[Gd2j(λ) 0]` is used instead `Gd2j(λ)`.

If `MDfreq = ω`, where `ω` is a vector of real frequency values, then each distance `dist[i,j]` 
is the maximum of pointwise distances evaluated for all frequencies in `ω` and 
the corresponding frequency `fpeak[i,j]`
is the value for which the maximum pointwise distance is achieved.

The stability boundary offset, `β`, to be used to assess the finite poles/zeros which belong to the
boundary of the stability domain can be specified via the keyword parameter `offset = β`.
Accordingly, for a continuous-time system, these are the finite poles/zeros having 
real parts within the interval `[-β,β]`, while for a discrete-time system, 
these are the finite pole/zeros having moduli within the interval `[1-β,1+β]`. 
The default value used for `β` is `sqrt(ϵ)`, where `ϵ` is the working machine precision. 

Pencil reduction algorithms are employed to compute range and coimage spaces involved in evaluating 
the ν-gap or the H∞-norm distances. These algorithms perform rank decisions based on rank 
revealing QR-decompositions with column pivoting 
if `fast = true` or the more reliable SVD-decompositions if `fast = false`.

The keyword arguments `atol1`, `atol2` and `rtol`, specify, respectively, 
the absolute tolerance for the nonzero elements of the state-space matrices 
   `A1i`, `A2j`, `B1i`, `B2j`, `C1i`, `C2j`, `D1i` and `D2j`,
the absolute tolerance for the nonzero elements of `E1i` and `E2j`,   
and the relative tolerance for the nonzero elements of all above matrices.  
The default relative tolerance is `nij*ϵ`, where `ϵ` is the working machine epsilon 
and `nij` is the maximum of the orders of the systems `sysm1.sys[i]` and `sysm2.sys[j]`. 
The keyword argument `atol` can be used to simultaneously set `atol1 = atol`, `atol2 = atol`. 

The keyword argument `rtolinf = tol` specifies the relative accuracy `tol` to be used 
to compute infinity norms. The default value used is `tol = 0.00001`.
   
_Method:_ The evaluation of ν-gap uses the definition proposed in [1],
extended to generalized LTI (descriptor) systems. 

_References:_

[1] G. Vinnicombe. Uncertainty and feedback: H∞ loop-shaping and the ν-gap metric. 
    Imperial College Press, London, 2001. 
"""
function mddist2c(sysm1::MDMModel, sysm2::MDMModel; MDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                cdinp::Bool = false, distance::AbstractString = "nugap", rtolinf::Real = 0.00001, 
                offset::Real = sqrt(eps(Float64)), atol::Real = zero(float(real(Float64))), 
                atol1::Real = atol, atol2::Real = atol, rtol::Real = 0, fast::Bool = true) 
    T = Float64
    if ismissing(MDfreq)   
       lfreq = 0
    else         
       isa(MDfreq,Vector) || (MDfreq = [MDfreq]) 
       lfreq = length(MDfreq);
    end
    N = length(sysm1.sys)
    M = length(sysm2.sys)
    p = size(sysm1.sys[1],1) 
    p == size(sysm2.sys[1],1) || error("the two multiple-models must have the same number of outputs")
    Ts = sysm1.sys[1].Ts  
    Ts == sysm2.sys[1].Ts || error("the two multiple-models must have the same sampling time")
    if isequal(distance,"nugap")
       jobdist = 0
    elseif isequal(distance,"Inf") || isequal(distance,"inf")
       jobdist = Inf
    elseif isequal(distance,"2")
       jobdist = 2
    else
       error("No such distance selection option")
    end
    mu = sysm1.mu
    mu == sysm2.mu || error("the two multiple-models must have the same number of control inputs")
    if cdinp 
        md1 = sysm1.md
        mdmax1 = maximum(md1)
        md2 = sysm2.md
        mdmax2 = maximum(md2)
        mdmax = max(mdmax1,mdmax2)
        mddif1 = mdmax .- md1
        mddif2 = mdmax .- md2
        m = mu+mdmax
        mdext = any(mddif1 .> 0) || any(mddif2 .> 0)
    else
        m = mu
        mdext = false
    end

    dist = similar(Matrix{T}, N, M)
    fpeak = similar(Matrix{T}, N, M)
    if jobdist == 0
       # compute the normalized right coprime factorizations only once
       R2 = similar(Vector{DescriptorStateSpace{T}}, N)
       for i = 1:N
          # R2 = [N2; M2]
          if mdext 
             R2[i] = grange([sysm1.sys[i][:,1:mu+md1[i]] zeros(T,p,mddif1[i]); eye(T,m)]; inner = true, 
                             atol1, atol2, rtol, fast, offset)[1]
          else
             R2[i] = grange([sysm1.sys[i][:,1:m]; eye(T,m)]; inner = true, 
                            atol1, atol2, rtol, fast, offset)[1]
          end
       end

       for j = 1:M
           if mdext 
              R1 = grange([sysm2.sys[j][:,1:mu+md2[j]] zeros(T,p,mddif2[j]); eye(T,m)]; inner = true, 
                           atol1, atol2, rtol, fast, offset)[1]
              L1 = gcrange([sysm2.sys[j][:,1:mu+md2[j]] zeros(T,p,mddif2[j]) eye(T,p)]; coinner = true, 
                            atol1, atol2, rtol, fast, offset)[1]
           else
              R1 = grange([sysm2.sys[j][:,1:m]; eye(T,m)]; inner = true, 
                           atol1, atol2, rtol, fast, offset)[1]
              L1 = gcrange([sysm2.sys[j][:,1:m] eye(T,p)]; coinner = true, 
                            atol1, atol2, rtol, fast, offset)[1]
           end
           for i = 1:N
               # check conditions on det(R2'*R1)
               syst = gir(R2[i]'*R1; atol1, atol2, rtol, fast)  
               infoz = gzeroinfo(syst; atol1, atol2, rtol, fast, offset)[2]
               # check invertibility and presence of zeros on the boundary of 
               # stability domain
               if infoz.nrank != order(syst)+m || infoz.nfsbz > 0
                  dist[i,j] = one(T); fpeak[i,j] = zero(T)
                  dist[j,i] = one(T); fpeak[j,i] = zero(T)
                  continue
               end
               # evaluate winding number 
               infop = gpoleinfo(syst; atol1, atol2, rtol, fast, offset)[2]
               wno = infoz.nfuz - infop.nfuev + infoz.niz - infop.nip
               # check condition on winding number 
               if wno != 0
                  # nonzero winding number
                  dist[i,j] = one(T); fpeak[i,j] = zero(T)
                  continue
               end
               # compute the underlying system to compute the ν-gap distance 
               # using the definition of Vinnicombe
               syst = L1*[zeros(T,m,p) -eye(T,m);eye(T,p,p+m)]*R2[i]
               if lfreq == 0
                  # compute the ν-gap using the definition of Vinnicombe
                  dist[i,j], fpeak[i,j] = ghinfnorm(syst; rtolinf, atol1, atol2, rtol, fast, offset)
               else
                  H = freqresp(syst,MDfreq)
                  tmax = opnorm(H[:,:,1]) 
                  fmax = MDfreq[1]
                  for ii = 2:lfreq 
                      temp = opnorm(H[:,:,ii])
                      temp > tmax && (tmax = temp; fmax = MDfreq[ii])
                  end 
                  dist[i,j] = tmax  
                  fpeak[i,j] = fmax                          
               end
            end
       end
       return dist, fpeak
    end

    # compute 2- or ∞-norm based distances
    for i = 1:N
        for j = 1:M
            if mdext 
               syst = [sysm1.sys[i][:,1:mu+md1[i]] zeros(T,p,mddif1[i])] - 
                      [sysm2.sys[j][:,1:mu+md2[j]] zeros(T,p,mddif2[j])]
            else
               syst = sysm1.sys[i][:,1:m] - sysm2.sys[j][:,1:m]
            end
            if lfreq == 0
                if jobdist == 2
                    dist[i,j] = gh2norm(syst; atolinf = rtolinf, atol1, atol2, rtol, fast, offset)
                    fpeak[i,j] = zero(T)
                else
                    dist[i,j], fpeak[i,j] = ghinfnorm(syst; rtolinf, atol1, atol2, rtol, fast, offset)
                end
            else
                H = freqresp(syst,MDfreq)
                tmax = opnorm(H[:,:,1]) 
                fmax = MDfreq[1]
                for i = 2:lfreq 
                    temp = opnorm(H[:,:,i])
                    temp > tmax && (tmax = temp; fmax = MDfreq[i])
                end 
                dist[i,j] = tmax  
                fpeak[i,j] = fmax                          
            end
        end
    end
    return dist, fpeak
end
mddist2c(sysm1::MDMModel, sysm2::Vector{<:MDModel}; kwargs...) = mddist2c(sysm1, mdmodset(sysm2); kwargs... )
mddist2c(sysm1::Vector{<:MDModel}, sysm2::MDMModel; kwargs...) = mddist2c(mdmodset(sysm1), sysm2; kwargs... )
mddist2c(sysm1::Vector{<:MDModel}, sysm2::Vector{<:MDModel}; kwargs...) = mddist2c(mdmodset(sysm1), mdmodset(sysm2); kwargs... )
"""
     mddist(sysm; MDfreq, cdinp = false, distance = "nugap",  rtolinf = 0.00001, 
            offset, atol, atol1 = atol, atol2 = atol, rtol = 0, fast = true) -> (dist, fpeak)

Compute the pairwise distances between the component models of `sysm::MDMModel`.  
If `sysm` contains `N` component models, then the resulting `dist` and `fpeak` are 
`N × N` real symmetric matrices, whose `[i,j]`-th elements contain the 
distance beetwen the systems `sysm.sys[i]` and `sysm.sys[j]`, and, 
respectively, the corresponding `peak frequency`
for which the value of the computed distance is achieved.  
`sysm` can be alternatively specified as a vector of component synthesis models. 

The `i`-th component system `sysm.sys[i]` has a descriptor realization of the form 
`sysm.sys[i] = (Ai-λEi,Bi,Ci,Di)` and the corresponding input-output form is

      yi = Gui(λ)*u + Gdi(λ)*di + Gwi(λ)*wi + Gvi(λ)*vi 

with the Laplace- or Z-transformed plant outputs `yi`, control inputs `u`, 
disturbance inputs `di`, noise inputs `wi`, and auxiliary inputs `vi`,  
and with `Gui(λ)`, `Gdi(λ)`, `Gwi(λ)` and `Gvi(λ)`, the corresponding transfer function matrices.  

The distance beetwen the systems `sysm.sys[i]` and `sysm.sys[j]` is evaluated as

      dist[i,j] = distf(G1i(λ),G2j(λ)),

where `distf(⋅,⋅)` is a distance function specified via the option parameter `distance` (see below)
and `G1i(λ)` and `G2j(λ)` are suitably defined transfer function matrices 
via the option parameter `cdinp` (see below). The corresponding peak frequency `fpeak[i,j]`
is the frequency value for which the distance is achieved.

`distance = job` specifies the distance function to be used as follows:

    job = "nugap"  - use the ν-gap distance ν(G1i(λ),G2j(λ)) defined in [1] (default)
    job = "inf"    - use the H∞-norm of the difference (i.e., ||G1i(λ)-G2j(λ)||_∞)
    job = "2"      - use the H2-norm of the difference (i.e., ||G1i(λ)-G2j(λ)||_2)

If `cdinp = false` (default), then `G1i(λ) = Gui(λ)` and `G2j(λ) = Guj(λ)`, while if `cdinp = true`
then `G1i(λ) = [Gui(λ) Gdi(λ)]` and `G2j(λ) = [Guj(λ) Gdj(λ)]`. 
If `Gdi(λ)` has less columns than `Gdj(λ)`, then `[Gdi(λ) 0]` (with suitably padded zero columns) 
is used instead `Gdi(λ)`, while 
if `Gdi(λ)` has more columns than `Gdj(λ)`, then `[Gdj(λ) 0]` is used instead `Gdj(λ)`.

If `MDfreq = ω`, where `ω` is a vector of real frequency values, then each distance `dist[i,j]` 
is the maximum of pointwise distances evaluated for all frequencies in `ω` and 
the corresponding frequency `fpeak[i,j]`
is the value for which the maximum pointwise distance is achieved.

The stability boundary offset, `β`, to be used to assess the finite poles/zeros which belong to the
boundary of the stability domain can be specified via the keyword parameter `offset = β`.
Accordingly, for a continuous-time system, these are the finite poles/zeros having 
real parts within the interval `[-β,β]`, while for a discrete-time system, 
these are the finite pole/zeros having moduli within the interval `[1-β,1+β]`. 
The default value used for `β` is `sqrt(ϵ)`, where `ϵ` is the working machine precision. 

Pencil reduction algorithms are employed to compute range and coimage spaces involved in evaluating 
the ν-gap or the H∞-norm distances. These algorithms perform rank decisions based on rank 
revealing QR-decompositions with column pivoting 
if `fast = true` or the more reliable SVD-decompositions if `fast = false`.

The keyword arguments `atol1`, `atol2` and `rtol`, specify, respectively, 
the absolute tolerance for the nonzero elements of the state-space matrices 
`Ai`, `Bi`, `Ci`,`Di`, `Aj`, `Bj`, `Cj` and `Dj`,
the absolute tolerance for the nonzero elements of `Ei` and `Ej`,   
and the relative tolerance for the nonzero elements of all above matrices.  
The default relative tolerance is `nij*ϵ`, where `ϵ` is the working machine epsilon 
and `nij` is the maximum of the orders of the systems `sysm.sys[i]` and `sysm.sys[j]`. 
The keyword argument `atol` can be used to simultaneously set `atol1 = atol`, `atol2 = atol`. 

The keyword argument `rtolinf = tol` specifies the relative accuracy `tol` to be used 
to compute infinity norms. The default value used is `tol = 0.00001`.
   
_Method:_ The evaluation of ν-gap uses the definition proposed in [1],
extended to generalized LTI (descriptor) systems. 

_References:_

[1] G. Vinnicombe. Uncertainty and feedback: H∞ loop-shaping and the ν-gap metric. 
    Imperial College Press, London, 2001. 
"""
function mddist(sysm::MDMModel; MDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                cdinp::Bool = false, distance::AbstractString = "nugap", rtolinf::Real = 0.00001, 
                offset::Real = sqrt(eps(Float64)), atol::Real = zero(float(real(Float64))), 
                atol1::Real = atol, atol2::Real = atol, rtol::Real = 0, fast::Bool = true) 
    T = Float64
    if ismissing(MDfreq)   
       lfreq = 0
    else         
       isa(MDfreq,Vector) || (MDfreq = [MDfreq]) 
       lfreq = length(MDfreq);
    end
    N = length(sysm.sys)
    p = size(sysm.sys[1],1) 
    if isequal(distance,"nugap")
       jobdist = 0
    elseif isequal(distance,"Inf") || isequal(distance,"inf")
       jobdist = Inf
    elseif isequal(distance,"2")
       jobdist = 2
    else
       error("No such distance selection option")
    end
    mu = sysm.mu
    if cdinp 
        md = sysm.md
        mdmax = maximum(md)
        mddif = mdmax .- md
        m = mu+mdmax
        mdext = any(mddif .> 0)
    else
        m = mu
        mdext = false
    end

    dist = similar(Matrix{T}, N, N)
    fpeak = similar(Matrix{T}, N, N)
    if jobdist == 0
       # compute the normalized right coprime factorizations only once
       R2 = similar(Vector{DescriptorStateSpace{T}}, N)
       for j = 1:N
          # R2 = [N2; M2]
          if mdext 
             R2[j] = grange([sysm.sys[j][:,1:mu+md[j]] zeros(T,p,mddif[j]); eye(T,m)]; inner = true, 
                             atol1, atol2, rtol, fast, offset)[1]
          else
             R2[j] = grange([sysm.sys[j][:,1:m]; eye(T,m)]; inner = true, 
                            atol1, atol2, rtol, fast, offset)[1]
          end
       end

       for i = 1:N
           if mdext 
              R1 = grange([sysm.sys[i][:,1:mu+md[i]] zeros(T,p,mddif[i]); eye(T,m)]; inner = true, 
                           atol1, atol2, rtol, fast, offset)[1]
              L1 = gcrange([sysm.sys[i][:,1:mu+md[i]] zeros(T,p,mddif[i]) eye(T,p)]; coinner = true, 
                            atol1, atol2, rtol, fast, offset)[1]
           else
              R1 = grange([sysm.sys[i][:,1:m]; eye(T,m)]; inner = true, 
                           atol1, atol2, rtol, fast, offset)[1]
              L1 = gcrange([sysm.sys[i][:,1:m] eye(T,p)]; coinner = true, 
                            atol1, atol2, rtol, fast, offset)[1]
           end
           dist[i,i] = zero(T); fpeak[i,i] = zero(T)
           for j = i+1:N
               # check conditions on det(R2'*R1)
               syst = gir(R2[j]'*R1; atol1, atol2, rtol, fast)  
               infoz = gzeroinfo(syst; atol1, atol2, rtol, fast, offset)[2]
               # check invertibility and presence of zeros on the boundary of 
               # stability domain
               if infoz.nrank != order(syst)+m || infoz.nfsbz > 0
                  dist[i,j] = one(T); fpeak[i,j] = zero(T)
                  dist[j,i] = one(T); fpeak[j,i] = zero(T)
                  continue
               end
               # evaluate winding number 
               infop = gpoleinfo(syst; atol1, atol2, rtol, fast, offset)[2]
               wno = infoz.nfuz - infop.nfuev + infoz.niz - infop.nip
               # check condition on winding number 
               if wno != 0
                  # nonzero winding number
                  dist[i,j] = one(T); fpeak[i,j] = zero(T)
                  dist[j,i] = one(T); fpeak[j,i] = zero(T)
                  continue
               end
               # compute the underlying system to compute the ν-gap distance 
               # using the definition of Vinnicombe
               syst = L1*[zeros(T,m,p) -eye(T,m);eye(T,p,p+m)]*R2[j]
               if lfreq == 0
                  # compute the ν-gap using the definition of Vinnicombe
                  dist[i,j], fpeak[i,j] = ghinfnorm(syst; rtolinf, atol1, atol2, rtol, fast, offset)
                  dist[j,i] = dist[i,j]
                  fpeak[j,i] = fpeak[i,j]  
               else
                  H = freqresp(syst,MDfreq)
                  tmax = opnorm(H[:,:,1]) 
                  fmax = MDfreq[1]
                  for i = 2:lfreq 
                      temp = opnorm(H[:,:,i])
                      temp > tmax && (tmax = temp; fmax = MDfreq[i])
                  end 
                  dist[i,j] = tmax  
                  dist[j,i] = tmax  
                  fpeak[i,j] = fmax                          
                  fpeak[j,i] = fmax                          
               end
            end
       end
       return dist, fpeak
    end

    # compute 2- or ∞-norm based distances
    for i = 1:N
        dist[i,i] = zero(T); fpeak[i,i] = zero(T)
        for j = i+1:N
            if mdext 
               syst = [sysm.sys[j][:,1:mu+md[j]] zeros(T,p,mddif[j])] - 
                      [sysm.sys[i][:,1:mu+md[i]] zeros(T,p,mddif[i])]
            else
               syst = sysm.sys[j][:,1:m] - sysm.sys[i][:,1:m]
            end
            if lfreq == 0
                if jobdist == 2
                    dist[i,j] = gh2norm(syst; atolinf = rtolinf, atol1, atol2, rtol, fast, offset)
                    fpeak[i,j] = zero(T)
                else
                    dist[i,j], fpeak[i,j] = ghinfnorm(syst; rtolinf, atol1, atol2, rtol, fast, offset)
                end
                dist[j,i] = dist[i,j]
                fpeak[j,i] = fpeak[i,j]  
            else
                H = freqresp(syst,MDfreq)
                tmax = opnorm(H[:,:,1]) 
                fmax = MDfreq[1]
                for i = 2:lfreq 
                    temp = opnorm(H[:,:,i])
                    temp > tmax && (tmax = temp; fmax = MDfreq[i])
                end 
                dist[i,j] = tmax  
                dist[j,i] = tmax  
                fpeak[i,j] = fmax                          
                fpeak[j,i] = fmax                          
            end
        end
    end
    return dist, fpeak
end
mddist(sysm::Vector{<:MDModel}; kwargs...) = mddist(mdmodset(sysm); kwargs... )


