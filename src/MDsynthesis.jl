"""
    amdsyn(sysm::MDMModel; rdim, MDSelect, nullspace = true, simple = false, minimal = true, 
                           emtest = false, normalize = false, fast = true, 
                           sdeg, smarg, poles, HDesign, MDtol, MDGainTol, MDfreq, 
                           epsreg = 0.1, sdegzer, nonstd = 1, rtolinf= 0.0001, 
                           tcond, offset, atol, atol1, atol2, atol3, rtol) 
                           -> (Q::MDFilter, R::MDFilterIF, info)

Solve the _approximate model detection problem_ (AMDP) for a given multiple synthesis model
`sysm::MDMModel`. The computed stable and proper filter objects `Q` and `R` contain the 
model detection filter, representing the solution of the AMDP, and its internal form, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.degs`, `info.MDperf`, 
`info.MDgap`, `info.HDesign` and `info.nonstandard`, 
contains additional synthesis related information (see below). 

For a multiple synthesis model `sysm` containing `N` stable component models, 
the `i`-th component system `sysm.sys[i]` has a descriptor realization of the form 
`sysm.sys[i] = (Ai-λEi,Bi,Ci,Di)` and the corresponding input-output form is

      yi = Gui(λ)*u + Gdi(λ)*di + Gwi(λ)*wi + Gvi(λ)*vi 

with the Laplace- or Z-transformed plant outputs `yi`, control inputs `u`, 
disturbance inputs `di`, noise inputs `wi`, and auxiliary inputs `vi`,  
and with `Gui(λ)`, `Gdi(λ)`, `Gwi(λ)` and 
`Gvi(λ)`, the corresponding transfer function matrices. 

The model detection filter object `Q`, contains in `Q.sys` the resulting bank of `N` filters. 
The `i`-th filter `Q.sys[i]` is in a standard state-space form and generates `r_i`, 
the `i`-th component (scalar or vector) of the overall residual vector `r := [r_1; r_2; ...; r_N]`. 
The corresponding input-output (implementation) form of the `i`-th filter is

            r_i = Qyi(λ)*y + Qui(λ)*u ,            

where `Qyi(λ)` and `Qui(λ)` are the transfer function matrices from the 
measured outputs and control inputs to the `i`-th residual component. 
The dimensions of measured outputs and control inputs are contained in the integers  
`Q.ny` and `Q.mu`, respectively.

The model detection filter internal form object `R`, contains `R.sys`, the resulting array of `N × N` 
internal form filters.  
The `(i,j)`-th component filter `R.sys[i,j]` is in a standard state-space form, 
generates the residual signal `r_ij`, and corresponds to the 
input-output form

       r_ij = Ruij(λ)*u + Rdij(λ)*dj + Rwij(λ)*wj + Rvij(λ)*vj ,

where 

       | Ruij(λ) Rdij(λ) Rwij(λ) Raij(λ) | := |Qyi(λ) Qui(λ)]*| Guj(λ) Gdj(λ) Gwj(λ) Gaj(λ) |. 
                                                              |   I    0      0      0      |

The solution of the AMDP ensures that for the `i`-th filter, `Ruii(λ) = 0`, `Rdii(λ) = 0`,  
`[Ruij(λ) Rdij(λ)]` ``\\neq`` `0` for `j` ``\\neq`` `i`, and, additionally, the maximization  of
the achieved gaps (see description of `info.MDgap`). 

Various user options can be specified via keyword arguments as follows:

`MDSelect = ifilt` specifies in the vector `ifilt` the indices of the desired filters to be designed
(default: `ifilt = 1:N`)

If `minimal = true` (default), least order filter synthesis is performed to determine 
each of the component filters
`Q.sys[i]` for `i = 1, ...,N` and `R.sys[i,j]` for `i, j = 1, ...,N`  while 
with `minimal = false` full order synthesis is performed.  

If `HDesign = H`, then `H` is an `N`-dimensional array of full row rank or empty design matrices `H = [H_1, ..., H_N]`,
where `H_i` is the design matrix employed for the synthesis of the `i`-th component filter (default: `HDesign = missing`)

`rdim = q` specifies the vector `q`, whose `i`-th component `q[i]` specifies 
the number of residual outputs for the `i`-th component filter `Q.sys[i]`. 
If `q` is a scalar, then a vector `rdim` with all components equal to `q` is assumed.
The default value of `q[i]` is chosen as follows: if `HDesign = missing` or `H_i` is empty then  
`q[i] = 1`, if `minimal = true`, or `q[i]` is the number of the nullspace basis 
vectors used for the synthesis of `Q.sys [i]`, if `minimal = false`; 
if  `H_i` specifies a full row rank design matrix, then `q[i]` is the row dimension of `H_i`. 

If `emdtest = false` (default), only the input channels are used for model detectability tests.
If `emdtest = true`, extended model detectability tests are performed using both control and
disturbance input channels.  

`MDfreq = ω` specifies a vector of real frequency values or a scalar real frequency value
for strong model detectability checks (default: `MDfreq = missing`).

If `nullspace = false` (default),  a full-order observer based nullspace basis is 
used for the synthesis of the `i`-th filter whenever the `i`-th component system has no disturbance inputs.
Otherwise, minimal proper nullspace bases are used. 
If `nullspace = true`, minimal proper nullspace bases are used for the synthesis of the
model detection filters. 

If `simple = true`, simple proper nullspace bases are emplyed for synthesis. 
The orders of the basis vectors employed for the synthesis of the `i`-th filter
are provided in `info.deg[i]`. 
If `simple = false` (default), then no simple bases are computed. 

`offset = β` specifies the boundary offset `β` to assess the stability of poles. 
Accordingly, for the stability of a continuous-time system all real parts of poles must be at most `-β`, 
while for the stability of a discrete-time system all moduli of poles must be at most `1-β`. 
The default value used for `β` is `sqrt(ϵ)`, where `ϵ` is the working machine precision. 

`smarg = α` specifies the stability margin which defines the stability 
domain `Cs` of poles, as follows: 
for a continuous-time system, `Cs` is the set of complex numbers 
with real parts at most `α`, 
while for a discrete-time system, `Cs` is the set of complex numbers with 
moduli at most `α < 1` (i.e., the interior of a disc of radius `α` centered in the origin). 
If `smarg` is missing, then the employed default values are `α = -β` 
for a continuous-time system and `α = 1-β` for a discrete-time system, 
where `β` is the boundary offset specified by the keyword argument `offset = β`. 

`sdeg = γ` is the prescribed stability degree for the poles of the filters `Q` and `R` 
(default: `γ = -0.05` for the real parts of poles for a continuous-time system and
`γ = 0.95` for the magnitudes of poles for a discrete-time system). 

`poles = v` specifies a complex vector `v` containing a complex conjugate set  
of desired poles within the stability domain `Cs` to be assigned for the filters `Q` and `R`
(default: `poles = missing`).

`tcond = tcmax` specifies the maximum alowed condition number `tcmax` 
of the employed non-orthogonal transformations (default: `tcmax = 1.e4`).

`MDtol = tol1` specifies the threshold `tol1` for model detectability checks
   (default: `tol1 = 0.0001`).

`MDGainTol = tol2` specifies the threshold `tol2` for strong model detectability checks
   (default: `tol2 = 0.01`). 

If `normalize = false` (default), the `i`-th component filter `Q.sys[i]` is scaled such that
the minimum gain of `R.sys[i,j]` for `j = 1, ..., N`,   `j` ``\\neq`` `i`, is equal to one. 
If `normalize = true`, the standard normalization of component filters is performed to ensure
equal gains for `R.sys[1,j]` and `R.sys[j,1]`. 
This option is only possible if `ifilt[1] = 1` (see `MDSelect`).

`epsreg = ϵ` specifies the value of the regularization parameter `ϵ` (default: `ϵ = 0.1`)

`sdegzer = δ` specifies the prescribed stability degree `δ` for zeros shifting
   (default: `δ = −0.05` for a continuous-time system `sysf.sys` and 
   `δ = 0.95` for a discrete-time system `sysf.sys`).

`nonstd = job` specifies the option to handle nonstandard optimization problems as follows:

      job = 1 – use the quasi-co-outer–co-inner factorization (default);
      job = 2 – use the modified co-outer–co-inner factorization with the
                regularization parameter `ϵ`;
      job = 3 – use the Wiener-Hopf type co-outer–co-inner factorization;
      job = 4 – use the Wiener-Hopf type co-outer-co-inner factorization with
                zero shifting of the non-minimum phase factor using the
                stabilization parameter `δ`;
      job = 5 – use the Wiener-Hopf type co-outer-co-inner factorization with
                the regularization of the non-minimum phase factor using the
                regularization parameter `ϵ`. 

The rank determinations in the performed reductions
are based on rank revealing QR-decompositions with column pivoting 
if `fast = true` (default) or the more reliable SVD-decompositions if `fast = false`.

The keyword arguments `atol1`, `atol2` and `rtol`, specify, respectively, 
the absolute tolerance for the nonzero elements of the state-space matrices 
`Ai`, `Bi`, `Ci` and `Di`, the absolute tolerance for the nonzero elements of `Ei`,   
and the relative tolerance for the nonzero elements of all above matrices.  
The default relative tolerance is `ni*ϵ`, where `ϵ` is the working machine epsilon 
and `ni` is the orders of the system `sysm.sys[i]`. 
The keyword argument `atol3` is an absolute tolerance for observability tests
   (default: internally determined value). 
The keyword argument `atol` can be used 
to simultaneously set `atol1 = atol`, `atol2 = atol` and `atol3 = atol`. 

The resulting named tuple `info` contains `(tcond, degs, MDperf, MDgap, HDesign, nonstandard) `, where:

`info.tcond` is the maximum of the condition numbers of the employed 
   non-orthogonal transformation matrices; a warning is issued if `info.tcond >= tcmax`;

`info.degs` is an `N`-dimensional vector, whose `i`-th component is an integer vector 
containing the degrees of the basis vectors of the employed simple
nullspace basis for the synthesis of the `i`-th filter component, if `simple = true`, 
and the degrees of the
basis vectors of an equivalent polynomial nullspace basis, if `simple = false`;

`info.MDperf` is an `N × N` array containing the resulting distance-mapping performance, 
representing  the peak gains of the associated internal representations of 
the model detection filters. If the `(i,j)`-th component filter `R.sys[i,j]` 
has the input-output form
 
     rij = Ruij(λ)*u + Rdij(λ)*dj + Rwij(λ)*wj + Rvij(λ)*vj ,

then, the `(i,j)`-th performance gain is evaluated as 
`info.MDperf[i,j] = ||Ruij(λ)||∞` if `emdtest = false` (default) or 
`info.MDperf[i,j] = ||[Ruij(λ) Rdij(λ)]||∞` if `emdtest = true`.  
If `MDfreq = ω`, then each gain `info.MDperf[i,j]` represents
the maximum of 2-norm pointwise gains evaluated for all frequencies in `ω`. 

`info.MDgap` is an `N`-dimensional vector of achieved noise gaps, where 
the `i`-th gap is computed as `info.MDgap[i] = β[i]/γ[i]`, where `β[i] = min(||Rij(λ)||∞` for `i` ``\\neq`` `j`) 
and `γ[i] = ||Rwii(λ)||∞`.  `Rij(λ)` is defined as `Rij(λ) = Ruij(λ)`  if `emtest = false` (default) 
or `Rij(λ) = [Ruij(λ) Rdij(λ)]` if `emtest = true`. If `MDfreq = ω`, where `ω` is a given vector of real frequency values, 
then each gain `β[i]` represents the minimum of the maximum of 2-norm pointwise gains evaluated for all frequencies in `ω`. 

`info.HDesign` is an `N`-dimensional vector, whose `i`-th component 
is the design matrix `H_i` employed for the synthesis of 
the `i`-th model detection filter.

`info.nonstandard` is an `N`-dimensional vector, whose `i`-th component is `true`
 if the synthesis of the `i`-th filter involved a nonstandard optimization problem and 
is `false` otherwise.
   
_Method:_ An extension of the Procedure AMD from [1] is implemented to solve 
the approximate model detection problem. This procedure relies on the nullspace-based
synthesis method of model detection filters proposed in [2].

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017; sec. 5.2.

[2] A. Varga, Least order fault and model detection using multi-models. 
       IEEE CDC'09, Shanghai, China, 2009.
"""
function amdsyn(sysm::MDMModel; rdim::Union{Vector{Int},Int,Missing} = missing, 
                      MDSelect::Vector{Int} = Vector(1:length(sysm.sys)), poles::Union{AbstractVector,Missing} = missing, 
                      sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing,  
                      nullspace::Bool = false, minimal::Bool = true, simple::Bool = false, emdtest::Bool = false,
                      MDtol::Real = 0.0001, MDGainTol::Real = 0.01, MDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                      tcond::Real = 1.e4, HDesign::Union{Vector{<:Matrix},Matrix,Missing} = missing, normalize::Bool = false,
                      epsreg::Real = 0.1, sdegzer::Union{Real,Missing} = missing, nonstd::Int = 1, rtolinf::Real = 0.0001,
                      offset::Real = sqrt(eps(Float64)), atol::Real = zero(Float64), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                      rtol::Real = (size(sysm.sys[1].A,1)+1)*eps(Float64)*iszero(max(atol1,atol2)), 
                      fast::Bool = true) 
   N = length(sysm.sys)   
   T1 = Float64
   # check stability of input model
   isstable(sysm) || error("the input multiple model must be stable")

   p = size(sysm.sys[1],1);       # number of measurable outputs              
   Ts = sysm.sys[1].Ts                  
   disc = (Ts != 0);  # system type (continuous- or discrete-time)
  
   # decode options

   M = length(MDSelect)
   any(MDSelect .<= 0) && error("the selected filter indices must be positive")   
   any(MDSelect .> N) && error("the selected filter indices must be at most $N")   
   allunique(MDSelect) || error("the selected filter indices must be unique")
   M > 0 && MDSelect[1] != 1 && (normalize = false)
  
   strongMD = !ismissing(MDfreq)
   strongMD && !isa(MDfreq,Vector) && (MDfreq = [MDfreq]) 
   lfreq = strongMD ? length(MDfreq) : 0
   
   # set default stability degree
   sdegdefault = disc ? 0.95 : -0.05
   
  
   # stability margin
   ismissing(smarg) && (smarg = disc ? 1-offset : -offset)  # set default stability margin
   
   poles_nomissing = !ismissing(poles)
    
   # sort desired poles
   if poles_nomissing 
      tempc = poles[imag.(poles) .> 0]
      if !isempty(tempc)
         tempc1 = conj(poles[imag.(poles) .< 0])
         isequal(tempc[sortperm(real(tempc))],tempc1[sortperm(real(tempc1))]) ||
                 error("poles must be a self-conjugated complex vector")
      end
      # check that all eigenvalues are inside of the stability region
      ( ((disc && any(abs.(poles) .> 1-offset) )  || (!disc && any(real.(poles) .> -offset)))  &&
            error("The elements of poles must lie in the stability region of interest") )
   end   
   
   ismissing(sdeg) && ismissing(poles) && (sdeg = sdegdefault)  # set desired stability degree to default value
  
   # desired stability degree for zero shifting
   ismissing(sdegzer) && (sdegzer = sdegdefault)

   # desired number of filter outputs   
   if !ismissing(rdim)
      if isa(rdim,Vector) 
         length(rdim) == M || error("dimension of rdim must be equal to the number of selected filters")
         minimum(rdim) > 0 || error("all components of rdim must be positive")
      else
         rdim > 0 || error("rdim must be positive")
         rdim = fill(rdim, M)
      end
   end

   # imposed design option to form linear combinations of basis vectors
   emptyHD = ismissing(HDesign)

   if !emptyHD
       if isa(HDesign,Matrix)
          size(HDesign,1) == rank(HDesign) || error("HDesign must have full row rank")
          HDesign = fill(HDesign,M)
       else
          size(HDesign,1) == M || error("number of HDesign components must be equal to the number of selected filters $M")
       end
       if !ismissing(rdim) 
          for i = 1:M
              if !isempty(HDesign[i])
                 mH = size(HDesign[i],1)
                 if mH > 0
                    mH == rdim[i] || error("row dimension of HDesign[$i] must be equal to rdim[$i]")
                    mH == rank(HDesign[i]) || error("HDesign[$i] must have full row rank")
                 end
              end
          end
       end
   else
      HDesign = fill(T1[],M)
   end

 
  
   # decode input information
   mu = sysm.mu  
   inpu = 1:mu 
   md = sysm.md 
   mw = sysm.mw 
   ma = sysm.ma 
   
   m = mu .+ (md+mw+ma)      # total number of inputs
    
   infotcond = similar(Vector{T1},M)
   infodegs = similar(Vector{Vector{Int}},M)
   HDesign1 = similar(Vector{Array{T1,2}},M)
   Qt = similar(Vector{DescriptorStateSpace{T1}},M)
   Rt = similar(Matrix{DescriptorStateSpace{T1}},M,N)
   distinf = fill(-1.,M,N) 
   ismissing(rdim) ? rdimtarget = missing : rdimtarget = copy(rdim)
   for i = 1:M
      isel = MDSelect[i]
      # solve the exact MD problem for the i-th system (with noise inputs)
      # 
      # compute a minimal left nullspace basis Q1i for the i-th system
      desc = (sysm.sys[isel].E != I)
      m2 = mw[isel]+ma[isel]
      sdegNS = strongMD ? sdegdefault : missing     
      if nullspace || simple || md[isel] > 0 || (desc && rcond(sysm.sys[isel].E) < 1.e-7 )
         # form [ Gu Gd Gw Ga; I 0 0 0] 
         syse = [sysm.sys[isel]; eye(mu,m[isel])];
         #
         # compute a left nullspace basis Q = Q1 of G1 = [Gu Gd; I 0] = 0 and
         # obtain QR = [ Q R ], where R = [ Rw Raux] = Q*[Gw Ga;0 0]
         qtemp, info1 = glnull(syse, m2; simple, atol1, atol2, rtol, fast, sdeg = sdegNS, offset) 
         degs = info1.degs
         tcond1 = info1.tcond
      elseif mu == 0 && md[isel] == 0
         qtemp = m2 == 0 ? dss(eye(T1,p); Ts) : [eye(T1,p) sysm.sys[isel]]
         degs = Int[]; tcond1 = 1.
      else
         # compute minimal basis as Q = Q1 = [ I -Gu] and set
         # QR = [ Q R ], where R = [ Gw Ga ]
         qtemp = [ eye(T1,p) dss(sysm.sys[isel].A, sysm.sys[isel].E, [-sysm.sys[isel].B[:,inpu] sysm.sys[isel].B[:,mu+1:mu+m2]],
                   sysm.sys[isel].C, [-sysm.sys[isel].D[:,inpu] sysm.sys[isel].D[:,mu+1:mu+m2]]; Ts)] 
         # perform stabilization if strong detectability has to be enforced
         degs = Int[]; tcond1 = 1.
      end
      infodegs[i] = degs
   
      nvec = size(qtemp,1);   # number of basis vectors
      # check solvability conditions
      nvec == 0 && error("emdsyn: empty nullspace basis - the EMDP is not solvable")

      nq = order(qtemp)             # order of the minimal basis
   
      # set H[i] for checking the solvability condition
      emptyHDi = isempty(HDesign[i]);   

      if emptyHDi
         S = trues(nvec,max(1,lfreq),N)
         Htemp = eye(T1,nvec)
         noutmax = nvec
      else
         mh, nh = size(HDesign[i])
         S = trues(mh,max(1,lfreq),N)
         degs = Int[];
         noutmax = mh
         if nh < nvec
            # pad with zeros: row rank is preserved
            Htemp = [ HDesign[i] zeros(T, mh, nvec-nh) ]
         else
            # remove trailing columns if necessary: rank may drop
            Htemp = HDesign[i][:,1:nvec];
            if nh > nvec && mh != rank(Htemp)
               error("The leading part of HDesign[$i] must have full row rank")
            end
         end
      end
      #strongMD && (qtemp = glcf(qtemp; atol1, atol2, rtol, fast, sdeg = sdegNS)[1])
      defer = !emptyHDi && minimal
      for j = 1:N
         # check j-th model detectability
         if isel != j
            temp = gir(Htemp*qtemp[:,1:p+mu]*[sysm.sys[j]; eye(T1,mu,m[j])]; atol)
            if strongMD
               if emdtest
                  St = fdisspec_(temp[:,1:mu+md[j]], MDfreq; FDGainTol = MDGainTol, atol1, atol2, rtol = 0, fast, stabilize = false)[1]
               else
                  St = fdisspec_(temp[:,1:mu], MDfreq; FDGainTol = MDGainTol, atol1, atol2, rtol = 0, fast, stabilize = false)[1]
               end
               Sj = trues(size(St,1),lfreq) 
               for k = 1:lfreq
                   Sj[:,k] = any(view(St,:,:,k),dims=2)
               end
               if any(maximum(Sj,dims=1)) 
                  S[:,:,j] = Sj
               else
                  if defer
                     @warn "Checking of strong solvability condition #$i for model #$j deferred"
                  else
                    emdtest ? error("Strong extended model detection not feasible for model #$j") : 
                              error("Strong model detection not feasible for model #$j ")
                  end
               end
            else
               if emdtest
                  Sj = fditspec_(temp[:,1:mu+md[j]]; atol1, atol2, rtol, FDtol = MDtol)
               else
                  Sj = fditspec_(temp[:,1:mu]; atol1, atol2, rtol, FDtol = MDtol)
               end
               if any(Sj) 
                  S[:,:,j] = maximum(Sj,dims=2)
               else
                  if defer
                     @warn "Checking of solvability condition #$i for model #$j deferred"
                  else
                    emdtest ? error("Extended model detection not feasible for model #$j") : 
                              error("Model detection not feasible for model #$j")
                  end
               end
            end
         end
      end
      # determine the rank of Rwii
      inpw = mu+md[isel]+1:m[isel]
      temp = gir(Htemp*qtemp[:,1:p]*sysm.sys[isel][:,inpw]; atol)
      rwi = gnrank(temp; atol)

      # setup the number of filter outputs
      if ismissing(rdimtarget)
         if emptyHDi 
            rdim = minimal ? 1 : nvec   
         else
            rdim = size(HDesign[i],1)
         end
      else
         rdim = min(rdimtarget[i],nvec) 
      end
      
      # adjust rmin to be at most the rank of Rwii
      rwi > 0 && (rdim = min(rdim,rwi)) 
       
      # Step 2): compute admissible Q2 to reduce the order of Q2*Q;  
      # update Q <- Q2*Q, R <- Q2*R 
   
      # reorder degs to correspond to the expected orders of basis vectors 
      # corresponding to the actual order of outputs of QR 
      reverse!(degs) 
      if rdim < nvec 
         # determine possible low order syntheses using nout >= rmin basis vectors
         # and the corresponding expected orders    
            
         finish = false    # set termination flag
         nout = rdim       # initialize number of selected basis vectors

         while !finish     
             # choose nout basis vectors, which potentially lead to a least order
             # filter with rdim outputs:
             # basesel(i,:) contains the indices of candidate basis vectors;
             # ordsel(i)    contains the presumably achievable least orders
             basesel, ordsel = emdbasesel(S, degs, rdim, nout, simple) 
             #
             # update the synthesis using the selections of candidate vector(s),
             # starting with the least (potentially) achievable order
             for ii = 1:size(basesel,1)
                 baseind = basesel[ii] # indices of current basis selection
                 if rdim == nout
                     hbase = eye(rdim)
                 else
                     hbase = rand(rdim,nout) 
                 end
                 ip = [baseind; setdiff(1:nvec,baseind)][:]
                 if simple && !isempty(degs)
                    if minimal
                       # determine Q2i such that Q2i*Q1i is admisible and 
                       # has least order
                       if emptyHDi 
                          # select vectors and elliminate unobservable dynamics  
                          noelim = falses(nq) 
                          #ell = sum(degs(1:basesel(i,1)-1)); 
                          ell = sum(degs[1:basesel[ii][1]-1]); 
                          for jj = 1:nout 
                              ellnext = sum(degs[1:baseind[jj]]);
                              noelim[ell+1:ellnext] .= true;
                              ell = ellnext
                          end
                       end
                       if rdim == nout
                          if emptyHDi 
                             ir = noelim
                             Ar, Er, Br, Cr, Dr = dssdata(qtemp[baseind,:])
                             qtest = dss(view(Ar,ir,ir), Er == I ? I : view(Er,ir,ir),view(Br,ir,:),view(Cr,:,ir), Dr; Ts)
                             h = Htemp[ip[1:rdim],:];
                          else
                             qtest = gir(Htemp*qtemp; atol1, atol2, rtol);
                          end
                       else
                          # this case is possible only if HDesign[i] is empty
                          # build rdim linear combinations of the first nout vectors 
                          ir = noelim
                          Ar, Er, Br, Cr, Dr = dssdata(qtest[baseind,:])
                          qtest = hbase*dss(view(Ar,ir,ir),Er == I ? I : view(Er,ir,ir),view(Br,ir,:),view(Cr,:,ir), Dr; Ts)
                          h = [ hbase zeros(T1, rdim, nvec-nout) ]
                          h = h[:,ip]  # permute columns to match unpermuted QR 
                       end
                    else
                       # determine Q2i = Hi such that Q2i*Q1i is admisible 
                       if rdim == nout
                          if emptyHDi 
                             h = Htemp[ip[1:rdim],:]
                             qtest = gir(qtemp[baseind,:]; atol1, atol2, rtol, infinite = false)
                          else
                             qtest = gir(Htemp*qtemp; atol1, atol2, rtol, infinite = false)
                          end
                       else
                          # this case is possible only if HDesign[i] is empty
                          # build rdim linear combinations of the first nout vectors 
                          qtest = gir(hbase*qtemp[baseind,:]; atol1, atol2, rtol, infinite = false) 
                          h = [ hbase zeros(T1, rdim, nvec-nout) ] 
                          h = h[:,ip]  # permute columns to match unpermuted QR 
                       end
                    end
                 else
                    if minimal
                       # determine Q2i such that Q2i*Q1i is admisible and 
                       # has least order
                       if rdim == nout
                          if emptyHDi 
                             qtest, _, info2 = glmcover1(qtemp[ip,:], rdim; atol1, atol2, rtol)
                             if !isempty(ordsel) && (order(qtest) != ordsel[ii])
                                @warn "amdsyn: expected reduced order not achieved"
                             end
                             h = Htemp[ip[1:rdim],:]
                          else
                             qtest, _, info2 = glmcover1([Htemp; eye(nvec)]*qtemp[ip,:], rdim; atol1, atol2, rtol)
                          end
                       else  
                          # this case is possible only if HDesign is empty
                          # build rdim linear combinations of the first nout vectors 
                          h = [ hbase zeros(T1, rdim, nvec-nout) ]; 
                          qtest, _, info2 = glmcover1([h; eye(nvec)]*qtemp[ip,:], rdim; atol1, atol2, rtol)
                          h = h[:,ip]  # permute columns to match unpermuted QR 
                        end
                    else
                       # determine Q2i = Hi such that Q2i*Q1i is admisible 
                       if rdim == nout
                          if emptyHDi
                             h = Htemp[ip[1:rdim],:]
                             qtest = gir(qtemp[baseind,:]; atol1, atol2, rtol, infinite = false) 
                          else
                             qtest = gir(Htemp*qtemp; atol1, atol2, rtol, infinite = false) 
                          end
                       else
                          qtest = gir(hbase*qtemp[baseind,:]; atol1, atol2, rtol, infinite = false) 
                          h = [ hbase zeros(T1, rdim, nvec-nout) ]; 
                          h = h[:,ip]  # permute columns to match unpermuted QR 
                       end
                    end
                 end
                 # check model detectability of the current design; 
                 if (rdim == nout && minimal) || rdim < nout
                    notOK = false
                    for j = 1:N
                        # check j-th model detectability
                        if isel != j
                           # dismiss design if check fails
                           temp = gir(qtest[:,1:p+mu]*[sysm.sys[j]; eye(T1,mu,m[j])]; atol1, atol2, rtol)
                           if strongMD
                              if emdtest
                                 St = fdisspec_(temp[:,1:mu+md[j]], MDfreq; FDGainTol = MDGainTol, atol1, atol2, rtol = 0, fast)[1]
                              else
                                 St = fdisspec_(temp[:,1:mu], MDfreq; FDGainTol = MDGainTol, atol1, atol2, rtol = 0, fast)[1]
                              end
                              Sj = trues(size(St,1),lfreq); 
                              for k = 1:lfreq
                                  Sj[:,k] = any(view(St,:,:,k),dims=2)
                              end
                              notOK = !any(maximum(Sj,dims=1)) 
                           else
                              if emdtest
                                 Sj = fditspec_(temp[:,1:mu+md[j]]; atol1, atol2, rtol, FDtol = MDtol)
                              else
                                 Sj = fditspec_(temp[:,1:mu]; atol1, atol2, rtol, FDtol = MDtol)
                              end
                              notOK = !any(Sj) 
                           end
                        end
                        notOK && break
                    end
                    if !notOK 
                       if !simple && minimal
                          # adjust condition number of employed transformations
                          tcond1 = max(tcond1, info2.fnorm, info2.tcond)
                          tcond1 > tcond && 
                             @warn "amdsyn: possible loss of numerical stability due to ill-conditioned transformations"
      #                     if ~emptyHD
      #                        info.HDesign = Htemp;
      #                     end
                       end
                       qtemp = qtest
                       finish = true
                       break
                    end
                 else
                    qtemp = qtest
                    finish = true
                    break
                 end
             end
             if !finish
                if emptyHDi
                   nout += 1
                   if nout > noutmax 
                      finish = true
                      @warn "model detectability not achieved with the chosen number of residuals"
                   end
                else
                  qtemp = gir(Htemp*qtemp; atol1, atol2, rtol);
                  minimal && 
                     (@warn "model detection with least order not feasible for the given tolerances")
                  finish = true;
                end
             end
         end
         # set Q[i] = Q2i*Q1i 
         Qt[i] = qtemp 
         if emptyHDi
            Htemp = h;
         end
      else
         hbase = eye(T1,rdim)
         baseind = 1:nvec;
         # set Q[i] = Q2i*Q1i 
         h = eye(T1,rdim)
         if emptyHDi
            Qt[i] = qtemp 
            Htemp = h
         else
            Qt[i] = Htemp*qtemp 
         end
      end

      # compute M such that M*Q[i] has a desired stability degree;  
      # update Q[i] <- M*Qi[i] 
      k = 1;
      if simple && isequal(hbase,I) && emptyHDi && rwi == 0
         # exploit the block diagonal structure of basis matrices al and cl
         # to compute block-diagonal Q3i
         al, el, bl, cl, dl, = dssdata(Qt[i])
         for ii = 1:length(baseind)
             blkord = degs[baseind[ii]]
             if blkord
                i1 = k:k+blkord-1; 
                Qi = glcf(dss(al[i1,i1],el[i1,i1],bl[i1,:],cl[ii:ii,i1],dl[ii:ii,:];Ts); 
                             atol1, atol2, atol3, sdeg, smarg, evals = poles)[1]
                al[i1,i1] = Qi.A; bl[i1,:] = Qi.B;  cl[ii,i1] = Qi.C; dl[ii,:] = Qi.D 
                el[i1,i1] = (Qi.e == I) ? eye(T1,blkord) : Qi.E  
                k += blkord
             end
         end
         Qt[i] = dss(al, el, bl, cl, dl; Ts)
      else
         rwi == 0 && (Qt[i] = glcf(Qt[i]; atol1, atol2, atol3, sdeg, smarg, evals = poles)[1])
      end
      
      for j = 1:N
         temp = gir(Qt[i][:,1:p+mu]*[sysm.sys[j]; eye(T1,mu,m[j])]; atol1, atol2, rtol)
         # compute gains
         if isel != j
            if strongMD
               gains = zero(T1)
               if emdtest
                  for k = 1:lfreq
                      gains = max(gains,norm(evalfr(temp[:,1:mu+md[j]]; fval = MDfreq[k])))
                  end
               else
                  for k = 1:lfreq
                      gains = max(gains,opnorm(evalfr(temp[:,1:mu]; fval = MDfreq[k])))
                  end
               end
               distinf[i,j] = gains 
            else
               if emdtest
                  distinf[i,j] = glinfnorm(temp[:,1:mu+md[j]]; rtolinf, atol1, atol2, rtol, fast, offset)[1]
               else
                  distinf[i,j] = glinfnorm(temp[:,1:mu]; rtolinf, atol1, atol2, rtol, fast, offset)[1]
               end
            end
         else
            distinf[i,isel] = 0
         end
         Rt[i,j] = temp
      end
      if !normalize 
         # scale to unit minimum gains
         scale = 1/minimum(view(distinf,i,[1:isel-1; isel+1:N]))
      else
         # scale to ensure distinf(i,1) = distinf(1,i)
         if isel == 1
            scale = 1/minimum(view(distinf,1,2:N))
         else
            scale = distinf[1,isel]/distinf[i,1]
         end
      end
      lmul!(scale,view(distinf,i,:))
      Qt[i]  = scale*gss2ss(Qt[i]; atol1, atol2, rtol)[1]
      for j = 1:N
          Rt[i,j] = scale*gss2ss(Rt[i,j]; atol1, atol2, rtol)[1]
      end

      infotcond[i] = tcond1
      HDesign1[i] = emptyHDi ? convert(Matrix{T1},h) : convert(Matrix{T1},Htemp)
   end

   gaps = fill(Inf,M) 
   nstd  =  falses(M)  
   if any(mw .> 0)
      for i = 1:M
          isel = MDSelect[i]
          mw[isel] == 0 && continue
         
          # determine the optimal factors Q3 to minimize the i-th gap 
  
          # compute the co-outer-co-inner factorizations of 
          # Rwii = [Rwoe 0]*Rwi = Rwoe*Rwi1
          inpw = mu+md[isel]+1:m[isel]
          Rwi, Rwo, info1 = goifac(Rt[i,i][:,inpw]; atol1, atol2, atol3, rtol)
          rw = size(Rwo,2);  # rank of Rwii
          if rw > 0
             nonstandard = info1.nfuz+info1.niuz > 0 
             nstd[i] = nonstandard
             # handle different cases
             if nonstandard && nonstd != 1
                # non-standard case: zeros on the boundary of the stability domain 
                if nonstd == 2
                   # use modified co-outer-co-inner factorization
                     Rwo = goifac([Rwo epsreg*eye(size(Rwo,1))]; atol1, atol2, atol3, rtol)[2]
                elseif nonstd == 3  
                   # use Wiener-Hopf type co-outer-co-inner factorization
                   # separate stable and strictly unstable zeros Rwo = Rwoe*Rwouz; 
                   Rwouz, Rwoe = gcrange(Rwo; atol1, atol2, rtol, zeros = "s-unstable")
                   Rwo = Rwoe*glinfnorm(Rwouz*Rwi[1:rw,:])[1] 
                elseif nonstd == 4  
                   # use modified Wiener-Hopf type co-outer-co-inner factorization
                   # with zero shifting of the non-minimum phase factor
                   # separate stable and strictly unstable zeros Rwo = Rwoe*Rwouz
                   # and update Rwo <- Rwoe
                   Rwouz, Rwoe = gcrange(Rwo; atol1, atol2, rtol, zeros = "s-unstable")
                   # set suitable bilinear transformation
                   if disc
                     sys1 = rtf('z')/sdegzer
                   else
                     s = rtf('s'); sys1 = rtf(1)
                     if info1.nfuz > 0
                        sys1 = sys1*(s-sdegzer);
                     end
                     if info1.niuz > 0
                        sys1 = sys1/(1-sdegzer*s);
                     end
                   end
                   # form shifted factor 
                   Rwo = Rwoe*gbilin(Rwouz, sys1; minimal = true, atol1, atol2, rtol)[1] 
                elseif nonstd == 5  
                   # use modified Wiener-Hopf type co-outer-co-inner factorization
                   # with regularization of the non-minimum phase factor
                   # separate stable and strictly unstable zeros Rwo = Rwoe*Rwouz 
                   # and update Rwo <- Rwoe
                   Rwouz, Rwoe = gcrange(Rwo; atol1, atol2, rtol, zeros = "s-unstable")
                   rw = size(Rwouz,1)
                   Rwo = Rwoe*goifac([Rwouz epsreg*eye(rw)]; atol1, atol2, atol3, rtol)[2]
                end
             end
             #   continue editing from here
             if rw == size(Rwo,1)
                # Q3 = inv(Rwo)
                Qt[i] = gminreal(Rwo\Qt[i]; atol1, atol2, rtol)
             else
                # regularization for non-invertible Rwo (this case should never occur)
                # Q3 = Rwoinv, where Rwoinv is a left inverse of Rwo
                # with stable spurious poles
                Rwoinv = glsol([Rwo;eye(rw)], rw; atol1, atol2, rtol, sdeg)[1]
                # form Qi= Rwoinv*[Rwo; I] 
                Qt[i] = gir(Rwoinv*Qt[i]; atol1, atol2, rtol, infinite = false)
             end
             if nonstandard && nonstd == 1
               # perform stabilization 
               # determine Q4 such that Q <- Q4*Q and R <- Q4*R are stable
               Qt[i] = glcf(Qt[i]; atol1, atol2, atol3, rtol, sdeg, smarg, evals = poles)[1]
            end
         
            for j = 1:N
               temp = gir(Qt[i][:,1:p+mu]*[sysm.sys[j]; eye(T1,mu,m[j])]; atol1, atol2, rtol)
               # compute gains
               if isel != j
                  if strongMD
                     gains = zero(T1)
                     if emdtest
                        for k = 1:lfreq
                            gains = max(gains,norm(evalfr(temp[:,1:mu+md[j]]; fval = MDfreq[k])))
                        end
                     else
                        for k = 1:lfreq
                            gains = max(gains,opnorm(evalfr(temp[:,1:mu]; fval = MDfreq[k])))
                        end
                     end
                     distinf[i,j] = gains 
                  else
                     if emdtest
                        distinf[i,j] = glinfnorm(temp[:,1:mu+md[j]]; rtolinf, atol1, atol2, rtol, fast, offset)[1]
                     else
                        distinf[i,j] = glinfnorm(temp[:,1:mu]; rtolinf, atol1, atol2, rtol, fast, offset)[1]
                     end
                  end
               else
                  distinf[i,isel] = 0
               end
               Rt[i,j] = temp
            end

             if !normalize
                # scale to unit minimum gains
                scale = 1/minimum(view(distinf,i,[1:isel-1; isel+1:N]))
             else
                # scale to ensure distinf(i,1) = distinf(1,i)
                if i == 1
                   scale = 1/minimum(view(distinf,1,2:N))
                else
                   scale = distinf[1,isel]/distinf[i,1]
                end
             end
             lmul!(scale,view(distinf,i,:))
   
             # scale and transform to standard state-space
             Qt[i]  = scale*gss2ss(Qt[i]; atol1, atol2, rtol)[1]
             for j = 1:N
                Rt[i,j] = scale*gss2ss(Rt[i,j]; atol1, atol2, rtol)[1]
             end
             beta = minimum(view(distinf, i, [1:i-1; i+1:N]))
             gaps[i] = beta/ghinfnorm(Rt[i,i][:,inpw]; rtolinf, atol1, atol2, rtol, fast, offset)[1]
          end
      end
   end    
   
   Q = MDFilter{T1}(Qt, p, mu)
   R = MDFilterIF{T1}(Rt, mu, md, mw, ma) 
   info = (tcond = infotcond, degs = infodegs, HDesign = HDesign1, MDperf = distinf, MDgap = gaps, nonstandard = nstd)

   return Q, R, info
   
   # end AMDSYN
end
"""
    emdsyn(sysm::MDMModel; rdim, MDSelect, nullspace = true, simple = false, minimal = true, 
                           emtest = false, normalize = false, fast = true, 
                           sdeg, smarg, poles, HDesign, MDtol, MDGainTol, MDfreq, 
                           tcond, offset, atol, atol1, atol2, atol3, rtol) 
                           -> (Q::MDFilter, R::MDFilterIF, info)

Solve the _exact model detection problem_ (EMDP) for a given multiple synthesis model
`sysm::MDMModel`. The computed stable and proper filter objects `Q` and `R` contain the 
model detection filter, representing the solution of the EMDP, and its internal form, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.degs`, `info.MDperf`, and `info.HDesign`, 
contains additional synthesis related information (see below). 

For a multiple synthesis model `sysm` containing `N` stable component models, 
the `i`-th component system `sysm.sys[i]` has a descriptor realization of the form 
`sysm.sys[i] = (Ai-λEi,Bi,Ci,Di)` and the corresponding input-output form is

      yi = Gui(λ)*u + Gdi(λ)*di + Gwi(λ)*wi + Gvi(λ)*vi 

with the Laplace- or Z-transformed plant outputs `yi`, control inputs `u`, 
disturbance inputs `di`, noise inputs `wi`, and auxiliary inputs `vi`,  
and with `Gui(λ)`, `Gdi(λ)`, `Gwi(λ)` and 
`Gvi(λ)`, the corresponding transfer function matrices. 

The model detection filter object `Q`, contains in `Q.sys` the resulting bank of `N` filters. 
The `i`-th filter `Q.sys[i]` is in a standard state-space form and generates `r_i`, 
the `i`-th component (scalar or vector) of the overall residual vector `r := [r_1; r_2; ...; r_N]`. 
The corresponding input-output (implementation) form of the `i`-th filter is

            r_i = Qyi(λ)*y + Qui(λ)*u   ,            

where `Qyi(λ)` and `Qui(λ)` are the transfer function matrices from the output and control inputs to the `i`-th residual component. 
The dimensions of output and control inputs are contained in the integers  
`Q.ny` and `Q.mu`, respectively.

The model detection filter internal form object `R`, contains `R.sys`, the resulting array of `N × N` 
internal form of the filters.  
The `(i,j)`-th component filter `R.sys[i,j]` is in a standard state-space form, 
generates the residual signal `r_ij`, and corresponds to the 
input-output form

       r_ij = Ruij(λ)*u + Rdij(λ)*dj + Rwij(λ)*wj + Rvij(λ)*vj ,

where 

       | Ruij(λ) Rdij(λ) Rwij(λ) Raij(λ) | := |Qyi(λ) Qui(λ)]*| Guj(λ) Gdj(λ) Gwj(λ) Gaj(λ) |. 
                                                              |   I    0      0      0      |

The solution of the EMDP ensures that for the `i`-th filter, `Ruii(λ) = 0`, `Rdii(λ) = 0`, and 
`[Ruij(λ) Rdij(λ)]` ``\\neq`` `0` for `j` ``\\neq`` `i`.

Various user options can be specified via keyword arguments as follows:

`MDSelect = ifilt` specifies in the vector `ifilt` the indices of the desired filters to be designed
(default: `ifilt = 1:N`)

If `minimal = true` (default), least order filter synthesis is performed to determine 
each of the component filters
`Q.sys[i]` for `i = 1, ...,N` and `R.sys[i,j]` for `i, j = 1, ...,N`  while 
with `minimal = false` full order synthesis is performed.  

If `HDesign = H`, then `H` is an `N`-dimensional array of full row rank or empty design matrices `H = [H_1, ..., H_N]`,
where `H_i` is the design matrix employed for the synthesis of the `i`-th component filter (default: `HDesign = missing`)

`rdim = q` specifies the vector `q`, whose `i`-th component `q[i]` specifies 
the number of residual outputs for the `i`-th component filter `Q.sys[i]`. 
If `q` is a scalar, then a vector `rdim` with all components equal to `q` is assumed.
The default value of `q[i]` is chosen as follows: if `HDesign = missing` or `H_i` is empty then  
`q[i] = 1`, if `minimal = true`, or `q[i]` is the number of the nullspace basis 
vectors used for the synthesis of `Q.sys [i]`, if `minimal = false`; 
if  `H_i` specifies a full row rank design matrix, then `q[i]` is the row dimension of `H_i`. 

If `emdtest = false` (default), only the input channels are used for model detectability tests.
If `emdtest = true`, extended model detectability tests are performed using both control and
disturbance input channels.  

`MDfreq = ω` specifies a vector of real frequency values or a scalar real frequency value
for strong model detectability checks (default: `MDfreq = missing`).

If `nullspace = false` (default),  a full-order observer based nullspace basis is 
used for the synthesis of the `i`-th filter whenever the `i`-th component system has no disturbance inputs.
Otherwise, minimal proper nullspace bases are used. 
If `nullspace = true`, minimal proper nullspace bases are used for the synthesis of the
model detection filters. 

If `simple = true`, simple proper nullspace bases are emplyed for synthesis. 
The orders of the basis vectors employed for the synthesis of the `i`-th filter
are provided in `info.deg[i]`. 
If `simple = false` (default), then no simple bases are computed. 

`offset = β` specifies the boundary offset `β` to assess the stability of poles. 
Accordingly, for the stability of a continuous-time system all real parts of poles must be at most `-β`, 
while for the stability of a discrete-time system all moduli of poles must be at most `1-β`. 
The default value used for `β` is `sqrt(ϵ)`, where `ϵ` is the working machine precision. 

`smarg = α` specifies the stability margin which defines the stability 
domain `Cs` of poles, as follows: 
for a continuous-time system, `Cs` is the set of complex numbers 
with real parts at most `α`, 
while for a discrete-time system, `Cs` is the set of complex numbers with 
moduli at most `α < 1` (i.e., the interior of a disc of radius `α` centered in the origin). 
If `smarg` is missing, then the employed default values are `α = -β` 
for a continuous-time system and `α = 1-β` for a discrete-time system, 
where `β` is the boundary offset specified by the keyword argument `offset = β`. 

`sdeg = γ` is the prescribed stability degree for the poles of the filters `Q` and `R` 
(default: `γ = -0.05` for the real parts of poles for a continuous-time system and
`γ = 0.95` for the magnitudes of poles for a discrete-time system). 

`poles = v` specifies a complex vector `v` containing a complex conjugate set  
of desired poles within the stability domain `Cs` to be assigned for the filters `Q` and `R`
(default: `poles = missing`).

`tcond = tcmax` specifies the maximum alowed condition number `tcmax` 
of the employed non-orthogonal transformations (default: `tcmax = 1.e4`).

`MDtol = tol1` specifies the threshold `tol1` for model detectability checks
   (default: `tol1 = 0.0001`).

`MDGainTol = tol2` specifies the threshold `tol2` for strong model detectability checks
   (default: `tol2 = 0.01`). 

If `normalize = false` (default), the `i`-th component filter `Q.sys[i]` is scaled such that
the minimum gain of `R.sys[i,j]` for `j = 1, ..., N`,   `j` ``\\neq`` `i`, is equal to one. 
If `normalize = true`, the standard normalization of component filters is performed to ensure
equal gains for `R.sys[1,j]` and `R.sys[j,1]`. 
This option is only possible if `ifilt[1] = 1` (see `MDSelect`).

The rank determinations in the performed reductions
are based on rank revealing QR-decompositions with column pivoting 
if `fast = true` (default) or the more reliable SVD-decompositions if `fast = false`.

The keyword arguments `atol1`, `atol2` and `rtol`, specify, respectively, 
the absolute tolerance for the nonzero elements of the state-space matrices 
`Ai`, `Bi`, `Ci` and `Di`, the absolute tolerance for the nonzero elements of `Ei`,   
and the relative tolerance for the nonzero elements of all above matrices.  
The default relative tolerance is `ni*ϵ`, where `ϵ` is the working machine epsilon 
and `ni` is the orders of the system `sysm.sys[i]`. 
The keyword argument `atol3` is an absolute tolerance for observability tests
   (default: internally determined value). 
The keyword argument `atol` can be used 
to simultaneously set `atol1 = atol`, `atol2 = atol` and `atol3 = atol`. 

The resulting named tuple `info` contains `(tcond, degs, MDperf, HDesign) `, where:

`info.tcond` is the maximum of the condition numbers of the employed 
   non-orthogonal transformation matrices; a warning is issued if `info.tcond >= tcmax`;

`info.degs` is an `N`-dimensional vector, whose `i`-th component is an integer vector 
containing the degrees of the basis vectors of the employed simple
nullspace basis for the synthesis of the `i`-th filter component, if `simple = true`, 
and the degrees of the
basis vectors of an equivalent polynomial nullspace basis, if `simple = false`;

`info.MDperf` is an `N × N` array containing the resulting distance-mapping performance, 
representing  the peak gains of the associated internal representations of 
the model detection filters. If the `(i,j)`-th component filter `R.sys[i,j]` 
has the input-output form
 
     rij = Ruij(λ)*u + Rdij(λ)*dj + Rwij(λ)*wj + Rvij(λ)*vj ,

then, the `(i,j)`-th performance gain is evaluated as 
`info.MDperf[i,j] = ||Ruij(λ)||∞` if `emdtest = false` (default) or 
`info.MDperf[i,j] = ||[Ruij(λ) Rdij(λ)]||∞` if `emdtest = true`.  
If `MDfreq = ω`, then each gain `info.MDperf[i,j]` represents
the maximum of 2-norm pointwise gains evaluated for all frequencies in `ω`.

`info.HDesign` is an `N`-dimensional vector, whose `i`-th component 
is the design matrix `H_i` employed for the synthesis of 
the `i`-th model detection filter.

   
_Method:_ The Procedure EMD from [1] is implemented to solve 
the exact model detection problem. For more details on 
the least order synthesis of model detection filters see [2].

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017; sec. 5.2.

[2] A. Varga, Least order fault and model detection using multi-models. 
       IEEE CDC'09, Shanghai, China, 2009.
"""
function emdsyn(sysm::MDMModel; rdim::Union{Vector{Int},Int,Missing} = missing, 
                      MDSelect::Vector{Int} = Vector(1:length(sysm.sys)), poles::Union{AbstractVector,Missing} = missing, 
                      sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing,  
                      nullspace::Bool = false, minimal::Bool = true, simple::Bool = false, emdtest::Bool = false,
                      MDtol::Real = 0.0001, MDGainTol::Real = 0.01, MDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                      tcond::Real = 1.e4, HDesign::Union{Vector{<:Matrix},Matrix,Missing} = missing, normalize::Bool = false,
                      offset::Real = sqrt(eps(Float64)), atol::Real = zero(Float64), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                      rtol::Real = (size(sysm.sys[1].A,1)+1)*eps(Float64)*iszero(max(atol1,atol2)), 
                      fast::Bool = true) 
   N = length(sysm.sys)   
   T1 = Float64
   # check stability of input model
   isstable(sysm) || error("the input multiple model must be stable")

   p = size(sysm.sys[1],1);       # number of measurable outputs              
   Ts = sysm.sys[1].Ts                  
   disc = (Ts != 0);  # system type (continuous- or discrete-time)
  
   # decode options

   M = length(MDSelect)
   any(MDSelect .<= 0) && error("the selected filter indices must be positive")   
   any(MDSelect .> N) && error("the selected filter indices must be at most $N")   
   M > 0 && MDSelect[1] != 1 && (normalize = false)

     
   strongMD = !ismissing(MDfreq)
   strongMD && !isa(MDfreq,Vector) && (MDfreq = [MDfreq]) 
   lfreq = strongMD ? length(MDfreq) : 0
   
   # set default stability degree
   sdegdefault = disc ? 0.95 : -0.05
   
  
   # stability margin
   ismissing(smarg) && (smarg = disc ? 1-offset : -offset)  # set default stability margin
   
   poles_nomissing = !ismissing(poles)
    
   # sort desired poles
   if poles_nomissing 
      tempc = poles[imag.(poles) .> 0]
      if !isempty(tempc)
         tempc1 = conj(poles[imag.(poles) .< 0])
         isequal(tempc[sortperm(real(tempc))],tempc1[sortperm(real(tempc1))]) ||
                 error("poles must be a self-conjugated complex vector")
      end
      # check that all eigenvalues are inside of the stability region
      ( ((disc && any(abs.(poles) .> 1-offset) )  || (!disc && any(real.(poles) .> -offset)))  &&
            error("The elements of poles must lie in the stability region of interest") )
   end   
   
   ismissing(sdeg) && ismissing(poles) && (sdeg = sdegdefault)  # set desired stability degree to default value
  
   # desired number of filter outputs   
   if !ismissing(rdim)
      if isa(rdim,Vector) 
         length(rdim) == M || error("dimension of rdim must be equal to the number of selected filters $M")
         minimum(rdim) > 0 || error("all components of rdim must be positive")
      else
         rdim > 0 || error("rdim must be positive")
         rdim = fill(rdim, M)
      end
   end

   # imposed design option to form linear combinations of basis vectors
   emptyHD = ismissing(HDesign)

   if !emptyHD
       if isa(HDesign,Matrix)
          size(HDesign,1) == rank(HDesign) || error("HDesign must have full row rank")
          HDesign = fill(HDesign,M)
       else
          size(HDesign,1) == M || error("number of HDesign components must be equal to the number of selected filters $M")
       end
       if !ismissing(rdim) 
          for i = 1:M
              if !isempty(HDesign[i])
                 mH = size(HDesign[i],1)
                 if mH > 0
                    mH == rdim[i] || error("row dimension of HDesign[$i] must be equal to rdim[$i]")
                    mH == rank(HDesign[i]) || error("HDesign[$i] must have full row rank")
                 end
              end
          end
       end
   else
      HDesign = fill(T1[],M)
   end

 
   # decode input information
   mu = sysm.mu  
   inpu = 1:mu 
   md = sysm.md 
   mw = sysm.mw 
   ma = sysm.ma 
   
   m = mu .+ (md+mw+ma)      # total number of inputs
   infotcond = similar(Vector{T1},M)
   infodegs = similar(Vector{Vector{Int}},M)
   HDesign1 = similar(Vector{Array{T1,2}},M)
   Qt = similar(Vector{DescriptorStateSpace{T1}},M)
   Rt = similar(Matrix{DescriptorStateSpace{T1}},M,N)
   distinf = fill(-1.,M,N) 
   ismissing(rdim) ? rdimtarget = missing : rdimtarget = copy(rdim)
   for i = 1:M
      isel = MDSelect[i]
      # solve the exact MD problem for the i-th system (with noise inputs)
      # 
      # compute a minimal left nullspace basis Q1i for the i-th system
      desc = (sysm.sys[isel].E != I)
      m2 = mw[isel]+ma[isel]
      sdegNS = strongMD ? sdegdefault : missing     
      if nullspace || simple || md[isel] > 0 || (desc && rcond(sysm.sys[isel].E) < 1.e-7 )
         # form [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
         syse = [sysm.sys[isel]; eye(mu,m[isel])];
         #
         # compute a left nullspace basis Q = Q1 of G1 = [Gu Gd; I 0] = 0 and
         # obtain QR = [ Q R ], where R = [ Rw Ra] = Q*[Gw Ga;0 0 0]
         qtemp, info1 = glnull(syse, m2; simple, atol1, atol2, rtol, fast, sdeg = sdegNS, offset) 
         degs = info1.degs
         tcond1 = info1.tcond
      elseif mu == 0 && md[isel] == 0
         # compute minimal basis as Q = Q1 = [ I ] and set
         # QR = [ Q R ], where R = [ Gw Ga ]
         qtemp = m2 == 0 ? dss(eye(T1,p); Ts) : [eye(T1,p) sysm.sys[isel]]
         degs = Int[]; tcond1 = 1.
      else
         # compute minimal basis as Q = Q1 = [ I -Gu] and set
         # QR = [ Q R ], where R = [ Gw Ga ]
         qtemp = [ eye(T1,p) dss(sysm.sys[isel].A, sysm.sys[isel].E, [-sysm.sys[isel].B[:,inpu] sysm.sys[isel].B[:,mu+1:mu+m2]],
                   sysm.sys[isel].C, [-sysm.sys[isel].D[:,inpu] sysm.sys[isel].D[:,mu+1:mu+m2]]; Ts)] 
         # perform stabilization if strong detectability has to be enforced
         degs = Int[]; tcond1 = 1.
      end
      infodegs[i] = degs
   
      nvec = size(qtemp,1);   # number of basis vectors
      # check solvability conditions
      nvec == 0 && error("emdsyn: empty nullspace basis - the EMDP is not solvable")

      nq = order(qtemp)             # order of the minimal basis
   
      # set H[i] for checking the solvability condition
      emptyHDi = isempty(HDesign[i]);   

      if emptyHDi
         S = trues(nvec,max(1,lfreq),N)
         Htemp = eye(nvec)
         noutmax = nvec
      else
         mh, nh = size(HDesign[i])
         S = trues(mh,max(1,lfreq),N)
         degs = Int[];
         noutmax = mh
         if nh < nvec
            # pad with zeros: row rank is preserved
            Htemp = [ HDesign[i] zeros(T, mh, nvec-nh) ]
         else
            # remove trailing columns if necessary: rank may drop
            Htemp = HDesign[i][:,1:nvec];
            if nh > nvec && mh != rank(Htemp)
               error("The leading part of HDesign[$i] must have full row rank")
            end
         end
      end

      #strongMD && (qtemp = glcf(qtemp; atol1, atol2, rtol, fast, sdeg = sdegNS)[1])
      defer = !emptyHDi && minimal
       for j = 1:N
         # check j-th model detectability
         if isel != j
            temp = gir(Htemp*qtemp[:,1:p+mu]*[sysm.sys[j]; eye(T1,mu,m[j])]; atol)
            if strongMD
               if emdtest
                  St = fdisspec_(temp[:,1:mu+md[j]], MDfreq; FDGainTol = MDGainTol, atol1, atol2, rtol = 0, fast, stabilize = false)[1]
               else
                  St = fdisspec_(temp[:,1:mu], MDfreq; FDGainTol = MDGainTol, atol1, atol2, rtol = 0, fast, stabilize = false)[1]
               end
               Sj = trues(size(St,1),lfreq) 
               for k = 1:lfreq
                   Sj[:,k] = any(view(St,:,:,k),dims=2)
               end
               if any(maximum(Sj,dims=1)) 
                  S[:,:,j] = Sj
               else
                  if defer
                     @warn "Checking of strong solvability condition #$i for model #$j deferred"
                  else
                    emdtest ? error("Strong extended model detection not feasible for model #$j") : 
                              error("Strong model detection not feasible for model #$j ")
                  end
               end
            else
               if emdtest
                  Sj = fditspec_(temp[:,1:mu+md[j]]; atol1, atol2, rtol, FDtol = MDtol)
               else
                  Sj = fditspec_(temp[:,1:mu]; atol1, atol2, rtol, FDtol = MDtol)
               end
               if any(Sj) 
                  S[:,:,j] = maximum(Sj,dims=2)
               else
                  if defer
                     @warn "Checking of solvability condition #$i for model #$j deferred"
                  else
                    emdtest ? error("Extended model detection not feasible for model #$j") : 
                              error("Model detection not feasible for model #$j")
                  end
               end
            end
         end
      end
      if ismissing(rdimtarget)
         if emptyHDi 
            rdim = minimal ? 1 : nvec   
         else
            rdim = size(HDesign[i],1)
         end
      else
         rdim = min(rdimtarget[i],nvec) 
      end
      
          
      # Step 2): compute admissible Q2 to reduce the order of Q2*Q;  
      # update Q <- Q2*Q, R <- Q2*R 
   
      # reorder degs to correspond to the expected orders of basis vectors 
      # corresponding to the actual order of outputs of QR 
      reverse!(degs) 
      if rdim < nvec 
         #if rdim < noutmax && mf > 0
            # determine possible low order syntheses using i >= rmin basis vectors
            # and the corresponding expected orders    
            
            finish = false    # set termination flag
            nout = rdim       # initialize number of selected basis vectors
            # if ~simple && minimal
            #    QR = xperm(QR,nq:-1:1);  # permute states to speedup glmcover1
            # end
            while !finish     
                # choose nout basis vectors, which potentially lead to a least order
                # filter with rdim outputs:
                # basesel(i,:) contains the indices of candidate basis vectors;
                # ordsel(i)    contains the presumably achievable least orders
                basesel, ordsel = emdbasesel(S, degs, rdim, nout, simple) 
                #
                # update the synthesis using the selections of candidate vector(s),
                # starting with the least (potentially) achievable order
                for ii = 1:size(basesel,1)
                    baseind = basesel[ii] # indices of current basis selection
                    if rdim == nout
                        hbase = eye(rdim)
                    else
                        hbase = rand(rdim,nout) 
                    end
                    ip = [baseind; setdiff(1:nvec,baseind)][:]
                    if simple && !isempty(degs)
                       if minimal
                          # determine Q2i such that Q2i*Q1i is admisible and 
                          # has least order
                          if emptyHDi 
                             # select vectors and elliminate unobservable dynamics  
                             noelim = falses(nq) 
                             #ell = sum(degs(1:basesel(i,1)-1)); 
                             ell = sum(degs[1:basesel[ii][1]-1]); 
                             for jj = 1:nout 
                                 ellnext = sum(degs[1:baseind[jj]]);
                                 noelim[ell+1:ellnext] .= true;
                                 ell = ellnext
                             end
                          end
                          if rdim == nout
                             if emptyHDi 
                                ir = noelim
                                Ar, Er, Br, Cr, Dr = dssdata(qtemp[baseind,:])
                                qtest = dss(view(Ar,ir,ir), Er == I ? I : view(Er,ir,ir),view(Br,ir,:),view(Cr,:,ir), Dr; Ts)
                                h = Htemp[ip[1:rdim],:];
                             else
                                qtest = gir(Htemp*qtemp; atol1, atol2, rtol);
                             end
                          else
                             # this case is possible only if HDesign[i] is empty
                             # build rdim linear combinations of the first nout vectors 
                             ir = noelim
                             Ar, Er, Br, Cr, Dr = dssdata(qtest[baseind,:])
                             qtest = hbase*dss(view(Ar,ir,ir),Er == I ? I : view(Er,ir,ir),view(Br,ir,:),view(Cr,:,ir), Dr; Ts)
                             h = [ hbase zeros(T1, rdim, nvec-nout) ]
                             h = h[:,ip]  # permute columns to match unpermuted QR 
                          end
                       else
                          # determine Q2i = Hi such that Q2i*Q1i is admisible 
                          if rdim == nout
                             if emptyHDi 
                                h = Htemp[ip[1:rdim],:]
                                qtest = gir(qtemp[baseind,:]; atol1, atol2, rtol, infinite = false)
                             else
                                qtest = gir(Htemp*qtemp; atol1, atol2, rtol, infinite = false)
                             end
                          else
                             # this case is possible only if HDesign[i] is empty
                             # build rdim linear combinations of the first nout vectors 
                             qtest = gir(hbase*qtemp[baseind,:]; atol1, atol2, rtol, infinite = false) 
                             h = [ hbase zeros(T1, rdim, nvec-nout) ] 
                             h = h[:,ip]  # permute columns to match unpermuted QR 
                          end
                       end
                    else
                       if minimal
                          # determine Q2i such that Q2i*Q1i is admisible and 
                          # has least order
                          if rdim == nout
                             if emptyHDi 
                                qtest, _, info2 = glmcover1(qtemp[ip,:], rdim; atol1, atol2, rtol)
                                if !isempty(ordsel) && (order(qtest) != ordsel[ii])
                                   @warn "emdsyn: expected reduced order not achieved"
                                end
                                h = Htemp[ip[1:rdim],:]
                             else
                                qtest, _, info2 = glmcover1([Htemp; eye(nvec)]*qtemp[ip,:], rdim; atol1, atol2, rtol)
                             end
                          else  
                             # this case is possible only if HDesign is empty
                             # build rdim linear combinations of the first nout vectors 
                             h = [ hbase zeros(T1, rdim, nvec-nout) ]; 
                             qtest, _, info2 = glmcover1([h; eye(nvec)]*qtemp[ip,:], rdim; atol1, atol2, rtol)
                             h = h[:,ip]  # permute columns to match unpermuted QR 
                           end
                       else
                          # determine Q2i = Hi such that Q2i*Q1i is admisible 
                          if rdim == nout
                             if emptyHDi
                                h = Htemp[ip[1:rdim],:]
                                qtest = gir(qtemp[baseind,:]; atol1, atol2, rtol, infinite = false) 
                             else
                                qtest = gir(Htemp*qtemp; atol1, atol2, rtol, infinite = false) 
                             end
                          else
                             qtest = gir(hbase*qtemp[baseind,:]; atol1, atol2, rtol, infinite = false) 
                             h = [ hbase zeros(T1, rdim, nvec-nout) ]; 
                             h = h[:,ip]  # permute columns to match unpermuted QR 
                          end
                       end
                    end
                    # check model detectability of the current design; 
                    if (rdim == nout && minimal) || rdim < nout
                       notOK = false
                       for j = 1:N
                           # check j-th model detectability
                           if isel != j
                              # dismiss design if check fails
                              temp = gir(qtest[:,1:p+mu]*[sysm.sys[j]; eye(T1,mu,m[j])]; atol1, atol2, rtol)
                              if strongMD
                                 if emdtest
                                    St = fdisspec_(temp[:,1:mu+md[j]], MDfreq; FDGainTol = MDGainTol, atol1, atol2, rtol = 0, fast)[1]
                                 else
                                    St = fdisspec_(temp[:,1:mu], MDfreq; FDGainTol = MDGainTol, atol1, atol2, rtol = 0, fast)[1]
                                 end
                                 Sj = trues(size(St,1),lfreq); 
                                 for k = 1:lfreq
                                     Sj[:,k] = any(view(St,:,:,k),dims=2)
                                 end
                                 notOK = !any(maximum(Sj,dims=1)) 
                              else
                                 if emdtest
                                    Sj = fditspec_(temp[:,1:mu+md[j]]; atol1, atol2, rtol, FDtol = MDtol)
                                 else
                                    Sj = fditspec_(temp[:,1:mu]; atol1, atol2, rtol, FDtol = MDtol)
                                 end
                                 notOK = !any(Sj) 
                              end
                           end
                           notOK && break
                       end
                       if !notOK 
                          if !simple && minimal
                             # adjust condition number of employed transformations
                             tcond1 = max(tcond1, info2.fnorm, info2.tcond)
                             tcond1 > tcond && 
                                @warn "emdsyn: possible loss of numerical stability due to ill-conditioned transformations"
         #                     if ~emptyHD
         #                        info.HDesign = Htemp;
         #                     end
                          end
                          qtemp = qtest
                          finish = true
                          break
                       end
                    else
                       qtemp = qtest
                       finish = true
                       break
                    end
                end
                if !finish
                   if emptyHDi
                      nout += 1
                      if nout > noutmax 
                         finish = true
                         @warn "model detectability not achieved with the chosen number of residuals"
                      end
                   else
                     qtemp = gir(Htemp*qtemp; atol1, atol2, rtol);
                     minimal && 
                        (@warn "model detection with least order not feasible for the given tolerances")
                     finish = true;
                   end
                end
            end
            # set Q[i] = Q2i*Q1i 
            Qt[i] = qtemp 
            if emptyHDi
               Htemp = h;
            end
         else
            hbase = eye(T1,rdim)
            baseind = 1:nvec;
            # set Q[i] = Q2i*Q1i 
            h = eye(T1,rdim)
            if emptyHDi
               Qt[i] = qtemp 
               Htemp = h
            else
               Qt[i] = Htemp*qtemp 
            end
         end
        
         # compute Q3i such that Q3i*Q[i] has a desired stability degree;  
         # update Q[i] <- Q3i*Qi[i] 
         k = 1;
         if simple && isequal(hbase,I) && emptyHDi
            # exploit the block diagonal structure of basis matrices al and cl
            # to compute block-diagonal Q3i
            al, el, bl, cl, dl, = dssdata(Qt[i])
            for ii = 1:length(baseind)
                blkord = degs[baseind[ii]]
                if blkord
                   i1 = k:k+blkord-1; 
                   Qi = glcf(dss(al[i1,i1],el[i1,i1],bl[i1,:],cl[ii:ii,i1],dl[ii:ii,:];Ts); 
                                atol1, atol2, atol3, sdeg, smarg, evals = poles)[1]
                   al[i1,i1] = Qi.A; bl[i1,:] = Qi.B;  cl[ii,i1] = Qi.C; dl[ii,:] = Qi.D 
                   el[i1,i1] = (Qi.e == I) ? eye(T1,blkord) : Qi.E  
                   k += blkord
                end
            end
            Qt[i] = dss(al, el, bl, cl, dl; Ts)
         else
            Qt[i] = glcf(Qt[i]; atol1, atol2, atol3, sdeg, smarg, evals = poles)[1]
         end
     
         for j = 1:N
           temp = gir(Qt[i][:,1:p+mu]*[sysm.sys[j]; eye(T1,mu,m[j])]; atol1, atol2, rtol)
           # compute gains
           if isel != j
              if strongMD
                 gains = zero(T1)
                 if emdtest
                    for k = 1:lfreq
                        gains = max(gains,norm(evalfr(temp[:,1:mu+md[j]]; fval = MDfreq[k])))
                    end
                 else
                    for k = 1:lfreq
                        gains = max(gains,opnorm(evalfr(temp[:,1:mu]; fval = MDfreq[k])))
                    end
                 end
                 distinf[i,j] = gains 
              else
                 if emdtest
                    distinf[i,j] = glinfnorm(temp[:,1:mu+md[j]])[1]
                 else
                    distinf[i,j] = glinfnorm(temp[:,1:mu])[1]
                 end
              end
           else
              distinf[i,isel] = 0
           end
           Rt[i,j] = temp
         end
     
       
         if !normalize 
            # scale to unit minimum gains
            scale = 1/minimum(view(distinf,i,[1:isel-1; isel+1:N]))
         else
            # scale to ensure distinf(i,1) = distinf(1,i)
            if isel == 1
               scale = 1/minimum(view(distinf,1,2:N))
            else
               scale = distinf[1,isel]/distinf[i,1]
            end
         end
         lmul!(scale,view(distinf,i,:))

         # scale and transform to standard state-space
         Qt[i]  = scale*gss2ss(Qt[i]; atol1, atol2, rtol)[1]
         for j = 1:N
            Rt[i,j] = scale*gss2ss(Rt[i,j]; atol1, atol2, rtol)[1]
         end
         infotcond[i] = tcond1
         HDesign1[i] = emptyHDi ? convert(Matrix{T1},h) : convert(Matrix{T1},Htemp)
   end
  
   Q = MDFilter{T1}(Qt, p, mu)
   R = MDFilterIF{T1}(Rt, mu, md, mw, ma) 
   info = (tcond = infotcond, degs = infodegs, HDesign = HDesign1, MDperf = distinf)

   return Q, R, info
   
   # end EMDSYN
end
function emdbasesel(S::BitArray, degs::Vector{Int}, rdim::Int, nout::Int, simple::Bool)
    #   emdbasesel(S, degs, rdim, nout, simple) -> (seli, selord)
    #
    #   Select admissible basis vectors for solving the exact model detection problem (EMDP) 
    #   using the binary structure matrix `S` corresponding to `nvec` basis vectors. 
    #   `S` is a `nvec x lfreq x N` binary 3-dimensional array or a `nvec x mf` 
    #   binary 2-dimensional array (reshaped internally as an `nvec x mf x 1` 
    #   binary 3-dimensional array).
    #   `seli` contains `nout`-touples of indices of basis vectors whose linear combination 
    #   is admissible, i.e. , `S[seli[i],k,j]` has all elements nonzero for all `k = 1:lfreq` 
    #   and `j = 1:N`. This ensures that the EMDP is solvable by using 
    #   model detection filters with `rdim <= nout` outputs. 
    #   If the associated `nvec` degrees contained in `degs` are provided, then
    #   `selord[i]` is the corresponding tentatively achievable least filter order.
    #   If `simple = true`, a simple basis is assumed, in which case, `degs[i]` is 
    #   also the order of the `i`-th basis vector. If `simple = false`, a minimum 
    #   rational basis is assumed. `selord` is empty if `degs` is empty. 
 
    #   Method: The selection approach is used in conjunction with the synthesis 
    #   Procedure EMD described in [1]. 
 
    #   References:
    #   [1] Varga A.
    #   Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017.
 
    ndim = ndims(S)
    if ndim == 3    
       nvec, lfreq, N = size(S)  # numbers of vectors, numbers of frequencies, numbers of models
    else
       nvec, N = size(S)   # numbers of vectors, number of models  
       lfreq = 1
       S = reshape(S,nvec,1,N)
    end
    nd = length(degs)
    nodegs = (nd == 0)
 
    nodegs || length(degs) == nvec || error("the dimension of degs must be equal to the number of rows of S")
    
    (rdim >=1 && rdim <= nvec) || error("ndim must have a positive value not exceeding $nvec")
    (nout >=rdim && nout <= nvec) || error("nout must have a value at least $rdim and at most $nvec")
    
    nvec == 1 && (return [[1]], nodegs ? Int[] : degs )
    
    # find rdim combinations of nout vectors which solve the EFDP 
    seli = collect(combinations(Vector(1:nvec),nout))
    ni = length(seli)
    selord = nodegs ? Int[] : fill(-1,ni) 
    nqmax = sum(degs)
    ii = trues(ni)
    for i = 1:ni
        indv = seli[i];
        # check admissibility
        #stest = all(view(S,indv,:,:))
        stest = true;
        for j = 1:N
            #any(view(S,indv,:,j)) || (stest = false; break)
            for k = 1:lfreq
               stest = stest & any(view(S,indv,k,j))
            end  
        end
        if stest
           if !nodegs
              # estimate orders 
              if simple || rdim == nout
                 # degree = the sums of degrees of selected vectors
                 selord[i] = sum(degs[indv])
              else
                 # degree = rdim times the maximum degree of selected vectors
                 selord[i] = min(nqmax,rdim*maximum(degs[indv]))
              end
           end
        else
           ii[i] = false
        end
   end
   
    seli = seli[ii]
    if !nodegs 
       selord = selord[ii];
       # sort row combinations to ensure increasing tentative orders  
       ii = sortperm(selord)
       seli = seli[ii]
       selord = selord[ii];
    end
    return seli, selord      
    # end EMDBASESEL
 end
 