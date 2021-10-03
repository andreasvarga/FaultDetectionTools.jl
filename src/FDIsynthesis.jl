"""
    efdsyn(sysf::FDIModel; rdim, nullspace = true, simple = false, minimal = true, 
                           sdeg, smarg, poles, HDesign, FDtol, FDGainTol, FDfreq, 
                           tcond, offset, atol, atol1, atol2, atol3, rtol, fast = true) 
                           -> (Q::FDFilter, R::FDFilterIF, info)

Solve the exact fault detection problem (EFDP) for a given synthesis model
`sysf` with additive faults. The computed stable and proper filter objects `Q` and `R` contain the 
fault detection filter, representing the solution of the EFDP, and its internal form, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.degs`, `info.S`, and `info.HDesign`, 
contains additional synthesis related information (see below). 

The continuous- or discete-time system `sysf.sys` is in a standard
or descriptor state-space form `sysf.sys = (A-λE,B,C,D)`, which corresponds to the input-output form  

       y = Gu(λ)*u + Gd(λ)*d + Gf(λ)*f + Gw(λ)*w + Ga(λ)*aux,

with the Laplace- or Z-transformed plant outputs `y`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Gu(λ)`, `Gd(λ)`, `Gf(λ)`, `Gw(λ)`, and `Ga(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysf.controls`, `sysf.disturbances`, `sysf.faults`, `sysf.noise` and `sysf.aux`, respectively.

The fault detection filter object `Q`, contains in `Q.sys` the resulting filter 
in a standard state-space form, which generates the residual signal `r`. 
The corresponding input-output (implementation) form is

            r = Qy(λ)*y + Qu(λ)*u               

where `Qy(λ)` and `Qu(λ)` are the transfer function matrices from the output and control inputs to the residual. 
The indices of output and control inputs are contained in the integer vectors 
`Q.outputs` and `Q.controls`, respectively.

The fault detection filter in internal form object `R`, contains `R.sys`, the resulting 
internal form of the filter 
in a standard state-space form, which generates the residual signal `r`, and corresponds to the 
input-output form

       r = Ru(λ)*u + Rd(λ)*d + Rf(λ)*f + Rw(λ)*w + Ra(λ)*aux ,

where 

       | Ru(λ) Rd(λ) Rf(λ) Rw(λ) Ra(λ) | = |Qy(λ) Qu(λ)|*| Gu(λ) Gd(λ) Gf(λ) Gw(λ) Ga(λ) |. 
                                                         |  I     0     0     0     0    |

The solution of the EFDP ensures that `Ru(λ) = 0`, `Rd(λ) = 0`, and `Rf(λ)` has all its columns nonzero. 
The indices of the inputs `u`, `d`, `f`, `w` and `aux` of the resulting filter `R.sys` are 
contained in the integer vectors `R.controls` (void), `R.disturbances` (void), 
`R.faults`, `R.noise` and `R.aux`, respectively.

The resulting filters `Q.sys` and `R.sys` have observable state-space realizations
`(AQ,BQ,CQ,DQ)` and `(AQ,BR,CQ,DR)`, respectively, and thus share the observable pairs `(AQ,CQ)`. 

Various user options can be specified via keyword arguments as follows:

If `minimal = true` (default), a least order filter synthesis is performed, while 
with `minimal = false` a full order synthesis is performed.  

If `HDesign = H`, a full row rank design matrix `H` is used to build `rdim = q` 
linear combinations of the left nullspace basis vectors (default: `HDesign = missing`)

`rdim = q` specifies the desired number `q` of residual outputs for `Q` and `R`. 
The default value of `q` is chosen as follows: if `HDesign = missing`, then  
`q = 1`, if `minimal = true`, or `q` is the number of the nullspace basis 
vectors used for the initial synthesis, if `minimal = false`; 
if `HDesign = H` specifies a full row rank design matrix `H`, 
then `q` is the row dimension of `H`. 

`FDFreq = freq` specifies a vector of real frequency values or a scalar real frequency value
for strong detectability checks (default: `FDFreq = missing`).

If `nullspace = true` (default), a minimal proper nullspace basis is used for the synthesis of the
fault detection filter. If `nullspace = false`, a full-order observer based nullspace basis is used. 
This option can be only used for a proper system without disturbance inputs. 

If `simple = true`, a simple proper nullspace basis is emplyed for synthesis. 
The orders of the basis vectors are provided in `info.deg`. 
If `simple = false` (default), then no simple basis is computed. 

To assess the stability of poles, a boundary offset  `β` can be specified via the keyword argument `offset = β`. 
Accordingly, for the stability of a continuous-time system all real parts of poles must be at most `-β`, 
while for the stability of a discete-time system all moduli of poles must be at most `1-β`. 
The default value used for `β` is `sqrt(ϵ)`, where `ϵ` is the working machine precision. 

The stability domain `Cs` of poles is defined by 
the keyword argument for the stability margin `smarg = α`, as follows: 
for a continuous-time system, `Cs` is the set of complex numbers 
with real parts at most `α`, 
while for a discete-time system, `Cs` is the set of complex numbers with 
moduli at most `α < 1` (i.e., the interior of a disc of radius `α` centered in the origin). 
If `smarg` is missing, then the employed default values are `α = -β` 
for a continuous-time system and `α = 1-β` for a discete-time system, 
where `β` is the boundary offset specified by the keyword argument `offset = β`. 

The poles of the filters `Q` and `R` can be specified using the keyword arguments `poles` and `sdeg`,
where `poles = v` specifies a complex vector `v` containing a complex conjugate set  
of desired poles within the stability domain `Cs` (default: `poles = missing`)
and `sdeg = γ` is the prescribed stability degree for the poles (default: `γ = -0.05` for the real parts of poles for a continuous-time system and
`γ = 0.95` for the magnitudes of poles for a discete-time system). 

The maximum alowed condition number `tcmax` of the employed non-orthogonal transformations 
can be specified as `tcond = tcmax` (Default: `tcmax = 1.e4`).

`FDtol = tol1` specifies the threshold `tol1` for fault detectability checks
   (default: `tol1 = 0.0001`).

`FDGainTol = tol2` specifies the threshold `tol2` for strong fault detectability checks
   (default: `tol2 = 0.01`). 

The rank determinations in the performed reductions
are based on rank revealing QR-decompositions with column pivoting 
if `fast = true` or the more reliable SVD-decompositions if `fast = false`.

The keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, 
the absolute tolerance for the nonzero elements of `A`, `B`, `C`, `D`,  
the absolute tolerance for the nonzero elements of `E`,  
and the relative tolerance for the nonzero elements of `A`, `B`, `C`, `D` and `E`.  
The default relative tolerance is `n*ϵ`, where `ϵ` is the working machine epsilon 
and `n` is the order of the system `sysf.sys`. 
The keyword argument `atol3` is an absolute tolerance for observability tests
(default: internally determined value). 
The keyword argument `atol` can be used 
to simultaneously set `atol1 = atol`, `atol2 = atol` and `atol3 = atol`. 

The resulting named tuple `info` contains `(tcond, degs, S, HDesign) `, where:

`info.tcond` is the maximum of the condition numbers of the employed 
   non-orthogonal transformation matrices; a warning is issued if `info.tcond >= cmax`;

`info.degs` is an integer vector containing the increasingly ordered degrees of a left minimal   
polynomial nullspace basis of `G(λ) := [ Gu(λ) Gd(λ); I 0]` (also the left Kronecker indices of `G(λ)`), if the 
state-space realization of `[Gu(λ) Gd(λ)]` is minimal;

`info.S` is the binary structure matrix corresponding to the computed left nullspace basis;

`info.HDesign` is the design matrix `H` employed for the synthesis of 
   the fault detection filter.
   
_Method:_ The Procedure EFD from [1] is implemented to solve 
the exact fault detection problem. For more details on 
the least order synthesis of fault detection filters see [2].

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017; sec. 5.2.

[2] A. Varga, On computing least order fault detectors using rational nullspace bases. 
IFAC SAFEPROCESS'03 Symposium, Washington DC, USA, 2003.
"""
function efdsyn(sysf::FDIModel{T}; rdim::Union{Int,Missing} = missing, poles::Union{AbstractVector,Missing} = missing, 
                      sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing, 
                      nullspace::Bool = true, minimal::Bool = true, simple::Bool = false,
                      FDtol::Real = 0.0001, FDGainTol::Real = 0.01, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                      tcond::Real = 1.e4, HDesign::Union{AbstractMatrix,Missing} = missing,
                      offset::Real = sqrt(eps(float(real(T)))), atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                      rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T)))))*iszero(max(atol1,atol2)), 
                      fast::Bool = true) where T
   Ts = sysf.sys.Ts                  
   disc = (Ts != 0);  # system type (continuous- or discete-time)
   
   # decode options
      
   strongFD = !ismissing(FDfreq)
   strongFD && !isa(FDfreq,Vector) && (FDfreq = [FDfreq]) 
   strongFD && (lfreq = length(FDfreq))

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
  
   # imposed design option to form linear combinations of basis vectors
   emptyHD = ismissing(HDesign)
   if !emptyHD
      !ismissing(rdim) && size(HDesign,1) != rdim && error("row dimension of HDesign must be equal to rdim")
      size(HDesign,1) == rank(HDesign) || error("HDesign must have full row rank")
   end
  
   # decode input information
   inpu = sysf.controls; mu = length(inpu)  
   inpd = sysf.disturbances; md = length(inpd) 
   inpf = sysf.faults; mf = length(inpf)  
   inpw = sysf.noise;  mw = length(inpw) 
   inpaux = sysf.aux;  maux = length(inpaux)  
   
   m = mu+md+mf+mw+maux;       # total number of inputs
   p = size(sysf.sys,1);       # number of measurable outputs
    
   if mf == 0 && minimal
      @warn "Minimal synthesis option not feasible in the case of no faults"
      minimal = false
   end
         
   # Step 1): nullspace based reduction
   #
   desc = (sysf.sys.E != I)
   m2 = mf+mw+maux
   sdegNS = strongFD ? sdegdefault : missing     
   if nullspace || simple || md > 0 || (desc && rcond(sysf.sys.e) < 1.e-7 )
      # form [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
      syse = [sysf.sys; eye(mu,m)];
      #
      # compute a left nullspace basis Q = Q1 of G1 = [Gu Gd; I 0] = 0 and
      # obtain QR = [ Q R ], where R = [ Rf Rw Raux] = Q*[Gf Gw Ga;0 0 0]
      QR, info1 = glnull(syse, m2; simple, atol1, atol2, rtol, fast, sdeg = sdegNS, offset) 
      degs = info1.degs
      tcond1 = info1.tcond
   elseif mu == 0 && md == 0
         QR = [eye(p) sysf.sys]
         degs = Int[]; tcond1 = 1.
   else
      # compute minimal basis as Q = Q1 = [ I -Gu] and set
      # QR = [ Q R ], where R = [ Gf Gw Ga ]
      QR = [ eye(p) dss(sysf.sys.A, sysf.sys.E, [-sysf.sys.B[:,inpu] sysf.sys.B[:,[inpf; inpw; inpaux]]],
             sysf.sys.C, [-sysf.sys.D[:,inpu] sysf.sys.D[:,[inpf; inpw; inpaux]]]; Ts)] 
      # perform stabilization if strong detectability has to be enforced
      strongFD  && (QR = glcf(QR; atol1, atol2, atol3, rtol, fast)[1]) 
      degs = Int[]; tcond1 = 1.
   end
   infodegs = degs
   
   nvec = size(QR,1);   # number of basis vectors
   # check solvability conditions
   nvec == 0 && error("empty nullspace basis: the EFDP is not solvable")

   nq = order(QR)             # order of the minimal basis
   
   # set H for checking the solvability condition
   if emptyHD
      Htemp = eye(nvec)
   else
      degs = Int[];
      rdim = size(HDesign,1); nh = size(HDesign,2)
      if nh < nvec
         # pad with zeros: row rank is preserved
         Htemp = [ HDesign zeros(T, rdim, nvec-nh) ]
      else
         # remove trailing columns if necessary: rank may drop
         Htemp = HDesign[:,1:nvec];
         if nh > nvec && rdim != rank(Htemp)
            error("The leading part of HDesign must have full row rank")
         end
      end
   end
   
   indf = (p+mu) .+ Vector(1:mf)   # input indices of Rf in QR
   if mf == 0
      # handle the case of no faults as a normal case
      S = falses(nvec,0);
   else
      if strongFD 
         S = fdisspec_(Htemp*QR[:,indf], FDfreq; FDGainTol, atol1, atol2, rtol = 0, fast)[1]
         # check strong detectability conditions 
         for ii = 1:lfreq
            all(maximum(S[:,:,ii],dims=1)) || error("strong detection of all faults not feasible")
         end
      else
         # check weak detectability conditions 
         S = fditspec_(Htemp*QR[:,indf]; atol1, atol2, rtol, FDtol)
         all(maximum(S,dims=1)) || error("detection of all faults not feasible")
      end
   end
   
   # setup the number of filter outputs
   if minimal
      #  least order design
      if ismissing(rdim)
         if emptyHD 
            rdim = 1   
         else
            rdim = size(HDesign,1)
         end
      else
         rdim = min(rdim,nvec)
      end
   else
      #  full order design
      if ismissing(rdim) 
         if emptyHD 
            rdim = nvec
         else
            rdim = min(size(HDesign,1),nvec);
         end
      else
         if mf == 0
            if rdim < nvec && emptyHD
               @warn "rdim reset to $nvec"
               rdim = nvec
            end
         else
            rdim = min(rdim,nvec); 
         end
      end
   end
          
   # Step 2): compute admissible Q2 to reduce the order of Q2*Q;  
   # update Q <- Q2*Q, R <- Q2*R 
   
   # reorder degs to correspond to the expected orders of basis vectors 
   # corresponding to the actual order of outputs of QR 
   reverse!(degs) 
   if rdim < nvec && mf > 0
      # determine possible low order syntheses using i >= rmin basis vectors
      # and the corresponding expected orders    
      
      finish = false    # set termination flag
      nout = rdim       # initialize number of selected basis vectors
      # if ~simple && minimal
      #    QR = xperm(QR,nq:-1:1);  # permute states to speedup glmcover1
      # end
      itry = 1; 
      while !finish     
          # choose nout basis vectors, which potentially lead to a least order
          # filter with rdim outputs:
          # basesel(i,:) contains the indices of candidate basis vectors;
          # ordsel(i)    contains the presumably achievable least orders
          basesel, ordsel = efdbasesel(S, degs, rdim, nout, simple) 
          #
          # update the synthesis using the selections of candidate vector(s),
          # starting with the least (potentially) achievable order
          for i = 1:size(basesel,1)
              baseind = basesel[i] # indices of current basis selection
              if rdim == nout
                  hbase = eye(rdim)
              else
                  hbase = rand(rdim,nout) 
              end
              ip = [baseind; setdiff(1:nvec,baseind)][:]
              if simple
                 if minimal
                    if emptyHD 
                       # select vectors and elliminate unobservable dynamics  
                       noelim = falses(nq) 
                       #ell = sum(degs(1:basesel(i,1)-1)); 
                       ell = sum(degs[1:basesel[i][1]-1]); 
                       for jj = 1:nout 
                           ellnext = sum(degs[1:baseind[jj]]);
                           noelim[ell+1:ellnext] .= true;
                           ell = ellnext
                       end
                    end
                    if rdim == nout
                       if emptyHD 
                          #QRfwtest = modred(QR[baseind,:],~noelim,'truncate');
                          ir = noelim
                          Ar, Er, Br, Cr, Dr = dssdata(QR[baseind,:])
                          QRfwtest = dss(view(Ar,ir,ir), Er == I ? I : view(Er,ir,ir),view(Br,ir,:),view(Cr,:,ir), Dr; Ts)
                          h = Htemp[ip[1:rdim],:];
                       else
                          QRfwtest = gir(Htemp*QR; atol1, atol2, rtol);
                       end
                    else
                       # this case is possible only if HDesign is empty
                       # build rdim linear combinations of the first nout vectors 
                       # QRfwtest = hbase*modred(QR(baseind,:),~noelim,'truncate');
                       ir = noelim
                       Ar, Er, Br, Cr, Dr = dssdata(QR[baseind,:])
                       QRfwtest = hbase*dss(view(Ar,ir,ir),Er == I ? I : view(Er,ir,ir),view(Br,ir,:),view(Cr,:,ir), Dr; Ts)
                       h = [ hbase zeros(T, rdim, nvec-nout) ]
                       h = h[:,ip]  # permute columns to match unpermuted QR 
                    end
                 else
                    if rdim == nout
                       if emptyHD 
                          h = Htemp[ip[1:rdim],:]
                          QRfwtest = gir(QR[baseind,:]; atol1, atol2, rtol, infinite = false)
                       else
                          QRfwtest = gir(Htemp*QR; atol1, atol2, rtol, infinite = false)
                       end
                    else
                       # this case is possible only if HDesign is empty
                       # build rdim linear combinations of the first nout vectors 
                       QRfwtest = gir(hbase*QR[baseind,:]; atol1, atol2, rtol, infinite = false) 
                       h = [ hbase zeros(T, rdim, nvec-nout) ] 
                       h = h[:,ip]  # permute columns to match unpermuted QR 
                    end
                 end
              else
                 if minimal
                    if rdim == nout
                       if emptyHD 
                          QRfwtest, _, info2 = glmcover1(QR[ip,:], rdim; atol1, atol2, rtol)
                          if !isempty(ordsel) && (order(QRfwtest) != ordsel[i])
                             @warn "efdsyyn: expected reduced order not achieved"
                          end
                          h = Htemp[ip[1:rdim],:]
                       else
                          QRfwtest, _, info2 = glmcover1([Htemp; eye(nvec)]*QR[ip,:], rdim; atol1, atol2, rtol)
                       end
                    else  
                       # this case is possible only if HDesign is empty
                       # build rdim linear combinations of the first nout vectors 
                       h = [ hbase zeros(T, rdim, nvec-nout) ]; 
                       QRfwtest, _, info2 = glmcover1([h; eye(nvec)]*QR[ip,:], rdim; atol1, atol2, rtol)
                       h = h[:,ip]  # permute columns to match unpermuted QR 
                    end
                 else
                    if rdim == nout
                       if emptyHD
                          h = Htemp[ip[1:rdim],:]
                          QRfwtest = gir(QR[baseind,:]; atol1, atol2, rtol, infinite = false) 
                       else
                          QRfwtest = gir(Htemp*QR; atol1, atol2, rtol, infinite = false) 
                       end
                    else
                       QRfwtest = gir(hbase*QR[baseind,:]; atol1, atol2, rtol, infinite = false) 
                       h = [ hbase zeros(T, rdim, nvec-nout) ]; 
                       h = h[:,ip]  # permute columns to match unpermuted QR 
                    end
                 end
              end
              # check complete fault detectability of the current design; 
              if (rdim == nout && minimal) || rdim < nout
                 # dismiss design if check fails
                 if strongFD 
                    Stest = fdisspec_(QRfwtest[:,indf], FDfreq; stabilize = true, block = true, 
                                              FDGainTol, atol1, atol2, atol3, rtol = 0, fast)[1]
                 else
                    Stest = fditspec_(QRfwtest[:,indf]; block = true, atol1, atol2, rtol, FDtol)
                 end
                 if all(Stest) 
                    if !simple && minimal
                       # adjust condition number of employed transformations
                       tcond1 = max(tcond1, info2.fnorm, info2.tcond)
                       tcond1 > tcond && 
                          @warn "efdsyn: possible loss of numerical stability due to ill-conditioned transformations"
   #                     if ~emptyHD
   #                        info.HDesign = Htemp;
   #                     end
                    end
                    QR = QRfwtest
                    finish = true
                    break
                 end
              else
                 QR = QRfwtest
                 finish = true
                 break
              end
          end
          nout += 1
          if nout > nvec
             if itry > 5
                finish = true
                @warn "fault detectability not achieved with the chosen number of residuals"
             else
                itry += 1
                nout -= 1
             end
          end
      end
      if emptyHD
         Htemp = h
      end
   else
      hbase = eye(rdim)
      if simple
         baseind = 1:rdim 
      else
         baseind = 1;
      end
      h = eye(rdim)
      if !emptyHD
         QR = Htemp*QR
      else
         # use full minimum basis 
         Htemp = h
      end
   end
   
   # Step 3): compute Q3 such that Q3*Q has a desired stability degree;  
   # update Q <- Q3*Q, R <- Q3*R 
   k = 1;
   if simple && isequal(hbase,I) && emptyHD
       # exploit the block diagonal structure of basis matrices al and cl
       # to compute block-diagonal Q3
       al, el, bl, cl, dl, = dssdata(QR)
       for i = 1:length(baseind)
           blkord = degs[baseind[i]]
           if blkord
              i1 = k:k+blkord-1; 
              QRfwi = glcf(dss(al[i1,i1],el[i1,i1],bl[i1,:],cl[i:i,i1],dl[i:i,:];Ts); 
                           atol1, atol2, atol3, sdeg, smarg, evals = poles)[1]
              al[i1,i1] = QRfwi.A; bl[i1,:] = QRfwi.B;  cl[i,i1] = QRfwi.C; dl[i,:] = QRfwi.D 
              el[i1,i1] = (QRfwi.e == I) ? eye(blkord) : QRfwi.E  
              k += blkord
           end
       end
       QR = dss(al, el, bl, cl, dl; Ts)
   else
       QR = glcf(QR; atol1, atol2, atol3, sdeg, smarg, evals = poles)[1]
   end
       
   # scale Rf to ensure unit minimum column gains
   if mf > 0
      if strongFD && minimum(FDfreq) == 0
         # compute minimum DC gains  
         dcg = dcgain(QR[:,indf]; atol1, atol2, rtol, fast) 
         y = maximum(abs.(dcg),dims=1)
         indj = sortperm(y[:])[1] 
         scale = y[indj]
         indi = sortperm(abs.(dcg[:,indj]))[end]
         sc = sign(dcg[indi,indj])/scale
      else
         # compute the minimum of H-inf norms of columns
         sc = 1/fdhinfminus(QR[:,indf])[1]
      end
      QR = sc*QR
   end
   
   # transform to standard state-space
   QR = gss2ss(QR; atol1, atol2, rtol)[1]
   
   Q = FDFilter(QR, p, mu)
   #R = fdRset(QR[:,p+mu+1:end],faults = Vector(1:mf), noise = mf .+ Vector(1:mw), aux = (mf+mw) .+ Vector(1:maux))
   R = FDFilterIF(QR,0,0,mf,mw,maux; moff = p+mu)
   info = (tcond = tcond1, degs = infodegs, S = S, HDesign = convert(Matrix{Float64},Htemp))

   return Q, R, info
   
   # end EFDSYN
end
"""
    efdsyn(sysf::FDIModel, S; rdim, nullspace = true, simple = false, minimal = true, 
                           sdeg, smarg, poles, HDesign, FDtol, FDGainTol, FDfreq, 
                           tcond, offset, atol, atol1, atol2, atol3, rtol, fast = true) 
                           -> (Q::FDFilter, R::FDFilterIF, info)

Solve the exact fault detection isolation problem (EFDIP) for a given synthesis model
`sysf` with additive faults and a given binary structure vector `S`. 
The computed stable and proper filter objects `Q` and `R` contain the 
fault detection filter, representing the solution of the EFDIP, and its internal form, respectively, and are determined such that
`R.sys[:,faults]` has its `j`-th column nonzero if `S[j] = 1` and the `j`-th column is zero if `S[j] = 0`. 
For the description of the keyword parameters see the function [`efdsyn`](@ref). 
"""
function efdsyn(sysf::FDIModel{T}, SFDI::Union{BitVector,Vector{Bool}}; kwargs...) where T
   mu = length(sysf.controls)  
   md = length(sysf.disturbances) 
   mf = length(sysf.faults)
   mw = length(sysf.noise) 
   maux = length(sysf.noise)  
   mf == length(SFDI) || error("number of faults must be equal to the length dimension of SFDI")
   indd = Vector(1:mf)[SFDI .== false] 
   indf = Vector(1:mf)[SFDI .== true] 
   sysc = fdimodset(sysf.sys, c = 1:mu, d = [Vector(mu .+ (1:md)); (mu+md) .+ indd], 
                    f = (mu+md) .+ indf, aux = (mu+md) .+ (1:mf+mw+maux));
   #
   Q, Rft, info = efdsyn(sysc; kwargs...)
   return Q, FDFilterIF(Rft.sys[:,Rft.aux],0,0,mf,mw,maux), info
end
function efdbasesel(S::BitArray, degs::Vector{Int}, rdim::Int, nout::Int, simple::Bool)
   #   efdbasesel(S, degs, rdim, nout, simple) -> (seli, selord)
   #
   #   Select admissible basis vectors for solving the exact fault detection problem (EFDP) 
   #   using the binary structure matrix `S` corresponding to `nvec` basis vectors. 
   #   `S` is a `nvec x mf x n` binary 3-dimensional array or a `nvec x mf` 
   #   binary 2-dimensional array (reshaped internally as an `nvec x mf x 1` 
   #   binary 3-dimensional array).
   #   `seli` contains `nout`-touples of indices of basis vectors whose linear combination 
   #   is admissible, i.e. , `S[seli[i],:,k]` has all columns nonzero for all `k = 1:n`. 
   #   This ensures that the EFDP is solvable by using 
   #   fault detection filters with `rdim <= nout` outputs. 
   #   If the associated `nvec` degrees contained in `degs` are provided, then
   #   `selord[i]` is the corresponding tentatively achievable least filter order.
   #   If `simple = true`, a simple basis is assumed, in which case, `degs[i]` is 
   #   also the order of the `i`-th basis vector. If `simple = false`, a minimum 
   #   rational basis is assumed. `selord` is empty if `degs` is empty. 

   #   Method: The selection approach is used in conjunction with the synthesis 
   #   Procedure EFD described in [1]. 

   #   References:
   #   [1] Varga A.
   #   Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017.


   ndim = ndims(S)
   if ndim == 3 
      nvec = size(S,1)  # numbers of vectors
      n = size(S,3)     # number of frequencies
   else
      nvec, mf = size(S)   # numbers of vectors, number of faults  
      n = 1
      S = reshape(S,nvec,mf,1)
   end
   nd = length(degs)
   nodegs = (nd == 0)

   nodegs || length(degs) == nvec || error("the dimension of degs must be equal to the number of rows of S")
   
   (rdim >=1 && rdim <= nvec) || error("ndim must have a positive value not exceeding $nvec")
   (nout >=rdim && nout <= nvec) || error("nout must have a value at least $rdim and at most $nvec")
   
   nvec == 1 && (return [1], nodegs ? Int[] : degs )
   
   # find rdim combinations of nout vectors which solve the EFDP 
   seli = collect(combinations(Vector(1:nvec),nout))
   ni = length(seli)
   selord = nodegs ? Int[] : fill(-1,ni) 
   nqmax = sum(degs)
   ii = trues(ni)
   for i = 1:ni
       indv = seli[i];
       # check admissibility
       stest = true;
       for k = 1:n
         all(maximum(view(S,indv,:,k),dims=1)) || (stest = false; break)
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
   # end EFDBASESEL
end
"""
    efdisyn(sysf::FDIModel, SFDI; rdim, nullspace = true, simple = false, minimal = true, separate = false,
                           sdeg, smarg, poles, HDesign, FDtol, FDGainTol, FDfreq, 
                           tcond, offset, atol, atol1, atol2, atol3, rtol, fast = true) 
                           -> (Q::FDFilter, R::FDFilterIF, info)

Solve the exact fault detection and isolation problem (EFDIP) for a given synthesis model
`sysf` with additive faults and a given binary structure matrix `SFDI` with `nb` rows (specifications). 
The computed stable and proper filter objects `Q` and `R` contain the 
fault detection and isolation filter, representing the solution of the EFDIP, and its internal form, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.degs` and `info.HDesign`, 
contains additional synthesis related information (see below). 

The continuous- or discete-time system `sysf.sys` is in a standard
or descriptor state-space form `sysf.sys = (A-λE,B,C,D)`, which corresponds to the input-output form  

       y = Gu(λ)*u + Gd(λ)*d + Gf(λ)*f + Gw(λ)*w + Ga(λ)*aux,

with the Laplace- or Z-transformed plant outputs `y`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Gu(λ)`, `Gd(λ)`, `Gf(λ)`, `Gw(λ)`, and `Ga(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysf.controls`, `sysf.disturbances`, `sysf.faults`, `sysf.noise` and `sysf.aux`, respectively.

The fault detection and isolation filter object `Q`, contains in `Q.sys` the resulting bank of `nb` filters. 
The i-th filter `Q.sys[i]` is in a standard state-space form and generates `r_i`, the `i`-th component (scalar or vector) 
of the overall residual vector `r := [r_1; r_2; ...; r_nb]`. 
The corresponding input-output (implementation) form of the `i`-th filter is

            r_i = Qyi(λ)*y + Qui(λ)*u   ,            

where `Qyi(λ)` and `Qui(λ)` are the transfer function matrices from the output and control inputs to the `i`-th residual component. 
The indices of output and control inputs are contained in the integer vectors 
`Q.outputs` and `Q.controls`, respectively.

The fault detection and isolation filter in internal form object `R`, contains `R.sys`, the resulting bank of `nb` 
internal form of the filters.  
The i-th filter `R.sys[i]` is in a standard state-space form, which generates the residual signal `r_i`, and corresponds to the 
input-output form

       r_i = Rui(λ)*u + Rdi(λ)*d + Rfi(λ)*f + Rwi(λ)*w + Rai(λ)*aux ,

where 

       | Rui(λ) Rdi(λ) Rfi(λ) Rwi(λ) Rai(λ) | := |Qyi(λ) Qui(λ)]*| Gu(λ) Gd(λ) Gf(λ) Gw(λ) Ga(λ) |. 
                                                                 |   I     0     0     0     0   |

The solution of the EFDIP ensures that for the `i`-th filter, `Rui(λ) = 0`, `Rdi(λ) = 0`, and 
`Rfi(λ)` has its `j`-th column nonzero if the `(i,j)`-th element of `SFDI` is nonzero. 
The indices of the inputs `u`, `d`, `f`, `w` and `aux` of the resulting filter `R.sys` are 
contained in the integer vectors `R.controls` (void), `R.disturbances` (void), 
`R.faults`, `R.noise` and `R.aux`, respectively.

The resulting component filters `Q.sys[i]` and `R.sys[i]` have observable state-space realizations
`(AQi,BQi,CQi,DQi)` and `(AQi,BRi,CQi,DRi)`, respectively, and thus share the observable pairs `(AQi,CQi)`. 

Various user options can be specified via keyword arguments as follows:

If `minimal = true` (default), least order filter synthesis is performed to determine each of the component filters
`Q.sys[i]` and `R.sys[i]` for `i = 1, ...,nb`, while 
with `minimal = false` full order synthesis is performed.  

If `HDesign = H`, then `H` is an `nb`-dimensional array of full row rank or empty design matrices `H = [H_1, ..., H_nb]`,
where `H_i` is the design matrix employed for the synthesis of the `i`-th component filter (default: `HDesign = missing`)

`rdim = q` specifies the vector `q`, whose `i`-th component `q[i]` specifies 
the number of residual outputs for the `i`-th component filters `Q.sys[i]` and `R.sys[i]`. 
If `q` is a scala, then a vector `rdim` with all components equal to `q` is assumed.
The default value of `q[i]` is chosen as follows: if `HDesign = missing` or `H_i` is empty then  
`q[i] = 1`, if `minimal = true`, or `q[i]` is the number of the nullspace basis 
vectors used for the synthesis of `Q.sys[i]` and `R.sys[i]`, if `minimal = false`; 
if  `H_i` specifies a full row rank design matrix, then `q[i]` is the row dimension of `H_i`. 

`FDFreq = freq` specifies a vector of real frequency values or a scalar real frequency value
for strong detectability checks (default: `FDFreq = missing`).

If `nullspace = true` (default), a minimal proper nullspace basis is used at the 
initial reduction step, if `separate = false`, 
or at all synthesis steps, if `separate = true`.
If `nullspace = false`, a full-order observer based nullspace basis is used at the 
initial reduction step, if `separate = false`, or at all synthesis steps, if `separate = true`.
This option can  only be used for a proper system without disturbance inputs. 

If `simple = true`, simple proper nullspace bases are emplyed for synthesis. 
The orders of the basis vectors employed for the synthesis of `i`-th filter
are provided in `info.deg[i]`. 
If `simple = false` (default), then no simple bases are computed. 

If `separate = false` (default), a two-step synthesis procedure is employed, 
where a minimal proper nullspace basis is used at the 
initial reduction step. 
If `separate = true`, the filter components are separately determined by solving
appropriately formulated fault detection problems. 

To assess the stability of poles, a boundary offset  `β` can be specified via the keyword argument `offset = β`. 
Accordingly, for the stability of a continuous-time system all real parts of poles must be at most `-β`, 
while for the stability of a discete-time system all moduli of poles must be at most `1-β`. 
The default value used for `β` is `sqrt(ϵ)`, where `ϵ` is the working machine precision. 

The stability domain `Cs` of poles is defined by 
the keyword argument for the stability margin `smarg = α`, as follows: 
for a continuous-time system, `Cs` is the set of complex numbers 
with real parts at most `α`, 
while for a discete-time system, `Cs` is the set of complex numbers with 
moduli at most `α < 1` (i.e., the interior of a disc of radius `α` centered in the origin). 
If `smarg` is missing, then the employed default values are `α = -β` 
for a continuous-time system and `α = 1-β` for a discete-time system, 
where `β` is the boundary offset specified by the keyword argument `offset = β`. 

The poles of the filters `Q` and `R` can be specified using the keyword arguments `poles` and `sdeg`,
where `poles = v` specifies a complex vector `v` containing a complex conjugate set  
of desired poles within the stability domain `Cs` (default: `poles = missing`)
and `sdeg = γ` is the prescribed stability degree for the poles (default: `γ = -0.05` for the real parts of poles for a continuous-time system and
`γ = 0.95` for the magnitudes of poles for a discete-time system). 

The maximum alowed condition number `tcmax` of the employed non-orthogonal transformations 
can be specified as `tcond = tcmax` (Default: `tcmax = 1.e4`).

`FDtol = tol1` specifies the threshold `tol1` for fault detectability checks
   (default: `tol1 = 0.0001`).

`FDGainTol = tol2` specifies the threshold `tol2` for strong fault detectability checks
   (default: `tol2 = 0.01`). 

The rank determinations in the performed reductions
are based on rank revealing QR-decompositions with column pivoting 
if `fast = true` or the more reliable SVD-decompositions if `fast = false`.

The keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, 
the absolute tolerance for the nonzero elements of `A`, `B`, `C`, `D`,  
the absolute tolerance for the nonzero elements of `E`,  
and the relative tolerance for the nonzero elements of `A`, `B`, `C`, `D` and `E`.  
The default relative tolerance is `n*ϵ`, where `ϵ` is the working machine epsilon 
and `n` is the order of the system `sysf.sys`. 
The keyword argument `atol3` is an absolute tolerance for observability tests
(default: internally determined value). 
The keyword argument `atol` can be used 
to simultaneously set `atol1 = atol`, `atol2 = atol` and `atol3 = atol`. 


The resulting named tuple `info` contains `(tcond, degs, HDesign) `, where:

`info.tcond` is an `nb`-dimensional vector, whose `i`-th component is the maximum of the condition numbers of the employed 
non-orthogonal transformation matrices employed for the synthesis of the `i`-th filter component; 
a warning is issued if any `info.tcond[i] >= cmax`;

`info.degs` is an `nb`-dimensional vector, whose `i`-th component is an integer vector 
containing the degrees of the basis vectors of the employed simple
nullspace basis for the synthesis of the i-th filter component, if `simple = true, 
and the degrees of the
basis vectors of an equivalent polynomial nullspace basis, if `simple = false`;

`info.HDesign` is an `nb`-dimensional vector, whose `i`-th component is 
is the design matrix `H_i` employed for the synthesis of 
the `i`-th fault detection filter.
   
_Method:_ The Procedure EFDI from [1] is implemented to solve 
the exact fault detection and isolation problem. 
This procedure relies on the nullspace-based synthesis method proposed in [2]. For more 
details on the least order synthesis of fault detection filters see [3].

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017; sec. 5.4.

[2] A. Varga,  On designing least order residual generators for fault detection 
    and isolation. 16th International Conference on Control Systems and  
    Computer Science, Bucharest, Romania, 2007.

[3] A. Varga, On computing least order fault detectors using rational nullspace bases. 
    IFAC SAFEPROCESS'03 Symposium, Washington DC, USA, 2003.
"""
function efdisyn(sysf::FDIModel{T}, SFDI::Union{BitMatrix,BitVector,Array{Bool,2},Array{Bool,1}} = trues(1,length(sysf.faults)); rdim::Union{Vector{Int},Int,Missing} = missing, poles::Union{AbstractVector,Missing} = missing, 
                      sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing, 
                      nullspace::Bool = true, minimal::Bool = true, simple::Bool = false, separate::Bool = false, 
                      FDtol::Real = 0.0001, FDGainTol::Real = 0.01, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                      tcond::Real = 1.e4, HDesign::Union{Vector{Matrix{T1}},Missing} = missing,
                      offset::Real = sqrt(eps(float(real(T)))), atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                      rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T)))))*iszero(max(atol1,atol2)), 
                      fast::Bool = true) where {T, T1 <: Real}
   Ts = sysf.sys.Ts                  
   disc = (Ts != 0);  # system type (continuous- or discete-time)
   
   # decode input information
   inpu = sysf.controls; mu = length(inpu)  
   inpd = sysf.disturbances; md = length(inpd) 
   inpf = sysf.faults; mf = length(inpf)  
   inpw = sysf.noise;  mw = length(inpw) 
   inpaux = sysf.aux;  maux = length(inpaux)  
   
   m = mu+md+mf+mw+maux;       # total number of inputs
   p = size(sysf.sys,1);       # number of measurable outputs
    
   if mf == 0 && minimal
      @warn "Minimal synthesis option not feasible in the case of no faults"
      minimal = false
   end
   
   # tolerance for rank tests 
   
   strongFD = !ismissing(FDfreq)
   strongFD && !isa(FDfreq,Vector) && (FDfreq = [FDfreq]) 
   strongFD && (lfreq = length(FDfreq))

   # set default stability degree
   if disc
      sdegdefault = 0.95;
   else
      sdegdefault = -0.05;
   end
   
   
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
   
   if isa(SFDI,Vector)
      nb = 1; mf2 = length(SFDI)
   else
      nb, mf2 = size(SFDI,1), size(SFDI,2) 
   end     # number of filters 
   mf == mf2 || error("number of faults must be equal to the column dimension of SFDI")

   if !ismissing(rdim)
      if isa(rdim,Vector) 
         length(rdim) == nb || error("dimension of rdim must be equal to the row dimension of SFDI")
         minimum(rdim) > 0 || error("all components of rdim must be positive")
      else
         rdim > 0 || error("rdim must be positive")
         rdim = fill(rdim, nb)
      end
   end
   
   # imposed design option to form linear combinations of basis vectors
   emptyHD = ismissing(HDesign)
   if !emptyHD
      size(HDesign,1) == nb || error("number of HDesign components must be equal to the row dimension of SFDI")
      if !ismissing(rdim) 
         for i = 1:nb
             if isempty(HDesign[i])
                mH = size(HDesign[i],1)
                if mH > 0
                   mH == rdim[i] || error("row dimension of HDesign[$i] must be equal to rdim[$i]")
                   mH == rank(HDesign[i]) || error("HDesign[$i] must have full row rank")
                end
             end
         end
      end
   end

   tcond1 = similar(Vector{T},nb)
   degs1 = similar(Vector{Vector{Int}},nb)
   HDesign1 = similar(Vector{Array{T,2}},nb)
   if separate
      Qt = similar(Vector{DescriptorStateSpace{T}},nb)
      Rt = similar(Vector{DescriptorStateSpace{T}},nb)
      for i = 1:nb
         indd = Vector(1:mf)[SFDI[i,:] .== false] 
         indf = Vector(1:mf)[SFDI[i,:] .== true] 
         # pack [Gu [Gd Gf1] Gf2 [Gf Gw Gaux]]
         sysc = fdimodset(sysf.sys, c = 1:mu, d = [Vector(mu .+ (1:md)); (mu+md) .+ indd], 
                          f = (mu+md) .+ indf, aux = (mu+md) .+ (1:mf+mw+maux));
         # solve the corresponding EFDP           
         Qti, Rti, infoi = try 
            efdsyn(sysc; rdim = rdim[i], HDesign = ismissing(HDesign) ? missing : HDesign[i], atol1, atol2, atol3, sdeg, smarg, poles, minimal,
                                      FDtol, FDfreq, FDGainTol, simple, tcond, offset); 
         catch err
            isnothing(findfirst("empty",string(err))) ? rethrow() : error("empty nullspace basis: the $i-th EFDIP is not solvable")
         end
         Qt[i] = Qti.sys
         Rt[i] = Rti.sys[:,Rti.aux]
         tcond1[i] = infoi.tcond
         degs1[i] = infoi.degs
         HDesign1[i] = infoi.HDesign
      end
      Q = FDIFilter(Qt, p, mu)
      R = FDIFilterIF(Rt,0,0,mf,mw,maux)
  else     
      # Step 1): nullspace based reduction
      #
      desc = (sysf.sys.E != I)
      m2 = mf+mw+maux
      sdegNS = strongFD ? sdegdefault : missing
      if nullspace || md > 0 || (desc && rcond(sysf.sys.e) < 1.e-7 )
         # form [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
         #syse = [sysf(:,[inpu inpd inpf inpw inpaux]); eye(mu,m)];
         syse = [sysf.sys; eye(mu,m)];
         #
         # compute a left nullspace basis Q = Q1 of G1 = [Gu Gd; I 0] = 0 and
         # obtain QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = Q*[Gf Gw Gaux;0 0 0]
         QR, info1 = glnull(syse, m2; atol1, atol2, rtol, fast, sdeg = sdegNS, offset) 
         tcond0 = info1.tcond
      elseif mu == 0 && md == 0
         # compute minimal basis as Q = Q1 = I  and set
         # QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = [ Gf Gw Gaux ]
         QR = [eye(p) sysf.sys]
         tcond0 = 1.
      else
         # compute minimal basis as Q = Q1 = [ I -Gu] and set
         # QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = [ Gf Gw Gaux ]
         QR = [ eye(p) dss(sysf.sys.A, sysf.sys.E, [-sysf.sys.B[:,inpu] sysf.sys.B[:,[inpf; inpw; inpaux]]],
                sysf.sys.C, [-sysf.sys.D[:,inpu] sysf.sys.D[:,[inpf; inpw; inpaux]]]; Ts)] 
         # perform stabilization if strong detectability has to be enforced
         strongFD  && (QR = glcf(QR; atol1, atol2, atol3, rtol, fast)[1]) 
         tcond0 = 1.
      end
      
      # check solvability conditions
      size(QR,1) == 0 && error("empty nullspace basis: the EFDIP is not solvable")
   
      # Step 2): of Procedure EFDI 
      # initialize overall filters Q and R
      QRt = similar(Vector{typeof(QR)},nb)
      for i = 1:nb
          indd = Vector(1:mf)[SFDI[i,:] .== false] 
          indf = Vector(1:mf)[SFDI[i,:] .== true] 
          # pack [Rfd1 Rff2 [Q1 Rf1 Rw1 Raux1]]
          sysc = fdimodset(QR, d = (p+mu) .+ indd, f = (p+mu) .+ indf, aux = Vector(1:p+mu+mf+mw+maux))
          # determine [Q1i*Rff2 [Q1i*Q1 Q1i*Rf1 Q1i*Rw1 Q1i*Raux1]]
          _, QRauxi, infoi = try 
             efdsyn(sysc; rdim = rdim[i], HDesign = ismissing(HDesign) ? missing : HDesign[i], atol1, atol2, atol3, sdeg, smarg, poles, minimal,
                                       FDtol, FDfreq, FDGainTol, simple, tcond, offset); 
          catch err
             isnothing(findfirst("empty",string(err))) ? rethrow() : error("empty nullspace basis: the $i-th reduced EFDP is not solvable")
          end
          # extract [Q1i*Q1 Q1i*Rf1 Q1i*Rw1 Q1i*Raux1 ]
          QRt[i] = QRauxi.sys[:,QRauxi.aux]
          # transform to standard state-space
          #Rt[i] = gss2ss(QRauxi.sys[:,QRauxi.aux]; atol1, atol2, rtol)[1]
          #QRt[i] = gss2ss(QRauxi.sys[:,[Vector(mf+1:mf+p+mu+mw+maux); Vector(1:mf)]]; atol1, atol2, rtol)[1]
          #QRt[i] = QRauxi.sys[:,[Vector(mf+1:mf+p+mu+mw+maux); Vector(1:mf)]]
          #QRi = QRauxi.sys[:,QRauxi.aux];   # extract [Q1i*Rf1 Q1i*Q1 Q1i*Rw1 Q1i*Raux1 ]
          #Qt[i] = QRi[:,1:p+mu]
          #Qt[i] = QRauxi[:,QRauxi.aux[1:p+mu]]  # extract Q1i*Q1
          #Rt[i] = QRi[:,(p+mu) .+ (1:mf+mw+maux)]
          #Rt[i] = QRauxi[:,QRauxi.aux[(p+mu) .+ (1:mf+mw+maux)]]  # extract [Q1i*Rf1 Q1i*Rw1 Q1i*Raux1 ]
          tcond1[i] = max(tcond0,infoi.tcond)
          degs1[i] = infoi.degs
          HDesign1[i] = infoi.HDesign
       end
       Q = FDIFilter(QRt, p, mu)
       R = FDIFilterIF(QRt,0,0,mf,mw,maux; moff = p+mu)
   end   
   info = (tcond = tcond1, degs = degs1, HDesign = HDesign1)

   return Q, R, info

   # end EFDISYN
end
function afdbasesel(S::BitArray, rwgain::Matrix, degs::Vector{Int}, rdim::Int, nout::Int, simple::Bool, atol::Real)
   #   afdbasesel(S, rwgain, degs, rdim, nout, simple, tol) -> (seli, selord)
   #
   #   Select admissible basis vectors for solving the approximate fault detection problem (AFDP)
   #   using the binary structure matrix `S`, corresponding to `nvec` basis vectors, and 
   #   the frequency gain matrix `rwgain`. 
   #   `S` is a `nvec x mf x n` binary 3-dimensional array or a `nvec x mf` 
   #   binary 2-dimensional array (reshaped internally as an `nvec x mf x 1` 
   #   binary 3-dimensional array).
   #   `seli` contains `nout`-touples of indices of basis vectors whose linear combination 
   #   is admissible, i.e. , `S[seli[i],:,k]` has all columns nonzero for all `k = 1:n` and,
   #   if `rwgain` is nonzero, then  `rwgain[seli[i,:],:] has full row rank.   
   #   This ensures that the AFDP is solvable by using fault detection filters with `rdim <= nout` outputs. 
   #   If the associated `nvec` degrees contained in `degs` are provided, then
   #   `selord[i]` is the corresponding tentatively achievable least filter order.
   #   If `simple = true`, a simple basis is assumed, in which case, `degs[i]` is 
   #   also the order of the `i`-th basis vector. If `simple = false`, a minimum 
   #   rational basis is assumed. `selord` is empty if `degs` is empty. 

   #   Method: The selection approach is used in conjunction with the synthesis 
   #   Procedure AFD described in [1]. 

   # References:
   # [1] Varga A.
   #     Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017.


   ndim = ndims(S)
   if ndim == 3 
      nvec = size(S,1)  # numbers of vectors
      n = size(S,3)     # number of frequencies
   else
      nvec, mf = size(S);   # numbers of vectors, number of faults  
      n = 1
      S = reshape(S,nvec,mf,1)
   end
   nd = length(degs)
   nodegs = (nd == 0)

   nodegs || length(degs) == nvec || error("the dimension of degs must be equal to the number of rows of S")
   
   (rdim >=1 && rdim <= nvec) || error("ndim must have a positive value not exceeding $nvec")
   (nout >=rdim && nout <= nvec) || error("nout must have a value at least $rdim and at most $nvec")
   
   nvec == 1 && (return [1], nodegs ? Int[] : degs )
   
   # determine rank of rwgain
   if size(rwgain,2) > 0 
      rw = atol > 0 ? rank(rwgain; atol) : rank(rwgain)
      rw > 0 && rdim > rw && error("rdim must not exceed the rank of rwgain")
   else
      rw = 0
   end

   # find rdim combinations of nout vectors which solve the AFDP 
   seli = collect(combinations(Vector(1:nvec),nout))
   ni = length(seli)
   selord = nodegs ? Int[] : fill(-1,ni) 
   nqmax = sum(degs)
   ii = trues(ni)
   for i = 1:ni
       indv = seli[i];
       # check admissibility
       stest = true;
       for k = 1:n
           all(maximum(view(S,indv,:,k),dims=1)) || (stest = false; break)
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
       # check full rank admissibility condition
       rw > 0 && ii[i] &&
          ((atol > 0 ? rank(view(rwgain,indv,:); atol) : rank(view(rwgain,indv,:))) < rdim) && (ii[i] = false)
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
   # end AFDBASESEL
end
"""
       afdredsyn(sysfred; kwargs...) -> (Qred, Rred, info) 

Solve the approximate fault detection problem (AFDP) for a reduced synthesis model
`sysfred`. The computed stable and proper filter objects `Qred` and `Rred` contain the 
fault detection filter, representing the solution of the AFDP, and its internal form, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.HDesign` and `info.gap`, 
contains additional synthesis related information (see below). 

The continuous- or discete-time system `sysfred.sys` is in a standard
or descriptor state-space form `sysfred.sys = (A-λE,B,C,D)`, which corresponds to the input-output form  

       y = Gf(λ)*f + Gw(λ)*w + Ga(λ)*aux,

with the Laplace- or Z-transformed plant outputs `y`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Gf(λ)`, `Gw(λ)`, and `Ga(λ)` the corresponding 
transfer-function matrices. It is assumed that `Gf(λ)` has all columns nonzero and `Gw(λ)` has full-row rank, a
The indices of fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysf.faults`, `sysf.noise` and `sysf.aux`, respectively. 
The indices `sysf.controls` and `sysf.disturbances` of control and disturbance inputs are assumed void. 

The fault detection filter object `Qred`, contains in `Qred.sys` the resulting filter 
in a standard state-space form, which generates the residual signal `r`. 
The corresponding input-output (implementation) form is

            r = Qy(λ)*y ,               

where `Qy(λ)` isthe transfer function matrix from the outputs to the residual. 
The indices of outputs are contained in the integer vector 
`Qred.outputs` and `Qred.controls` is void.

The fault detection filter in internal form object `Rred`, contains `Rred.sys`, the resulting 
internal form of the filter 
in a standard state-space form, which generates the residual signal `r`, and corresponds to the 
input-output form

       r = Rf(λ)*f + Rw(λ)*w + Ra(λ)*aux ,

where 

       | Rf(λ) Rw(λ) Ra(λ) | = Qy(λ)*| Gf(λ) Gw(λ) Ga(λ) |. 

The solution of the AFDP ensures that `Rf(λ)` has all its columns nonzero
and the H∞-norm of `Rw(λ)` satisfies `||Rw(λ)||∞ < γ`, where the bound `γ` is 
specified via the keyword argument `gamma`.
The indices of the inputs `f`, `w` and `aux` of the resulting filter `Rre.sys` are 
contained in the integer vectors  
`R.faults`, `R.noise` and `R.aux`, respectively, while the indices 
of the inputs `u` and `d`, i.e., `R.controls` and `R.disturbances` are void.

_Note:_ This function is intended to performs the second synthesis step of the
AFD Procedure in [1] and `Ga(λ) = Q1(λ)`, where `Q1(λ)` is the transfer function of the 
filter resulted from a nullspace basis computation 
at the first synthesis step. The resulting `Ra(λ)` is the transfer function matrix of the 
updated filter `Qy(λ)*Q1(λ)`. Information on `Q1(λ)`, which are relevant for a least 
order synthesis are provided in the keyword arguments `simple`, `degs` and `S` 
"""
function afdredsyn(sysfred::Union{FDFilterIF{T},FDIModel{T}}; rdim::Union{Int,Missing} = missing, poles::Union{AbstractVector,Missing} = missing, 
   sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing, 
   minimal::Bool = true, simple::Bool = false,
   FDtol::Real = 0.0001, FDGainTol::Real = 0.01, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
   tcond::Real = 1.e4, HDesign::Union{AbstractMatrix,Missing} = missing,
   offset::Real = sqrt(eps(float(real(T)))), atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
   rtol::Real = ((size(sysfred.sys.A,1)+1)*eps(real(float(one(T)))))*iszero(max(atol1,atol2)), 
   fast::Bool = true, exact::Bool = false, S::Union{BitMatrix,BitVector,Array{Bool,2},Array{Bool,1},Missing} = missing, 
   degs::Union{Vector{Int},Missing} = missing, gamma::Real = 1, epsreg::Real = 0.1, 
   sdegzer::Union{Real,Missing} = missing, nonstd::Int = 1, freq::Real = rand()) where T
   #
   #  The OPTIONS structure allows to specify various user options, as
   #  follows:
   #  OPTIONS.tol       - relative tolerance for rank computations
   #                      (Default: internally determined value)
   #  OPTIONS.tolmin    - absolute tolerance for observability tests
   #                      (Default: internally determined value)
   #  OPTIONS.FDTol     - threshold for fault detectability checks
   #                      (Default: 0.0001))
   #  OPTIONS.FDGainTol - threshold for strong fault detectability checks
   #                      (Default: 0.01)
   #  OPTIONS.rdim      - desired number of residual outputs for Qred
   #                      (Default: number of outputs of the reduced system)  
   #  OPTIONS.FDFreq    - vector of real frequency values for strong  
   #                      detectability checks (Default: [])
   #  OPTIONS.smarg     - stability margin for the poles of filters Q and R
   #                      (Default: -sqrt(eps) for a continuous-time system 
   #                               and 1-sqrt(eps) for a discete-time system) 
   #  OPTIONS.sdeg      - prescribed stability degree for the poles of the 
   #                      filters Q and R
   #                      (Default:  -0.05 for a continuous-time system and
   #                                  0.95  for a discete-time system) 
   #  OPTIONS.poles     - complex vector containing a complex conjugate set  
   #                      of desired poles (within the stability domain) 
   #                      to be assigned for the filters Q and R (Default: []) 
   #  OPTIONS.simple    - option to employ a simple proper basis for synthesis 
   #                      true  - use a simple basis; the orders of the  
   #                              basis vectors are provided in INFO.deg
   #                      false - no simple basis computed (default) 
   #  OPTIONS.minimal   - option to perform a least order filter synthesis 
   #                      true  - perform least order synthesis (default)
   #                      false - perform full order synthesis   
   #  OPTIONS.exact     - option to perform exact filter synthesis 
   #                      true  - perform exact synthesis 
   #                      false - perform approximate synthesis (default)  
   #  OPTIONS.tcond     - maximum alowed condition number of the employed 
   #                      non-orthogonal transformations (Default: 1.e4).
   #  OPTIONS.HDesign   - full row rank design matrix H to build OPTIONS.rdim 
   #                      linear combinations of the left nullspace basis 
   #                      vectors (Default: [])
   #  OPTIONS.gamma     - upper bound on ||Rw||_inf  (Default: 1).
   #  OPTIONS.epsreg    - regularization parameter   (Default: 0.1).
   #  OPTIONS.sdegzer   - prescribed stability degree for zros shifting  
   #                      (Default:  -0.05 for a continuous-time system and
   #                                  0.95  for a discete-time system) 
   #  OPTIONS.nonstd    - option to handle nonstandard optimization problems
   #                      1 - use the quasi-co-outer-co-inner factorization
   #                          (Default: [])
   #                      2 - use the modified co-outer-co-inner factorization
   #                          with the regularization parameter OPTIONS.epsreg
   #                      3 - use the Wiener-Hopf type co-outer-co-inner 
   #                          factorization
   #                      4 - use the Wiener-Hopf type co-outer-co-inner 
   #                          factorization with zero shifting of the 
   #                          non-minimum phase factor using the
   #                          stabilization parameter OPTIONS.sdeg 
   #                      5 - use the Wiener-Hopf type co-outer-co-inner 
   #                          factorization with the regularization of the 
   #                          non-minimum phase factor using the
   #                          regularization parameter OPTIONS.epsreg
   #  OPTIONS.degs      - increasingly ordered degrees of a left minimal   
   #                      polynomial nullspace basis of G := [ Gu Gd; I 0] 
   #                      (also the left Kronecker indices of G), if the 
   #                      state-space realization of [Gu Gd ] is minimal                      
   #  OPTIONS.S         - binary structure matrix corresponding to the  
   #                      computed left nullspace basis.
   #
   #  INFO is a structure containing additional information, as follows: 
   #  INFO.tcond        - maximum of the condition numbers of the employed 
   #                      non-orthogonal transformation matrices; a warning is 
   #                      issued if INFO.tcond >= OPTIONS.tcond.
   #  INFO.HDesign      - the design matrix H employed for the synthesis of 
   #                      the fault detection filter.
   #  INFO.gap          - achieved gap min_j ||Rf_j||_inf/||Rw||_inf.
   
   #  Copyright 2018 A. Varga
   #  Author:    A. Varga, 04-05-2018.
   #  Revisions: A. Varga, 25-08-2018,11-06-2019.
   #
   #  Method: The Procedure AFD from [1] is implemented, which is based 
   #  on the synthesis method proposed in [2]. For the regularization approach 
   #  based on the modified co-outer-co-inner factorization see [3]. For the
   #  exact synthesis, the Procedure EFD from [1] is performed. 
   #
   #  References:
   #  [1] A. Varga
   #      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
   #      Springer Verlag, 2017; sec. 5.3.
   #  [2] A. Varga
   #      General computational approach for optimal fault detection. 
   #      Proc. IFAC Symposium SAFEPROCESS, Barcelona, Spain, pp. 107–112, 2009.
   #  [3] K. Glover and A. Varga
   #      On solving non-standard H-/H_2/inf fault detection problems. 
   #      Proc. IEEE CDC, Orlando, FL, USA, pp. 891–896, 2011.

   nvec, m = size(sysfred.sys)        
   # nvec - number of reduced outputs
   # m    - total number of inputs
   
   # decode input information 
   length(sysfred.controls) == 0 || error("the reduced system must not have control inputs")
   length(sysfred.disturbances) == 0 || error("the reduced system must not have disturbance inputs")
   inpf = sysfred.faults; mf = length(inpf)
   inpw = sysfred.noise; mw = length(inpw)
   inpaux = sysfred.aux;  maux = length(inpaux)  
   indf = inpf  
   Ts = sysfred.sys.Ts                  
   disc = (Ts != 0);  # system type (continuous- or discete-time)

   strongFD = !ismissing(FDfreq)
   strongFD && !isa(FDfreq,Vector) && (FDfreq = [FDfreq]) 
   ismissing(S) && (S = strongFD ? 
             fdisspec(sysfred, FDfreq; stabilize = true, FDGainTol, atol1, atol2, rtol, fast) :
             fditspec(sysfred; block = false, FDtol, atol1, atol2, rtol, fast) )


   # set default stability degree
   sdegdefault = disc ? 0.95 : -0.05
   
   ismissing(sdeg) && (sdeg = sdegdefault)  # set desired stability degree to default value
   
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

   # imposed design option to form linear combinations of basis vectors
   emptyHD = ismissing(HDesign)
   if !emptyHD
      !ismissing(rdim) && size(HDesign,1) != rdim && error("row dimension of HDesign must be equal to rdim")
      size(HDesign,1) == rank(HDesign) || error("HDesign must have full row rank")
   end   
      
   # information on nullspace basis vector degrees
   ismissing(degs) || length(degs) == nvec ||
        error("dimension of degs must be equal to the number of reduced outputs")
   ismissing(degs) && (degs = Int[])          
   # structure matrix corresponding to employed basis
   ismissing(S) || size(S,2) == mf ||
         error("column dimension of S must be equal to the number of faults")
          
   # desired stability degree for zero shifting
   ismissing(sdegzer) && (sdegzer = sdegdefault)
     
   
   # # set options for LCF-based stabilization to be used for admisibility checks
   # opts_glcf_default = struct('tol',tol,'tolmin',tolmin);
   # # set options for LCF-based stabilization to be used for final synthesis
   # opts_glcf = struct('tol',tol,'tolmin',tolmin, ...
   #                    'sdeg',sdeg,'smarg',smarg,'poles',poles);
   
   tcond1 = 1                 # condition number of employed transformations
   nq = order(sysfred.sys)     # order of sysfred
   
   # set H for checking the solvability condition
   if emptyHD
      Htemp = eye(nvec)
   else
      degs = Int[];
      rdim = size(HDesign,1); nh = size(HDesign,2)
      if nh < nvec
         # pad with zeros: row rank is preserved
         Htemp = [ HDesign zeros(T, rdim, nvec-nh) ]
      else
         # remove trailing columns if necessary: rank may drop
         Htemp = HDesign[:,1:nvec];
         if nh > nvec && rdim != rank(Htemp)
            error("The leading part of HDesign must have full row rank")
         end
      end
   end
   
   # # setup the number of filter outputs
   # if isempty(rdim)
   #    if emptyHD 
   #       rdim = 1;   
   #    else
   #       rdim = size(HDesign,1);
   #    end
   # else
   #    rdim = min(rdim,nvec); 
   # end
   # setup the number of filter outputs
   if minimal
      #  least order design
      if ismissing(rdim)
         if emptyHD 
            rdim = 1   
         else
            rdim = size(HDesign,1)
         end
      else
         rdim = min(rdim,nvec)
      end
   else
      #  full order design
      if ismissing(rdim) 
         if emptyHD 
            rdim = nvec
         else
            rdim = min(size(HDesign,1),nvec);
         end
      else
         if mf == 0
            if rdim < nvec && emptyHD
               @warn "rdim reset to $nvec"
               rdim = nvec
            end
         else
            rdim = min(rdim,nvec); 
         end
      end
   end

   
   # no check of the solvability condition needed
   rw = (mw == 0) ? 0 : gnrank(Htemp*sysfred.sys[:,inpw]; atol, rtol) 
          
   # Compute admissible Q2 to reduce the order of Q2*Q;  
   # update Q <- Q2*Q, Rf = Q2*Gf 
   
   sysfredupd = [sysfred.sys I]
   #InpG = sysfredupd.InputGroup; 
   
   # reorder degs to correspond to the expected orders of basis vectors 
   # corresponding to the actual order of outputs of QR 
   ismissing(degs) || reverse!(degs) 
   if rdim < nvec 
      # determine possible low order syntheses using i >= rmin basis vectors
      # and the corresponding expected orders    
      
      finish = false;    # set termination flag
      nout = rdim;       # initialize number of selected basis vectors
      if mw > 0 && !exact
         rwgain = evalfr(Htemp*sysfredupd[:,inpw],freq);
      else
         rwgain = zeros(size(Htemp,1),0)  # set rwgain an empty matrix
      end
      itry = 1; 
      while !finish     
          # choose nout basis vectors, which potentially lead to a least order
          # filter with rdim outputs:
          # basesel(i,:) contains the indices of candidate basis vectors;
          # ordsel(i)    contains the presumably achievable least orders
          basesel, ordsel = afdbasesel(S,rwgain,degs,rdim,nout,simple,atol)
          #
          # update the synthesis using the selections of candidate vector(s),
          # starting with the least (potentially) achievable order
          for i = 1:size(basesel,1);
              baseind = basesel[i] # indices of current basis selection
              if rdim == nout
                  hbase = eye(rdim)
              else
                  hbase = rand(rdim,nout) 
              end
              ip = [baseind; setdiff(1:nvec,baseind)][:]
              if simple
                 if minimal
                    if emptyHD 
                       # select vectors and elliminate unobservable dynamics  
                       noelim = falses(nq) 
                       #ell = sum(degs(1:basesel(i,1)-1)); 
                       ell = sum(degs[1:basesel[i][1]-1]); 
                       for jj = 1:nout 
                           ellnext = sum(degs[1:baseind[jj]]);
                           noelim[ell+1:ellnext] .= true;
                           ell = ellnext
                       end
                    end
                    if rdim == nout
                       if emptyHD 
                          #QRfwtest = modred(sysfredupd(baseind,:),~noelim,'truncate');
                          #h = Htemp(ip(1:rdim),:);
                          #QRfwtest = modred(QR[baseind,:],~noelim,'truncate');
                          ir = noelim
                          Ar, Er, Br, Cr, Dr = dssdata(sysfredupd[baseind,:])
                          QRfwtest = dss(view(Ar,ir,ir), Er == I ? I : view(Er,ir,ir),view(Br,ir,:),view(Cr,:,ir), Dr; Ts)
                          h = Htemp[ip[1:rdim],:]
                       else
                          QRfwtest = gir(Htemp*sysfredupd; atol1, atol2, rtol)
                       end
                    else
                       # this case is possible only if HDesign is empty
                       # build rdim linear combinations of the first nout vectors 
                       ir = noelim
                       Ar, Er, Br, Cr, Dr = dssdata(sysfredupd[baseind,:])
                       QRfwtest = hbase*dss(view(Ar,ir,ir),Er == I ? I : view(Er,ir,ir),view(Br,ir,:),view(Cr,:,ir), Dr; Ts)
                       h = [ hbase zeros(T, rdim, nvec-nout) ]
                       h = h[:,ip]  # permute columns to match unpermuted QR 
                    end
                 else
                    if rdim == nout
                       if emptyHD 
                          h = Htemp[ip[1:rdim],:]
                          QRfwtest = gir(sysfredupd[baseind,:]; atol1, atol2, rtol, infinite = false)
                       else
                          QRfwtest = gir(Htemp*sysfredupd; atol1, atol2, rtol, infinite = false)
                       end
                    else
                       QRfwtest = gir(hbase*sysfredupd[baseind,:]; atol1, atol2, rtol, infinite = false) 
                       h = [ hbase zeros(T, rdim, nvec-nout) ] 
                       h = h[:,ip]  # permute columns to match unpermuted QR 
                     end
                 end
              else
                 if minimal
                    # build output permutation vector for glmcover1  
                    if rdim == nout
                       if emptyHD 
                          QRfwtest, _, info2 = glmcover1(sysfredupd[ip,:], rdim; atol1, atol2, rtol)
                          if !isempty(ordsel) && (order(QRfwtest) != ordsel[i])
                             @warn "afdsyyn: expected reduced order not achieved"
                          end
                          h = Htemp[ip[1:rdim],:]
                       else
                          QRfwtest, _, info2 = glmcover1([Htemp; eye(nvec)]*sysfredupd[ip,:], rdim; atol1, atol2, rtol)
                       end
                    else  
                       # this case is possible only if HDesign is empty
                       # build rdim linear combinations of the first nout vectors 
                       h = [ hbase zeros(T, rdim, nvec-nout) ]; 
                       QRfwtest, _, info2 = glmcover1([h; eye(nvec)]*sysfredupd[ip,:], rdim; atol1, atol2, rtol)
                       h = h[:,ip]  # permute columns to match unpermuted QR 
                    end
                 else
                    if rdim == nout
                       if emptyHD
                          h = Htemp[ip[1:rdim],:]
                          QRfwtest = gir(sysfredupd[baseind,:]; atol1, atol2, rtol, infinite = false) 
                       else
                          QRfwtest = gir(Htemp*QR; atol1, atol2, rtol, infinite = false) 
                       end
                    else
                       QRfwtest = gir(hbase*sysfredupd[baseind,:]; atol1, atol2, rtol, infinite = false) 
                       h = [ hbase zeros(T, rdim, nvec-nout) ]; 
                       h = h[:,ip]  # permute columns to match unpermuted sysfredupd 
                    end
                 end
              end
              # check admissibility of the current design; 
              if (rdim == nout && minimal) || rdim < nout
                 # dismiss design if check fails
                 if strongFD 
                    Stest = fdisspec_(QRfwtest[:,indf], FDfreq; stabilize = true, block = true, 
                                              FDGainTol, atol1, atol2, atol3, rtol = 0, fast)[1]
                 else
                    Stest = fditspec_(QRfwtest[:,indf]; block = true, atol1, atol2, rtol, FDtol)
                 end
                 rwtest = rw > 0 ? rank(evalfr(QRfwtest[:,inpw],freq; atol, rtol)) : rdim
                 # check complete fault detectability of the current design 
                 # and full row rank condition  
                 if all(Stest) && rdim == rwtest
                    if !simple && minimal
                       # adjust condition number of employed transformations
                       tcond1 = max(tcond1, info2.fnorm, info2.tcond)
                       tcond1 > tcond && 
                          @warn "efdsyn: possible loss of numerical stability due to ill-conditioned transformations"
   #                     if ~emptyHD
   #                        info.HDesign = Htemp;
   #                     end
                    end
                    sysfredupd = QRfwtest
                    finish = true
                    break
                 end
              else
                 sysfredupd = QRfwtest
                 finish = true
                 break
              end
          end
          nout += 1
          if nout > nvec
             if itry > 5
                finish = true
                @warn "fault detectability not achieved with the chosen number of residuals"
             else
                itry += 1
                nout -= 1
             end
          end
      end
      if emptyHD
         Htemp = h
      end
   else
      hbase = eye(rdim)
      if simple
         baseind = 1:rdim 
      else
         baseind = 1;
      end
      h = eye(rdim)
      if !emptyHD
         sysfredupd = Htemp*sysfredupd
      else
         # use full minimum basis 
         Htemp = h
      end
   end
   
   # compute M such that M*Q has a desired stability degree;  
   # update Q <- M*Q and R <- M*R 
   # this operation is performed only if rank is null or for exact synthesis
   if rw == 0 || exact
      k = 1;
      if simple && isequal(hbase,I) && emptyHD 
         # exploit the block diagonal structure of basis matrices al and cl
         # to compute block-diagonal M
         al, el, bl, cl, dl, = dssdata(sysfredupd)
         for i = 1:length(baseind)
            blkord = degs[baseind[i]]
            if blkord
               i1 = k:k+blkord-1; 
               QRfwi = glcf(dss(al[i1,i1],el[i1,i1],bl[i1,:],cl[i:i,i1],dl[i:i,:];Ts); 
                            atol1, atol2, atol3, sdeg, smarg, evals = poles)[1]
               al[i1,i1] = QRfwi.A; bl[i1,:] = QRfwi.B;  cl[i,i1] = QRfwi.C; dl[i,:] = QRfwi.D 
               el[i1,i1] = (QRfwi.e == I) ? eye(blkord) : QRfwi.E  
               k += blkord
            end
         end
         sysfredupd = dss(al, el, bl, cl, dl; Ts)
      else
         sysfredupd = glcf(sysfredupd; atol1, atol2, atol3, sdeg, smarg, evals = poles)[1]
      end
   end
       
   # finish if no noise input or if all noise inputs are decoupled or
   # exact synthesis is performed
   if mw == 0 || rw == 0 || exact
      # scale Rf to ensure unit minimum column gains
      if mf > 0
         if strongFD && minimum(FDfreq) == 0
            # compute minimum DC gains  
            dcg = dcgain(sysfredupd[:,indf]; atol1, atol2, rtol, fast) 
            y = maximum(abs.(dcg),dims=1)
            indj = sortperm(y[:])[1] 
            scale = y[indj]
            indi = sortperm(abs.(dcg[:,indj]))[end]
            sc = sign(dcg[indi,indj])/scale
         else
            # compute the minimum of H-inf norms of columns
            sc = 1/fdhinfminus(sysfredupd[:,indf])[1]
         end
         sysfredupd = sc*sysfredupd
      end
      
      # transform to standard state-space
      sysfredupd = gss2ss(sysfredupd; atol1, atol2, rtol)[1]
      # set(sysfredupd,'InputGroup',InpG)

      Qred = FDFilter(sysfredupd[:,end-nvec+1:end], nvec, 0)
      Rred = FDFilterIF(sysfredupd,0,0,mf,mw,maux)
      if rw > 0
         beta = fdhinfminus(sysfredupd[:,inpf],FDfreq)[1]
         gap = beta/ghinfnorm(sysfredupd[:,inpw])[1]
      else
         gap = Inf
      end
   
      info = (tcond = tcond1, HDesign = convert(Matrix{Float64},Htemp), gap = gap)

      return Qred, Rred, info
   end
   
   # determine the optimal factor Q3 to minimize the gap
   # and update Q <- Q3*Q and R <- Q3*R 
     
   # compute the extended quasi-co-outer-co-inner factorization  
   # Rw = [Rwoe 0]*Rwi = Rwoe*Rwi1
   Rwi, Rwo, info1 = goifac(sysfredupd[:,inpw]; atol1, atol2, atol3, rtol)
   rw = size(Rwo,2);  # rank of Rw
   nonstandard = info1.nfuz+info1.niuz > 0; 
   
   # handle different cases
   if nonstandard
      # non-standard case: zeros on the boundary of the stability domain 
      if nonstd == 2
         # use modified co-outer-co-inner factorization
            _, Rwo = goifac([Rwo epsreg*eye(size(Rwo,1))]; atol1, atol2, atol3, rtol)
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
         Rwouz, Rwo = gcrange(Rwo; atol1, atol2, rtol, zeros = "s-unstable")
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
         Rwouz = gbilin(Rwouz,sys1; minimal = true, atol1, atol2, rtol)[1] 
      elseif nonstd == 5  
         # use modified Wiener-Hopf type co-outer-co-inner factorization
         # with regularization of the non-minimum phase factor
         # separate stable and strictly unstable zeros Rwo = Rwoe*Rwouz 
         # and update Rwo <- Rwoe
         Rwouz, Rwo = gcrange(Rwo; atol1, atol2, rtol, zeros = "s-unstable")
         rw = size(Rwouz,1)
         _, Rwouz = goifac([Rwouz epsreg*eye(rw)]; atol1, atol2, atol3, rtol)
      end
   end
  
   if rw == size(Rwo,1)
      # Q3 = inv(Rwo)
      # extract descriptor state-space data
      aQR, eQR, bQR, cQR, dQR = dssdata(sysfredupd)
      # form [Rwo Rf Rw Raux Q] 
      RwoQR = dss(aQR, eQR, [Rwo.B bQR],cQR,[Rwo.D dQR];Ts)
      # form QR = inv(Rwo)*[ Rf Rw Raux Q] 
      sysfredupd = grsol(RwoQR, m+nvec; atol1, atol2, rtol)[1]
      if nonstandard && (nonstd == 4 || nonstd == 5)
         sysfredupd = gminreal(Rwouz\sysfredupd; atol1, atol2, rtol)
      end
   else
     # regularization for non-invertible Rwo (this case should never occur)
      # Q3 = Rwoinv, where Rwoinv is a left inverse of Rwo
      # with stable spurious poles
      Rwoinv = glsol([Rwo;eye(rw)], rw; atol1, atol2, rtol, sdeg)[1]
      # form QR = Rwoinv*[Rf Rw Raux Q] 
      sysfredupd = gir(Rwoinv*sysfredupd; atol1, atol2, rtol, infinite = false)
      if nonstandard && (nonstd == 4 || nonstd == 5)
         sysfredupd = gminreal(Rwouz\sysfredupd; atol1, atol2, rtol)
      end
   end
   if nonstandard && nonstd == 1
      # perform stabilization 
      # determine Q4 such that Q <- Q4*Q and R <- Q4*R are stable
      sysfredupd = glcf(sysfredupd; atol1, atol2, atol3, rtol, sdeg, smarg, evals = poles)[1]
   end
   
   # scale to enforce ||Rw||_inf = gamma
   scale = gamma/ghinfnorm(sysfredupd[:,inpw])[1]
   sysfredupd = gss2ss(scale*sysfredupd)[1]
  
   Qred = FDFilter(sysfredupd[:,end-nvec+1:end], nvec, 0)
   Rred = FDFilterIF(sysfredupd,0,0,mf,mw,maux)
   beta = fdhinfminus(sysfredupd[:,inpf],FDfreq)[1]
   info = (tcond = tcond1, HDesign = convert(Matrix{Float64},Htemp), gap = beta/gamma)
   
   return Qred, Rred, info
   # end AFDREDSYN
   end
   
   
   