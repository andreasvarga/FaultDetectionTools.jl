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

The continuous- or discrete-time system `sysf.sys` is in a standard
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

`FDfreq = freq` specifies a vector of real frequency values or a scalar real frequency value
for strong detectability checks (default: `FDfreq = missing`).

If `nullspace = true` (default), a minimal proper nullspace basis is used for the synthesis of the
fault detection filter. If `nullspace = false`, a full-order observer based nullspace basis is used. 
This option can be only used for a proper system without disturbance inputs. 

If `simple = true`, a simple proper nullspace basis is emplyed for synthesis. 
The orders of the basis vectors are provided in `info.deg`. 
If `simple = false` (default), then no simple basis is computed. 

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
   non-orthogonal transformation matrices; a warning is issued if `info.tcond >= tcmax`;

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
   disc = (Ts != 0);  # system type (continuous- or discrete-time)
   
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
   nvec == 0 && error("efdsyn: empty nullspace basis - the EFDP is not solvable")

   nq = order(QR)             # order of the minimal basis
   
   # set H for checking the solvability condition
   if emptyHD
      Htemp = eye(nvec)
      noutmax = nvec
   else
      degs = Int[];
      rdim = size(HDesign,1); nh = size(HDesign,2)
      noutmax = rdim
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
            all(maximum(S[:,:,ii],dims=1)) || error("efdsyn: strong detection of all faults not feasible")
         end
      else
         # check weak detectability conditions 
         S = fditspec_(Htemp*QR[:,indf]; atol1, atol2, rtol, FDtol)
         all(maximum(S,dims=1)) || error("efdsyn: detection of all faults not feasible")
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
   #if rdim < noutmax && mf > 0
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
                             @warn "efdsyn: expected reduced order not achieved"
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
          #if nout > nvec
          if nout > noutmax
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

The continuous- or discrete-time system `sysf.sys` is in a standard
or descriptor state-space form `sysf.sys = (A-λE,B,C,D)`, which corresponds to the input-output form  

       y = Gu(λ)*u + Gd(λ)*d + Gf(λ)*f + Gw(λ)*w + Ga(λ)*aux,

with the Laplace- or Z-transformed plant outputs `y`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Gu(λ)`, `Gd(λ)`, `Gf(λ)`, `Gw(λ)`, and `Ga(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysf.controls`, `sysf.disturbances`, `sysf.faults`, `sysf.noise` and `sysf.aux`, respectively.

The fault detection and isolation filter object `Q`, contains in `Q.sys` the resulting bank of `nb` filters. 
The `i`-th filter `Q.sys[i]` is in a standard state-space form and generates `r_i`, the `i`-th component (scalar or vector) 
of the overall residual vector `r := [r_1; r_2; ...; r_nb]`. 
The corresponding input-output (implementation) form of the `i`-th filter is

            r_i = Qyi(λ)*y + Qui(λ)*u   ,            

where `Qyi(λ)` and `Qui(λ)` are the transfer function matrices from the output and control inputs to the `i`-th residual component. 
The indices of output and control inputs are contained in the integer vectors 
`Q.outputs` and `Q.controls`, respectively.

The fault detection and isolation filter in internal form object `R`, contains `R.sys`, the resulting bank of `nb` 
internal form of the filters.  
The `i`-th filter `R.sys[i]` is in a standard state-space form, which generates the residual signal `r_i`, and corresponds to the 
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
If `q` is a scalar, then a vector `rdim` with all components equal to `q` is assumed.
The default value of `q[i]` is chosen as follows: if `HDesign = missing` or `H_i` is empty then  
`q[i] = 1`, if `minimal = true`, or `q[i]` is the number of the nullspace basis 
vectors used for the synthesis of `Q.sys[i]` and `R.sys[i]`, if `minimal = false`; 
if  `H_i` specifies a full row rank design matrix, then `q[i]` is the row dimension of `H_i`. 

`FDfreq = freq` specifies a vector of real frequency values or a scalar real frequency value
for strong detectability checks (default: `FDfreq = missing`).

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
a warning is issued if any `info.tcond[i] >= tcmax`;

`info.degs` is an `nb`-dimensional vector, whose `i`-th component is an integer vector 
containing the degrees of the basis vectors of the employed simple
nullspace basis for the synthesis of the `i`-th filter component, if `simple = true, 
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
function efdisyn(sysf::FDIModel{T}, SFDI::Union{BitMatrix,BitVector,Array{Bool,2},Array{Bool,1}} = trues(1,length(sysf.faults)); 
                      rdim::Union{Vector{Int},Int,Missing} = missing, poles::Union{AbstractVector,Missing} = missing, 
                      sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing, 
                      nullspace::Bool = true, minimal::Bool = true, simple::Bool = false, separate::Bool = false, 
                      FDtol::Real = 0.0001, FDGainTol::Real = 0.01, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                      tcond::Real = 1.e4, HDesign::Union{Vector{Matrix{T1}},Missing} = missing,
                      offset::Real = sqrt(eps(float(real(T)))), atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                      rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T)))))*iszero(max(atol1,atol2)), 
                      fast::Bool = true) where {T, T1 <: Real}
   Ts = sysf.sys.Ts                  
   disc = (Ts != 0);  # system type (continuous- or discrete-time)
   
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
            findfirst("empty",string(err)) === nothing &&   
            findfirst("detection",string(err)) === nothing && 
            findfirst("condition",string(err)) === nothing &&  error("$err")  
            t1 = (mod(i,10) == 1 ? "$i-st" : "")
            t2 = (mod(i,10) == 2 ? "$i-nd" : "")
            t3 = (mod(i,10) == 3 ? "$i-rd" : "")
            t4 = (mod(i,10) > 3 ? "$i-th" : "")
            error("the $t1$t2$t3$t4 EFDIP is not solvable")
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
             efdsyn(sysc; rdim = ismissing(rdim) ? missing : rdim[i], HDesign = ismissing(HDesign) ? missing : HDesign[i], atol1, atol2, atol3, sdeg, smarg, poles, minimal,
                                       FDtol, FDfreq, FDGainTol, simple, tcond, offset); 
          catch err
             findfirst("empty",string(err)) === nothing &&   
             findfirst("detection",string(err)) === nothing && 
             findfirst("condition",string(err)) === nothing &&  error("$err")  
             t1 = (mod(i,10) == 1 ? "$i-st" : "")
             t2 = (mod(i,10) == 2 ? "$i-nd" : "")
             t3 = (mod(i,10) == 3 ? "$i-rd" : "")
             t4 = (mod(i,10) > 3 ? "$i-th" : "")
             error("the $t1$t2$t3$t4  EFDIP is not solvable")
          end
          # extract [Q1i*Q1 Q1i*Rf1 Q1i*Rw1 Q1i*Raux1 ]
          QRt[i] = QRauxi.sys[:,QRauxi.aux]
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
function afdredsyn(sysfred::Union{FDFilterIF{T},FDIModel{T}}; rdim::Union{Int,Missing} = missing, poles::Union{AbstractVector,Missing} = missing, 
   sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing, 
   minimal::Bool = true, simple::Bool = false,
   FDtol::Real = 0.0001, FDGainTol::Real = 0.01, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
   tcond::Real = 1.e4, HDesign::Union{AbstractMatrix,Missing} = missing,
   offset::Real = sqrt(eps(float(real(T)))), atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
   rtol::Real = ((size(sysfred.sys.A,1)+1)*eps(real(float(one(T)))))*iszero(max(atol1,atol2)), 
   fast::Bool = true, exact::Bool = false, S::Union{BitArray{3},BitMatrix,BitVector,Array{Bool,2},Array{Bool,1},Missing} = missing, 
   degs::Union{Vector{Int},Missing} = missing, gamma::Real = 1, epsreg::Real = 0.1, 
   sdegzer::Union{Real,Missing} = missing, nonstd::Int = 1, freq::Real = rand()) where T
   
"""
   afdredsyn(sysfred::FDIModel; rdim, simple = false, minimal = true, exact = false, 
   gamma = 1, epsreg = 0.1, sdegzer, nonstd = 1, freq, sdeg, smarg, poles, HDesign, 
   FDtol, FDGainTol, FDfreq, tcond, offset, atol, atol1, atol2, atol3, rtol, fast = true) 
      -> (Qred::FDFilter, Rred::FDFilterIF, info) 

Solve the approximate fault detection problem (AFDP) for a reduced synthesis model
`sysfred`. The computed stable and proper filter objects `Qred` and `Rred` contain the 
fault detection filter, representing the solution of the AFDP, and its internal form, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.HDesign`, `info.freq` and  `info.gap`, 
contains additional synthesis related information (see below). 

The continuous- or discrete-time system `sysfred.sys` is in a standard
or descriptor state-space form `sysfred.sys = (A-λE,B,C,D)`, which corresponds to the input-output form  

   y = Gf(λ)*f + Gw(λ)*w + Ga(λ)*aux,

with the Laplace- or Z-transformed plant outputs `y`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Gf(λ)`, `Gw(λ)`, and `Ga(λ)` the corresponding 
transfer-function matrices. It is assumed that `Gf(λ)` has all columns nonzero.
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
The indices of the inputs `f`, `w` and `aux` of the resulting filter `Rred.sys` are 
contained in the integer vectors  
`Rred.faults`, `Rred.noise` and `Rred.aux`, respectively, while the indices 
of the inputs `u` and `d`, i.e., `Rred.controls` and `Rred.disturbances` are void.

The resulting filters `Qred.sys` and `Rred.sys` have observable state-space 
realizations `(AQ,BQ,CQ,DQ)` and `(AQ,BR,CQ,DR)`, respectively, and thus  
share the observable pairs `(AQ,CQ)`. 

Various user options can be specified via keyword arguments as follows:

If `minimal = true` (default), a least order filter synthesis is performed, while 
with `minimal = false` a full order synthesis is performed.  

If `exact = true`, an exact synthesis (without optimization) is performed, while 
with `exact = false` (default), an approximate synthesis is performed.  

If `HDesign = H`, a design matrix `H` of full row rank `q` is used to build `q` 
linear combinations of the left nullspace basis vectors of 
`G(λ) := [ Gu(λ) Gd(λ); I 0]`.

`rdim = q` specifies the desired number `q` of residual outputs for `Qred` and `Rred`.

If `simple = true`, a simple proper nullspace basis is emplyed for synthesis. 
The orders of the basis vectors employed for the synthesis
are provided in `info.deg` (see below). 
If `simple = false` (default), then no simple basis is computed. 

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

`freq = val` specifies the values of a test frequency to be employed to check the 
full row rank admissibility condition (default: randomly generated in the interval `(0,1)`).

`tcond = tcmax` specifies the maximum alowed condition number `tcmax` of 
the employed non-orthogonal transformations (default: `tcmax = 1.e4`).

`FDtol = tol1` specifies the threshold `tol1` for fault detectability checks
   (default: `tol1 = 0.0001`).

`FDGainTol = tol2` specifies the threshold `tol2` for strong fault detectability checks
   (default: `tol2 = 0.01`). 

`gamma = γ` specifies the allowed upper bound for `∥Rw(λ)∥∞` (default: `γ = 1`).

`epsreg = ϵ` specifies the value of the regularization parameter `ϵ` (default: `ϵ = 0.1`)

`sdegzer = δ` specifies the prescribed stability degree `δ` for zeros shifting
   (default: `δ = −0.05` for a continuous-time system `sysf.sys` and 
   `δ = 0.95` for a discrete-time system `sysf.sys`).

`nonstd = job` specifies the option to handle nonstandard optimization problems, as follows:
      job = 1 – use the quasi-co-outer–co-inner factorization (default);
      job = 2 – use the modified co-outer–co-inner factorization with the
                regularization parameter ϵ;
      job = 3 – use the Wiener-Hopf type co-outer–co-inner factorization;
      job = 4 – use the Wiener-Hopf type co-outer-co-inner factorization with
                zero shifting of the non-minimum phase factor using the
                stabilization parameter δ;
      job = 5 – use the Wiener-Hopf type co-outer-co-inner factorization with
                the regularization of the non-minimum phase factor using the
                regularization parameter ϵ. 

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

The resulting named tuple `info` contains `(tcond, HDesign, freq, gap)`, where:

`info.tcond` is the maximum of the condition numbers of the employed 
   non-orthogonal transformation matrices; a warning is issued if `info.tcond >= tcmax`;

`info.HDesign` is the design matrix `H` employed for the synthesis of 
   of the fault detection filter;

`info.freq` is the frequency value employed to check the full 
row rank admissibility condition;

`info.gap` is the achieved gap `∥Rf(λ)∥∞−/∥Rw(λ)∥∞`, where the H−minus index is computed
over the whole frequency range, if `FDfreq = missing`, or over
the frequency values contained in `freq` if `FDfreq = freq`.

_Method:_ The Procedures EFD and AFD from [1] are implemented to solve 
the exact and approximate fault detection problems. 
For the details of the regularization techniques employed in 
the non-standard case see [2] and [3].

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017.

[2] K. Glover and A. Varga, On solving non-standard H-/H_2/inf fault detection problems. 
              Proc. IEEE CDC, Orlando, FL, USA, pp. 891–896, 2011.

[3] A. Varga, Fault Detection and Isolation Tools (FDITOOLS) User's Guide, 
              [arXiv:1703.08480](https://arxiv.org/pdf/1703.08480). 
"""

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
   disc = (Ts != 0);  # system type (continuous- or discrete-time)

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
   ismissing(degs) && (degs = Int[])  
   isempty(degs) || length(degs) == nvec ||
        error("dimension of degs must be equal to the number of reduced outputs")
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
         rwgain = evalfr(Htemp*sysfredupd[:,inpw],freq; atol1, atol2, rtol)
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
                             @warn "afdsyn: expected reduced order not achieved"
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
                          QRfwtest = gir(Htemp*sysfredupd; atol1, atol2, rtol, infinite = false) 
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
                 rwtest = rw > 0 ? rank(evalfr(QRfwtest[:,inpw], freq; atol1, atol2, rtol)) : rdim
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
   
      info = (tcond = tcond1, HDesign = convert(Matrix{Float64},Htemp), freq = freq, gap = gap)

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
   if nonstandard && nonstd != 1
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
   info = (tcond = tcond1, HDesign = convert(Matrix{Float64},Htemp), freq = freq, gap = beta/gamma)
   
   return Qred, Rred, info
   # end AFDREDSYN
end
"""
    afdsyn(sysf::FDIModel; rdim, nullspace = true, simple = false, minimal = true, exact = false, 
                           gamma = 1, epsreg = 0.1, sdegzer, nonstd = 1, freq, sdeg, smarg, poles, 
                           HDesign, HDesign2, scale2, FDtol, FDGainTol, FDfreq, 
                           tcond, offset, atol, atol1, atol2, atol3, rtol, fast = true) 
                              -> (Q::FDFilter, R::FDFilterIF, info)

Solve the approximate fault detection problem (AFDP) for a given synthesis model
`sysf` with additive faults. The computed stable and proper filter objects `Q` and `R` contain the 
fault detection filter, representing the solution of the AFDP, and its internal form, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.degs`, `info.degs2`, 
`info.S`, `info.S2`, `info.HDesign`, `info.HDesign2`, `info.freq` and `info.gap`
contains additional synthesis related information (see below). 

The continuous- or discrete-time system `sysf.sys` is in a standard
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

            r = Qy(λ)*y + Qu(λ)*u  ,             

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

The solution of the AFDP ensures that `Ru(λ) = 0`, `Rd(λ) = 0`, `Rf(λ)` has all its columns nonzero
and the H∞-norm of `Rw(λ)` satisfies `||Rw(λ)||∞ < γ`, where the bound `γ` is 
specified via the keyword argument `gamma`.
The indices of the inputs `u`, `d`, `f`, `w` and `aux` of the resulting filter `R.sys` are 
contained in the integer vectors `R.controls` (void), `R.disturbances` (void), 
`R.faults`, `R.noise` and `R.aux`, respectively.

The transfer function matrices `Q(λ) = [ Qy(λ) Qu(λ) ]` and `R(λ) = [ Ru(λ) Rd(λ) Rf(λ) Rw(λ) Ra(λ) ]` 
of the resulting filters `Q.sys` and `R.sys`, respectively,
have, in general, the partitioned forms

     Q(λ) = [ Q1(λ) ] ,   R(λ) = [ R1(λ) ] ,                      (1)
            [ Q2(λ) ]            [ R2(λ) ]

where the filters `Q1(λ)` and `R1(λ)` with `q1` outputs are the solution of an 
AFDP, while the filters `Q2(λ)` and `R2(λ)` with q2 outputs are the solution of an  
exact fault detection problem formulated for a reduced system obtained  
by decoupling the control and disturbance inputs from the residuals (see [4]).  
The overall resulting filters `Q.sys` and `R.sys` have observable state-space 
realizations `(AQ,BQ,CQ,DQ)` and `(AQ,BR,CQ,DR)`, respectively, and thus  
share the observable pairs `(AQ,CQ)`. 

Various user options can be specified via keyword arguments as follows:

If `minimal = true` (default), a least order filter synthesis is performed, while 
with `minimal = false` a full order synthesis is performed.  

If `exact = true`, an exact synthesis (without optimization) is performed, while 
with `exact = false` (default), an approximate synthesis is performed.  

If `HDesign = H1`, a design matrix `H1` of full row rank `q1` is used to build `q1` 
linear combinations of the left nullspace basis vectors of 
`G1(λ) := [ Gu(λ) Gd(λ); I 0]`. `H1` is used in the synthesis of the components 
`Q1(λ)` and `R1(λ)` in `(1)` (default: `HDesign = missing`).

If `HDesign2 = H2`, a design matrix `H2` of full row rank `q2` is used to build `q2` 
linear combinations of the left nullspace basis vectors of 
`G2(λ) := [ Gu(λ) Gd(λ) Gw(λ); I 0 0]`. `H2` is used in the synthesis of the components 
`Q2(λ)` and `R2(λ)` in `(1)` (default: `HDesign2 = missing`)

`rdim = q` specifies the desired number `q` of residual outputs for `Q` and `R`. 
If `rdim = missing`, the default value of `q` is chosen as `q = q1 + q2`, where 
the default values of `q1` and `q2` are chosen taking into account the rank `rw` 
of `Rw(λ)` in the reduced system (see [4]), as follows: 
if `HDesign = missing`, then  
`q1 = min(1,rw)`, if `minimal = true`, or `q1 = rw`, if `minimal = false`; 
if `HDesign = H`, then `q1` is the row dimension of the design matrix `H2`.
if `HDesign2 = missing`, then  
`q2 = 1-min(1,rw)`, if `minimal = true`, or `q2 = nvec-rw`, if `minimal = false`,
where `nvec` is the number of the nullspace basis 
vectors used for the initial synthesis (see [1]); 
if `HDesign2 = H2`, then `q2` is the row dimension of the design matrix `H2`.

`FDfreq = freq` specifies a vector of real frequency values or a scalar real frequency value
for strong detectability checks (default: `FDfreq = missing`).

If `nullspace = true` (default), a minimal proper nullspace basis is used for 
the initial synthesis of the fault detection filter. 
If `nullspace = false`, a full-order observer based nullspace basis is used. 
This option can be only used for a proper system without disturbance inputs. 

If `simple = true`, a simple proper nullspace basis is emplyed for synthesis. 
The orders of the basis vectors employed for the synthesis
are provided in `info.deg` and `info.deg2` (see below). 
If `simple = false` (default), then no simple basis is computed. 

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

`scale2 = σ2` specifies the scaling factor `σ2` to be employed for the components 
`Q2(λ)` and `R2(λ)` in `(1)`, i.e., 
use  `σ2*Q2(λ)` and `σ2*R2(λ)` instead of `Q2(λ)` and `R2(λ)`. 
(default: `σ2` is chosen to ensure the minimum gap provided by `Q1(λ)`) 

`freq = val` specifies the values of a test frequency to be employed to check the 
full row rank admissibility condition (default: randomly generated in the interval `(0,1)`).

`tcond = tcmax` specifies the maximum alowed condition number `tcmax` of 
the employed non-orthogonal transformations (default: `tcmax = 1.e4`).

`FDtol = tol1` specifies the threshold `tol1` for fault detectability checks
   (default: `tol1 = 0.0001`).

`FDGainTol = tol2` specifies the threshold `tol2` for strong fault detectability checks
   (default: `tol2 = 0.01`). 

`gamma = γ` specifies the allowed upper bound for `∥Rw(λ)∥∞` (default: `γ = 1`).

`epsreg = ϵ` specifies the value of the regularization parameter `ϵ` (default: `ϵ = 0.1`)

`sdegzer = δ` specifies the prescribed stability degree `δ` for zeros shifting
   (default: `δ = −0.05` for a continuous-time system `sysf.sys` and 
   `δ = 0.95` for a discrete-time system `sysf.sys`).

`nonstd = job` specifies the option to handle nonstandard optimization problems, as follows:

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

The resulting named tuple `info` contains `(tcond, degs, degs2, S, S2, HDesign, HDesign2, freq, gap)`, where:

`info.tcond` is the maximum of the condition numbers of the employed 
   non-orthogonal transformation matrices; a warning is issued if `info.tcond >= tcmax`;

`info.degs` is an integer vector containing the increasingly ordered degrees of a left minimal   
polynomial nullspace basis of `G1(λ) := [ Gu(λ) Gd(λ); I 0]` (also the left Kronecker indices of `G1(λ)`), if the 
state-space realization of `[Gu(λ) Gd(λ)]` is minimal. 
This information has been used in the case 
`minimal = true` to determine the least order of components `Q1(λ)` and `R1(λ)` in `(1)`.

`info.degs2` is an integer vector containing the increasingly ordered degrees of a left minimal   
polynomial nullspace basis of `G2(λ) := [ Gu(λ) Gd(λ) Gw(λ); I 0 0]` (also the left Kronecker indices of `G2(λ)`), if the 
state-space realization of `[Gu(λ) Gd(λ) Gw(λ)]` is minimal. 
This information has been used in the case 
`minimal = true` to determine the least order of components `Q2(λ)` and `R2(λ)` in `(1)`.

`info.S` is the binary structure matrix of the reduced system 
corresponding to the computed left nullspace basis of `G1(λ) := [ Gu(λ) Gd(λ); I 0]`;

`info.S2` is the binary structure matrix of the reduced system 
corresponding to the computed left nullspace basis of `G2(λ) := [ Gu(λ) Gd(λ)  Gw(λ); I 0 0]`;

`info.HDesign` is the design matrix `H1` employed for the synthesis of 
   the components `Q1(λ)` and `R1(λ)` in `(1)` of the fault detection filter;

`info.HDesign2` is the design matrix `H2` employed for the synthesis of 
the components `Q2(λ)` and `R2(λ)` in `(1)` of the fault detection filter;

`info.freq` is the frequency value employed to check the full 
row rank admissibility condition;

`info.gap` is the achieved gap `∥Rf(λ)∥∞−/∥Rw(λ)∥∞`, where the H−minus index is computed
over the whole frequency range, if `FDfreq = missing`, or over
the frequency values contained in `freq` if `FDfreq = freq`.

_Method:_ An extension of the Procedure AFD from [1] is implemented to solve 
the approximate fault detection problem (see also [2] and Remark 5.10 of [1]). 
The employed regularization approach, based on the modified co-outer-co-inner 
factorization, is discussed in [3], see also Remark 5.8 of [1]. 
For the details of the implemented method, see the documentation of the _afdsyn_ function in [4]. 

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017; sec. 5.3.

[2] A. Varga, General computational approach for optimal fault detection. 
              Proc. IFAC Symposium SAFEPROCESS, Barcelona, Spain, pp. 107–112, 2009.

[3] K. Glover and A. Varga, On solving non-standard H-/H_2/inf fault detection problems. 
              Proc. IEEE CDC, Orlando, FL, USA, pp. 891–896, 2011.

[4] A. Varga, Fault Detection and Isolation Tools (FDITOOLS) User's Guide, 
              [arXiv:1703.08480](https://arxiv.org/pdf/1703.08480). 


"""
function afdsyn(sysf::FDIModel{T}; rdim::Union{Int,Missing} = missing, poles::Union{AbstractVector,Missing} = missing, 
   sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing, 
   minimal::Bool = true, nullspace::Bool = true, simple::Bool = false, 
   FDtol::Real = 0.0001, FDGainTol::Real = 0.01, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
   tcond::Real = 1.e4, HDesign::Union{AbstractMatrix,Missing} = missing, HDesign2::Union{AbstractMatrix,Missing} = missing,
   offset::Real = sqrt(eps(float(real(T)))), atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
   rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T)))))*iszero(max(atol1,atol2)), 
   fast::Bool = true, exact::Bool = false, gamma::Real = 1, epsreg::Real = 0.1, 
   sdegzer::Union{Real,Missing} = missing, nonstd::Int = 1, freq::Real = rand(), scale2::Union{Real,Missing} = missing) where T

   Ts = sysf.sys.Ts                  
   disc = (Ts != 0);  # system type (continuous- or discrete-time)
   
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
  
   # imposed design matrix H1 to form linear combinations of basis vectors
   emptyHD1 = ismissing(HDesign)
   if !emptyHD1
      rdim1 = size(HDesign,1)
      !ismissing(rdim) && rdim1 > rdim && error("row dimension of HDesign must not exceed rdim")
      rdim1 == rank(HDesign) || error("HDesign must have full row rank")
   end

   # imposed design matrix H2 to form linear combinations of basis vectors
   emptyHD2 = ismissing(HDesign2)
   if !emptyHD2
      rdim2 = size(HDesign2,1)
      !ismissing(rdim) && rdim2 > rdim && error("row dimension of HDesign2 must not exceed rdim")
      rdim2 == rank(HDesign2) || error("HDesign2 must have full row rank")
   end

   if !ismissing(rdim) && !emptyHD1 && !emptyHD2 && rdim1+rdim2 != rdim
      error("the sum of row dimensions of HDesign and HDesign2 must be equal to rdim")
   end

   !ismissing(scale2) && iszero(scale2) && error("scale2 must be nonzero") 
   
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
      syse = [sysf.sys; eye(mu,m)]
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
   nvec = size(QR,1);   # number of basis vectors
   # check solvability conditions
   nvec == 0 && error("afdsyn: empty nullspace basis - the AFDP is not solvable")

   indf = (p+mu) .+ Vector(1:mf)             # input indices of Rf in QR
   indw = (p+mu+mf) .+ Vector(1:mw)          # input indices of Rw in QR
   indaux = (p+mu+mf+mw) .+ Vector(1:maux)   # input indices of Ra in QR

   # compute rank of Rw
   rw = gnrank(QR[:,indw]; atol1, atol2, rtol)

   nq = order(QR)             # order of the minimal basis
   
   # set H1 for checking the solvability condition
   if emptyHD1
      Htemp1 = eye(nvec)
   else
      degs = Int[];
      nh = size(HDesign,2)
      if nh < nvec
         # pad with zeros: row rank is preserved
         Htemp1 = [ HDesign zeros(T, rdim1, nvec-nh) ]
      else
         # remove trailing columns if necessary: rank may drop
         Htemp1 = HDesign[:,1:nvec];
         if nh > nvec && rdim1 != rank(Htemp1)
            error("The leading $rdim1 × $nvec part of HDesign must have full row rank")
         end
      end
      min(rw, rdim1) == gnrank(Htemp1*QR[:,indw]; atol1, atol2, rtol) ||
         error("afdsyn: the admissibility condition for HDesign is not fulfilled")
   end

   # set H2 for checking the solvability condition
   nvec2 = nvec-rw
   if emptyHD2
      Htemp2 = eye(nvec2)
   else
      degs2 = Int[]
      nh = size(HDesign2,2)
      if nh < nvec2
         # pad with zeros: row rank is preserved
         Htemp2 = [ HDesign2 zeros(T, rdim2, nvec2-nh) ]
      else
         # remove trailing columns if necessary: rank may drop
         Htemp2 = HDesign2[:,1:nvec2]
         if nh > nvec2 && rdim2 != rank(Htemp2)
            error("The leading $rdim2 × $nvec2 part of HDesign2 must have full row rank")
         end
      end
   end

   # adjust rdim, rdim1 and rdim2
   !ismissing(rdim) && (rdim = min(rdim,nvec))
   !ismissing(HDesign) && (rdim1 = min(rdim1,rw))
   !ismissing(HDesign2) && (rdim2 = min(rdim2,nvec2))
   
   # setup the number of outputs of filters Q1 and Q2
   if minimal
      #  least order design
      if ismissing(rdim)
         # set default output dimensions 
         emptyHD1 && (rdim1 = rw > 0 ? 1 : 0 )  
         emptyHD2 && (rdim2 = rw > 0 ? 0 : 1 )  
      else
         # set output dimensions for a given rdim 
         if emptyHD1 && emptyHD2
            if rdim <= rw
               rdim1 = rdim; rdim2 = 0
            else
               rdim1 = rw; rdim2 = rdim-rw
            end
         elseif emptyHD1 
            rdim1 = rdim-rdim2
         else
            rdim2 = rdim-rdim1
         end
      end
   else
      #  full order design
      if ismissing(rdim) 
         if emptyHD1 && emptyHD2
            rdim1 = rw 
            rdim2 = nvec-rw
         elseif emptyHD1  
            rdim1 = rw 
         else
            rdim2 = nvec-rw
         end
         # emptyHD1 && (rdim1 = rw )  
         # emptyHD2 && (rdim2 = nvec-rw) 
      else
         # set output dimensions for a given rdim 
         if emptyHD1 && emptyHD2
            if rdim <= rw
               rdim1 = rdim; rdim2 = 0
            else
               rdim1 = rw; rdim2 = rdim-rw
            end
         elseif emptyHD1 
            rdim1 = rdim-rdim2
         else
            rdim2 = rdim-rdim1
         end
      end
   end

   if rdim1 > 0
      # determine the structure matrix S1 underlying the synthesis of Q1
      if strongFD 
         S1 = fdisspec_(Htemp1*QR[:,indf], FDfreq; stabilize = true, atol1, atol2, atol3, rtol, FDGainTol)[1]
      else
         S1 = fditspec_(Htemp1*QR[:,indf]; atol1, atol2, rtol, FDtol)
      end
   else
      if strongFD 
         S1 = falses(0,mf,lfreq)
      else
         S1 = falses(0,mf)
      end
   end
   
   if rdim2 > 0
      # the case rw < nvec
      # compute a left nullspace basis Qt such that Qt*Rw = 0 and
      # obtain QRt = [ Qt Rt ], where Rt = [ Q Rft Rwt Rat] = Qt*[Q1 Rf Rw Ra]
      m2 = p+mu+mf+mw+maux
      QRt, infot = glnull(QR[:,[indw; Vector(1:p+mu); indf; indw; indaux]], m2; simple, atol1, atol2, rtol, fast, sdeg = sdegNS, offset) 
      QR2 = QRt[:,nvec .+ [Vector(1:p+mu); indf; indw; indaux]]
      nvec2 == size(QR2,1) || error("something wrong: try to adapt the rank decision tolerance")
      degs2 = infot.degs
      # determine the structure matrix S2 underlying the synthesis of Q2
      if strongFD 
         S2 = fdisspec_(Htemp2*QR2[:,indf], FDfreq; atol1, atol2, atol3, rtol, FDGainTol)[1]
      else
         S2 = fditspec_(Htemp2*QR2[:,indf]; atol1, atol2, rtol, FDtol)
      end
   else
      if strongFD 
         S2 = falses(0,mf,lfreq)
      else
         S2 = falses(0,mf)
      end
      degs2 = Int[]
   end

   if mf == 0
      # handle the case of no faults as a normal case
      S1 = falses(nvec,0)
      S2 = falses(nvec2,0)
      foff = 0   # set offset value
   else
      S = [S1;S2]
      if strongFD 
         # check strong detectability conditions 
         for ii = 1:lfreq
            all(maximum(S[:,:,ii],dims=1)) || error("afdsyn: strong detection of all faults not feasible")
         end
      else
         # check weak detectability conditions 
         all(maximum(S,dims=1)) || error("afdsyn: detection of all faults not feasible")
      end
      foff = indf[1]-1
   end
   if rdim1 > 0
      # synthesis of Q1 and R1
      #QR1 = QR;
      if strongFD 
         Sc = maximum(S1[:,:,1],dims=1)
         for ii = 2:lfreq
             Sc = maximum([Sc;S1[:,:,ii]],dims=1)
         end
         indf1 = Vector(1:mf)[Sc[:]]
      else
         indf1 = Vector(1:mf)[maximum(S1,dims=1)[:]]
      end
      Sc = strongFD ? S1[:,indf1,:] : S1[:,indf1]
      # build sysc = [Rf1 Rw [Q1 Rf Rw Ra]] 
      sysc = fdimodset(QR, f = foff .+ indf1, n = indw, aux = [Vector(1:p+mu); indf; indw; indaux])
      # R1 = Q2*[Rf1 Rw [Q1 Rf Rw Ra]]
      _, R1, info1 = afdredsyn(sysc; rdim = rdim1, simple, minimal, S = Sc, degs,   
                               gamma, epsreg, sdegzer, nonstd, freq,
                               sdeg, smarg, poles, HDesign, FDtol, FDGainTol, FDfreq, 
                               tcond, offset, exact, atol1, atol2, atol3, rtol, fast)
      tcond1 = info1.tcond
      HDesign1 = info1.HDesign
      gap1 = info1.gap
      QR1 = R1.sys[:,R1.aux]
      freq = info1.freq
   else
      QR1 = dss(zeros(0,size(QR,2)); Ts)
      tcond1 = 1
      gap1 = Inf
      HDesign1 = zeros(0,nvec)
   end
   
   if rdim2 > 0
      # synthesis of Q2 and R2
      if strongFD 
         Sc = maximum(S2[:,:,1],dims=1)
         for ii = 2:lfreq
             Sc = maximum([Sc;S2[:,:,ii]],dims=1)
         end
         indf2 = Vector(1:mf)[Sc[:]]
      else
         indf2 = Vector(1:mf)[maximum(S2,dims=1)[:]]
      end
      Sc = strongFD ? S2[:,indf2,:] : S2[:,indf2]
      # build sysc = [Rf2 [Q1 Rf Rw Ra]] 
      sysc = fdimodset(QR2, f = foff .+ indf2, aux = [Vector(1:p+mu); indf; indw; indaux])
      # R1 = Q2*[Rf2 [Q1 Rf Rw Ra]]
      _, R2, info2 = afdredsyn(sysc; rdim = rdim2, simple, minimal, S = Sc, degs = degs2,   
                               gamma, epsreg, sdegzer, nonstd, freq,
                               sdeg, smarg, poles, HDesign = HDesign2, FDtol, FDGainTol, FDfreq, 
                               tcond, offset, atol, atol1, atol2, atol3, rtol, fast)
      tcond2 = info2.tcond
      HDesign2 = info2.HDesign
      gap2 = info2.gap
      QR2 = R2.sys[:,R2.aux]
      # scale QR2 to ensure gap1 as minimum gap 
      if rdim1 > 0 
         ismissing(scale2) && (scale2 = gap1/ghinfnorm(QR2[:,indf])[1])
         QR2 = scale2*QR2
      end        
   else
      QR2 = dss(zeros(0,size(QR,2)); Ts)
      tcond2 = 1
      gap2 = Inf
      HDesign2 = zeros(0,nvec)
   end
   QR = [QR1;QR2] 
    
   Q = FDFilter(QR, p, mu)
   R = FDFilterIF(QR,0,0,mf,mw,maux; moff = p+mu)
   info = (tcond = max(tcond1,tcond2), degs = degs, degs2 = degs2, S = S1, S2 = S2, 
           HDesign = HDesign1, HDesign2 = HDesign2, freq = freq, 
           gap = fdif2ngap(R,FDfreq)[1])

   return Q, R, info
   
   # end AFDSYN
end
"""
    afdsyn(sysf::FDIModel, SFDI; rdim, nullspace = true, simple = false, minimal = true, exact = false, 
                           gamma = 1, epsreg = 0.1, sdegzer, nonstd = 1, freq, sdeg, smarg, poles, 
                           HDesign, HDesign2, scale2, FDtol, FDGainTol, FDfreq, 
                           tcond, offset, atol, atol1, atol2, atol3, rtol, fast = true)  
                           -> (Q::FDFilter, R::FDFilterIF, info)

Solve the approximate fault detection and isolation problem (AFDIP) for a given synthesis model
`sysf` with additive faults and a given binary structure vector `SFDI`. 
The computed stable and proper filter objects `Q` and `R` contain the 
fault detection filter, representing the solution of the AFDIP, and its internal form, respectively, 
and are determined such that the transfer function matrix of 
`R.sys[:,faults]` has its `j`-th column nonzero if `SFDI[j] = 1`. 
If the solution of a strong AFDIP is feasible, then the `j`-th column is zero if `SFDI[j] = 0`. 
If only a the solution of a  weak AFDIP is feasible, then the `j`-th column may be nonzero if `SFDI[j] = 0`. 
For the description of the keyword parameters see the function [`afdsyn`](@ref). 
"""
function afdsyn(sysf::FDIModel{T}, SFDI::Union{BitVector,AbstractVector{Bool}}; kwargs...) where T
   mu = length(sysf.controls)  
   md = length(sysf.disturbances) 
   mf = length(sysf.faults)
   mw = length(sysf.noise) 
   maux = length(sysf.aux)  
   mf == length(SFDI) || error("number of faults must be equal to the length dimension of SFDI")
   indd = Vector(1:mf)[SFDI .== false] 
   indf = Vector(1:mf)[SFDI .== true] 
   sysc = fdimodset(sysf.sys, c = 1:mu, d = [Vector(mu .+ (1:md)); (mu+md) .+ indd], 
                    f = (mu+md) .+ indf, n = (mu+md+mf) .+ (1:mw), aux = (mu+md) .+ (1:mf+mw+maux));
   #
   Q, Rft, info = try 
      afdsyn(sysc; kwargs...)
   catch err
      findfirst("empty nullspace basis",string(err)) === nothing &&   
      findfirst("detection of all faults not feasible",string(err)) === nothing && 
      findfirst("the admissibility condition",string(err)) === nothing &&  error("$err")  
      @warn "afdsyn: solution of strong AFDIP failed: trying to solve a weak AFDIP"           
      sysc = fdimodset(sysf.sys, c = 1:mu, d = mu .+ (1:md), 
                         f = (mu+md) .+ indf, n = [(mu+md) .+ indd; (mu+md+mf) .+ (1:mw)], 
                         aux = (mu+md) .+ (1:mf+mw+maux));
      Q, Rft, info = afdsyn(sysc; kwargs...)
      return Q, FDFilterIF(Rft.sys[:,Rft.aux],0,0,mf,mw,maux), info
   end
   return Q, FDFilterIF(Rft.sys[:,Rft.aux],0,0,mf,mw,maux), info
end
"""
    afdisyn(sysf::FDIModel, SFDI; rdim, nullspace = true, simple = false, minimal = true, separate = false,
                           gamma = 1, epsreg = 0.1, sdegzer, nonstd = 1, freq, sdeg, smarg, poles, 
                           HDesign, HDesign2, scale2, FDtol, FDGainTol, FDfreq, 
                           tcond, offset, atol, atol1, atol2, atol3, rtol, fast = true) 
                           -> (Q::FDFilter, R::FDFilterIF, info)

Solve the approximate fault detection and isolation problem (AFDIP) for a given synthesis model
`sysf` with additive faults and a given binary structure matrix `SFDI` with `nb` rows (specifications). 
The computed stable and proper filter objects `Q` and `R` contain the 
fault detection and isolation filter, representing the solution of the AFDIP, and its internal form, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.degs`, `info.degs2`, `info.HDesign`, 
`info.HDesign2` and`info.gap` contains additional synthesis related information (see below). 

The continuous- or discrete-time system `sysf.sys` is in a standard
or descriptor state-space form `sysf.sys = (A-λE,B,C,D)`, which corresponds to the input-output form  

       y = Gu(λ)*u + Gd(λ)*d + Gf(λ)*f + Gw(λ)*w + Ga(λ)*aux,

with the Laplace- or Z-transformed plant outputs `y`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Gu(λ)`, `Gd(λ)`, `Gf(λ)`, `Gw(λ)`, and `Ga(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysf.controls`, `sysf.disturbances`, `sysf.faults`, `sysf.noise` and `sysf.aux`, respectively.

The fault detection and isolation filter object `Q`, contains in `Q.sys` the resulting bank of `nb` filters. 
The `i`-th filter `Q.sys[i]` is in a standard state-space form and generates `r_i`, the `i`-th component (scalar or vector) 
of the overall residual vector `r := [r_1; r_2; ...; r_nb]`. 
The corresponding input-output (implementation) form of the `i`-th filter is

            r_i = Qyi(λ)*y + Qui(λ)*u   ,            

where `Qyi(λ)` and `Qui(λ)` are the transfer function matrices from the output and control inputs to the `i`-th residual component. 
The indices of output and control inputs are contained in the integer vectors 
`Q.outputs` and `Q.controls`, respectively.

The fault detection and isolation filter in internal form object `R`, contains `R.sys`, the resulting bank of `nb` 
internal form of the filters.  
The `i`-th filter `R.sys[i]` is in a standard state-space form, which generates the residual signal `r_i`, and corresponds to the 
input-output form

       r_i = Rui(λ)*u + Rdi(λ)*d + Rfi(λ)*f + Rwi(λ)*w + Rai(λ)*aux ,

where 

       | Rui(λ) Rdi(λ) Rfi(λ) Rwi(λ) Rai(λ) | := |Qyi(λ) Qui(λ)]*| Gu(λ) Gd(λ) Gf(λ) Gw(λ) Ga(λ) |. 
                                                                 |   I     0     0     0     0   |

The solution of the AFDIP ensures that for the `i`-th filter, `Rui(λ) = 0`, `Rdi(λ) = 0`, 
`Rfi(λ)` has its `j`-th column nonzero if the `(i,j)`-th element of `SFDI` is nonzero, 
and the H∞-norm of `Rwi(λ)` satisfies `||Rwi(λ)||∞ < γ`, where the bound `γ` is 
specified via the keyword argument `gamma`.
The indices of the inputs `u`, `d`, `f`, `w` and `aux` of the resulting filter `R.sys` are 
contained in the integer vectors `R.controls` (void), `R.disturbances` (void), 
`R.faults`, `R.noise` and `R.aux`, respectively.

The transfer function matrices `Qi(λ) := [ Qyi(λ) Qui(λ) ]` and `Ri(λ) := [ Rui(λ) Rdi(λ) Rfi(λ) Rwi(λ) Rai(λ) ]` 
of the `i`-th components of the resulting filters `Q.sys` and `R.sys`, respectively,
have, in general, the partitioned forms

     Qi(λ) = [ Q1i(λ) ] ,   Ri(λ) = [ R1i(λ) ] ,                      (1)
             [ Q2i(λ) ]             [ R2i(λ) ]

where the filters `Q1i(λ)` and `R1i(λ)` with `q1i` outputs are the solution of an 
AFDP, while the filters `Q2i(λ)` and `R2i(λ)` with q2i outputs are the solution of an  
exact fault detection problem formulated for a reduced system obtained  
by decoupling the control and disturbance inputs from the residuals (see [4]).  
The overall resulting component filters `Q.sys[i]` and `R.sys[i]` have observable state-space 
realizations `(AQi,BQi,CQi,DQi)` and `(AQi,BRi,CQi,DRi)`, respectively, and thus  
share the observable pairs `(AQi,CQi)`. 

Various user options can be specified via keyword arguments as follows:

If `minimal = true` (default), least order filter synthesis is performed to determine each of the component filters
`Q.sys[i]` and `R.sys[i]` for `i = 1, ...,nb`, while 
with `minimal = false` full order synthesis is performed.  

If `exact = true`, an exact synthesis (without optimization) is performed, while 
with `exact = false` (default), an approximate synthesis is performed.  

If `HDesign = H1`, then `H1` is an `nb`-dimensional array of full row rank or empty design matrices,
where `H1[i]` is the design matrix employed for the synthesis  
of the components `Q1i(λ)` and `R1i(λ)` in `(1)` of the `i`-th filter 
(default: `HDesign = missing`)

If `HDesign2 = H2`, then `H2` is an `nb`-dimensional array of full row rank or empty design matrices, 
where `H2[i]` is the design matrix employed for the synthesis 
of the components `Q2i(λ)` and `R2i(λ)` in `(1)` of the `i`-th filter  
(default: `HDesign2 = missing`)

`rdim = q` specifies the vector `q`, whose `i`-th component `q[i]` specifies 
the number of residual outputs for the `i`-th component filters `Q.sys[i]` and `R.sys[i]`. 
If `q` is a scalar, then a vector `rdim` with all components equal to `q` is assumed.
If `rdim = missing`, the default value of `q[i]` is chosen as `q[i] = q1i + q2i`, where 
the default values of `q1i` and `q2i` are chosen taking into account the rank `rwi` 
of `Rwi(λ)` in the reduced system (see [2]),  as follows: 
if `HDesign = missing`, then  
`q1i = min(1,rwi)`, if `minimal = true`, or `q1i = rwi`, if `minimal = false`; 
if `HDesign = H1`, then `q1i` is the row dimension of the nonemty design matrix `H1[i]`, or
if `H1[i]` is empty, the above choice for `HDesign = missing` is employed;
if `HDesign2 = missing`, then  
`q2i = 1-min(1,rwi)`, if `minimal = true`, or `q2i` is set to its maximum achievable value,
 if `minimal = false` (see [1]); 
if `HDesign2 = H2`, then `q2i` is the row dimension of the nonemty design matrix `H2[i]`, or 
if `H2[i]` is empty, the above choice for `HDesign2 = missing` is employed.

`FDfreq = freq` specifies a vector of real frequency values or a scalar real frequency value
for strong detectability checks (default: `FDfreq = missing`).

If `nullspace = true` (default), a minimal proper nullspace basis is used at the 
initial reduction step, if `separate = false`, 
or at all synthesis steps, if `separate = true`.
If `nullspace = false`, a full-order observer based nullspace basis is used at the 
initial reduction step, if `separate = false`, or at all synthesis steps, if `separate = true`.
This option can  only be used for a proper system without disturbance inputs. 

If `simple = true`, simple proper nullspace bases are emplyed for synthesis. 
If `simple = false` (default), then no simple bases are computed. 

If `separate = false` (default), a two-step synthesis procedure is employed, 
where a minimal proper nullspace basis is used at the 
initial reduction step. 
If `separate = true`, the filter components are separately determined by solving
appropriately formulated fault detection problems. 

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

`FDtol = tol1` specifies the threshold `tol1` for fault detectability checks
   (default: `tol1 = 0.0001`).

`FDGainTol = tol2` specifies the threshold `tol2` for strong fault detectability checks
   (default: `tol2 = 0.01`). 

`scale2 = σ2` specifies the vector of scaling factors `σ2` to be employed for the components 
`Q2i(λ)` and `R2i(λ)` in `(1)`, i.e., 
use  `σ2[i]*Q2i(λ)` and `σ2[i]*R2i(λ)` instead of `Q2i(λ)` and `R2i(λ)`. 
(default: For `scale2 = missing`, each `σ2[i]` is chosen to ensure the minimum gap provided by `Q1i(λ)`) 

`gamma = γ` specifies the allowed upper bound for the resulting `∥Rwi(λ)∥∞` (default: `γ = 1`).

`epsreg = ϵ` specifies the value of the regularization parameter `ϵ` used in
[`afdsyn`](@ref) (default: `ϵ = 0.1`)

`sdegzer = δ` specifies the prescribed stability degree `δ` for zeros shifting
   (default: `δ = −0.05` for a continuous-time system `sysf.sys` and 
   `δ = 0.95` for a discrete-time system `sysf.sys`).

`nonstd = job` specifies the option to handle nonstandard optimization problems
in [`afdsyn`](@ref), as follows:

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

The resulting named tuple `info` contains `(tcond, HDesign, HDesign2, freq, gap)`, where:

`info.tcond` is an `nb`-dimensional vector, whose `i`-th component is the maximum of the condition numbers of the employed 
non-orthogonal transformation matrices employed for the synthesis of the `i`-th filter component; 
a warning is issued if any `info.tcond[i] >= tcmax`;

`info.HDesign = H1` is an `nb`-dimensional vector of design matrices, 
whose `i`-th component `H1[i]` is the design matrix to be employed for the synthesis 
of the components `Q1i(λ)` and `R1i(λ)` in `(1)` of 
the `i`-th fault detection filter.

`info.HDesign2 = H2` is an `nb`-dimensional vector of design matrices, 
whose `i`-th component `H2[i]` is the design matrix to be employed for the synthesis 
of the components `Q2i(λ)` and `R2i(λ)` in `(1)` of  
the `i`-th fault detection filter.

`info.freq` is the frequency value employed to check the full 
row rank admissibility condition.

`info.gap` is an `nb`-dimensional vector, whose `i`-th component is the 
achieved gap for the synthesis of the `i`-th filter component.
   
_Method:_ The Procedure AFDI from [1] is implemented to solve 
the approximate fault detection and isolation problem. For implementation details, see [2].

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017; sec. 5.5.

[2] A. Varga, Fault Detection and Isolation Tools (FDITOOLS) User's Guide, 
              [arXiv:1703.08480](https://arxiv.org/pdf/1703.08480). 
"""
function afdisyn(sysf::FDIModel{T}, SFDI::Union{BitMatrix,BitVector,Array{Bool,2},Array{Bool,1}} = trues(1,length(sysf.faults)); 
                      rdim::Union{Vector{Int},Int,Missing} = missing, poles::Union{AbstractVector,Missing} = missing, 
                      sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing, 
                      nullspace::Bool = true, minimal::Bool = true, simple::Bool = false, separate::Bool = false, 
                      FDtol::Real = 0.0001, FDGainTol::Real = 0.01, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                      tcond::Real = 1.e4, HDesign::Union{Vector{Matrix{T1}},Missing} = missing, HDesign2::Union{Vector{Matrix{T2}},Missing} = missing,
                      offset::Real = sqrt(eps(float(real(T)))), atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                      rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T)))))*iszero(max(atol1,atol2)), 
                      fast::Bool = true, exact::Bool = false, gamma::Real = 1, epsreg::Real = 0.1, 
                      sdegzer::Union{Real,Missing} = missing, nonstd::Int = 1, freq::Real = rand(), scale2::Union{Vector{Real},Missing} = missing) where {T, T1 <: Real, T2 <: Real}
   Ts = sysf.sys.Ts                  
   disc = (Ts != 0);  # system type (continuous- or discrete-time)
   
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
   emptyHD1 = ismissing(HDesign)
   if !emptyHD1
      size(HDesign,1) == nb || error("number of HDesign components must be equal to the row dimension of SFDI")
      rdim1 = size.(HDesign,1)
      if !ismissing(rdim) 
         any(rdim1 .> rdim) && error("row dimensions of HDesign must not exceed rdim")
         for i = 1:nb
             if !isempty(HDesign[i])
                if rdim1[i] > 0
                   rdim1[i] == rank(HDesign[i]) || error("HDesign[$i] must have full row rank")
                end
             end
         end
      end
   end

   emptyHD2 = ismissing(HDesign2)
   if !emptyHD2
      size(HDesign2,1) == nb || error("number of HDesign2 components must be equal to the row dimension of SFDI")
      rdim2 = size.(HDesign2,1)
      if !ismissing(rdim) 
         any(rdim2 .> rdim) && error("row dimensions of HDesign2 must not exceed rdim")
         for i = 1:nb
             if !isempty(HDesign2[i])
                if rdim2[i] > 0
                   rdim2[i] == rank(HDesign2[i]) || error("HDesign2[$i] must have full row rank")
                end
             end
         end
      end
   end

   if !ismissing(rdim) && !emptyHD1 && !emptyHD2 && (rdim1 + rdim2) != rdim
      error("the sum of row dimensions of HDesign and HDesign2 must be equal to rdim")
   end

   if !ismissing(scale2) 
      any(iszero.(scale2)) && error("scale2 must be nonzero") 
      nb == length(scale2) || error("scale2 must be an $nb-dimensional vector with nonzero components")
   end

   tcond1 = similar(Vector{T},nb)
   degs1 = similar(Vector{Vector{Int}},nb)
   degs2 = similar(Vector{Vector{Int}},nb)
   HDesign1s = similar(Vector{Array{T,2}},nb)
   HDesign2s = similar(Vector{Array{T,2}},nb)
   gap = similar(Vector{T},nb)
   if separate
      Qt = similar(Vector{DescriptorStateSpace{T}},nb)
      Rt = similar(Vector{DescriptorStateSpace{T}},nb)
      for i = 1:nb
          # solve the corresponding AFDP           
          Qti, Rti, infoi = try 
                 afdsyn(sysf, view(SFDI,i,:); rdim = ismissing(rdim) ? missing : rdim[i], 
                        HDesign = ismissing(HDesign) ? missing : HDesign[i], 
                        HDesign2 = ismissing(HDesign2) ? missing : HDesign2[i], 
                        atol1, atol2, atol3, rtol, sdeg, smarg, poles, minimal,
                        FDtol, FDfreq, FDGainTol, simple, exact, tcond, offset, gamma, epsreg, 
                        sdegzer, nonstd, freq, scale2 = ismissing(scale2) ? missing : scale2[i]); 
          catch err
             findfirst("empty",string(err)) === nothing &&   
             findfirst("detection",string(err)) === nothing && 
             findfirst("condition",string(err)) === nothing &&  error("$err")  
             t1 = (mod(i,10) == 1 ? "$i-st" : "")
             t2 = (mod(i,10) == 2 ? "$i-nd" : "")
             t3 = (mod(i,10) == 3 ? "$i-rd" : "")
             t4 = (mod(i,10) > 3 ? "$i-th" : "")
             error("the $t1$t2$t3$t4  AFDIP is not solvable")
          end
          Qt[i] = Qti.sys
          Rt[i] = Rti.sys
          tcond1[i] = infoi.tcond
          degs1[i] = infoi.degs
          degs2[i] = infoi.degs2
          HDesign1s[i] = infoi.HDesign
          HDesign2s[i] = infoi.HDesign2
          gap[i] = infoi.gap
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
      size(QR,1) == 0 && error("empty nullspace basis: the AFDIP is not solvable")
   
      # Step 2): of Procedure EFDI 
      # initialize overall filters Q and R
      QRt = similar(Vector{typeof(QR)},nb)
      for i = 1:nb
          indd = Vector(1:mf)[SFDI[i,:] .== false] 
          indf = Vector(1:mf)[SFDI[i,:] .== true] 
          # pack [Rfd1 Rff2 [Q1 Rf1 Rw1 Raux1]]
          sysc = fdimodset(QR, d = (p+mu) .+ indd, f = (p+mu) .+ indf, n = (p+mu+mf) .+ (1:mw), aux = Vector(1:p+mu+mf+mw+maux))
          # determine [Q1i*Rff2 [Q1i*Q1 Q1i*Rf1 Q1i*Rw1 Q1i*Raux1]]
          _, QRauxi, infoi = try 
             afdsyn(sysc; rdim = ismissing(rdim) ? missing : rdim[i],  HDesign = ismissing(HDesign) ? missing : HDesign[i], 
             HDesign2 = ismissing(HDesign2) ? missing : HDesign2[i], atol1, atol2, atol3, rtol, sdeg, smarg, poles, minimal,
                                       FDtol, FDfreq, FDGainTol, simple, exact, tcond, offset, gamma, epsreg, 
                                       sdegzer, nonstd, freq, scale2 = ismissing(scale2) ? missing : scale2[i]); 
          catch err
             findfirst("empty nullspace basis",string(err)) === nothing &&   
             findfirst("detection of all faults not feasible",string(err)) === nothing && 
             findfirst("the admissibility condition",string(err)) === nothing &&  error("$err")  
             t1 = (mod(i,10) == 1 ? "$i-st" : "")
             t2 = (mod(i,10) == 2 ? "$i-nd" : "")
             t3 = (mod(i,10) == 3 ? "$i-rd" : "")
             t4 = (mod(i,10) > 3 ? "$i-th" : "")
             @warn "afdsyn: solution of strong AFDIP failed for the $t1$t2$t3$t4 reduced AFDP: trying to solve a weak AFDIP"           
             sysc = fdimodset(QR, f = (p+mu) .+ indf, n = [(p+mu) .+ indd; (p+mu+mf) .+ (1:mw)], aux = Vector(1:p+mu+mf+mw+maux))
             _, QRauxi, infoi = try 
               afdsyn(sysc; rdim = ismissing(rdim) ? missing : rdim[i],  HDesign = ismissing(HDesign) ? missing : HDesign[i], 
               HDesign2 = ismissing(HDesign2) ? missing : HDesign2[i], atol1, atol2, atol3, rtol, sdeg, smarg, poles, minimal,
                                         FDtol, FDfreq, FDGainTol, simple, exact, tcond, offset, gamma, epsreg, 
                                         sdegzer, nonstd, freq, scale2 = ismissing(scale2) ? missing : scale2[i]); 
            catch err
               (isnothing(findfirst("empty",string(err))) && isnothing(findfirst("detection",string(err))) && isnothing(findfirst("condition",string(err)))) ? rethrow() : error("the $t1$t2$t3 reduced AFDP is not solvable")
            end
          end
          # extract [Q1i*Q1 Q1i*Rf1 Q1i*Rw1 Q1i*Raux1 ]
          QRt[i] = QRauxi.sys[:,QRauxi.aux]
          tcond1[i] = max(tcond0,infoi.tcond)
          degs1[i] = infoi.degs
          degs2[i] = infoi.degs2
          HDesign1s[i] = infoi.HDesign
          HDesign2s[i] = infoi.HDesign2
          gap[i] = infoi.gap
       end
       Q = FDIFilter(QRt, p, mu)
       R = FDIFilterIF(QRt,0,0,mf,mw,maux; moff = p+mu)
   end   
   info = (tcond = tcond1, degs = degs1, degs2 = degs2, HDesign = HDesign1s, HDesign2 = HDesign2s, 
           freq = freq, gap = fdif2ngap(R, SFDI, FDfreq)[1])

   return Q, R, info

   # end AFDISYN
end
function emmbasesel(rgain::Matrix, degs::Vector{Int}, nout::Int, simple::Bool, atol::Real)
   #   emmbasesel(rgain, degs, nout, simple, atol) -> (seli, selord)
   #
   #   Select admissible basis vectors for solving the strong fault detection and isolation problem (strong EFDIP)
   #   using the `nvec × mf` full column rank frequency gain matrix `rgain`. 
   #   `seli` contains `nout`-touples (`mf ≤ nout`) of indices of basis vectors whose linear combination 
   #   is admissible, i.e. , the strong EFDIP is solvable by using fault detection filters with `mf` outputs. 
   #   If the associated `nvec` degrees contained in `degs` are provided, then
   #   `selord[i]` is the corresponding tentatively achievable least filter order.
   #   If `simple = true`, a simple basis is assumed, in which case, `degs[i]` is 
   #   also the order of the `i`-th basis vector. If `simple = false`, a minimum 
   #   rational basis is assumed. `selord` is empty if `degs` is empty. 
   #   `atol` is an aboslute tolerance for rank determinations. 


   #   Method: The selection approach is used in conjunction with the synthesis 
   #   Procedure EMMS described in [1]. 

   # References:
   # [1] Varga A.
   #     Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017.

   nvec, mf = size(rgain); 
   nd = length(degs)
   nodegs = (nd == 0)
   rdim = mf

   nodegs || length(degs) == nvec || error("the dimension of degs must be equal to the number of rows of rgain")
   
   (rdim >=1 && rdim <= nvec) || error("mf must have a positive value not exceeding $nvec")
   (nout >=rdim && nout <= nvec) || error("nout must have a value at least $rdim and at most $nvec")
   
   nvec == 1 && (return [1], nodegs ? Int[] : degs )
   

   # find rdim combinations of nout vectors which solve the AFDP 
   seli = collect(combinations(Vector(1:nvec),nout))
   ni = length(seli)
   selord = nodegs ? Int[] : fill(-1,ni) 
   nqmax = sum(degs)
   ii = trues(ni)
   for i = 1:ni
       indv = seli[i];
       # check admissibility
       if  rank(view(rgain,indv,:); atol) == mf
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
   # end EMMBASESEL
end
"""
    emmsyn(sysf::FDIModel, sysr::FDFilterIF; nullspace = true, simple = false, minimal = true, regmin = true, normalize = "gain", 
                           sdeg, smarg, poles, freq, HDesign, tcond, offset, 
                           atol, atol1, atol2, atol3, rtol, fast = true) 
                           -> (Q::FDFilter, R::FDFilterIF, info)

Solve the _exact model-matching problem_ (EMMP) for a given synthesis model `sysf::FDIModel` with additive faults 
and a given stable reference filter `sysr::FDFilterIF`. 
The computed stable and proper filter objects `Q` and `R` contain the 
fault detection filter, representing the solution of the EMMP, and its internal form, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.degs`, `info.M`, `info.freq` 
and `info.HDesign`, contains additional synthesis related information (see below). 

The continuous- or discrete-time system `sysf.sys` is in a standard
or descriptor state-space form `sysf.sys = (A-λE,B,C,D)`, which corresponds to the input-output form  

       y = Gu(λ)*u + Gd(λ)*d + Gf(λ)*f + Gw(λ)*w + Ga(λ)*aux,

with the Laplace- or Z-transformed plant outputs `y`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Gu(λ)`, `Gd(λ)`, `Gf(λ)`, `Gw(λ)`, and `Ga(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysf.controls`, `sysf.disturbances`, `sysf.faults`, `sysf.noise` and `sysf.aux`, respectively.

The continuous- or discrete-time reference filter `sysr.sys` is in a standard
or descriptor state-space form `sysr.sys = (Ar-λEr,Br,Cr,Dr)`, which corresponds to the input-output form  

       yr = Mru(λ)*u + Mrd(λ)*d + Mrf(λ)*f + Mrw(λ)*w + Mra(λ)*aux,

with the Laplace- or Z-transformed reference filter outputs `yr`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Mru(λ)`, `Mrd(λ)`, `Mrf(λ)`, `Mrw(λ)`, and `Mra(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysr.controls`, `sysr.disturbances`, `sysr.faults`, `sysr.noise` and `sysr.aux`, respectively.
If any of the above vectors is void, then the corresponding transfer function matrix is considered null. 

The fault detection filter object `Q`, contains in `Q.sys` the resulting filter 
in a standard state-space form, which generates the residual signal `r`. 
The corresponding input-output (implementation) form is

            r = Qy(λ)*y + Qu(λ)*u               

where `Qy(λ)` and `Qu(λ)` are the transfer function matrices from the output and control inputs to the residual. 
The indices of output and control inputs are contained in the integer vectors 
`Q.outputs` and `Q.controls`, respectively.

The fault detection filter internal form object `R`, contains `R.sys`, the resulting 
internal form of the filter 
in a standard state-space form, which generates the residual signal `r`, and corresponds to the 
input-output form

       r = Ru(λ)*u + Rd(λ)*d + Rf(λ)*f + Rw(λ)*w + Ra(λ)*aux ,

where 

       | Ru(λ) Rd(λ) Rf(λ) Rw(λ) Ra(λ) | = |Qy(λ) Qu(λ)|*| Gu(λ) Gd(λ) Gf(λ) Gw(λ) Ga(λ) |. 
                                                         |  I     0     0     0     0    |

The solution of the _standard_ EMMP is computed if `sysr.noise` and `sysr.aux` are void and
ensures that `Ru(λ) = M(λ)*Mru(λ)`, `Rd(λ) = M(λ)*Mrd(λ)` and `Rf(λ) = M(λ)*Mrf(λ)`, where `M(λ)` is 
the transfer function matrix of a stable, diagonal and invertible updating filter returned in `info.M`. 
This filter is determined to guarantee the stability of resulting filters `Q` and `R`.  
If `sysr.noise` and `sysr.aux` are not both void, then 
the _extended_ EMMP is solved which additionally ensures `Rw(λ) = M(λ)*Mrw(λ)` and `Ra(λ) = M(λ)*Mra(λ)`. 
The indices of the inputs `u`, `d`, `f`, `w` and `aux` of the resulting filter `R.sys` are 
contained in the integer vectors `R.controls`, `R.disturbances`, `R.faults`, `R.noise` and `R.aux`, respectively.

Various user options can be specified via keyword arguments as follows:

If `minimal = true` (default), a least order filter synthesis is performed, while 
with `minimal = false` no least order synthesis is performed.  

If `regmin = true` (default), the regularization (see [1]) is performed for the case when `sysr.controls` and
`sysr.disturbances` are void with the selection of 
a least order left annihilator `Nl(λ)` such that  `Nl(λ)*[Gu(λ) Gd(λ); I 0 ]`. 
If `regmin = false`, the regularization is performed by choosing 
`Nl(λ)` a minimal left nullspace basis of `G(λ) = [Gu(λ) Gd(λ); I 0 ]`.  

If `HDesign = H` is a full row rank design matrix, then `H*Nl(λ)` is used 
instead `Nl(λ)` (default: `HDesign = missing`).

An initial reduction step is performed using the nullspace-based approach (see [1]) 
if `sysr.controls`, `sysr.disturbances`, `sysr.noise` and `sysr.aux` are void and
`minimal = false`. In this case,  
if `nullspace = true` (default),
a minimal proper nullspace basis is used at the initial reduction step, while,
if `nullspace = false`, a full-order observer based nullspace basis is used at the 
initial reduction step.
This later option can  only be used for a proper system without disturbance inputs. 
The `nullspace` option is ignored if any of `sysr.controls`, `sysr.disturbances`, `sysr.noise` or
`sysr.aux` is non-void or if `minimal = true` 

If `simple = true`, a simple proper nullspace basis `Nl(λ)` 
is emplyed as left annihilator for synthesis. 
The orders of the basis vectors are provided in `info.deg`. 
If `simple = false` (default), then a minimal proper nullspace basis is computed. 

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

`freq = val` specifies the value of a test frequency to be employed to 
check the full column rank (i.e., left-invertibility) solvability condition 
(default: randomly generated in the interval `(0,1)`). 
The employed value of `freq` is returned in `info.freq`.

`normalize = job` specifies the option for the normalization  
of the diagonal elements of the updating matrix `M(λ)` as follows:

      job = "gain"    – scale with the gains of the zero-pole-gain representation (default);
      job = "dcgain"  – scale with the DC-gains;
      job = "infnorm" – scale with the values of infinity-norms. 

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

The resulting named tuple `info` contains `(tcond, degs, M, freq, HDesign) `, where:

`info.tcond` is the maximum of the condition numbers of the employed 
   non-orthogonal transformation matrices; a warning is issued if `info.tcond >= tcmax`;

`info.degs` is an integer vector containing the increasingly ordered degrees of a left minimal   
polynomial nullspace basis of `G(λ) := [ Gu(λ) Gd(λ); I 0]` (also the left Kronecker indices of `G(λ)`), if the 
state-space realization of `[Gu(λ) Gd(λ)]` is minimal;

`info.M` is the employed stable and invertible updating filter used to solve the EMMP, 
with a diagonal transfer function matrix `M(λ)`; 

`info.freq` is the employed frequency used to check left invertibility 
(set to `missing` if no frequency-based left invertibility check was performed)

`info.HDesign` is the design matrix `H` employed for the synthesis of 
   the fault detection filter `Q`; `H = missing` if no design matrix was involved.
   
_Method:_ The synthesis Procedures EMM and EMMS from [1] are implemented.
  Procedure EMM relies on the model-matching synthesis method proposed in 
  [2], while Procedure EMMS uses the inversion-based method proposed in [3]. 
  Procedure EMM is generally employed, unless a strong exact fault 
  detection and isolation problem (strong EFDIP) is solved, in
  which case Procedure EMMS is used. 

  The strong EFDIP corresponds to the choice of the reference filter `sysr` such that
  `Mru(λ) = 0`, `Mrd(λ) = 0`, `Mrf(λ)` is invertible, `Mrw(λ) = 0` and `Mra(λ) = 0`. 
  In this case, only the indices of fault inputs `sysr.faults` must be specified
  and the indices of the rest of inputs must be void. 
  The solution of a fault estimation problem can be targeted by
  choosing `Mrf(λ) = I` and checking that the resulting `info.M = I`. 

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017; sec. 5.6.

[2] A. Varga, New computational approach for the design of fault
      detection and isolation filters. 
      In M. Voicu (Ed.), "Advances in Automatic Control", vol. 754 of 
      The Kluwer International Series in Engineering and Computer Science, 
      Kluwer Academic Publishers, pp. 367-381, 2003.

[3] A. Varga. New computational paradigms in solving fault detection 
      and isolation problems. Annual Reviews in Control, 37:25–42, 2013. 
"""
function emmsyn(sysf::FDIModel{T1}, sysr::Union{FDFilterIF{T2},FDIModel{T2}}; poles::Union{AbstractVector,Missing} = missing, 
                      sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing, 
                      nullspace::Bool = true, minimal::Bool = true, simple::Bool = false, regmin::Bool = true, 
                      normalize::AbstractString = "gain", freq::Real = rand(), 
                      #FDtol::Real = 0.0001, FDGainTol::Real = 0.01, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                      tcond::Real = 1.e4, HDesign::Union{AbstractMatrix,Missing} = missing,
                      offset::Real = sqrt(eps(float(real(T1)))), atol::Real = zero(float(real(T1))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                      rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T1)))))*iszero(max(atol1,atol2)), 
                      fast::Bool = true) where {T1,T2}

   Ts = DescriptorSystems.promote_Ts(sysf.sys.Ts,sysr.sys.Ts)
   disc = (Ts != 0);  # system type (continuous- or discrete-time)

   rdim = size(sysr.sys,1)    
 
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
   
   # set default stability degree
   sdegdefault = disc ? 0.95 : -0.05
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
   mru = length(sysr.controls); 
   mru == 0 || mru == mu ||  error("Incompatible control input groups in sysf and sysr") 
   mrd = length(sysr.disturbances); 
   mrd == 0 || mrd == md ||  error("Incompatible disturbance input groups in sysf and sysr") 
   mrf = length(sysr.faults); 
   mrf == 0 || mrf == mf ||  error("Incompatible fault input groups in sysf and sysr") 
   mrw = length(sysr.noise); 
   mrw == 0 || mrw == mw ||  error("Incompatible noise input groups in sysf and sysr") 
   mra = length(sysr.aux); 
   mra > 0 || mra == maux ||  error("Incompatible auxiliary input groups in sysf and sysr") 
       
   m = mu+md+mf+mw+maux;       # total number of inputs
   p = size(sysf.sys,1);       # number of measurable outputs
    
   if mf == 0 && minimal
      @warn "Minimal synthesis option not feasible in the case of no faults"
      minimal = false
   end

   if mru+mrd+mrw+mra == 0 && !minimal
      # perform either Procedure EMM or EMMS if SYSR has no 'controls', 
      # 'disturbances', 'noise' and 'aux' input groups and no least order option is selected 
                   
      # form [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
      syse = [sysf.sys; eye(mu,m)];
      m2 = mf+mw+maux
      desc = (sysf.sys.E != I)
      #
      # compute a left nullspace basis Q = Q1 of G1 = [Gu Gd; I 0] = 0 and
      # obtain QR = [ Q1 R1 ], where R1 = [ Rf1 Rw1 Raux1] = Q1*[Gf Gw Ga;0 0 0]
      # QR, info1 = glnull(syse, m2; simple, atol1, atol2, rtol, fast) 
      if nullspace || md > 0 || (desc && rcond(sysf.sys.e) < 1.e-7 )
         # form [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
         #syse = [sysf(:,[inpu inpd inpf inpw inpaux]); eye(mu,m)];
         syse = [sysf.sys; eye(mu,m)];
         #
         # compute a left nullspace basis Q = Q1 of G1 = [Gu Gd; I 0] = 0 and
         # obtain QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = Q*[Gf Gw Gaux;0 0 0]
         QR, info1 = glnull(syse, m2; simple, atol1, atol2, rtol, fast, sdeg = sdegdefault, offset) 
         tcond1 = info1.tcond
         degs = info1.degs
      elseif mu == 0 && md == 0
         # compute minimal basis as Q = Q1 = I  and set
         # QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = [ Gf Gw Gaux ]
         QR = [eye(p) sysf.sys]
         tcond1 = 1.
         degs = Int[]
      else
         # compute minimal basis as Q = Q1 = [ I -Gu] and set
         # QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = [ Gf Gw Gaux ]
         QR = [ eye(p) dss(sysf.sys.A, sysf.sys.E, [-sysf.sys.B[:,inpu] sysf.sys.B[:,[inpf; inpw; inpaux]]],
                sysf.sys.C, [-sysf.sys.D[:,inpu] sysf.sys.D[:,[inpf; inpw; inpaux]]]; Ts)] 
         # perform stabilization if strong detectability has to be enforced
         tcond1 = 1.
         degs = Int[]
      end
   
      nvec = size(QR,1);   # number of basis vectors
      # check solvability conditions
      nvec == 0 && error("emmsyn: empty nullspace basis - the EMMP is not solvable")
      # degs = info1.degs
      # tcond1 = info1.tcond
      infodegs = degs

      indf = (p+mu) .+ Vector(1:mf)            # input indices of Rf1 in QR
      indfwa = (p+mu) .+ Vector(1:mf+mw+maux)  # input indices of [ Rf1 Rw1 Raux1] in QR
      Rfwa = QR[:,indfwa]
   
      Rftest = evalfr(QR[:,indf],freq; atol1, atol2, rtol)
   
      if mf+mw+maux > 0 
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
         # this is the usual case with nonempty [Rf1 Rw1 Raux1]
         if rdim == mf && nvec >= mf && 
            rank(evalfr(sysr.sys[:,sysr.faults],freq; atol1, atol2, rtol);atol) == mf &&
            rank(Htemp*Rftest;atol) == mf
            # perform Procedure EMMS if Gw = 0, Gaux = 0, Rf is invertible,
            # and the reduced Rf1 is left invertible
            # flip degrees of a minimal polynomial basis
            reverse!(degs)
            finish = (nvec == mf)  # set termination flag
            nq = order(QR)
        
            if !finish
               nout = mf       # initialize number of selected basis vectors
               # if !simple && regmin
               #    # permute states to speedup glmcover1 
               #    QR = xperm(QR,nq:-1:1);  
               # end
            else
               if emptyHD 
                  h = eye(nvec)
               end
            end
      
            while !finish     
                # choose nout basis vectors, which potentially lead to a least order
                # filter with rdim outputs:
                # basesel(i,:) contains the indices of candidate basis vectors;
                # ordsel(i)    contains the presumably achievable least orders
                basesel, ordsel = emmbasesel(Htemp*Rftest, degs, nout, simple, atol) 
                #
                # update the synthesis using the selections of candidate vector(s),
                # starting with the least (potentially) achievable order
                for i = 1:size(basesel,1)
                    baseind = basesel[i] # indices of current basis selection
                    if nout == mf
                       hbase = eye(mf);
                    else
                       hbase = rand(mf,nout); 
                    end
                    ip = [baseind; setdiff(1:nvec,baseind)][:]
                    if simple
                       # handle simple basis
                       # here only the case rdim = nout can happen
                       if regmin
                          if emptyHD 
                             # select vectors and elliminate unobservable dynamics  
                             noelim = falses(nq) 
                             ell = sum(degs[1:basesel[i][1]-1]); 
                             for jj = 1:nout 
                                 ellnext = sum(degs[1:baseind[jj]]);
                                 noelim[ell+1:ellnext] .= true;
                                 ell = ellnext
                             end
                             ir = noelim
                             Ar, Er, Br, Cr, Dr = dssdata(QR[baseind,:])
                             QRfwtest = dss(view(Ar,ir,ir), Er == I ? I : view(Er,ir,ir),view(Br,ir,:),view(Cr,:,ir), Dr; Ts)
                             h = Htemp[ip[1:rdim],:]
                          else
                             QRfwtest = gir(Htemp*QR[ip,:]; atol1, atol2, rtol)
                          end
                       else
                          if emptyHD 
                              h = Htemp[ip[1:rdim],:]
                              QRfwtest = gir(QR[baseind,:]; atol1, atol2, rtol, infinite = false)
                           else
                              QRfwtest = gir(Htemp*QR; atol1, atol2, rtol, infinite = false)
                           end
                       end
                    else
                       # handle minimal basis
                       # use the above output permutation vector ip for glmcover1  
                       if regmin
                          if rdim == nout
                             if emptyHD 
                                QRfwtest, _, info2 = glmcover1(QR[ip,:], rdim; atol1, atol2, rtol)
                                if !isempty(ordsel) && (order(QRfwtest) != ordsel[i])
                                   @warn "emmsyn: expected reduced order not achieved"
                                end
                                h = Htemp[ip[1:rdim],:]
                             else
                                QRfwtest, _, info2 = glmcover1([Htemp; eye(nvec)]*QR[ip,:], rdim; atol1, atol2, rtol)
                             end
                          else
                             # the case rdim < nout can only happen if no
                             # HDesign is explicitly provided
                             Htemp = blockdiag(hbase,eye(nvec-nout)) 
                             QRfwtest, _, info2 = glmcover1([Htemp; eye(nvec)]*QR[ip,:], rdim; atol1, atol2, rtol)
                          end
                       else
                          # here only the case rdim = nout can happen
                          if emptyHD
                             h = Htemp[ip[1:rdim],:]
                             QRfwtest = gir(QR[baseind,:]; atol1, atol2, rtol, infinite = false) 
                          else
                             QRfwtest = gir(Htemp*QR; atol1, atol2, rtol, infinite = false) 
                          end
                       end
                    end
                    # check invertibility of compressed Rf1; 
                    if !simple && regmin
                       # dismiss minimal design if the check fails
                       Rtest1 = evalfr(QRfwtest[:,indf],freq; atol1, atol2, rtol) 
                       if rank(Rtest1; atol) == mf
                          # adjust condition number of employed transformations
                          tcond1 = max(tcond1, info2.fnorm, info2.tcond)
                          tcond1 > tcond && 
                             @warn "emmsyn: possible loss of numerical stability due to ill-conditioned transformations"
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
                nout > nvec && !finish &&
                        error("something wrong: try perhaps with another test frequency")
            end
            emptyHD && (Htemp = h)
            # compute the irreducible realization of Qtilde = Rf*(inv(Rf1)*Q1)  
            # by first solving the linear rational matrix equation Rf*Q2 = Q1
            Q2 = grsol(QR[:,[indf; 1:p+mu]],p+mu; atol1, atol2, rtol)[1]
            Qtilde = gminreal(sysr.sys[:,sysr.faults]*Q2; atol1, atol2, rtol)  
         else     
            # perform Procedure EMM in the general case
            # 2) Solve Q2*Rfwa = [Rf Rw Raux] and form Q = Q2*Q1 .
            Htemp = missing
            freq = missing
            if mrw+mra == 0 
               Q2, info2, = glsol(Rfwa[:,1:mf], sysr.sys; atol1, atol2, rtol)
            else
               Q2, info2, = glsol(Rfwa,sysr.sys; atol1, atol2, rtol)
            end
            Qtilde = Q2*QR[:,1:p+mu]
            tcond1 = max(tcond1, info2.fnorm, info2.tcond)
         end
         # 3) compute diagonal M such that M*Q has a desired stability degree;  
         #    update Q <- M*Q
         # LCF is performed after removing unobservable eigenvalues
         Qt = dss(zeros(0,p+mu); Ts); M = dss(zeros(0,0); Ts)
         for i = 1:rdim
             Qi, Mi = glcf(gir(Qtilde[i,:]; atol1, atol2, rtol, contr = false); atol1, atol2, atol3, rtol, sdeg, smarg, evals = poles, 
                                                          mininf = true, mindeg = true)
             if isequal(normalize,"gain") 
               sc = zpk(dss2rm(Mi)[1,1])[3]
             elseif isequal(normalize,"dcgain")
               sc = dcgain(Mi)[1,1]
             else
               sc = ghinfnorm(Mi)[1]
             end
             Qt = [Qt; Qi/sc]; M = append(M,Mi/sc)
         end
      else
         # this is the case with [Rf1 Rw1 Ra1] empty; M must not be diagonal
         # 3) compute M such that M*Q has a desired stability degree;  
         #    update Q <- M*Q
         Qt, M = glcf(QR; atol1, atol2, atol3, rtol, sdeg, smarg, evals = poles, 
                         mininf = true, mindeg = true)
         Htemp = missing
      end
      Rt = M*sysr.sys     # R = M*SYSR    
      if mw+maux > 0 && mrw+mra == 0
         Rt = gir([Rt Qt*[sysf.sys[:,[inpw;inpaux]]; zeros(mu,mw+maux)]]; atol1, atol2, rtol) 
      end
   else
      # apply the two-step procedure to solve the EMMP if the least order
      # synthesis option has been selected or SYSR has either 
      # 'controls' or 'disturbances' input groups, or both
      if mrw+mra == 0
         # apply the two-step procedure for the case
         # form Ge = [ Gu Gd Gf; I 0 0 ] 
         m1 = m-mw-maux;
         syse = [sysf.sys[:,1:m1]; eye(mu,m1)]
         rinp = zeros(0,m1)
         mru == 0 || (rinp = [rinp; eye(mu,m1)])
         mrd == 0 || (rinp = [rinp; zeros(md,mu) eye(md,m1-mu)])
         mrf == 0 || (rinp = [rinp; zeros(mf,mu+md) eye(mf,m1-mu-md)])
         # form explicitly Rref = [ Mru_i Mrd_i Mrf_i ] 
         Rref = sysr.sys*rinp;
      else
         # apply the two-step procedure for the general case
         # form Ge = [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
         syse = [sysf.sys; eye(mu,m)];
         rinp = zeros(0,m);
         mru == 0 || (rinp = [rinp; eye(mu,m)])
         mrd == 0 || (rinp = [rinp; zeros(md,mu) eye(md,m-mu)])
         mrf == 0 || (rinp = [rinp; zeros(mf,mu+md) eye(mf,m-mu-md)])
         mrw == 0 || (rinp = [rinp; zeros(mw,mu+md+mf) eye(mf,m-mu-md-mf)])
         mra == 0 || (rinp = [rinp; zeros(maux,m-maux) eye(maux)])
         # form explicitly Rref = [Mru_i Mrd_i Mrf_i Mrw_i ] 
         Rref = sysr.sys*rinp;
      end
      Qtilde, info2,  = glsol(syse, Rref; atol1, atol2, rtol, mindeg = true)

      tcond1 = max(info2.fnorm,info2.tcond)
      infodegs = Int[]
      Htemp = missing
      freq = missing
      Qt = dss(zeros(0,p+mu); Ts); M = dss(zeros(0,0); Ts)
      for i = 1:rdim
          Qi, Mi = glcf(gminreal(Qtilde[i,:]; atol1, atol2, rtol); atol1, atol2, atol3, rtol, sdeg, smarg, evals = poles, 
                                                       mininf = true, mindeg = true)
          if isequal(normalize,"gain") 
             sc = zpk(dss2rm(Mi; atol1, atol2, rtol)[1,1])[3]
          elseif isequal(normalize,"dcgain")
             sc = dcgain(Mi)[1,1]
          else
             sc = ghinfnorm(Mi)[1]
          end
          Qt = [Qt; Qi/sc]; M = append(M,Mi/sc)
      end
      Rt = M*sysr.sys     # R = M*SYSR   
 
      if mw+maux > 0 && mrw+mra == 0
         Rt = gir([Rt Qt*[sysf[:,[inpw;inpaux]]; zeros(mu,mw+maux)]]; atol1, atol2, rtol) 
      end
   end
   # transform to standard state-space  
   Q = FDFilter(gss2ss(Qt; atol1, atol2, rtol)[1], p, mu)

   R = FDFilterIF(gss2ss(Rt; atol1, atol2, rtol)[1], mru, mrd, mrf, mw > 0 ? mw : mrw, maux > 0 ? maux : mra)
   info = (tcond = tcond1, degs = infodegs, M = M, freq = freq, HDesign = ismissing(Htemp) ? missing : convert(Matrix{Float64},Htemp))
   return Q, R, info
   
   # end EMMSYN
end
"""
    emmsyn(sysf::FDIModel, sysr::FDIFilterIF; nullspace = true, simple = false, minimal = true, regmin = true, normalize = "gain", 
                           sdeg, smarg, poles, freq, HDesign, tcond, offset, 
                           atol, atol1, atol2, atol3, rtol, fast = true) 
                           -> (Q::FDIFilter, R::FDIFilterIF, info)

Solve the _exact model-matching problem_ (EMMP) for a given synthesis model `sysf::FDIModel` with additive faults 
and a given bank of stable reference filters `sysr::FDIFilterIF`. 
The computed stable and proper filter objects `Q` and `R` contain the 
bank of fault detection filters, 
representing the component-wise solution of the EMMP, and 
their internal forms, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.degs`, `info.M`, `info.freq` 
and `info.HDesign`, contains additional synthesis related information (see below). 

The continuous- or discrete-time system `sysf.sys` is in a standard
or descriptor state-space form `sysf.sys = (A-λE,B,C,D)`, which corresponds to the input-output form  

       y = Gu(λ)*u + Gd(λ)*d + Gf(λ)*f + Gw(λ)*w + Ga(λ)*aux,

with the Laplace- or Z-transformed plant outputs `y`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Gu(λ)`, `Gd(λ)`, `Gf(λ)`, `Gw(λ)`, and `Ga(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysf.controls`, `sysf.disturbances`, `sysf.faults`, `sysf.noise` and `sysf.aux`, respectively.

The continuous- or discrete-time reference filters packed in `sysr` 
are in a standard
or descriptor state-space form, where the `i`-th filter 
`sysr.sys[i] = (Ari-λEri,Bri,Cri,Dri)` corresponds to the input-output form  

       yri = Mrui(λ)*u + Mrdi(λ)*d + Mrfi(λ)*f + Mrwi(λ)*w + Mrai(λ)*aux,

with the Laplace- or Z-transformed reference filter outputs `yri`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Mrui(λ)`, `Mrdi(λ)`, `Mrfi(λ)`, `Mrwi(λ)`, and `Mrai(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysr.controls`, `sysr.disturbances`, `sysr.faults`, `sysr.noise` and `sysr.aux`, respectively.
If any of the above vectors is void, then the corresponding transfer function matrices are considered null. 

The fault detection and isolation filter object `Q`, contains in its `i`-th 
component `Q.sys[i]` the resulting filter 
in a standard state-space form, which generates the `i`-th component  `ri` of the residual signal. 
The corresponding input-output (implementation) form is

            ri = Qyi(λ)*y + Qui(λ)*u               

where `Qyi(λ)` and `Qui(λ)` are the transfer function matrices from the output and control inputs to the residual. 
The indices of output and control inputs are contained in the integer vectors 
`Q.outputs` and `Q.controls`, respectively.

The fault detection and isolation filter internal form object `R`, 
contains in its `i`-th component `R.sys[i]`, the resulting 
internal form of the filter 
in a standard state-space form, 
which generates the `i`-th component `ri` of residual signal , and corresponds to the 
input-output form

       ri = Rui(λ)*u + Rdi(λ)*d + Rfi(λ)*f + Rwi(λ)*w + Rai(λ)*aux ,

where 

       | Rui(λ) Rdi(λ) Rfi(λ) Rwi(λ) Rai(λ) | = |Qyi(λ) Qui(λ)|*| Gu(λ) Gd(λ) Gf(λ) Gw(λ) Ga(λ) |. 
                                                                |  I     0     0     0     0    |

The component-wise solution of the _standard_ EMMP is computed if `sysr.noise` and `sysr.aux` are void and
ensures that `Rui(λ) = Mi(λ)*Mrui(λ)`, `Rdi(λ) = Mi(λ)*Mrdi(λ)` and `Rfi(λ) = Mi(λ)*Mrfi(λ)`, where `Mi(λ)` is 
the transfer function matrix of a stable, diagonal and invertible updating filter returned in the `i`-th component of
the vector `info.M`. 
This filter is determined to guarantee the stability of the `i`-th components of resulting filters `Q` and `R`.  
If `sysr.noise` and `sysr.aux` are not both void, then 
the _extended_ EMMP is component-wise solved which additionally ensures `Rwi(λ) = Mi(λ)*Mrwi(λ)` and `Rai(λ) = Mi(λ)*Mrai(λ)`. 
The indices of the inputs `u`, `d`, `f`, `w` and `aux` of the resulting filter `R.sys` are 
contained in the integer vectors `R.controls`, `R.disturbances`, `R.faults`, `R.noise` and `R.aux`, respectively.

Various user options can be specified via keyword arguments as follows:

If `minimal = true` (default), least order filter syntheses are performed, while 
with `minimal = false` no least order synthesis are performed.  

If `regmin = true` (default), the regularization (see [1]) is performed for the case when `sysr.controls` and
`sysr.disturbances` are void with the selection of a
least order left annihilator `Nl(λ)` such that  `Nl(λ)*[Gu(λ) Gd(λ); I 0 ]`. 
If `regmin = false`, the regularization is performed by choosing 
`Nl(λ)` a minimal left nullspace basis of `G(λ) = [Gu(λ) Gd(λ); I 0 ]`.  

If `HDesign = H` is a vector of full row rank design matrices, 
then `H[i]*Nl(λ)` is used 
instead `Nl(λ)` for the synthesis of the `i`-th filter (default: `HDesign = missing`).

An initial reduction step is performed using the nullspace-based approach (see [1]) 
if `sysr.controls`, `sysr.disturbances`, `sysr.noise` and `sysr.aux` are void and
`minimal = false`. In this case,  
if `nullspace = true` (default),
a minimal proper nullspace basis is used at the initial reduction step, while,
if `nullspace = false`, a full-order observer based nullspace basis is used at the 
initial reduction step.
This later option can  only be used for a proper system without disturbance inputs. 
The `nullspace` option is ignored if any of `sysr.controls`, `sysr.disturbances`, `sysr.noise` or
`sysr.aux` is non-void or if `minimal = true` 

If `simple = true`, a simple proper nullspace basis `Nl(λ)` 
is emplyed as left annihilator for synthesis. 
The orders of the basis vectors are provided in `info.deg`. 
If `simple = false` (default), then a minimal proper nullspace basis is computed. 

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

`freq = val` specifies the value of a test frequency to be employed to 
check the full column rank (i.e., left-invertibility) solvability condition 
(default: randomly generated in the interval `(0,1)`). 
The employed value of `freq` is returned in `info.freq`.

`normalize = job` specifies the option for the normalization  
of the diagonal elements of the updating matrices `Mi(λ)` as follows:

      job = "gain"    – scale with the gains of the zero-pole-gain representation (default);
      job = "dcgain"  – scale with the DC-gains;
      job = "infnorm" – scale with the values of infinity-norms. 

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

The resulting named tuple `info` contains `(tcond, degs, M, freq, HDesign) `, where:

`info.tcond` is the maximum of the condition numbers of the employed 
   non-orthogonal transformation matrices; a warning is issued if `info.tcond >= tcmax`;

`info.degs` is an integer vector containing the increasingly ordered degrees of a left minimal   
polynomial nullspace basis of `G(λ) := [ Gu(λ) Gd(λ); I 0]` (also the left Kronecker indices of `G(λ)`), if the 
state-space realization of `[Gu(λ) Gd(λ)]` is minimal;

`info.M` is a vector of descriptor systems, whose `i`-th system `info.M[i]` 
contains the employed stable and invertible updating filter used to solve 
the `i`-th EMMP, with a diagonal transfer function matrix `Mi(λ)`; 

`info.freq` is the employed frequency used to check left invertibility 
(set to `missing` if no frequency-based left invertibility check was performed)

`info.HDesign` is a vector of design matrices `H`, where `H[i]` is the design matrix 
employed for the synthesis of the `i`-th component of the fault detection filter `Q`; 
`H[i]` is an empty matrix if no design matrix was involved.
   
_Method:_ The synthesis Procedures EMM and EMMS from [1] are used 
to determine the component filters.
  Procedure EMM relies on the model-matching synthesis method proposed in 
  [2], while Procedure EMMS uses the inversion-based method proposed in [3]. 
  Procedure EMM is generally employed, unless a strong exact fault 
  detection and isolation problem (strong EFDIP) is solved, in
  which case Procedure EMMS is used. 

  The strong EFDIP corresponds to the choice of each component of the bank of 
  reference filters `sysr` such that
  `Mrui(λ) = 0`, `Mrdi(λ) = 0`, `Mrfi(λ)` is invertible, `Mrwi(λ) = 0` and `Mrai(λ) = 0`. 
  In this case, only the indices of fault inputs `sysr.faults` must be specified
  and the indices of the rest of inputs must be void. 
  The solution of a fault estimation problem can be targeted by
  choosing `Mrfi(λ) = ei`, where `ei` is the `i`-th row of the appropriate identity matrix,  and checking that the resulting `info.M = 1`. 

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017; sec. 5.6.

[2] A. Varga, New computational approach for the design of fault
      detection and isolation filters. 
      In M. Voicu (Ed.), "Advances in Automatic Control", vol. 754 of 
      The Kluwer International Series in Engineering and Computer Science, 
      Kluwer Academic Publishers, pp. 367-381, 2003.

[3] A. Varga. New computational paradigms in solving fault detection 
      and isolation problems. Annual Reviews in Control, 37:25–42, 2013. 
"""
function emmsyn(sysf::FDIModel{T1}, sysref::FDIFilterIF{T2}; 
                poles::Union{AbstractVector,Missing} = missing, 
                sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing, 
                nullspace::Bool = true, minimal::Bool = true, simple::Bool = false, regmin::Bool = true, 
                normalize::AbstractString = "gain", freq::Real = rand(), 
                tcond::Real = 1.e4, HDesign::Union{AbstractMatrix,Missing} = missing,
                offset::Real = sqrt(eps(float(real(T1)))), atol::Real = zero(float(real(T1))), 
                atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T1)))))*iszero(max(atol1,atol2)), 
                fast::Bool = true) where {T1,T2} 

   N = length(sysref.sys)
   DescriptorSystems.promote_Ts(sysf.sys.Ts, sysref.sys[1].Ts)  # check sampling time

   emptyHD = ismissing(HDesign)
   T = promote_type(T1,T2)
   
   HDesign = similar(Array{Matrix{T},1},N)
   Q = similar(sysref.sys,N) 
   R = similar(sysref.sys,N) 
   M = similar(sysref.sys,N) 
   tcond1 = 1
   infodegs = 0
   inpu, inpd, inpf, inpw, inpa = 0,0,0,0,0
   for i = 1:N
       Mri = FDFilterIF(sysref.sys[i], sysref.controls, sysref.disturbances, 
                        sysref.faults, sysref.noise, sysref.aux)
       Qi, Ri, infoi = emmsyn(sysf, Mri; sdeg, smarg, poles, nullspace, minimal, simple, regmin, normalize, freq, 
                                         tcond, HDesign = emptyHD ? missing : HDesign[i], offset, 
                                         atol1, atol2, atol3, rtol, fast)
       Q[i] = Qi.sys
       R[i] = Ri.sys
       M[i] = infoi.M
       HDesign[i] = ismissing(infoi.HDesign) ? zeros(T,0,0) : infoi.HDesign
       tcond1 = max(tcond1, infoi.tcond)
       i == 1 && (infodegs = infoi.degs; inpu = Ri.controls; inpd = Ri.disturbances; inpf = Ri.faults; inpw = Ri.noise; inpa = Ri.aux)
   end
   p = size(sysf.sys,1)
   mu = length(sysf.controls)
   info = (tcond = tcond1, degs = infodegs, M = M, freq = freq, HDesign = all(ismissing.(HDesign)) ? missing : HDesign)
   return FDIFilter(Q, p, mu), 
          FDIFilterIF(R, controls = inpu, disturbances = inpd, faults = inpf, noise = inpw, aux = inpa), info
end
function ammbasesel(rgain::Matrix, degs::Vector{Int}, nout::Int, simple::Bool, atol::Real, strongfdi::Bool )
   #   ammbasesel(rgain, degs, nout, simple, atol, strongfdi) -> (seli, selord)
   #
   #   Select admissible basis vectors for solving the approximate model-matching problem (AMMP)
   #   using the `nvec × mr` frequency gain matrix `rgain`. 
   #   `seli` contains `nout`-touples (`mr ≤ nout`) of indices of basis vectors whose linear combination 
   #   is admissible, i.e. , the strong fault detection and isolation problem is solvable if `strongfdi = true`,
   #   or the fault detection problem is solvable if `strongfdi = false`. 
   #   The number of selected vectors must satisfy `nout ≤ min(nvec,mr)`. 
   #   Each row seli[i,:] contains nout indices of basis vectors, such that rgain[seli[i,:],:] 
   #   has full row rank nout if `strongfdi = true`, or 
   #   rgain[seli[i,:],:] has all columns nonzero if `strongfdi = false`.
   #   If the associated `nvec` degrees contained in `degs` are provided, then
   #   `selord[i]` is the corresponding tentatively achievable least filter order.
   #   If `simple = true`, a simple basis is assumed, in which case, `degs[i]` is 
   #   also the order of the `i`-th basis vector. If `simple = false`, a minimum 
   #   rational basis is assumed. `selord` is empty if `degs` is empty. 
   #   `atol` is an aboslute tolerance for rank determinations. 


   #   Method: The selection approach is used in conjunction with the synthesis 
   #   Procedure AMMS described in [1]. 

   # References:
   # [1] Varga A.
   #     Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017.

   nvec, mr = size(rgain); 
   nd = length(degs)
   nodegs = (nd == 0)
   
   nodegs || length(degs) == nvec || error("the dimension of degs must be equal to the number of rows of rgain")
   nout <= min(mr,nvec) || error("nout must have a value at most $(min(mr,nvec))")
   
   nvec == 1 && (return [1], nodegs ? Int[] : degs )
   

   # find nout combinations of nvec vectors which solve the AMMP 
   seli = collect(combinations(Vector(1:nvec),nout))
   ni = length(seli)
   selord = nodegs ? Int[] : fill(-1,ni) 
   nqmax = sum(degs)
   ii = falses(ni)
   for i = 1:ni
       indv = seli[i];
       # check admissibility
       if strongfdi 
          if rank(view(rgain,indv,:); atol) == nout
             ii[i] = true
             if !nodegs
                # estimate orders 
                if simple 
                   # degree = the sums of degrees of selected vectors
                   selord[i] = sum(degs[indv])
                else
                   # degree = rdim times the maximum degree of selected vectors
                   selord[i] = min(nqmax,nout*maximum(degs[indv]))
                end
             end
          end
       else
          # evaluate minimum column norm
          beta = norm(view(rgain,indv,1)) 
          for j = 2:mr
              beta  = min(beta,norm(view(rgain,indv,j)))
          end  
          if (atol > 0 && beta > atol) || (atol == 0 && beta > nvec*mr*eps(norm(rgain)))
             ii[i] = true
             if !nodegs
                # estimate orders 
                if simple 
                   # degree = the sums of degrees of selected vectors
                   selord[i] = sum(degs[indv])
                else
                   # degree = rdim times the maximum degree of selected vectors
                   selord[i] = min(nqmax,nout*maximum(degs[indv]))
                end
             end
           end
       end
   end
  
   seli = seli[ii]

   if !nodegs 
      selord = selord[ii]
      # sort row combinations to ensure increasing tentative orders  
      ii = sortperm(selord)
      seli = seli[ii]
      selord = selord[ii]
   end
   return seli, selord      
   # end AMMBASESEL
end
"""
    ammsyn(sysf::FDIModel, sysr::FDFilterIF; nullspace = true, simple = false, mindeg = false, 
                           regmin = true, normalize = "infnorm", H2syn = false, reltol = 1.e-4, 
                           sdeg, smarg, poles, freq, HDesign, tcond, offset, 
                           atol, atol1, atol2, atol3, rtol, fast = true) 
                           -> (Q::FDFilter, R::FDFilterIF, info)

Solve the _approximate model-matching problem_ (AMMP) for a given synthesis model `sysf::FDIModel` with additive faults 
and a given stable reference filter `sysr::FDFilterIF`. 
The computed stable and proper filter objects `Q` and `R` contain the 
fault detection filter, representing a solution of the AMMP, and its internal form, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.degs`, `info.M`, `info.freq`, 
`info.HDesign`, `info.nonstandard`, `info.gammaopt0`, `info.gammaopt` and `info.gammasub`,  
contains additional synthesis related information (see below). 

The continuous- or discrete-time system `sysf.sys` is in a standard
or descriptor state-space form `sysf.sys = (A-λE,B,C,D)`, which corresponds to the input-output form  

       y = Gu(λ)*u + Gd(λ)*d + Gf(λ)*f + Gw(λ)*w + Ga(λ)*aux,

with the Laplace- or Z-transformed plant outputs `y`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Gu(λ)`, `Gd(λ)`, `Gf(λ)`, `Gw(λ)`, and `Ga(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysf.controls`, `sysf.disturbances`, `sysf.faults`, `sysf.noise` and `sysf.aux`, respectively.

The continuous- or discrete-time reference filter `sysr.sys` is in a standard
or descriptor state-space form `sysr.sys = (Ar-λEr,Br,Cr,Dr)`, which corresponds to the input-output form  

       yr = Mru(λ)*u + Mrd(λ)*d + Mrf(λ)*f + Mrw(λ)*w + Mra(λ)*aux,

with the Laplace- or Z-transformed reference filter outputs `yr`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Mru(λ)`, `Mrd(λ)`, `Mrf(λ)`, `Mrw(λ)`, and `Mra(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysr.controls`, `sysr.disturbances`, `sysr.faults`, `sysr.noise` and `sysr.aux`, respectively.
If any of the above vectors is void, then the corresponding transfer function matrix is considered null. 

The fault detection filter object `Q`, contains in `Q.sys` the resulting filter 
in a standard state-space form, which generates the residual signal `r`. 
The corresponding input-output (implementation) form is

            r = Qy(λ)*y + Qu(λ)*u               

where `Qy(λ)` and `Qu(λ)` are the transfer function matrices from the output and control inputs to the residual. 
The indices of output and control inputs are contained in the integer vectors 
`Q.outputs` and `Q.controls`, respectively.

Let define 

      Ge(λ) = | Gu(λ) Gd(λ) Gf(λ) Gw(λ) Ga(λ) |,  Mr(λ) = | Mru(λ) Mrd(λ) Mrf(λ) Mrw(λ) Mra(λ) | . 
              |  I     0     0     0     0    |

In the standard case, `Ge(λ)` has no zeros on the boundary of the 
stability domain, and the resulting stable filter `Q(λ) := |Qy(λ) Qu(λ)|` is 
`Q(λ) = Q0(λ)`, where `Q0(λ)` is
the optimal solution of the H∞- or H2-norm error minimization problem

    gammaopt0 = ||Q0(λ)*Ge(λ)-M(λ)*Mr(λ)|| = min,         (1)

where `M(λ) = M0(λ)` is an updating factor chosen as follows:
`M0(λ) = I` in the case of emplyoing the 
H∞ norm, while in the case of employing the H2 norm, `M0(λ) = I` for
a discrete-time system or, for a continuous-time system, 
`M0(λ)` is determined a stable, diagonal, and invertible transfer function 
matrix, which ensures the existence of a finite H2-norm.

In the non-standard case, `Ge(λ)`  has zeros on the boundary of the 
stability domain, and the resulting optimal filter `Q0(λ)`, which solves 
the H∞- or H2-norm error minimization problem (1) is a possibly 
unstable or improper. A second updating factor `M1(λ)` is determined, with   
the same properties as `M0(λ)`, which ensures that the computed stable and  
proper filter `Q(λ) := M1(λ)*Q0(λ)` represents a suboptimal solution of an 
updated H∞- or H2-norm error minimization problem, for which the   
achieved suboptimal model-matching performance is

    gammasub = ||Q(λ)*Ge(λ)-M(λ)*Mr(λ)|| ,            (2)

where `M(λ) := M1(λ)*M0(λ)`. The _optimal_ solution `Qt(λ)` of the 
updated H∞- or H2-norm error minimization problem 

    gammaopt = ||Qt(λ)*Ge(λ)-M(λ)*Mr(λ)|| = min ,     (3)  

is still possibly unstable or improper. The values of `gammaopt0`, `gammaopt` and `gammasub`
are returned in `info.gammaopt0`, `info.gammaopt` and `info.gammasub`, respectively.

The fault detection filter internal form object `R`, contains `R.sys`, the resulting 
internal form of the filter 
in a standard state-space form, which generates the residual signal `r`, and corresponds to the 
input-output form

       r = Ru(λ)*u + Rd(λ)*d + Rf(λ)*f + Rw(λ)*w + Ra(λ)*aux ,

where 

       | Ru(λ) Rd(λ) Rf(λ) Rw(λ) Ra(λ) | = Q(λ)*Ge(λ). 

The indices of the inputs `u`, `d`, `f`, `w` and `aux` of the resulting filter `R.sys` are 
contained in the integer vectors `R.controls`, `R.disturbances`, `R.faults`, `R.noise` and `R.aux`, respectively.
The state-space realization of the resulting `M(λ)` is returned in `info.M`.

Various user options can be specified via keyword arguments as follows:

If `H2syn = false` (default), a H∞-norm based synthesis is performed, while 
if `H2syn = true`, a H2-norm based synthesis is performed. 

`reltol = tol` specifies the relative tolerance `tol` for the desired 
accuracy of γ-iteration (default:  `tol = 1.e-4`).
   
If `mindeg = true`, a least order filter synthesis is performed, if possible, while 
with `minimal = false` (default) no least order synthesis is performed.  

If `regmin = true` (default), the regularization (see [1]) is performed for the case 
when `sysr.controls` and/or `sysr.disturbances` are void with the selection of 
a least order left annihilator `Nl(λ)` of `G(λ) = [Gu(λ) Gd(λ); I 0 ]`. 
If `regmin = false`, the regularization is performed by choosing a left annihilator
`Nl(λ)` as a minimal left nullspace basis of `G(λ)`.  

If `HDesign = H` is a full row rank design matrix, then `H*Nl(λ)` is used as left annihilator
instead `Nl(λ)` (default: `HDesign = missing`).

If `nullspace = true` (default) and `sysr.controls` and/or `sysr.disturbances` are void, 
a minimal proper nullspace basis is used at the initial reduction step.
If `nullspace = false` and `sysr.controls` and/or `sysr.disturbances` are void, 
a full-order observer based nullspace basis is used at the 
initial reduction step.
This option can  only be used for a proper system without disturbance inputs. 
The `nullspace` option is ignored if both `sysr.controls` and `sysr.disturbances` are non-void. 

If `simple = true`, a simple proper nullspace basis  
is emplyed as left annihilator for synthesis. 
The orders of the basis vectors are provided in `info.deg`. 
If `simple = false` (default), then a minimal proper nullspace basis is computed. 

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

`freq = val` specifies the value of a test frequency to be employed to 
check the full column rank (i.e., left-invertibility) solvability condition 
(default: randomly generated in the interval `(0,1)`). 
The employed value of `freq` is returned in `info.freq`.

`normalize = job` specifies the option for the normalization  
of the diagonal elements of the updating matrix `M(λ)` as follows:

      job = "gain"    – scale with the gains of the zero-pole-gain representation;
      job = "dcgain"  – scale with the DC-gains;
      job = "infnorm" – scale with the values of infinity-norms (default). 

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

The resulting named tuple `info` contains `(tcond, degs, M, freq, HDesign, gammaopt0, gammaopt, gammasub, nonstandard) `, where:

`info.tcond` is the maximum of the condition numbers of the employed 
   non-orthogonal transformation matrices; a warning is issued if `info.tcond >= tcmax`;

`info.degs` is an integer vector containing the increasingly ordered degrees of a left minimal   
polynomial nullspace basis of `G(λ) := [ Gu(λ) Gd(λ); I 0]` (also the left Kronecker indices of `G(λ)`), if the 
state-space realization of `[Gu(λ) Gd(λ)]` is minimal;

`info.M` is the employed stable and invertible updating filter used to solve the AMMP, 
with a stable, diagonal and invertible transfer function matrix `M(λ)`; 

`info.freq` is the employed frequency used to check left invertibility 
(set to `missing` if no frequency-based left invertibility check was performed)

`info.HDesign` is the design matrix `H` employed for the synthesis of 
   the fault detection filter `Q`; `H = missing` if no design matrix was involved;
   
`info.gammaopt0` is the optimal performance `gammaopt0` for the original problem (1); 
   
`info.gammaopt` is the optimal performance `gammaopt` for the updated problem (3); 

`info.gammasub` is the suboptimal performance `gammasub` in (2); 

`info.nonstandard` is set to `true` for a non-standard problem 
   (i.e., `Ge(λ)` has zeros on the boundary of the stability domain), and set to 
   `false` for a standard problem 
      (i.e., `Ge(λ)` has no zeros on the boundary of the stability domain). 

_Method:_ The synthesis Procedure AMMS from [1] is implemented. The 
Procedure AMMS relies on the approximate model-matching synthesis method 
proposed in [2]. For more details on computational aspects see [3].  

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017; sec. 5.6.

[2] A. Varga, Integrated computational algorithm for solving 
    H_inf-optimal FDI problems. In Proc. of the IFAC World Congress, 
    Milano, Italy, pp. 10187–10192, 2011.

[3] A. Varga. Descriptor system techniques in solving H_2/H-Inf-optimal
    fault detection and isolation problems". In L. T. Biegler,  
    S. L. Campbell, and V. Mehrmann (Eds.), Control and Optimization 
    with Differential-Algebraic Constraints, vol. 23 of Advances in 
    Design and Control, pp. 105–125. SIAM, 2012. 
"""
function ammsyn(sysf::FDIModel{T1}, sysr::Union{FDFilterIF{T2},FDIModel{T2}}; poles::Union{AbstractVector,Missing} = missing, 
                      sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing, 
                      nullspace::Bool = true, mindeg::Bool = false, simple::Bool = false, regmin::Bool = true, 
                      normalize::AbstractString = "infnorm", freq::Real = rand(), reltol::Real = 0.0001, H2syn::Bool = false, 
                      tcond::Real = 1.e4, HDesign::Union{AbstractMatrix,Missing} = missing,
                      offset::Real = sqrt(eps(float(real(T1)))), atol::Real = zero(float(real(T1))), 
                      atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                      rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T1)))))*iszero(max(atol1,atol2)), 
                      fast::Bool = true) where {T1,T2}

   Ts = DescriptorSystems.promote_Ts(sysf.sys.Ts,sysr.sys.Ts)
   disc = (Ts != 0);  # system type (continuous- or discrete-time)
   T = promote_type(T1,T2)

   rdim = size(sysr.sys,1)    
 
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
   
   # set default stability degree
   sdegdefault = disc ? 0.95 : -0.05
   ismissing(sdeg) && ismissing(poles) && (sdeg = sdegdefault)  # set desired stability degree to default value
  
   # imposed design option to form linear combinations of basis vectors
   emptyHD = ismissing(HDesign)
   if !emptyHD
      !ismissing(rdim) && size(HDesign,1) != rdim && error("row dimension of HDesign must be equal to rdim")
      size(HDesign,1) == rank(HDesign) || error("HDesign must have full row rank")
   end
   normflag = H2syn ? 2 : Inf
  
   # decode input information
   inpu = sysf.controls; mu = length(inpu)  
   inpd = sysf.disturbances; md = length(inpd) 
   inpf = sysf.faults; mf = length(inpf)  
   inpw = sysf.noise;  mw = length(inpw) 
   inpaux = sysf.aux;  maux = length(inpaux)  
   mru = length(sysr.controls); 
   mru == 0 || mru == mu ||  error("Incompatible control input groups in sysf and sysr") 
   mrd = length(sysr.disturbances); 
   mrd == 0 || mrd == md ||  error("Incompatible disturbance input groups in sysf and sysr") 
   mrf = length(sysr.faults); 
   mrf == 0 || mrf == mf ||  error("Incompatible fault input groups in sysf and sysr") 
   mrw = length(sysr.noise); 
   mrw == 0 || mrw == mw ||  error("Incompatible noise input groups in sysf and sysr") 
   mra = length(sysr.aux); 
   mra > 0 || mra == maux ||  error("Incompatible auxiliary input groups in sysf and sysr") 
       
   m = mu+md+mf+mw+maux;       # total number of inputs
   p = size(sysf.sys,1);       # number of measurable outputs
    
   # if mf == 0 && minimal
   #    @warn "Minimal synthesis option not feasible in the case of no faults"
   #    minimal = false
   # end

   rd = gnrank(sysf.sys[:,inpd]; atol1, atol2, rtol)
   SMr = fditspec(sysr; atol1, atol2, rtol)
   
   if mru+mrd+mrw+mra == 0 
      weakfdi = false
      rref = gnrank(sysr.sys[:,sysr.faults]; atol1, atol2, rtol)
      if rdim <= p-rd
         strongfdi = (rdim == mf && rref == mf)
         strongfdi || (weakfdi = !all(maximum(SMr,dims=1)))
      else
         weakfdi = true
         strongfdi = false
      end
   else
      weakfdi = true
      strongfdi = false
   end
   
   if mru+mrd+mrw+mra == 0 && !weakfdi
      # perform either Procedure EMM or EMMS if SYSR has no 'controls', 
      # 'disturbances', 'noise' and 'aux' input groups and no least order option is selected 
      # Step 1): nullspace based reduction
      #
      desc = (sysf.sys.E != I)
      m2 = mf+mw+maux
      if nullspace || md > 0 || (desc && rcond(sysf.sys.e) < 1.e-7 )
         # form [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
         #syse = [sysf(:,[inpu inpd inpf inpw inpaux]); eye(mu,m)];
         syse = [sysf.sys; eye(mu,m)];
         #
         # compute a left nullspace basis Q = Q1 of G1 = [Gu Gd; I 0] = 0 and
         # obtain QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = Q*[Gf Gw Gaux;0 0 0]
         QR, info1 = glnull(syse, m2; simple, atol1, atol2, rtol, fast, sdeg = sdegdefault, offset) 
         tcond1 = info1.tcond
         degs = info1.degs
      elseif mu == 0 && md == 0
         # compute minimal basis as Q = Q1 = I  and set
         # QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = [ Gf Gw Gaux ]
         QR = [eye(p) sysf.sys]
         tcond1 = 1.
         degs = Int[]
      else
         # compute minimal basis as Q = Q1 = [ I -Gu] and set
         # QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = [ Gf Gw Gaux ]
         QR = [ eye(p) dss(sysf.sys.A, sysf.sys.E, [-sysf.sys.B[:,inpu] sysf.sys.B[:,[inpf; inpw; inpaux]]],
                sysf.sys.C, [-sysf.sys.D[:,inpu] sysf.sys.D[:,[inpf; inpw; inpaux]]]; Ts)] 
         # perform stabilization if strong detectability has to be enforced
         tcond1 = 1.
         degs = Int[]
      end
                     
      nvec = size(QR,1);   # number of basis vectors
      # check solvability conditions
      nvec == 0 && error("ammsyn: empty nullspace basis - the nullspace-based approach to solve the AMMP is not applicable")
      infodegs = degs

      indf = (p+mu) .+ Vector(1:mf)            # input indices of Rf1 in QR
      indfwa = (p+mu) .+ Vector(1:mf+mw+maux)  # input indices of [ Rf1 Rw1 Raux1] in QR
      mfw = mf+mw
      mfwa = mfw+maux
      Rfwa = QR[:,indfwa]
      indR = (p+mu) .+ (1:mfw);     # input indices of R = [Rf1 Rw1] in QR
   
      Rftest = evalfr(QR[:,indf], freq; atol1, atol2, rtol)
   
      if  !weakfdi
         # set H for checking the solvability condition
         if emptyHD
            Htemp = eye(nvec)
         else
            degs = Int[];
            nh = size(HDesign,2)
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
         if strongfdi
            # check strong isolability
            rf = rank(Htemp*Rftest; atol, rtol)
            mf == rf || @warn "The system sysf is not strongly isolable"
         else
            # check complete detectability
            S = fditspec_(Htemp*QR[:,indf]; atol1, atol2, rtol)
            all(maximum(S,dims=1)) || error("The system sysf is not completely detectable")
         end
         if nvec > rdim
            # 2) Determine Q2 such that Q2*Rf1 has maximal full row rank
            #    and Q2*Q1 has least McMillan degree
            # flip degrees of a minimal polynomial basis
            reverse!(degs)
            finish = false
            nq = order(QR)
            emptyHD && (h = eye(nvec))
            nout = rdim
             
            while !finish     
                # choose nout basis vectors, which potentially lead to a least order
                # filter with rdim outputs:
                # basesel(i,:) contains the indices of candidate basis vectors;
                # ordsel(i)    contains the presumably achievable least orders
                basesel, ordsel = ammbasesel(Htemp*Rftest, degs, nout, simple, atol,strongfdi) 
                #
                # update the synthesis using the selections of candidate vector(s),
                # starting with the least (potentially) achievable order
                for i = 1:size(basesel,1)
                    baseind = basesel[i] # indices of current basis selection
                    if nout == mf
                       hbase = eye(mf);
                    else
                       hbase = rand(mf,nout); 
                    end
                    ip = [baseind; setdiff(1:nvec,baseind)][:]
                    if simple
                       # handle simple basis
                       # here only the case rdim = nout can happen
                       if regmin
                          if emptyHD 
                             # select vectors and elliminate unobservable dynamics  
                             noelim = falses(nq) 
                             ell = sum(degs[1:basesel[i][1]-1]); 
                             for jj = 1:nout 
                                 ellnext = sum(degs[1:baseind[jj]]);
                                 noelim[ell+1:ellnext] .= true;
                                 ell = ellnext
                             end
                             ir = noelim
                             Ar, Er, Br, Cr, Dr = dssdata(QR[baseind,:])
                             QRfwtest = dss(view(Ar,ir,ir), Er == I ? I : view(Er,ir,ir),view(Br,ir,:),view(Cr,:,ir), Dr; Ts)
                             h = Htemp[ip[1:rdim],:]
                          else
                             QRfwtest = gir(Htemp*QR[ip,:]; atol1, atol2, rtol)
                          end
                       else
                          if emptyHD 
                              h = Htemp[ip[1:rdim],:]
                              QRfwtest = gir(QR[baseind,:]; atol1, atol2, rtol, infinite = false)
                           else
                              QRfwtest = gir(Htemp*QR; atol1, atol2, rtol, infinite = false)
                           end
                       end
                    else
                       # handle minimal basis
                       # use the above output permutation vector ip for glmcover1  
                       if regmin
                          if rdim == nout
                             if emptyHD 
                                QRfwtest, _, info2 = glmcover1(QR[ip,:], rdim; atol1, atol2, rtol)
                                if !isempty(ordsel) && (order(QRfwtest) != ordsel[i])
                                   @warn "ammsyn: expected reduced order not achieved"
                                end
                                h = Htemp[ip[1:rdim],:]
                             else
                                QRfwtest, _, info2 = glmcover1([Htemp; eye(nvec)]*QR[ip,:], rdim; atol1, atol2, rtol)
                             end
                          else
                             # the case rdim < nout can only happen if no
                             # HDesign is explicitly provided
                             Htemp = blockdiag(hbase,eye(nvec-nout)) 
                             QRfwtest, _, info2 = glmcover1([Htemp; eye(nvec)]*QR[ip,:], rdim; atol1, atol2, rtol)
                          end
                       else
                          # here only the case rdim = nout can happen
                          if emptyHD
                             h = Htemp[ip[1:rdim],:]
                             QRfwtest = gir(QR[baseind,:]; atol1, atol2, rtol, infinite = false) 
                          else
                             QRfwtest = gir(Htemp*QR; atol1, atol2, rtol, infinite = false) 
                          end
                       end
                    end
                    # check admissibility of compressed Rf1; 
                    if !simple && regmin
                       # dismiss minimal design if the check fails
                       Rftest1 = evalfr(QRfwtest[:,indf],freq; atol1, atol2, rtol) 
                       if strongfdi
                          rank(Rftest1; atol) == nout && (finish = true)
                       else
                          beta = norm(view(Rftest1,:,1)) 
                          for j = 2:mf
                              beta  = min(beta,norm(view(Rftest1,:,j)))
                          end
                          ((atol > 0 && beta > atol) || (atol == 0 && beta > nvec*mf*eps(norm(Rftest1)))) && (finish = true)
                        end
                        if finish
                           # adjust condition number of employed transformations
                           tcond1 = max(tcond1, info2.fnorm, info2.tcond)
                           tcond1 > tcond && 
                               @warn "ammsyn: possible loss of numerical stability due to ill-conditioned transformations"
                           QR = QRfwtest
                           break
                        end
                    end
                end
                nout += 1
                nout > nvec && !finish &&
                        error("something wrong: try perhaps with another test frequency")
            end
            emptyHD && (Htemp = h)
         end
         # 3): determine Q3 = inv(Ro), from the extended quasi-co-outer-co-inner
         #     factorization R = [Ro 0]Ri;
         #     update Q <- Q3*Q and R <- Q3*R
         Ri, Ro, info1 = goifac(QR[:,indR]; atol1, atol2, rtol)
         lo, ro = size(Ro)
       
         # detect nonstandard problem
         nonstandard = (info1.nfuz+info1.niuz > 0)

         if lo == ro && (order(Ro) == order(QR)) 
            # extract descriptor state-space data
            ae, ee, be, ce, de = dssdata(QR)
            # form [Ro Q R] (recall that Ro and QR share the same A, E and C matrices)
            # check concatanation formulas
            RoQR = dss(ae, ee, [Ro.B be], ce, [Ro.D de]; Ts)
            # form QR = inv(Ro)*[Q R]
            QR = grsol(RoQR,p+mu+length(indR); atol1, atol2, rtol)[1]
         else
            # compute a left inverse of Ro with stable free poles
            Roinv = glsol([Ro;eye(ro)],ro; atol1, atol2, rtol, sdeg)[1] 
            # update QR <- Roinv*QR 
            QR = gir(Roinv*QR; atol1, atol2, rtol, infinite = false)
         end
         
         # form explicitly Rref = [ Mrf 0 ]
         Rref = [ sysr.sys[:,sysr.faults] zeros(rdim,mw)];
       
         # compute F = [ F1 F2 ] =  Rref*Gi'
         F = Rref*Ri';

         # 4): Compute the solution Q4 of the LDP min||[F1-Q4 F2]||
         #     and update Q <- Q4*Q and R <- Q4*R
         if H2syn
            if atol1 > 0
               told = atol1;
            else
               told = 1.e4*eps;
            end
            if !disc && norm(F.D[:,ro+1:end],Inf) < told
               Mib = dss(eye(T,rdim); Ts) 
               F.D[:,ro+1:end] = zeros(size(F,1),mfw-ro)
            else
               # update Mr to achieve finite L2-norm
               # Mib = dss(sdeg, 1, -sdeg, 0; Ts)*eye(rdim);
               # mib = dss(sdeg, 1, -sdeg, 0; Ts)
               # Mib = mib
               # for i = 2:rdim
               #    Mib = append(Mib,mib)
               # end
               # Rref = Mib*Rref
               Mib = dsdiag(dss(sdeg, 1, -sdeg, 0; Ts), rdim)
               Rref =  Mib * Rref
               F = Rref*Ri'
            end
            # solve the H2-LDP min||[ F1-Q4 F2 ] ||_2
            Qt, Qtu = gsdec(F[:,1:ro]; atol1, atol2, rtol, job = "stable")
            QR = Qt*QR;  #  update [ Q R ] <- Q4*[ Q R ]
            gopt0 = gl2norm(gir([Qtu F[:,ro+1:end]]; atol1, atol2, rtol); offset, atol, atol1, atol2, atol3, rtol, fast )
         else
            # solve the H_inf-LDP min ||[ F1-Q4 F2 ] ||_inf
            Qt, gopt0 = glinfldp(F,mfw-ro; reltol, offset, atol1, atol2, rtol, fast )
            QR = Qt*QR;  #  update [ Q R ] <- Q4*[ Q R ]
            Mib = dss(eye(T,rdim); Ts)
         end

         if nonstandard
            # 5) compute diagonal Q5 such that Q5*Q has a desired stability degree;
            #    update Q <- Q5*Q and R <- Q5*R
            Qt = dss(zeros(0,size(QR,2)); Ts); M = dss(zeros(0,0); Ts)
            for i = 1:rdim
                Qti, Mi = glcf(gir(QR[i,:]; atol1, atol2, rtol, contr = false); atol1, atol2, atol3, rtol, sdeg, smarg, evals = poles, 
                                                          mininf = true, mindeg = true)
                if isequal(normalize,"gain") 
                   sc = zpk(dss2rm(Mi)[1,1])[3]
                elseif isequal(normalize,"dcgain")
                   sc = dcgain(Mi)[1,1]
                else
                   sc = ghinfnorm(Mi)[1]
                end
                Qt = [Qt; Qti/sc]; M = append(M,Mi/sc)
            end
            Mib = Mib*M
            QR = Qt
            # compute the optimal distance for the updated problem 
            if H2syn
               _, Qtu = gsdec(Mib*F[:,1:ro]; atol1, atol2, rtol, job = "stable")
               gopt = gl2norm(gir([Qtu Mib*F[:,ro+1:end]]; atol1, atol2, rtol); offset, atol, atol1, atol2, atol3, rtol, fast )
            else
               _, gopt = glinfldp(Mib*F,mfw-ro; offset, atol1, atol2, rtol, fast )
            end
         else
            gopt = gopt0
         end
         # transform to standard state-space
         QR = gss2ss(QR; atol1, atol2, rtol)[1]  
         Q = FDFilter(QR, p, mu)
         R = FDFilterIF(QR, mru, mrd, mrf, mw > 0 ? mw : mrw, maux > 0 ? maux : mra; moff = p+mu)
         M = Mib     
      end
   end
   if weakfdi
      # apply the three-step procedure to solve the AMMP if SYSR has  
      # 'controls', 'disturbances' or 'noise' input groups
      degs = []; 
      finish = false;
      if mru+mrd == 0
         # perform nullspace-based synthesis if SYSR has no 'controls' and
         # 'disturbances' input groups
         # 1) Compute Q1, the left nullspace of [Gu Gd;I 0], and  
         #    R = [Rf1 Rw1] = Q1*[Gf Gw;0 0]
         #    QR contains Q1 in QR(:,1:p+mu) and R in QR(:,p+mu+1:end),
         #    where Rf1 is in QR(:,p+mu+1:p+mu+mf) and 
         #    Rw1 is in QR(:,p+mu+mf+1:end)
   
         desc = (sysf.sys.E != I)
         mr = mf+mw+maux
         if nullspace || md > 0 || (desc && rcond(sysf.sys.e) < 1.e-7 )
            # form [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
            #syse = [sysf(:,[inpu inpd inpf inpw inpaux]); eye(mu,m)];
            syse = [sysf.sys; eye(mu,m)];
            #
            # compute a left nullspace basis Q = Q1 of G1 = [Gu Gd; I 0] = 0 and
            # obtain QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = Q*[Gf Gw Gaux;0 0 0]
            QR, info1 = glnull(syse, mr; simple, atol1, atol2, rtol, fast, sdeg = sdegdefault, offset) 
            tcond1 = info1.tcond
        elseif mu == 0 && md == 0
            # compute minimal basis as Q = Q1 = I  and set
            # QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = [ Gf Gw Gaux ]
            QR = [eye(p) sysf.sys]
            tcond1 = 1.
         else
            # compute minimal basis as Q = Q1 = [ I -Gu] and set
            # QR = [ Q1 R1 ], where R1 = [Rf1 Rw1 Raux1] = [ Gf Gw Gaux ]
            QR = [ eye(p) dss(sysf.sys.A, sysf.sys.E, [-sysf.sys.B[:,inpu] sysf.sys.B[:,[inpf; inpw; inpaux]]],
                   sysf.sys.C, [-sysf.sys.D[:,inpu] sysf.sys.D[:,[inpf; inpw; inpaux]]]; Ts)] 
            # perform stabilization if strong detectability has to be enforced
            tcond1 = 1.
         end
         if size(QR,1) > 0
            # non-empty nullspace
            finish = true
            # separate implementation and internal forms
            Q1 = gir(QR[:,1:p+mu]; atol1, atol2, rtol )
            R1 = QR[:,(p+mu) .+ (1:mr)]
         end
      end
      if !finish && mru == 0 
         # 1) Compute Q1, the left nullspace of [Gu;I], and  
         #    R = [Rd1 Rf1 Rw1] = Q1*[Gd Gf Gw;0 0 0]
         #    QR contains Q1 in QR(:,1:p+mu) and R in QR(:,p+mu+1:end),
         #    where Rd1 is in QR(:,p+mu+1:p+mu+md), Rf1 is in 
         #    QR(:,p+mu+md+1:p+mu+md+mf) and Rw1 is in QR(:,p+mu+md+mf+1:end)
      
         # set options for nullspace computation
         mr = md+mf+mw+maux;
         if nullspace || (desc && rcond(sysf.sys.e) < 1.e-7 )
            # form [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
            syse = [sysf.sys; eye(mu,m)];
            #
            # compute a left nullspace basis Q = Q1 of G1 = [Gu; I] = 0 and
            # obtain QR = [ Q1 R1 ], where R1 = [Rd1 Rf1 Rw1 Raux1] = Q*[Gd Gf Gw Gaux;0 0 0 0]
            QR, info1 = glnull(syse, mr; simple, atol1, atol2, rtol, fast, sdeg = sdegdefault, offset) 
            tcond1 = info1.tcond
        elseif mu == 0 
            # compute minimal basis as Q = Q1 = I  and set
            # QR = [ Q1 R1 ], where R1 = [Rd1 Rf1 Rw1 Raux1] = [ GD Gf Gw Gaux ]
            QR = [eye(p) sysf.sys]
            tcond1 = 1.
         else
            # compute minimal basis as Q = Q1 = [ I -Gu] and set
            # QR = [ Q1 R1 ], where R1 = [ Rd1 Rf1 Rw1 Raux1] = [ Gd Gf Gw Gaux ]
            QR = [ eye(p) dss(sysf.sys.A, sysf.sys.E, [-sysf.sys.B[:,inpu] sysf.sys.B[:,[inpd; inpf; inpw; inpaux]]],
                   sysf.sys.C, [-sysf.sys.D[:,inpu] sysf.sys.D[:,[inpd; inpf; inpw; inpaux]]]; Ts)] 
            # perform stabilization if strong detectability has to be enforced
            tcond1 = 1.
         end
         if size(QR,1) > 0
            # non-empty nullspace
            finish = true
            # separate implementation and internal forms
            Q1 = gir(QR[:,1:p+mu]; atol1, atol2, rtol )
            R1 = QR[:,(p+mu) .+ (1:mr)]
         end
      end
      if !finish && mrd == 0 
         # 1) Compute Q1, the left nullspace of [Gd;0], and  
         #    R = [Ru1 Rf1 Rw1] = Q1*[Gu Gf Gw;I 0 0]
         #    QR contains Q1 in QR(:,1:p+mu) and R in QR(:,p+mu+1:end),
         #    where Ru1 is in QR(:,p+mu+1:p+mu+mu), Rf1 is in 
         #    QR(:,p+mu+mu+1:p+mu+mu+mf) and Rw1 is in QR(:,p+mu+mu+mf+1:end)
         # set options for nullspace computation
         mr = mu+mf+mw+maux;
         # form [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
         syse = [sysf.sys; eye(mu,m)];
         #
         # compute a left nullspace basis Q = Q1 of G1 = [Gu; I] = 0 and
         # obtain QR = [ Q1 R1 ], where R1 = [Rd1 Rf1 Rw1 Raux1] = Q*[Gd Gf Gw Gaux;0 0 0 0]
         QR, info1 = glnull(syse, mr; simple, atol1, atol2, rtol, fast, sdeg = sdegdefault, offset) 
         tcond1 = info1.tcond
         if size(QR,1) > 0
            # non-empty nullspace
            finish = true
            # separate implementation and internal forms
            Q1 = gir(QR[:,1:p+mu]; atol1, atol2, rtol )
            R1 = QR[:,(p+mu) .+ (1:mr)]
         end
      end
      if !finish
         # 1) Set Q1 = I and R = [Ru1 Rd1 Rf1 Rw1] = Q1*[Gu Gd Gf Gw;I 0 0 0]
         #    QR contains Q1 in QR(:,1:p+mu) and R in QR(:,p+mu+1:end),
         #    where Ru1 is in QR(:,p+mu+1:p+mu+mu), Rd1 is in 
         #    QR(:,p+mu+1:p+mu+mu+md), Rf1 is in QR(:,p+mu+mu+md+1:p+mu+mu+md+mf)
         #    and Rw1 is in QR(:,p+mu+mu+md+mf+1:end)
         mr = m;
         Q1 = dss(eye(p+mu);Ts)
         R1 = [sysf; eye(mu,m)]
         syse = R1;
         tcond1 = 1;
      end
           
      if mra == 0
         # apply the two-step procedure for the case
         # form Ge = [ Gu Gd Gf Gw; I 0 0 0] 
         m1 = m-maux;
         syse = [sysf.sys[:,1:m1]; eye(mu,m1)]
         rinp = zeros(0,m1)
         mru == 0 || (rinp = [rinp; eye(mu,m1)])
         mrd == 0 || (rinp = [rinp; zeros(md,mu) eye(md,m1-mu)])
         mrf == 0 || (rinp = [rinp; zeros(mf,mu+md) eye(mf,m1-mu-md)])
         mrw == 0 || (rinp = [rinp; zeros(mw,mu+md+mf) eye(mw,m-mu-md-mf)])
         # form explicitly Rref = [ Mru_i Mrd_i Mrf_i ] 
         Rref = sysr.sys*rinp;
      else
         # apply the two-step procedure for the general case
         # form Ge = [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
         syse = [sysf.sys; eye(mu,m)];
         rinp = zeros(0,m);
         mru == 0 || (rinp = [rinp; eye(mu,m)])
         mrd == 0 || (rinp = [rinp; zeros(md,mu) eye(md,m-mu)])
         mrf == 0 || (rinp = [rinp; zeros(mf,mu+md) eye(mf,m-mu-md)])
         mrw == 0 || (rinp = [rinp; zeros(mw,mu+md+mf) eye(mw,m-mu-md-mf)])
         mra == 0 || (rinp = [rinp; zeros(maux,m-maux) eye(maux)])
         # form explicitly Rref = [Mru_i Mrd_i Mrf_i Mrw_i ] 
         Rref = sysr.sys*rinp;
      end
      # 2) Compute the approximate solution of Q2*Ge = Rref .
      Q2, info1 = glasol(R1, Rref[:,m-mr+1:end]; reltol, atol1, atol2, rtol, mindeg, L2sol = H2syn, sdeg, poles, fast, offset) 
      Qib = gir(Q2*Q1; atol1, atol2, rtol)
      gopt0 = info1.mindist
      nonstandard = info1.nonstandard
      if nonstandard
         # 3) compute diagonal M such that M*Q has a desired stability degree;  
         #    update Q <- M*Q
         Qt = dss(zeros(0,p+mu); Ts); Mib = dss(zeros(0,0); Ts)
         for i = 1:rdim
             Qti, Mi = glcf(gminreal(Qib[i,:]; atol1, atol2, rtol); atol1, atol2, atol3, rtol, sdeg, smarg, evals = poles, 
                                                          mininf = true, mindeg = true)
             if isequal(normalize,"gain") 
                sc = zpk(dss2rm(Mi; atol1, atol2, rtol)[1,1])[3]
             elseif isequal(normalize,"dcgain")
                sc = dcgain(Mi)[1,1]
             else
                sc = ghinfnorm(Mi)[1]
             end
             Qt = [Qt; Qti/sc]; Mib = append(Mib,Mi/sc)
         end
         Qib = Qt
         # compute optimal performance
         _, info1 = glasol(R1, Mib*Rref[:,m-mr+1:end]; reltol, atol1, atol2, rtol, mindeg, L2sol = H2syn, sdeg, poles, fast, offset) 
         gopt = info1.mindist 
         M = Mib  
      else
         gopt = gopt0
         M = dss(eye(T,rdim);Ts)
      end
      Q = FDFilter(gss2ss(Qib; atol1, atol2, rtol)[1], p, mu)

      R = FDFilterIF(gir(Q.sys*syse[:,m-mr+1:end]; atol1, atol2, rtol), mru, mrd, mrf, mw > 0 ? mw : mrw, maux > 0 ? maux : mra)
      infodegs = Int[]
      Htemp = missing
      freq = missing
   end   

   gsub = nonstandard ? fdimmperf(R,M*sysr,normflag; atolinf = atol1) : gopt
   # transform to standard state-space  
   #Q = FDFilter(gss2ss(Qt; atol1, atol2, rtol)[1], p, mu)

   #R = FDFilterIF(gss2ss(Rt; atol1, atol2, rtol)[1], mru, mrd, mrf, mw > 0 ? mw : mrw, maux > 0 ? maux : mra)
   info = (tcond = tcond1, degs = infodegs, M = M, freq = freq, HDesign = ismissing(Htemp) ? missing : convert(Matrix{Float64},Htemp),
           gammaopt0 = gopt0, gammaopt = gopt, gammasub = gsub, nonstandard = nonstandard)
   return Q, R, info
   
   # end AMMSYN
end
"""
    ammsyn(sysf::FDIModel, sysr::FDIFilterIF; nullspace = true, simple = false, mindeg = false, 
                           regmin = true, normalize = "infnorm", H2syn = false, reltol = 1.e-4, 
                           sdeg, smarg, poles, freq, HDesign, tcond, offset, 
                           atol, atol1, atol2, atol3, rtol, fast = true) 
                           -> (Q::FDFilter, R::FDFilterIF, info)

Solve the _approximate model-matching problem_ (AMMP) for a given synthesis model `sysf::FDIModel` with additive faults 
and a given bank of stable reference filters `sysr::FDIFilterIF`.
The computed stable and proper filter objects `Q` and `R` contain the 
bank of fault detection filters, representing the component-wise solution of the AMMP, and 
their internal forms, respectively.

The returned named tuple `info`, with the components `info.tcond`, `info.degs`, `info.M`, `info.freq`, 
`info.HDesign`, `info.nonstandard`, `info.gammaopt0`, `info.gammaopt` and `info.gammasub`,  
contains additional synthesis related information (see below). 

The continuous- or discrete-time system `sysf.sys` is in a standard
or descriptor state-space form `sysf.sys = (A-λE,B,C,D)`, which corresponds to the input-output form  

       y = Gu(λ)*u + Gd(λ)*d + Gf(λ)*f + Gw(λ)*w + Ga(λ)*aux,

with the Laplace- or Z-transformed plant outputs `y`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Gu(λ)`, `Gd(λ)`, `Gf(λ)`, `Gw(λ)`, and `Ga(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysf.controls`, `sysf.disturbances`, `sysf.faults`, `sysf.noise` and `sysf.aux`, respectively.

The continuous- or discrete-time reference filters packed in `sysr` 
are in a standard or descriptor state-space form, where the `i`-th filter 
`sysr.sys[i] = (Ari-λEri,Bri,Cri,Dri)` corresponds to the input-output form  

       yri = Mrui(λ)*u + Mrdi(λ)*d + Mrfi(λ)*f + Mrwi(λ)*w + Mrai(λ)*aux,

with the Laplace- or Z-transformed reference filter outputs `yri`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Mrui(λ)`, `Mrdi(λ)`, `Mrfi(λ)`, `Mrwi(λ)`, and `Mrai(λ)` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysr.controls`, `sysr.disturbances`, `sysr.faults`, `sysr.noise` and `sysr.aux`, respectively.
If any of the above vectors is void, then the corresponding transfer function matrices are considered null. 

The fault detection and isolation filter object `Q`, contains in its `i`-th 
component `Q.sys[i]` the resulting filter 
in a standard state-space form, which generates the `i`-th component  `ri` of the residual signal. 
The corresponding input-output (implementation) form is

            ri = Qyi(λ)*y + Qui(λ)*u               

where `Qyi(λ)` and `Qui(λ)` are the transfer function matrices from the output and control inputs to the residual. 
The indices of output and control inputs are contained in the integer vectors 
`Q.outputs` and `Q.controls`, respectively.

Let define 

      Ge(λ) = | Gu(λ) Gd(λ) Gf(λ) Gw(λ) Ga(λ) |,  Mri(λ) = | Mrui(λ) Mrdi(λ) Mrfi(λ) Mrwi(λ) Mrai(λ) | . 
              |  I     0     0     0     0    |

In the standard case, `Ge(λ)` has no zeros on the boundary of the 
stability domain, and each resulting component stable filter `Qi(λ) := |Qyi(λ) Qui(λ)|` is 
`Qi(λ) = Q0i(λ)`, where `Q0i(λ)` is
the optimal solution of the H∞- or H2-norm error minimization problem

    gammaopt0i = ||Q0i(λ)*Ge(λ)-Mi(λ)*Mri(λ)|| = min,         (1)

where `Mi(λ) = M0i(λ)` is an updating factor chosen as follows:
`M0i(λ) = I` in the case of emplyoing the 
H∞ norm, while in the case of employing the H2 norm, `M0i(λ) = I` for
a discrete-time system or, for a continuous-time system, 
`M0i(λ)` is determined a stable, diagonal, and invertible transfer function 
matrix, which ensures the existence of a finite H2-norm.

In the non-standard case, `Ge(λ)`  has zeros on the boundary of the 
stability domain, and each resulting optimal component filter `Q0i(λ)`, which solves 
the H∞- or H2-norm error minimization problem (1) is a possibly 
unstable or improper. A second updating factor `M1i(λ)` is determined, with   
the same properties as `M0i(λ)`, which ensures that the computed stable and  
proper filter `Qi(λ) := M1i(λ)*Q0i(λ)` represents a suboptimal solution of an 
updated H∞- or H2-norm error minimization problem, for which the   
achieved suboptimal model-matching performance is

    gammasubi = ||Qi(λ)*Ge(λ)-Mi(λ)*Mri(λ)|| ,            (2)

where `Mi(λ) := M1i(λ)*M0i(λ)`. The _optimal_ solutions `Qti(λ)` of the 
updated H∞- or H2-norm error minimization problem 

    gammaopti = ||Qti(λ)*Ge(λ)-Mi(λ)*Mri(λ)|| = min ,     (3)  

is still possibly unstable or improper. The values of `gammaopt0i`, `gammaopti` and `gammasubi`
are returned in the `i`-th components of the vectors `info.gammaopt0`, `info.gammaopt` and `info.gammasub`, respectively.

The fault detection and isolation filter internal form object `R`, 
contains in its `i`-th component `R.sys[i]`, the resulting 
internal form of the filter 
in a standard state-space form, 
which generates the `i`-th component `ri` of residual signal , and corresponds to the 
input-output form

       ri = Rui(λ)*u + Rdi(λ)*d + Rfi(λ)*f + Rwi(λ)*w + Rai(λ)*aux ,

where 

       | Rui(λ) Rdi(λ) Rfi(λ) Rwi(λ) Rai(λ) | = Qi(λ)*Ge(λ). 

The indices of the inputs `u`, `d`, `f`, `w` and `aux` of the resulting filter `R.sys` are 
contained in the integer vectors `R.controls`, `R.disturbances`, `R.faults`, `R.noise` and `R.aux`, respectively.
The state-space realization of the resulting `Mi(λ)` is returned in the `i`-th component of the
vector `info.M`. 

Various user options can be specified via keyword arguments as follows:

If `H2syn = false` (default), a H∞-norm based synthesis is performed, while 
if `H2syn = true`, a H2-norm based synthesis is performed. 

`reltol = tol` specifies the relative tolerance `tol` for the desired 
accuracy of γ-iteration (default:  `tol = 1.e-4`).
   
If `mindeg = true`, least order filter syntheses are performed, if possible, while 
with `minimal = false` (default) no least order synthesis are performed.  

If `regmin = true` (default), the regularization (see [1]) is performed for the case 
when `sysr.controls` and/or `sysr.disturbances` are void with the selection of 
a least order left annihilator `Nl(λ)` of `G(λ) = [Gu(λ) Gd(λ); I 0 ]`. 
If `regmin = false`, the regularization is performed by choosing a left annihilator
`Nl(λ)` as a minimal left nullspace basis of `G(λ)`.  

If `HDesign = H` is a vector of full row rank design matrices, 
then `H[i]*Nl(λ)` is used as left annihilator
instead `Nl(λ)` for the synthesis of the `i`-th filter (default: `HDesign = missing`).

If `nullspace = true` (default) and `sysr.controls` and/or `sysr.disturbances` are void, 
a minimal proper nullspace basis is used at the initial reduction step.
If `nullspace = false` and `sysr.controls` and/or `sysr.disturbances` are void, 
a full-order observer based nullspace basis is used at the 
initial reduction step.
This option can  only be used for a proper system without disturbance inputs. 
The `nullspace` option is ignored if both `sysr.controls` and `sysr.disturbances` are non-void. 

If `simple = true`, a simple proper nullspace basis  
is emplyed as left annihilator for synthesis. 
The orders of the basis vectors are provided in `info.deg`. 
If `simple = false` (default), then a minimal proper nullspace basis is computed. 

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

`freq = val` specifies the value of a test frequency to be employed to 
check the full column rank (i.e., left-invertibility) solvability condition 
(default: randomly generated in the interval `(0,1)`). 
The employed value of `freq` is returned in `info.freq`.

`normalize = job` specifies the option for the normalization  
of the diagonal elements of the updating matrices `Mi(λ)` as follows:

      job = "gain"    – scale with the gains of the zero-pole-gain representation;
      job = "dcgain"  – scale with the DC-gains;
      job = "infnorm" – scale with the values of infinity-norms (default). 

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

The resulting named tuple `info` contains `(tcond, degs, M, freq, HDesign, gammaopt0, gammaopt, gammasub, nonstandard) `, where:

`info.tcond` is the maximum of the condition numbers of the employed 
   non-orthogonal transformation matrices; a warning is issued if `info.tcond >= tcmax`;

`info.degs` is an integer vector containing the increasingly ordered degrees of a left minimal   
polynomial nullspace basis of `G(λ) := [ Gu(λ) Gd(λ); I 0]` (also the left Kronecker indices of `G(λ)`), if the 
state-space realization of `[Gu(λ) Gd(λ)]` is minimal;

`info.M` is a vector of descriptor systems, whose `i`-th system `info.M[i]` 
contains the employed stable and invertible updating filter used to solve 
the `i`-th AMMP, with a diagonal transfer function matrix `Mi(λ)`; 

`info.freq` is the employed frequency used to check left invertibility 
(set to `missing` if no frequency-based left invertibility check was performed)

`info.HDesign` is a vector of design matrices `H`, where `H[i]` is the design matrix 
employed for the synthesis of the `i`-th component of the fault detection filter `Q`; 
`H[i]` is an empty matrix if no design matrix was involved.
   
`info.gammaopt0` is a vector whose `i`-th component is the optimal performance `gammaopt0i` for the `i`-th original problem (1); 
   
`info.gammaopt` is a vector whose `i`-th component is the optimal performance `gammaopti` for the `i`-th updated problem (3); 

`info.gammasub` is a vector whose `i`-th component is the suboptimal performance `gammasubi` in (2); 

`info.nonstandard` is set to `true` for a non-standard problem 
   (i.e., `Ge(λ)` has zeros on the boundary of the stability domain), and set to 
   `false` for a standard problem 
      (i.e., `Ge(λ)` has no zeros on the boundary of the stability domain). 

_Method:_ The synthesis Procedure AMMS from [1] is used 
to determine the component filters. The 
Procedure AMMS relies on the approximate model-matching synthesis method 
proposed in [2]. For more details on computational aspects see [3].  

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
              Springer Verlag, 2017; sec. 5.6.

[2] A. Varga, Integrated computational algorithm for solving 
    H_inf-optimal FDI problems. In Proc. of the IFAC World Congress, 
    Milano, Italy, pp. 10187–10192, 2011.

[3] A. Varga. Descriptor system techniques in solving H_2/H-Inf-optimal
    fault detection and isolation problems". In L. T. Biegler,  
    S. L. Campbell, and V. Mehrmann (Eds.), Control and Optimization 
    with Differential-Algebraic Constraints, vol. 23 of Advances in 
    Design and Control, pp. 105–125. SIAM, 2012. 
"""
function ammsyn(sysf::FDIModel{T1}, sysref::FDIFilterIF{T2}; 
                poles::Union{AbstractVector,Missing} = missing, 
                sdeg::Union{Real,Missing} = missing, smarg::Union{Real,Missing} = missing, 
                nullspace::Bool = true, mindeg::Bool = false, simple::Bool = false, regmin::Bool = true, 
                normalize::AbstractString = "infnorm", freq::Real = rand(), reltol::Real = 0.0001, H2syn::Bool = false, 
                tcond::Real = 1.e4, HDesign::Union{AbstractMatrix,Missing} = missing,
                offset::Real = sqrt(eps(float(real(T1)))), atol::Real = zero(float(real(T1))), 
                atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T1)))))*iszero(max(atol1,atol2)), 
                fast::Bool = true) where {T1,T2} 

   N = length(sysref.sys)
   DescriptorSystems.promote_Ts(sysf.sys.Ts, sysref.sys[1].Ts)  # check sampling time

   emptyHD = ismissing(HDesign)
   T = promote_type(T1,T2)
   
   HDesign = similar(Array{Matrix{T},1},N)
   Q = similar(sysref.sys,N) 
   R = similar(sysref.sys,N) 
   M = similar(sysref.sys,N) 
   gopt0 = similar(Vector{T},N)
   gopt = similar(Vector{T},N)
   gsub = similar(Vector{T},N)
   tcond1 = 1
   infodegs = 0
   inpu, inpd, inpf, inpw, inpa = 0,0,0,0,0
   nonstandard = true
   for i = 1:N
       Mri = FDFilterIF(sysref.sys[i], sysref.controls, sysref.disturbances, 
                        sysref.faults, sysref.noise, sysref.aux)
       Qi, Ri, infoi = ammsyn(sysf, Mri; sdeg, smarg, poles, nullspace, H2syn, reltol, 
                                         mindeg, simple, regmin, normalize, freq, 
                                         tcond, HDesign = emptyHD ? missing : HDesign[i], offset, 
                                         atol1, atol2, atol3, rtol, fast)
       Q[i] = Qi.sys
       R[i] = Ri.sys
       M[i] = infoi.M
       HDesign[i] = ismissing(infoi.HDesign) ? zeros(T,0,0) : infoi.HDesign
       tcond1 = max(tcond1, infoi.tcond)
       gopt0[i] = infoi.gammaopt0
       gopt[i] = infoi.gammaopt
       gsub[i] = infoi.gammasub
       i == 1 && (infodegs = infoi.degs; inpu = Ri.controls; inpd = Ri.disturbances; inpf = Ri.faults; inpw = Ri.noise; inpa = Ri.aux; nonstandard = infoi.nonstandard)
   end
   p = size(sysf.sys,1)
   mu = length(sysf.controls)
   info = (tcond = tcond1, degs = infodegs, M = M, freq = freq, HDesign = all(ismissing.(HDesign)) ? missing : HDesign,
           gammaopt0 = gopt0, gammaopt = gopt, gammasub = gsub, nonstandard = nonstandard)
   return FDIFilter(Q, p, mu), 
          FDIFilterIF(R, controls = inpu, disturbances = inpd, faults = inpf, noise = inpw, aux = inpa), info
end
