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

       y = Gu*u + Gd*d + Gf*f + Gw*w + Ga*aux,

with the Laplace- or Z-transformed plant outputs `y`, control inputs `u`, 
disturbance inputs `d`, fault inputs `f`, noise inputs `w` and auxiliary 
inputs `aux`, and with `Gu`, `Gd`, `Gf`, `Gw`, and `Ga` the corresponding 
transfer-function matrices.
The indices of control, disturbance, fault, noise and auxiliary inputs are contained in the associated integer vectors 
`sysf.controls`, `sysf.disturbances`, `sysf.faults`, `sysf.noise` and `sysf.aux`, respectively.

The fault detection filter object `Q`, contains in `Q.sys` the resulting filter 
in a standard state-space form, which generates the residual signal `r`. 
The corresponding input-output (implementation) form is

            r = Qy*y + Qu*u               

and the indices of output and control inputs are contained in the integer vectors 
`Q.outputs` and `Q.controls`, respectively.

The fault detection filter in internal form object `R`, contains `R.sys`, the resulting 
internal form of the filter 
in a standard state-space form, which generates the residual signal `r`, and corresponds to the 
input-output form

       r = Ru*u + Rd*d + Rf*f + Rw*w + Ra*aux ,

where 

       | Ru Rd Rf Rw Ra | = |Qy Qu]*| Gu Gd Gf Gw Ga |. 
                                    | I  0  0  0  0  |

The solution of the EFDP ensures that `Ru = 0`, `Rd = 0`, and `Rf` has all its columns nonzero. 
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

`FDFreq = freq` specifies a vector of real frequency values or a scalar real frquency value
for strong detectability checks (default: `FDFreq = missing`).

If `nullspace = true` (default), a minimal proper nullspace basis is used for the synthesis of the
fault detection filter. If `nullspace = false`, a full-order observer based nullspace basis is used. 
This option can be only used for a proper system without disturbance inputs. 

If `simple = true`, a simple proper nullspace basis is emplyed for synthesis. 
The orders of the basis vectors are provided in `info.deg`. 
If `simple = false` (default), then no simple basis is computed. 

To assess the stability of poles, a boundary offset  `β` can be specified via the keyword argument `offset = β`. 
Accordingly, for the stability of a continuous-time system all real parts of poles must be at most `-β`, 
while for the stability of a discrete-time system all moduli of poles must be at most `1-β`. 
The default value used for `β` is `sqrt(ϵ)`, where `ϵ` is the working machine precision. 

The stability domain `Cs` of poles is defined by 
the keyword argument for the stability margin `smarg = α`, as follows: 
for a continuous-time system, `Cs` is the set of complex numbers 
with real parts at most `α`, 
while for a discrete-time system, `Cs` is the set of complex numbers with 
moduli at most `α < 1` (i.e., the interior of a disc of radius `α` centered in the origin). 
If `smarg` is missing, then the employed default values are `α = -β` 
for a continuous-time system and `α = 1-β` for a discrete-time system, 
where `β` is the boundary offset specified by the keyword argument `offset = β`. 

The poles of the filters `Q` and `R` can be specified using the keyword arguments `poles` and `sdeg`,
where `poles = v` specifies a complex vector `v` containing a complex conjugate set  
of desired poles within the stability domain `Cs` (default: `poles = missing`)
and `sdeg = γ` is the prescribed stability degree for the poles (default: `γ = -0.05` for the real parts of poles for a continuous-time system and
`γ = 0.95` for the magnitudes of poles for a discrete-time system). 

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
polynomial nullspace basis of `G := [ Gu Gd; I 0]` (also the left Kronecker indices of `G`), if the 
state-space realization of `[Gu Gd]` is minimal;

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
                      tcond::Real = 1.e4, HDesign::Union{AbstractVecOrMat,Missing} = missing,
                      offset::Real = sqrt(eps(float(real(T)))), atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                      rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T)))))*iszero(max(atol1,atol2)), 
                      fast::Bool = true) where T
   Ts = sysf.sys.Ts                  
   disc = (Ts != 0);  # system type (continuous- or discete-time)
   if disc
      smax = 1-offset; smin = 0;
   else
      smax = -offset;  smin = -Inf;
   end
   
   # decode options
   
   # tolerance for rank tests
   
   tolmin = atol3 # tolerance for observability tests
   
   
   strongFD = !ismissing(FDfreq)
   strongFD && !isa(FDfreq,Vector) && (FDfreq = [FDfreq]) 
   strongFD && (lfreq = length(FDfreq))

      # set default stability degree
   if disc
      sdegdefault = 0.95;
   else
      sdegdefault = -0.05;
   end
   
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
   
   # set options for nullspace computation
   # if strongFD
   #    # set options for nullspace computation with stabilization
   #    opts_glnull = struct('tol',tol,'m2',mf+mw+maux,'simple',simple,'sdeg',sdegdefault);
   # else 
   #    opts_glnull = struct('tol',tol,'m2',mf+mw+maux,'simple',simple);
   # end
   # set options for LCF-based stabilization to be used for solvability checks
   # opts_glcf_default = struct('tol',tol,'tolmin',tolmin);
   # set options for LCF-based stabilization to be used for final synthesis
   # opts_glcf = struct('tol',tol,'tolmin',tolmin, ...
   #                    'sdeg',sdeg,'smarg',smarg,'poles',poles);
       
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
      # obtain QR = [ Q R ], where R = [ Rf Rw Raux] = Q*[Gf Gw Ga;0 0 0]
      QR, info1 = glnull(syse, m2; atol1, atol2, rtol, fast, sdeg = sdegNS, offset) 
      degs = info1.degs
      tcond1 = info1.tcond
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
   
   indf = (p+mu) .+ Vector(1:mf);   # input indices of Rf in QR
   if mf == 0
      # handle the case of no faults as a normal case
      S = falses(nvec,0);
   else
      if strongFD 
         S = fdisspec(Htemp*QR[:,indf], FDfreq; FDGainTol, atol1, atol2, rtol = 0, fast)[1]
         # check strong detectability conditions 
         for ii = 1:lfreq
            all(maximum(S[:,:,ii],dims=1)) || error("strong detection of all faults not feasible")
         end
      else
         # check weak detectability conditions 
         S = fditspec(Htemp*QR[:,indf]; atol1, atol2, rtol, FDtol)
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
                    for ii = 1:lfreq
                        Stest = fdisspec(glcf(QRfwtest[:,indf]; atol1, atol2, atol3, rtol)[1], FDfreq[ii]; 
                                     FDGainTol, atol1, atol2, rtol = 0, fast)[1]
                        all(maximum(Stest[:,:,1],dims=1)) || break
                    end
                 else
                    Stest = fditspec(QRfwtest[:,indf]; atol1, atol2, rtol, FDtol)
                 end
                 if all(maximum(Stest,dims=1)) 
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
          nout = nout+1
          if nout > nvec
             if itry > 5
                finish = true
                @warn "fault detectability not achieved with the chosen number of residuals"
             else
                itry = itry+1
                nout = nout-1
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
       QR = dss(al, el, bl, cl, dl; Ts)[1]
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
function efdbasesel(S::BitArray, degs::Vector{Int}, rdim::Int, nout::Int, simple::Bool)
   #   efdbasesel(S, degs, rdim, nout, simple) -> seli, selord
   #
   #   Select admissible basis vectors for solving the EFDP using the binary 
   #   structure matrix S corresponding to nvec basis vectors. 
   #   seli contains nout-touples of indices of basis vectors whose linear combination 
   #   is admissible, i.e. , S(seli[i],:) has all columns nonzero. This ensures that
   #   the exact fault detection problem (EFDP) is solvable by using 
   #   fault detection filters with rdim <= nout outputs. 
   #   If the associated nvec degrees contained in degs are provided, then
   #   selord[i] is the corresponding tentatively achievable least filter order.
   #   If simple = true, a simple basis is assumed, in which case, degs[i] is 
   #   also the order of the i-th basis vector. If simple = false, a minimum 
   #   rational basis is assumed. selord is empty if degs is empty. 

   #   Method: The selection approach is used in conjunction with the synthesis 
   #   Procedure EFD described in [1]. 

   # References:
   # [1] Varga A.
   #     Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017.


   ndim = ndims(S)
   if ndim == 3 
      nvec = size(S,1)
      n = size(S,3)
   else
      nvec, mf = size(S);   # numbers of vectors, number of faults  
      n = 1
      S = reshape(S,nvec,mf,1)
   end
   nd = length(degs)
   nodegs = (nd == 0)

   nodegs || length(degs) == nvec || error("the dimension of degs must must be equal to the number of rows of S")
   
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
           stest = stest & all(maximum(S[indv,:,k],dims=1))
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
