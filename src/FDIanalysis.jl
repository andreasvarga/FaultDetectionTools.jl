"""
    fdichkspec(sysf::FDIModel, SFDI::BitMatrix; sdeg, FDtol, FDGainTol, FDfreq, 
                 atol, atol1, atol2, atol3, rtol, fast = true, minimal = false) -> (rdims, orders, leastorders)

Check for a given synthesis model `sysf::FDIModel` the feasibility of a set of fault detection and isolation specifications `SFDI`. 
If `SFDI` has `N` rows (i.e., contains `N` specifications), then the `N`-dimensional integer vectors `rdims`, `orders`, `leastorders` 
are returned and contain information related to the synthesis of FDI filters to achieve the feasible specifications. 
For the `i`-th specification contained in `SFDI[i,:]`, `rdims[i]` contains the number of residual 
outputs of a minimal nullspace basis based FDI filter which can be used to achieve this specification. 
If `rdims[i] = 0`, then the `i`-th specification is not feasible. For a feasible `i`-th specification, `orders[i]` 
contains the order of the minimal nullspace basis based FDI filter which can be used to achieve this specification. 
If the `i`-th specification is not feasible, then `orders[i]` is set to `-1`.
If `minimal = true`, `leastorders[i]` contains the least achievable order for a scalar output FDI filter which can be used  
to achieve the `i`-th specification. If `minimal = false` or if the `i`-th specification is not 
feasible, then `leastorders[i]` is set to `-1`.

`FDFreq = freq` specifies a vector of real frequency values or a scalar real frquency value
for strong detectability checks (default: `FDFreq = missing`).

`FDtol = tol1` specifies the threshold `tol1` for assessing weak specifications
                      (see also function [`fditspec`](@ref)) (default: `tol1 = 0.0001`).

`FDGainTol = tol2` specifies the threshold `tol2` for assessing strong specifications,  
i.e., the threshold for nonzero frequency responce gains for all frequency values specified in
`freq` (see also function [`fdisspec`](@ref)) (default: `tol2 = 0.01`). 

The keyword argument `sdeg = β` specifies a prescribed stability degree `β` for the poles of the internally 
generated candidate filters, such that the real parts of
filters poles must be less than or equal to `β`, in the continuous-time case, and 
the magnitudes of filter poles must be less than or
equal to `β`, in the discrete-time case. If `sdeg = missing` then no then no stabilization is performed if and `FDFreq = missing`.
If `sdeg = missing` and `FDFreq = freq`, then the fllowing default values are employed : `β = -0.05`, in continuous-time case, and  `β = 0.95`, 
in discrete-time case. 

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

_Method:_ The nullspace method of [1] is successively employed to 
determine FDI filters as minimal left nullspace bases which solve 
suitably formulated fault detection problems. 

_References:_

[1] A. Varga, On computing nullspace bases – a fault detection perspective. 
      Proc. IFAC 2008 World Congress, Seoul, Korea, pages 6295–6300, 2008.
"""
function fdichkspec(sysf::FDIModel{T}, SFDI::Union{BitMatrix,Matrix{Bool}}; sdeg::Union{Real,Missing} = missing, minimal::Bool = false,  
                    FDtol::Real = 0.0001, FDGainTol::Real = 0.01, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                    atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                    rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T)))))*iszero(max(atol1,atol2)), 
                    fast::Bool = true) where T

   mu = length(sysf.controls)  
   md = length(sysf.disturbances)  
   mud = mu+md
   mf = length(sysf.faults)  
   m = mu+md+mf
   strongFD = !ismissing(FDfreq)
   strongFD && !isa(FDfreq,Vector) && (FDfreq = [FDfreq]) 
   strongFD && ismissing(sdeg) && (sdeg = sysf.sys.Ts != 0 ? 0.95 : -0.05)   

   size(SFDI,2) == mf || error("Structure matrix incompatible with the number of faults")
  
   N = size(SFDI,1)
   p = size(sysf.sys,1)
   strongFD && (lfreq = length(FDfreq))     # number of frequency values
   rdims = fill(0,N)
   orders = fill(-1,N)
   leastorders = fill(-1,N)

   N == 0 && (return rdims, orders, leastorders)
   
   for i = 1:N
       # set input groups for first design
       inddf = Vector(1:mf)[SFDI[i,:] .== false] 
       indff = Vector(1:mf)[SFDI[i,:] .== true] 
       isempty(indff) && (rdims[i] = 1; orders[i] = 0; leastorders[i] = minimal ? 0 : -1; continue)
       syse = [sysf.sys[:,[1:mud; mud .+ inddf; mud .+ indff]]; eye(mu,m)]
       m2f = length(indff)
       #
       # compute a left nullspace basis Q of G = [Gu Gd Gdf; I 0 0] = 0 
       # obtain QR = [ Q Rf ], where Rf = Q*[Gf;0] 
       QR, info = glnull(syse, m2f; sdeg, atol1, atol2, rtol, fast) 
       nvec = size(QR,1);          # number of basis vectors
       # check solvability conditions
       nvec == 0 && continue
      
       nq = order(QR)               # order of minimal basis
       degs = reverse(info.degs)    # degrees of a minimal polynomial basis
       m2f == 0 && (rdims[i] = nvec; orders[i] = nq; leastorders[i] = minimum(degs); continue)
      
       indf = (p+mu) .+ Vector(1:m2f)
       feasible = true;
       if strongFD 
          t = gir(QR[:,indf]; atol = atol1)
          S = fdisspec_(t, FDfreq; stabilize = true, block = !minimal, FDGainTol, atol1, atol2, rtol = 0, fast)[1]
          # check strong detectability conditions 
          for ii = 1:lfreq
             all(maximum(S[:,:,ii],dims=1)) || (feasible = false; break)
          end
       else
          # check weak detectability conditions 
          S = fditspec_(QR[:,indf]; block = !minimal, atol1, atol2, rtol, FDtol)
          all(maximum(S,dims=1)) || (feasible = false)
       end
       if feasible
          rdims[i] = nvec
          orders[i] = nq 
          if minimal
             if nvec > 1
                 finish = false    # set termination flag
                 nout = 1          # initialize number of selected basis vectors
                 #QR = xperm(QR,nq:-1:1);  # permute states to speedup glmcover1
                 while !finish     
                    # choose nout basis vectors, which potentially lead to a least order
                    # filter with rdim outputs:
                    # basesel[i,:] contains the indices of candidate basis vectors;
                    # ordsel[i]    contains the presumably achievable least orders
                    basesel, ordsel = efdbasesel(S, degs, 1, nout, false) 
                    #
                    # update the synthesis using the selections of candidate vector(s),
                    # starting with the least (potentially) achievable order
                    for ibas = 1:size(basesel,1)
                        baseind = basesel[ibas] # indices of current basis selection
                        ip = [baseind; setdiff(1:nvec,baseind)][:]
                        if nout == 1
                           QRfwtest = glmcover1(QR[ip,:], 1; atol1, atol2, rtol)[1]
                           if !isempty(ordsel) && (order(QRfwtest) != ordsel[ibas])
                              @warn "fdichkspec: expected reduced order not achieved"
                           end
                        else  
                           # build a linear combination of the first nout vectors 
                           nrest = nvec - nout
                           h = [ rand(1,nout) zeros(T, 1, nrest) ] 
                           QRfwtest = glmcover1([h; zeros(T,nrest,nout) eye(nrest)]*QR[ip,:], 1; atol1, atol2, rtol)[1]
                        end
                        # check complete fault detectability of the current design; 
                        # dismiss design if check fails
                        t = gminreal(QRfwtest[:,indf]; atol1, atol2, rtol, fast)
                        if strongFD 
                           Stest = fdisspec_(t, FDfreq; stabilize = true, block = true,
                                             FDGainTol, atol1, atol2, atol3, rtol = 0, fast)[1]
                        else
                           Stest = fditspec_(t; atol1, atol2, rtol, FDtol, block = true)
                        end
                        if all(Stest) 
                           finish = true
                           leastorders[i] = order(QRfwtest)
                           break
                        end
                     end
                     nout += 1
                     nout > nvec && (finish = true; leastorders[i] = nq)  # least order is 
                  end
             else
                  leastorders[i] = nq  #least order for nvec = 1
             end
          end
       end
   end
   return rdims, orders, leastorders   
   # end FDICHKSPEC  
end
"""
    S = fdigenspec(sysf::FDIModel; sdeg, FDtol, FDGainTol, FDfreq, atol, atol1, atol2, atol3, rtol, fast = true) 

Generate all achievable specifications `S` for a given synthesis model `sysf` with additive faults. 
Each row of the resulting binary matrix `S` contains a nonzero specification (or fault signature) which can
be achieved using a linear fault detection filter (e.g., as obtainable with the help of function [`efdisyn`](@ref)).

`FDFreq = freq` specifies a vector of real frequency values or a scalar real frquency value
for strong detectability checks (default: `FDFreq = missing`).

`FDtol = tol1` specifies the threshold `tol1` for assessing weak specifications
                      (see also function [`fditspec`](@ref)) (default: `tol1 = 0.0001`).

`FDGainTol = tol2` specifies the threshold `tol2` for assessing strong specifications,  
i.e., the threshold for nonzero frequency responce gains for all frequency values specified in
`freq` (see also function [`fdisspec`](@ref)) (default: `tol2 = 0.01`). 

The keyword argument `sdeg = β` specifies a prescribed stability degree `β` for the poles of the internally 
generated candidate filters, such that the real parts of
filters poles must be less than or equal to `β`, in the continuous-time case, and 
the magnitudes of filter poles must be less than or
equal to `β`, in the discrete-time case. If `sdeg = missing` then no stabilization is performed if `FDFreq = missing`.
If `sdeg = missing` and `FDFreq = freq`, then the following default values are employed : `β = -0.05`, in continuous-time case, and  `β = 0.95`, 
in discrete-time case. 

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

_Method:_ The Procedure GENSPEC from [1] is implemented. 
The nullspace method [2] is recursively employed to generate candidate fault detection and isolation filters,
whose internal forms provide the structure matrices corresponding to the achievable weak specifications, if `freq = missing`, 
or strong specifications for the frequencies conatined in `freq`. The generation method is also described in [3].

_References:_

[1] A. Varga, Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
     Springer Verlag, 2017; sec. 5.4.
   
[2] A. Varga, On computing nullspace bases – a fault detection perspective. 
      Proc. IFAC 2008 World Congress, Seoul, Korea, pages 6295–6300, 2008.

[3] A. Varga, On computing achievable fault signatures. Proc. SAFEPROCESS'2009, Barcelona, Spain. 
"""
function fdigenspec(sysf::FDIModel{T}; sdeg::Union{Real,Missing} = missing, 
                    FDtol::Real = 0.0001, FDGainTol::Real = 0.01, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                    atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                    rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T)))))*iszero(max(atol1,atol2)), 
                    fast::Bool = true) where T

   mu = length(sysf.controls)  
   md = length(sysf.disturbances)  
   mf = length(sysf.faults)  
   m = mu+md+mf
   strongFD = !ismissing(FDfreq)
   strongFD && !isa(FDfreq,Vector) && (FDfreq = [FDfreq]) 
   strongFD && ismissing(sdeg) && (sdeg = sysf.sys.Ts != 0 ? 0.95 : -0.05)   
   
   return  genspec_in([sysf.sys[:,1:m]; eye(T,mu,m)]; m1 = mu+md, sdeg, FDtol, FDGainTol, FDfreq, 
                       atol1, atol2, atol3, rtol, fast) 
   
   # end FDIGENSPEC  
end
function genspec_in(sys::DescriptorStateSpace{T}; inp1::Vector{Int} = Int[], inp2::Vector{Int} = Int[], m1::Int = 0,
    sdeg::Union{Real,Missing} = missing, FDtol::Real = 0.0001, FDGainTol::Real = 0.01, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
    atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
    rtol::Real = ((size(sysf.sys.A,1)+1)*eps(real(float(one(T)))))*iszero(max(atol1,atol2)), 
    fast::Bool = true) where T
#  Recursive function used to generate the achievable specifications.
#
#     S = genspec_in(sys::DescriptorStateSpace; inp1, inp2, md, FDfreq, FDtol, FDGainTol, atol1, atol2, atol3, rtol, fast) 
#
#  determines, in the rows of the binary matrix S, all achievable fault detection specifications for the
#  partitioned state-space system sys = [sys[:,inp1] sys[:,inp2]], with the
#  corresponding input-output representation of the form,
#
#                y = Gd*d + Gf*f ,
#
#  where y is the output, d and f are the disturbance and fault
#  inputs, respectively, and Gd and Gf are the transfer-function
#  matrices of sys[:,inp1] and sys[:,inp2], respectively.
#  If inp1 and inp2 are both empty (e.g., on the first call), then 
#  the value of m1 specifies the number of disturbance inputs and 
#  inp1 and inp2 are set to inp1 = 1:m1 and inp2 = m1+1:m, 
#  where m is the number of inputs of sys.
#  Each row of the resulting S contains the an achievable specification,
#  obtainable by using a certain scalar output fault detection and isolation filter
#  with the input-output form
#
#                r = Q*y ,
#
#  which ensures that all disturbance inputs d are decoupled
#  from the residual r (i.e., Q*Gd = 0) and the residual r is
#  sensitive to a subset of fault inputs (i.e., Q*Gf != 0).
#  If FDfreq = missing, then the i-th row S[i,:] is the weak structure matrix of Q*Gf, 
#  whose element S[i,j] = true if Q*Gf[:,j] is nonzero and S[i,j] = false if Q*Gf[:,j] = 0.
#  The check for nonzero Q*Gf[:,j] is performed by using the function
#  fditspec to evaluate the corresponding weak structure matrix. 
#  If FDfreq = freq specifies the frequency values in freq, then
#  S[i,j] = true if |Q*Gf[:,j]| is larger than or equal to the threshold FDGainTol
#  for all specified frequency values and S[i,j] = false otherwise. 
#  For checking purposes, the function fdisspec is used to
#  evaluate the corresponding strong structure matrix. 
#  For the description of keyword parameters atol1, atol2, atol3, rtol and fast
#  see the documentation of the function fdispecgen. 

   # set persistent variables used for recursive calls: level, ksave, nf
   global level 
   global ksave 
   global nf
   
   # determine input information to define SYS1 and SYS2 (relevant for recursive usage)
   md = isempty(inp1) ? 0 : length(inp1)
   mf = isempty(inp2) ? 0 : length(inp2)
          
   # initialize recursion level
   md == 0 && mf == 0 && (level = 0)   # iteration level is set to zero on first call
   
   if level == 0
      # at first call the number of disturbance inputs are specified via m1
      # and the rest are fault inputs
      md = m1
      mf = size(sys,2) - md  # the number of fault inputs
    
      # finish if there are no faults
      mf == 0 && (return falses(0,0))
   
      # set nf and ksave
      nf = ismissing(FDfreq) ? 0 : length(FDfreq);
      ksave = 1;  # parameter used to reduce the total number of iterations    
   end
   
   # compute SYSN*SYS2, where SYSN is a left nullspace of SYS1 
   md > 0 && (sys = glnull(sys, mf; sdeg, atol1, atol2, rtol, fast)[1][:,end-mf+1:end])
   # compute weak specification for the current level
   S = fditspec_(sys; atol1, atol2, rtol, FDtol) 

   if nf > 0
      # determine strong structure matrix
      # enforce that strong specifications are a subset of weak specifications
      #S2 = fdisspec_(sys, FDfreq; stabilize = true, atol1, atol2, rtol, FDGainTol)[1]; 
      t = gminreal(sys; atol1, atol2, rtol)
      Sc = maximum(S,dims=1) # include the cumulated specification
      S2 = fdisspec_(t, FDfreq; stabilize = true, atol1, atol2, rtol, FDGainTol)[1]; 
      ns = size(S,1)
      tv = trues(ns)
      tvc = true
      for i = 1:nf
          for j = 1:ns
              tv[j] && (tv[j] = isequal(view(S,j,:),view(S2,j,:,i)))
          end
          tvc && (tvc = isequal(Sc,maximum(view(S2,:,:,i),dims=1)))
      end
      S = S[tv,:]
      tvc && (S = [S; Sc] )
   end  
   # add cummulated specification 
   size(S,1) > 1 && (S = [S; maximum(S,dims=1)])

   # exit level if SYS has only one output, or if S = [] or S = 0
   if size(sys,1) == 1 || isempty(S) || all(.!S)   
      level > 0 && (level -= 1)
      # eliminate duplicate or zero specifications
      size(S,1) > 1 && (S = unique(S,dims=1); S = S[any(S,dims=2)[:],:])
      return S
   end
   #   loop over all fault inputs
   for k1 = 1:mf
      k1 < ksave && continue
      level += 1
      # redefine disturbance and fault inputs    
      # inp1 = 1; 
      # sys.InputGroup.inp2 = 1:mf-1; 
      inp1 = collect(1:1)
      inp2 = collect(1:mf-1)
      ksave = k1
      #s1 = genspec_in(sys(:,[k1, 1:k1-1,k1+1:mf])); 
      s1 = genspec_in(sys[:,[[k1]; collect(1:k1-1); collect(k1+1:mf)]]; 
                      inp1, inp2, sdeg, FDtol, FDGainTol, FDfreq, atol1, atol2, atol3, rtol, fast)
      ksave = k1
      #  add specifications from the current level
      S = [S; s1[:,1:k1-1] falses(size(s1,1),1) s1[:,k1:mf-1]]
   end 
      
   # eliminate duplicate or zero specifications
   size(S,1) > 1 && (S = [S; maximum(S,dims=1)])
   S = unique(S,dims=1)
   S = S[any(S,dims=2)[:],:]
   level -= 1
   return S
       
   # end GENSPEC_IN  
end
fdichkspec(sysf::DescriptorStateSpace{T},SFDI::Union{BitMatrix,Matrix{Bool}}; kwargs...) where T = fdichkspec(fdimodset(sysf,faults = 1:size(sysf,2)),SFDI; kwargs...)
fdigenspec(sysf::DescriptorStateSpace{T}; kwargs...) where T = fdigenspec(fdimodset(sysf,faults = 1:size(sysf,2)); kwargs...)


