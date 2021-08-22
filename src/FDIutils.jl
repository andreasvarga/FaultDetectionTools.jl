"""
    fdIFeval(sysQ::FDFilter, sysf::FDIModel; minimal = false, atol, atol1 = atol, atol2 = atol, rtol, fast = true) -> sysR::FDFilterIF

Compute the internal form `sysR` of the fault detection filter `sysQ` applied to the synthesis model `sysf`. 
If `sysf` has the partitioned transfer function matrix `G(λ) = [ Gu(λ)  Gd(λ) Gf(λ) Gw(λ) Gv(λ) ]` in accordance with
the partitioned system inputs as `controls`, `disturbances`, `faults`, `noise` and `auxiliary` inputs, respectively,
and `Q(λ) = [ Qy(λ) Qu(λ) ]` is the partitioned transfer function matrix of the fault detection filter `sysQ` 
in accordance with the partitioned filter inputs as `outputs` and `controls`, then 
the transfer function matrix `R(λ)` of the resulting internal form `sysR` is given by
     
     R(λ) = | Qy(λ)  Qu(λ) | * | Gu(λ)  Gd(λ) Gf(λ) Gw(λ) Gv(λ) |
                               |  I       0     0     0     0   |

A minimal descriptor realization is computed if `minimal = true` and a possibly non-minimal 
realization is determined if `minimal = false` (default). 

The minimal realization computation relies on pencil manipulation algorithms which 
employ rank determinations based on either the use of 
rank revealing QR-decomposition with column pivoting, if `fast = true`, or the SVD-decomposition.
The rank decision based on the SVD-decomposition is generally more reliable, but the involved computational effort is higher.

If `(Ar-λEr,Br,Cr,Dr)` is the full order descriptor realization of `sysR.sys`, then 
the keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Ar`, `Br`, `Cr`, `Dr`, the absolute tolerance for the nonzero elements of `Er`,  
and the relative tolerance for the nonzero elements of `Ar`, `Br`, `Cr`, `Dr` and `Er`.
The default relative tolerance is `nϵ`, where `ϵ` is the working _machine epsilon_ 
and `n` is the order of the system `sysR`.  
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 

"""
function fdIFeval(Q::FDFilter, sysf::FDIModel; minimal::Bool = false, 
                  atol::Real = 0., atol1::Real = atol, atol2::Real = atol, rtol::Real = 0., fast::Bool = true)
   # compute the internal form as
   #   R = Q * [ Gu Gd Gf Gw Gaux]
   #           [ I  0  0  0  0   ]
   p, m = size(sysf.sys)
   mu = length(sysf.controls)
   (length(Q.outputs) == p && length(Q.controls) == mu) || error("filter Q is incompatible with the given system")
   # return fdRset(gminreal(Q.sys*[sysf.sys;  I zeros(eltype(sysf.sys),mu,m-mu)]; atol1, atol2, rtol); 
   #               sysf.controls, sysf.disturbances, sysf.faults, sysf.noise, sysf.aux) 
   Rsys = Q.sys*[sysf.sys;  I zeros(eltype(sysf.sys),mu,m-mu)]
   return FDFilterIF(minimal ? gminreal(Rsys; atol1, atol2, rtol) : Rsys, 
                 sysf.controls, sysf.disturbances, sysf.faults, sysf.noise, sysf.aux) 
   
end
"""
     fdhinfminus(sys,freq) -> (β, ind, fr)

Compute for a stable descriptor system `sys = (A-λE,B,C,D)` the `H∞-` index `β` of its
transfer function matrix `G(λ)`. If `freq = missing` (default), then `β` is the 
minimum `H∞-norm` of the columns of `G`, `ind` is the index of the minimum-norm column and `fr` is 
the frequency where the minimum `H∞-norm` of the columns is achieved. If `freq` is a real value or 
a real vector of frequency values, then `β` is the minimum of the 2-norms of the columns of the 
frequency responses of `G` evaluated for all values contained in `freq`, `ind` is the index of column 
for which the minimum is achieved and `fr` is the corresponding frequency. 
"""
function fdhinfminus(sys::DescriptorStateSpace{T}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
   p, m = size(sys)

   m == 0 && (return T[], Int[], T[])
   
   isstable(sys) || error("the system is unstable")
   β = Inf; ind = 0; f = -1;
   if ismissing(freq) 
      # evaluate β as the minimum of H-infinity norms of the columns of G
      for j = 1:m
         temp, fr = ghinfnorm(sys[:,j])
         β <= temp || (β = temp; ind = j; f = fr)
      end
   else
      # evaluate β as the minimum of the norms of columns of the frequency 
      # responses of G evaluated over all frequencies contained in FREQ 
      if !isa(freq, Vector) 
         # use evalfr if only one frequency is present
         H = evalfr(sys; fval = freq) 
         for j = 1:m
             temp = norm(view(H,:,j))
             β <= temp || (β = temp; ind = j)
         end
         f = freq
      else
         T1 = T <: BlasFloat ? T : promote_type(Float64,T) 
         ONE = one(T1)
         a, e, b, c, d = dssdata(T1,sys)
         Ts = abs(sys.Ts)
         disc = !iszero(Ts)
         sw = disc ? im*Ts : im
         desc = !(e == I)
         # Determine the complex Hessenberg-form system to be used for efficient
         # frequency response computation.
         ac, ec, bc, cc, dc = chess(a, e, b, c, d)
         H = similar(dc, eltype(dc), p, m)
         bct = similar(bc) 
         for i = 1:length(freq)
             if isinf(freq[i])
                # exceptional call to evalfr
                H = evalfr(sys, fval=Inf)
             else
                copyto!(H, dc)
                copyto!(bct, bc)
                w = disc ? -exp(sw*freq[i]) : -sw*freq[i]
                desc ? ldiv!(UpperHessenberg(ac+w*ec),bct) : ldiv!(ac,bct,shift = w)
                mul!(H, cc, bct, -ONE, ONE)
             end
             for j = 1:m
                 temp = norm(view(H,:,j))
                 β <= temp || (β = temp; ind = j; f = freq[i])
             end
         end
      end
   end
   return β, ind, f
end
"""
     fdhinfmax(sys,freq) -> (γ, ind, fr)

Compute for a descriptor system `sys = (A-λE,B,C,D)`, `γ` - the maximum norm of the columns of its
transfer function matrix `G(λ)`. If `freq = missing` (default), then `γ` is the 
maximum `H∞-norm` of the columns of `G`, `ind` is the index of the maximum-norm column and `fr` is 
the frequency where the maximum `H∞-norm` of the columns is achieved. If `freq` is a real value or 
a real vector of frequency values, then `γ` is the maximum of the 2-norms of the columns of the 
frequency responses of `G` evaluated for all values contained in `freq`, `ind` is the index of column 
for which the maximum is achieved and `fr` is the corresponding frequency. 
"""
function fdhinfmax(sys::DescriptorStateSpace{T}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
   p, m = size(sys)

   m == 0 && (return T[], Int[], T[])
   
   γ = 0; ind = 0; f = -1;
   if ismissing(freq) 
      # evaluate γ as the maximum of H-infinity norms of the columns of G
      for j = 1:m
         temp, fr = ghinfnorm(sys[:,j])
         γ >= temp || (γ = temp; ind = j; f = fr)
      end
   else
      # evaluate γ as the maximum of the norms of columns of the frequency 
      # responses of G evaluated over all frequencies contained in FREQ 
      if !isa(freq, Vector) 
         # use evalfr if only one frequency is present
         H = evalfr(sys; fval = freq) 
         for j = 1:m
             temp = norm(view(H,:,j))
             γ >= temp || (γ = temp; ind = j)
         end
         f = freq
      else
         T1 = T <: BlasFloat ? T : promote_type(Float64,T) 
         ONE = one(T1)
         a, e, b, c, d = dssdata(T1,sys)
         Ts = abs(sys.Ts)
         disc = !iszero(Ts)
         sw = disc ? im*Ts : im
         desc = !(e == I)
         # Determine the complex Hessenberg-form system to be used for efficient
         # frequency response computation.
         ac, ec, bc, cc, dc = chess(a, e, b, c, d)
         H = similar(dc, eltype(dc), p, m)
         bct = similar(bc) 
         for i = 1:length(freq)
             if isinf(freq[i])
                # exceptional call to evalfr
                H = evalfr(sys, fval=Inf)
             else
                copyto!(H, dc)
                copyto!(bct, bc)
                w = disc ? -exp(sw*freq[i]) : -sw*freq[i]
                desc ? ldiv!(UpperHessenberg(ac+w*ec),bct) : ldiv!(ac,bct,shift = w)
                mul!(H, cc, bct, -ONE, ONE)
             end
             for j = 1:m
                 temp = norm(view(H,:,j))
                 γ >= temp || (γ = temp; ind = j; f = freq[i])
             end
         end
      end
   end
   return γ, ind, f
end
function chess(a, e, b, c, d)
   # reduce the descriptor system (A-λE,B,C,D) to an equivalent complex descriptor system 
   # (Ac-λEc,Bc,Cc,Dc) such that Ac-λEc is in an upper Hessenberg form 
   # Note: This function is used in the context of efficient frequence response computations. 
   T = eltype(a)
   desc = !(e == I)
   if desc
       # Reduce (A,E) to (generalized) upper-Hessenberg form for
       # fast frequency response computation 
       at = copy(a)
       et = copy(e)
       bt = copy(b)
       # first reduce E to upper triangular form 
       _, tau = LinearAlgebra.LAPACK.geqrf!(et)
       T <: Complex ? tran = 'C' : tran = 'T'
       LinearAlgebra.LAPACK.ormqr!('L', tran, et, tau, at)
       LinearAlgebra.LAPACK.ormqr!('L', tran, et, tau, bt)
       # reduce A to Hessenberg form and keep E upper triangular
       _, _, Q, Z = MatrixPencils.gghrd!('I', 'I',  1, size(at,1), at, et, similar(at),similar(et))
       if T <: Complex 
          bc = Q'*bt; cc = c*Z; ac = at; ec = et; dc = d;
       else
          bc = complex(Q'*bt); cc = complex(c*Z); ac = complex(at); ec = complex(et); dc = complex(d);
       end
   else
       # Reduce A to upper Hessenberg form for
       # fast frequency response computation 
       Ha = hessenberg(a)
       ac = Ha.H
       if T <: Complex 
          bc = Ha.Q'*b; cc = c*Ha.Q; dc = d; 
       else
          bc = complex(Ha.Q'*b); cc = complex(c*Ha.Q); dc = complex(d);
       end
       ec = e
   end
   return ac, ec, bc, cc, dc
end