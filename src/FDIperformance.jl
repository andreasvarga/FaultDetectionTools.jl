"""
    S = fditspec(sysrf; FDfreq = missing, block = false, FDtol, FDStol, 
                        atol = 0, atol1 = atol, atol2 = atol, rtol, fast = true) 

Compute the weak or strong binary structure matrix `S` of the transfer function matrix of a 
linear time-invariant system `sysrf` 
(typically representing the transfer channel from the fault inputs to residuals).
`sysrf` is either a descriptor system representation `sysrf = (Af-lambda*Ef,Bf,Cf,Df)` 
with a  `q x mf` transfer function matrix `Rf(λ)` or
can be a fault detection filter internal form object `sysrf::FDFilterIF`, in which case 
only the fault inputs channel `sysrf.sys[:,sysrf.faults] := (Af-lambda*Ef,Bf,Cf,Df)` is selected. 

If `FDfreq = missing` (default), then `S` contains the weak structure matrix of `Rf(λ)`. 
For `block = false`, `S` is determined as a `q x mf` 
binary matrix (BitMatrix), whose `(i,j)`-th element is `S[i,j] = 1`, if the `(i,j)`-th element of `Rf(λ)` is 
nonzero, and otherwise, `S[i,j] = 0`. 
For `block = true`, `S` is determined as a `1 x mf` binary matrix, whose `(1,j)`-th element is `S[1,j] = 1`, 
if the `j`-th column of `Rf(λ)` is nonzero, and otherwise, `S[1,j] = 0`. 

If `FDfreq = freq` specifies a vector `freq` of `nf` real frequencies 
which characterize the classes of persistent fault signals, then `S` contains the strong 
structure matrix of `Rf(λ)` with respect to a set of `nf` complex frequencies `Ω`, defined as follows: 
if `f` is a real frequency in `freq`, then the corresponding complex frequency in `Ω` 
is `λ := im*f`, for a continuous-time system `sysrf`,
or `λ := exp(im*f*abs(Ts))`, for a discrete-time system `sysfr` with sampling-time `Ts`.  
For `block = false`, `S` is determined as a `q x mf x nf` 
binary matrix, whose `(i,j,k)`-th element is `S[i,j,k] = 1`, if the `(i,j)`-th element of `Rf(λ)` is 
nonzero for the `k`-th frequency in `freq`, and otherwise, `S[i,j,k] = 0`. 
For `block = true`, `S` is determined as a `1 x mf x nf` binary matrix, whose `(1,j,k)`-th element is `S[1,j,k] = 1`, 
if the `j`-th column of `Rf(λ)` is nonzero for the `k`-th frequency in `FDfreq`, and otherwise, `S[1,j,k] = 0`. 

The keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Af`, `Bf`, `Cf`, `Df`, the absolute tolerance for the nonzero elements of `Ef`,  
and the relative tolerance for the nonzero elements of `Af`, `Bf`, `Cf`, `Df` and `Ef`. 
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 

`FDtol = tol1` specifies an absolute threshold `tol1` for the magnitudes of nonzero elements in the system matrices 
`Bf` and `Df` and is used to determine the weak structure matrix. 
Its default value is `tol1 = 0.0001*max(1, norm(Bf,1), norm(Df,1))`. 

`FDStol = tol2` specifies an absolute  threshold `tol2` for the magnitudes of nonzero elements in the system matrices 
`Af`, `Ef`, `Bf`, `Cf` and `Df` and is used to determine the strong structure matrix. 
Its default value is 
`tol2 = epsm*max(1, norm(Ef,1), norm(Af,1), norm(Bf,1), norm(Cf,Inf), norm(Df,1)))`, 
where `epsm` is the working machine precision.

_Method:_ For the definition of the structure matrix, see [1]. For the
determination of the weak structure matrix, minimal realizations
are determined for each column of `Rf(λ)` if `block = true` or for 
each element of `Rf(λ)` if `block = false` and the nonzero columns or 
elements in each column are identified  (see Corollary 7.1 of [1]).
For the determination of the strong structure matrix, minimal realizations
are determined for each column of `Rf(λ)` if `block = true` or for 
each element of `Rf(λ)` if `block = false` and  the full rank of the
corresponding system matrix is checked for all frequencies in `FDfreq`
(see Corollary 7.2 in [1]) (i.e., the lack of zeros in all frequencies).

_References:_

[1] Varga A. Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017; sec.3.4.
"""
function fditspec(sysrf::DescriptorStateSpace{T}; FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                  block::Bool = false, FDtol::Real = 0., FDStol::Real = 0., 
                  atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, 
                  rtol::Real =  (size(sysrf.A,1)+1)*eps(float(one(real(T))))*iszero(max(atol1,atol2)), fast::Bool = true) where T
   p, mf = size(sysrf) 
   n = order(sysrf)
   q = block ? 1 : p  # number of rows of S
   Sstrong = !ismissing(FDfreq)
   if Sstrong
      isa(FDfreq,Vector) || (FDfreq = [FDfreq]) 
      lfreq = length(FDfreq);
      w = im*FDfreq;                   # w = j*freq
      Ts = abs(sysrf.Ts);
      Ts > 0 && ( w = exp(Ts*w))     # w = exp(j*Ts*freq)
      S = falses(q,mf,lfreq)
   else
      S = falses(q,mf)
   end   
   a, e, b, c, d = dssdata(sysrf)
   standard = isequal(e,I)
   if n == 0
      S1 =  abs.(d) .> FDtol 
      S = block ? maximum(S1, dims=1) : S1
      return Sstrong ? repeat(S,1,1,lfreq) : S
   end
   if Sstrong
      #FDStol <= 0. && (FDStol = 0.0001*max(1., norm(a,1), norm(b,1), norm(c,Inf), norm(d,1), standard ? 0 : norm(e,1))) 
      FDStol <= 0. && (FDStol = eps(max(1., norm(a,1), norm(b,1), norm(c,Inf), norm(d,1), standard ? 0 : norm(e,1)))) 
      S = trues(q, mf, lfreq)
      # employ structural analysis to compute weak/strong structure matrix  
      for j = 1:mf 
         # elliminate uncontrollable eigenvalues for the i-th column of B
         a1, e1, b1, c1, d1 = lsminreal(a, e, view(b,:,j), c, view(d,:,j); obs = false, noseig = false, atol1, atol2, rtol, fast)
         if block
            for k = 1:lfreq
                # check if freq(k) is a zero of the j-th column 
                s = isinf(FDfreq[k]) ? svdvals!([e1 b1; c1 d1]) : svdvals!([a1-w[k]*e1 b1; c1 d1])
                S[1, j, k] = (s[end] > FDStol)
            end
         else
            for i = 1:p 
                # elliminate unobservable and non-dynamic eigenvalues for the i-th row of C
                a2, e2, b2, c2, d2 = lsminreal(a1, e1, b1, view(c1,i:i,:),view(d1,i:i,:); contr = false, noseig = true, atol1, atol2, rtol, fast)
                for k = 1:lfreq
                    # check if freq(k) is a zero of the (i,j)-th element 
                    s = isinf(FDfreq[k]) ? svdvals!([e2 b2; c2 d2]) : svdvals!([a2-w[k]*e2 b2; c2 d2])
                    S[i,j,k] = s[end] > FDStol
                end
            end
         end
      end
   else
      FDtol <= 0. && (FDtol = 0.0001*max(1., norm(sysrf.B,1), norm(sysrf.D,1)))
      standard && all(abs.(d) .> FDtol) && (return trues(q,mf))
      S = falses(q, mf)
      if block
         for j = 1:mf
             # compute minimal realization of the j-th column of Rf
             _, _, b1, _, d1 = lsminreal(a, e, view(b,:,j:j), c, view(d,:,j:j); atol1, atol2, rtol, fast)
             S[1,j] |=  (any(abs.(view(d1,:,1)) .> FDtol) || any(abs.(view(b1,:,1)) .> FDtol))
         end
      else
         for j = 1:mf 
             # elliminate uncontrollable eigenvalues for the j-th column of B
             a1, e1, b1, c1, d1 = lsminreal(a, e, view(b,:,j), c, view(d,:,j); obs = false, noseig = false, atol1, atol2, rtol, fast)
             for i = 1:p 
                # elliminate unobservable and non-dynamic eigenvalues for the (i,j)-th element
                _, _, b2, _, d2 = lsminreal(a1, e1, b1, view(c1,i:i,:),view(d1,i:i,:); contr = false, noseig = true, atol1, atol2, rtol, fast)
                S[i,j] |=  (abs(d2[1,1]) > FDtol || any(abs.(view(b2,:,1)) .> FDtol))
             end
         end
      end
   end
   return S
end
fditspec(sysr::FDFilterIF{T}; kwargs...) where T = fditspec(sysr.sys[:,sysr.faults]; kwargs...) 
"""
     fdisspec(sysrf, freq; block = false, FDGainTol = 0.01, 
                     atol = 0, atol1 = 0, atol2 = 0, rtol = 0, fast = true) -> (S, gains)

Compute for the linear time-invariant system `sysrf` and 
a given set of real frequencies `freq`,  
the strong binary structure matrix `S` of the transfer function matrix of `sysrf` 
and the corresponding frequency response gains `gains`. 
`sysrf` typically represents the transfer channel from the fault inputs to residuals
and is either a descriptor system representation `sysrf = (Af-lambda*Ef,Bf,Cf,Df)` 
with a  `q x mf` transfer function matrix `Rf(λ)` or
is a fault detection filter internal form object `sysrf::FDFilterIF`, in which case 
only the fault inputs channel `sysrf.sys[:,sysrf.faults] := (Af-lambda*Ef,Bf,Cf,Df)` is selected. 

`freq` must contain a real frequency value or a vector of `nf` real frequencies 
which characterize the classes of persistent fault signals 
(default: `freq = 0`, i.e., characterizing constant faults).
`S` contains the strong 
structure matrix of `Rf(λ)` with respect to a set of `nf` complex frequencies `Ω`, defined as follows: 
if `f` is a real frequency in `freq`, then the corresponding complex frequency in `Ω` 
is `λ := im*f`, for a continuous-time system `sysrf`,
or `λ := exp(im*f*abs(Ts))`, for a discrete-time system `sysfr` with sampling-time `Ts`.  
If any of the frequency values in `freq` is a pole of `sysf`, then `sysf` is replaced by the stable numerator 
of a left coprime factorization of `sysf`. 

`FDGainTol = tol` specifies an absolute  threshold `tol` for the nonzero magnitudes of 
the frequency response gains (default: `tol = 0.01`). 


For `block = false`, `gains` is a `p x mf x nf` matrix, whose `(i,j,k)`-th element is the magnitude of the 
`(i,j)`-th element of `Rf(λ)` evaluated for the `k`-th frequency in `freq`. 
`S` is determined as a `q x mf x nf` 
binary matrix, whose `(i,j,k)`-th element is `S[i,j,k] = 1`, if the `(i,j,k)`-th element of `gains` is 
is larger then or equal to `tol`, and otherwise, `S[i,j,k] = 0`. 
For `block = true`, `gains` is an `1 x mf x nf` matrix, whose `(1,j,k)`-th element is the norm of the 
`j`-th column of `Rf(λ)` evaluated for the `k`-th frequency in `freq`.
`S` is determined as an `1 x mf x nf` binary matrix, whose `(1,j,k)`-th element is `S[1,j,k] = 1`, 
if the `(1,j,k)`-th element of `gains` is larger then or equal to `tol` 
and otherwise, `S[1,j,k] = 0`. 

The keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Af`, `Bf`, `Cf`, `Df`, the absolute tolerance for the nonzero elements of `Ef`,  
and the relative tolerance for the nonzero elements of `Af`, `Bf`, `Cf`, `Df` and `Ef`. 
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 

The computation of minimal realizations of individual input-output channels relies on pencil manipulation algorithms,
which employ rank determinations based on either the use of 
rank revealing QR-decomposition with column pivoting, if `fast = true`, or the SVD-decomposition.
The rank decision based on the SVD-decomposition is generally more reliable, but the involved computational effort is higher.

_Method:_ `S` is evaluated using the definition of the strong structure matrix in [1]. 

_References:_

[1] Varga A. Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017; sec. 3.4.

"""
function fdisspec(sysrf::DescriptorStateSpace{T}, freq::Union{AbstractVector{<:Real},Real} = 0; FDGainTol::Real = 0.01, block::Bool = false, 
                  atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, 
                  rtol::Real =  ((max(size(sysrf.A)...))+1)*eps(float(one(real(T))))*iszero(max(atol1,atol2)), fast::Bool = true) where T
   p, mf = size(sysrf) 
   isa(freq,Vector) || (freq = [freq]) 
   lfreq = length(freq);
   w = im*freq;                   # w = j*freq
   Ts = abs(sysrf.Ts);
   Ts > 0 && ( w = exp(Ts*w))     # w = exp(j*Ts*freq)

   try
      if block
         gs = evalfr(sysrf, w[1]; atol1, atol2, rtol, fast) 
         any(isinf.(gs)) && error("fdisspec:pole - the frequency $(w[1]) is a system pole")
         smat = falses(1, mf, lfreq)
         gains = zeros(T, 1, mf, lfreq)
         for j = 1:mf 
             gsj = norm(view(gs,:,j))
             gains[1,j,1] = gsj
             smat[1,j,1] = (gsj .> FDGainTol)
         end
         for i = 2:lfreq
             gs = evalfr(sysrf, w[i]; atol1, atol2, rtol, fast) 
             any(isinf.(gs)) && error("fdisspec:pole - the frequency $(w[i]) is a system pole")
             for j = 1:mf 
                 gsj = norm(view(gs,:,j))
                 gains[1,j,i] = gsj
                 smat[1,j,i] = (gsj .> FDGainTol)
             end
         end
      else
         gs = abs.(evalfr(sysrf, w[1]; atol1, atol2, rtol, fast)) 
         any(isinf.(gs)) && error("fdisspec:pole - the frequency $(w[1]) is a system pole")   
         smat = falses(p, mf, lfreq)
         gains = zeros(T, p, mf, lfreq)
         gains[:,:,1] = gs
         smat[:,:,1] = (gs .> FDGainTol)
         for i = 2:lfreq
            gs = abs.(evalfr(sysrf,w[i]; atol1, atol2, rtol, fast))
            any(isinf.(gs)) && error("fdisspec:pole - the frequency $(w[i]) is a system pole")
            gains[:,:,i] = gs
            smat[:,:,i] = (gs .> FDGainTol)
         end
      end
      return smat, gains
   catch err   
      findfirst("fdisspec:pole", string(err)) === nothing && rethrow() 
      return fdisspec(glcf(sysrf; atol1, atol2, rtol)[1], freq; FDGainTol, block, atol1, atol2, rtol, fast)
   end
end
fdisspec(sysr::FDFilterIF{T}, freq::Union{AbstractVector{<:Real},Real} = 0; kwargs...) where T = fdisspec(sysr.sys[:,sysr.faults], freq; kwargs...) 
"""
     fdscond(sysf,freq) -> (scond, β, γ)

Compute for a stable descriptor system `sysf = (A-λE,B,C,D)` with the transfer function matrix `Rf(λ)`, 
`β` - the H∞- index of `Rf(λ)`, `γ` - the maximum of the columns norms of `Rf(λ)` and 
the fault detection sensitivity condition `scond` evaluated as `scond := β/γ`. 
If `freq` is a vector of real frequency values, then `β` and `γ`
are evaluated over the frequencies contained in `freq`. 
If `sysf` is a fault detection filter internal form object `sysf::FDFilterIF` then `scond`, `β` and `γ`
are evaluated only for the fault inputs channel `sysf.sys[:,sysf.faults]`. 
"""
function fdscond(sysf::DescriptorStateSpace{T}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
   p, m = size(sysf)

   m == 0 && (return T[], T[], T[])

   isstable(sysf) || error("the system is unstable")

   β = Inf
   γ = 0 
   if ismissing(freq) 
      # evaluate β and γ as the minimum and maximum of H-infinity norms of the columns of G
      for j = 1:m
         temp = ghinfnorm(sysf[:,j])[1]
         γ >= temp || (γ = temp)
         β <= temp || (β = temp)
      end
   else
      # evaluateβ and γ as the minimum and maximum of the norms of columns of the frequency 
      # responses of G evaluated over all frequencies contained in FREQ 
      if !isa(freq, Vector) 
         # use evalfr if only one frequency is present
         H = evalfr(sysf; fval = freq) 
         for j = 1:m
             temp = norm(view(H,:,j))
             γ >= temp || (γ = temp)
             β <= temp || (β = temp)
         end
      else
         T1 = T <: BlasFloat ? T : promote_type(Float64,T) 
         ONE = one(T1)
         a, e, b, c, d = dssdata(T1,sysf)
         Ts = abs(sysf.Ts)
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
                 γ >= temp || (γ = temp)
                 β <= temp || (β = temp)
             end
         end
      end
   end
   return β/γ, β, γ
end
fdscond(sysr::FDFilterIF{T}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T = fdscond(sysr.sys[:,sysr.faults], freq) 

