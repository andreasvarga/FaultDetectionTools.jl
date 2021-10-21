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
function binmat2dec(S::Union{BitArray,VecOrMat{<:Number}})
# Convert a binary matrix with zero/non-zero elements to an equivalent decimal number integer vector
   len = size(S,2)
   n = Int.(S[:,1] .!= 0)
   for i = 2:len
       n .= 2*n .+ Int.(S[:,i] .!= 0)
   end
   return n
end
function dec2binmat(n::Union{Int,Vector{Int}}, numbits::Int = 1)
# Convert a decimal non-negative integer or integer vector with decimal non-negative values 
# to their binary matrix representation with `numbits` bits.
   minimum(n) < 0 && error("only non-negative values allowed")
   maxn = maximum(n)
   len = maxn == 0 ? max(1,numbits) : max(Int(floor(log2(maxn)+1)),numbits) 
   ns = size(n,1)
   v = isa(n,Int) ? string.([n], base = 2, pad = len) : string.(n, base = 2, pad = len) 
   S = falses(ns,len)  
   for i = 1:ns
       for j = 1:len
           S[i,j] = parse(Int,v[i][j]) == 1
       end
   end
   return S
end
"""
    S = fditspec_(sysrf::DescriptorStateSpace; FDfreq = missing, block = false, poleshift = false, 
                 FDtol, FDStol, atol = 0, atol1 = atol, atol2 = atol, rtol, fast = true) 

Compute the weak or strong binary structure matrix `S` of the transfer function matrix of a 
linear time-invariant system `sysrf` 
(typically representing the transfer channel from the fault inputs to residuals).
`sysrf` has a descriptor system realization of the form `sysrf = (Af-lambda*Ef,Bf,Cf,Df)` 
with a  `q x mf` transfer function matrix `Rf(λ)`. 
For the description of keyword parameters see the documentation of [`fditspec`](@ref). 
"""
function fditspec_(sysrf::DescriptorStateSpace{T}; FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                  block::Bool = false, poleshift::Bool = false, FDtol::Real = 0., FDStol::Real = 0., 
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
      S1 =  abs.(d) .> (FDtol <= 0 ? 0.0001*max(1., norm(d,1)) : FDtol) 
      S = block ? maximum(S1, dims=1) : S1
      return Sstrong ? repeat(S,1,1,lfreq) : S
   end
   if Sstrong
      #FDStol <= 0. && (FDStol = 0.0001*max(1., norm(a,1), norm(b,1), norm(c,Inf), norm(d,1), standard ? 0 : norm(e,1))) 
      FDStol <= 0. && (FDStol = eps(max(1., norm(a,1), norm(b,1), norm(c,Inf), norm(d,1), standard ? 0 : norm(e,1)))) 
      S = trues(q, mf, lfreq)
      p1, n1 = size(c,1), size(c,2)
      # employ structural analysis to compute weak/strong structure matrix  
      ispole = false
      for k = 1:lfreq
         # check if freq(k) is a pole 
         ispole = isinf(FDfreq[k]) ? (e == I ? false : rank(e) < n1) : rank(a-w[k]*e) < n1 
         ispole || break
      end
      ispole && !poleshift && error("sysf has poles in Ω")
      if block
         ispole && (f = rand(T,n1,p1); mul!(a,f,c,1,1); mul!(b,f,d,1,1) ) 
         for j = 1:mf 
             # elliminate uncontrollable and unobservable eigenvalues for the j-th column of B
             a1, e1, b1, c1, d1 = lsminreal(a, e, view(b,:,j), c, view(d,:,j); noseig = false, atol1, atol2, rtol, fast)
             for k = 1:lfreq
                # check if freq(k) is a zero of the j-th column 
                s = isinf(FDfreq[k]) ? svdvals!([e1 b1; c1 d1]) : svdvals!([a1-w[k]*e1 b1; c1 d1])
                S[1, j, k] = (s[end] > FDStol)
            end
         end
      else
         for i = 1:p 
             # elliminate unobservable eigenvalues for the i-th row of C
             a1, e1, b1, c1, d1 = lsminreal(a, e, b, view(c,i:i,:), view(d,i:i,:); contr = false, noseig = false, atol1, atol2, rtol, fast)
             n1 = size(a1,1)
             ispole && (f1 = rand(T,n1,1); mul!(a1,f1,c1,1,1); mul!(b1,f1,d1,1,1) ) 
             for j = 1:mf 
                 # elliminate uncontrollable and non-dynamic eigenvalues for the j-th column of B
                 a2, e2, b2, c2, d2 = lsminreal(a1, e1, view(b1,:,j:j), c1, view(d1,:,j:j); obs = false, noseig = true, atol1, atol2, rtol, fast)
                 for k = 1:lfreq
                    # check if freq(k) is a zero of the (i,j)-th element 
                    s = isinf(FDfreq[k]) ? svdvals!([e2 b2; c2 d2]) : svdvals!([a2-w[k]*e2 b2; c2 d2])
                    S[i,j,k] = s[end] > FDStol
                 end
             end
         end
      end
   else
      FDtol <= 0. && (FDtol = 0.0001*max(1., norm(b,1), norm(d,1)))
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
"""
     fdisspec_(sysrf::DescriptorStateSpace, freq; block = false, stabilize = false, FDGainTol = 0.01, 
                     atol, atol1, atol2, atol3, rtol, fast = true) -> (S, gains)

Compute the strong binary structure matrix `S` of the transfer function matrix of a 
linear time-invariant system `sysrf` 
(typically representing the transfer channel from the fault inputs to residuals).
`sysrf` has a descriptor system realization of the form `sysrf = (Af-lambda*Ef,Bf,Cf,Df)` 
with a  `q x mf` transfer function matrix `Rf(λ)`. 
For the description of keyword parameters see the documentation of [`fdisspec`](@ref). 
"""
function fdisspec_(sysrf::DescriptorStateSpace{T}, freq::Union{AbstractVector{<:Real},Real} = 0; stabilize::Bool = false, FDGainTol::Real = 0.01, block::Bool = false, 
                   atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, atol3::Real = atol, 
                   rtol::Real =  ((max(size(sysrf.A)...))+1)*eps(float(one(real(T))))*iszero(max(atol1,atol2)), fast::Bool = true) where T
   p, mf = size(sysrf) 
   isa(freq,Vector) || (freq = [freq]) 
   lfreq = length(freq);
   w = im*freq;                   # w = j*freq
   Ts = abs(sysrf.Ts);
   Ts > 0 && ( w = exp(Ts*w))     # w = exp(j*Ts*freq)

   ispole = false
   n = size(sysrf.A,1)
   for k = 1:lfreq
      # check if freq(k) is a pole 
      ispole = isinf(freq[k]) ? (sysrf.e == I ? false : rank(sysrf.E) < n) : rank(sysrf.A-w[k]*sysrf.E) < n 
      ispole || break
   end
   ispole && !stabilize && error("sysf has poles in Ω")
   if block
      stabilize && (sysrf = glcf(sysrf; atol1, atol2, atol3, rtol)[1])
      # gs = evalfr(sysrf, w[1]; atol1, atol2, rtol, fast) 
      # any(isinf.(gs)) && error("fdisspec:pole - the frequency $(w[1]) is a system pole")
      smat = falses(1, mf, lfreq)
      gains = zeros(T, 1, mf, lfreq)
      for i = 1:lfreq
          gs = evalfr(sysrf, w[i]; atol1, atol2, rtol, fast) 
          any(isinf.(gs)) && error("fdisspec_:pole - the frequency $(w[i]) is a system pole")
          for j = 1:mf 
              gsj = norm(view(gs,:,j))
              gains[1,j,i] = gsj
              smat[1,j,i] = (gsj .> FDGainTol)
          end
      end
   else
      smat = falses(p, mf, lfreq)
      gains = zeros(T, p, mf, lfreq)
      for k = 1:p
         t = gir(sysrf[k,:]; contr = false, atol1, atol2, rtol)
         for i = 1:lfreq
             gs = stabilize ? abs.(evalfr(glcf(t; atol1, atol2, atol3, rtol)[1], w[i]; atol1, atol2, rtol, fast)) :
                              abs.(evalfr(t, w[i]; atol1, atol2, rtol, fast))
             any(isinf.(gs)) && error("fdisspec_:pole - the frequency $(w[i]) is a system pole")
             gains[k,:,i] = gs
             smat[k,:,i] = (gs .> FDGainTol)
         end
      end
   end
   return smat, gains
end
"""
     fdiscond_(sysrf::DescriptorStateSpace, freq) -> (scond, β, γ)

Compute for a stable descriptor system `sysrf = (A-λE,B,C,D)` with the transfer function matrix `Rf(λ)`, 
`β` - the H∞- index of `Rf(λ)`, `γ` - the maximum of the columns norms of `Rf(λ)` and 
`scond` - the column-gains sensitivity condition evaluated as `scond := β/γ`. 
If `freq` is a vector of real frequency values, then `β` and `γ`
are evaluated over the frequencies contained in `freq`. 
"""
function fdiscond_(sysrf::DescriptorStateSpace{T}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
   p, m = size(sysrf)

   m == 0 && (return T[], T[], T[])

   isstable(sysrf) || error("the system is unstable")

   β = Inf
   γ = 0 
   if ismissing(freq) 
      # evaluate β and γ as the minimum and maximum of H-infinity norms of the columns of G
      for j = 1:m
         temp = ghinfnorm(sysrf[:,j])[1]
         γ >= temp || (γ = temp)
         β <= temp || (β = temp)
      end
   else
      # evaluate β and γ as the minimum and maximum of the norms of columns of the frequency 
      # responses of G evaluated over all frequencies contained in FREQ 
      if !isa(freq, Vector) 
         # use evalfr if only one frequency is present
         H = evalfr(sysrf; fval = freq) 
         for j = 1:m
             temp = norm(view(H,:,j))
             γ >= temp || (γ = temp)
             β <= temp || (β = temp)
         end
      else
         T1 = T <: BlasFloat ? T : promote_type(Float64,T) 
         ONE = one(T1)
         a, e, b, c, d = dssdata(T1,sysrf)
         Ts = abs(sysrf.Ts)
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
