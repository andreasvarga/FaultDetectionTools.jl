"""
    S = fditspec(sysr::FDFilterIF; FDfreq = missing, block = false, poleshift = false, 
                 FDtol, FDStol, atol = 0, atol1 = atol, atol2 = atol, rtol, fast = true) 

Compute the weak or strong binary structure matrix `S` of the transfer function matrix `Rf(λ)` of 
the transfer channel from the fault inputs to residuals of 
a fault detection filter internal form object `sysr::FDFilterIF`.  
For a filter `sysr` with `q` residual outputs and `mf` fault inputs, 
`Rf(λ)` is the `q x mf` transfer function matrix of the fault inputs channel with the descriptor system representation
`sysr.sys[:,sysr.faults] := (Af-lambda*Ef,Bf,Cf,Df)`. 

If `FDfreq = missing` (default), then `S` contains the weak structure matrix of `Rf(λ)`. 
For `block = false`, `S` is determined as a `q x mf` 
binary matrix (BitMatrix), whose `(i,j)`-th element is `S[i,j] = 1`, if the `(i,j)`-th element of `Rf(λ)` is 
nonzero, and otherwise, `S[i,j] = 0`. 
For `block = true`, `S` is determined as a `1 x mf` binary matrix, whose `(1,j)`-th element is `S[1,j] = 1`, 
if the `j`-th column of `Rf(λ)` is nonzero, and otherwise, `S[1,j] = 0`. 

If `FDfreq = freq` specifies a vector `freq` of `nf` real frequencies 
which characterize the classes of persistent fault signals, then 
for a suitable proper and invertible `M(λ)` (see below),  
`S` contains the strong structure matrix of `M(λ)*Rf(λ)` with respect to a set of `nf` complex frequencies `Ω`, defined as follows: 
if `f` is a real frequency in `freq`, then the corresponding complex frequency in `Ω` 
is `λ := im*f`, for a continuous-time system,
or `λ := exp(im*f*abs(Ts))`, for a discrete-time system with sampling-time `Ts`. 

`FDtol = tol1` specifies an absolute threshold `tol1` for the magnitudes of nonzero elements in the system matrices 
`Bf` and `Df` and is used to determine the weak structure matrix. 
Its default value is `tol1 = 0.0001*max(1, norm(Bf,1), norm(Df,1))`. 

`FDStol = tol2` specifies an absolute  threshold `tol2` for the magnitudes of nonzero elements in the system matrices 
`Af`, `Ef`, `Bf`, `Cf` and `Df` and is used to determine the strong structure matrix. 
Its default value is 
`tol2 = epsm*max(1, norm(Ef,1), norm(Af,1), norm(Bf,1), norm(Cf,Inf), norm(Df,1)))`, 
where `epsm` is the working machine precision.

For `block = false`, then, if `poleshift = true`, `M(λ)` is chosen diagonal such that `M(λ)*Rf(λ)`
has no poles in `Ω` and if `poleshift = false` (default), `M(λ) = I` is used and 
an error is issued if `Rf(λ)` has poles in `Ω`. 
`S` is determined as a `q x mf` binary matrix, whose `(i,j)`-th element is `S[i,j] = 1`, 
if the `(i,j)`-th element of `M(λ)*Rf(λ)` 
evaluated for all frequencies in `Ω` is nonzero, and otherwise, `S[i,j] = 0`.  

For `block = true`, then, if `poleshift = true`, `M(λ)` is chosen such that `M(λ)*Rf(λ)` 
as no poles in `Ω` and if `poleshift = false` (default), `M(λ) = I` is used and 
an error is issued if `Rf(λ)` has poles in `Ω`. 
`S` is determined as an `1 x mf` binary matrix, whose `(1,j)`-th element is `S[1,j] = 1`, 
if the `j`-th column of `M(λ)*Rf(λ)` evaluated for all frequencies in `Ω`
is nonzero and otherwise `S[1,j] = 0`. 

The keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Af`, `Bf`, `Cf`, `Df`, the absolute tolerance for the nonzero elements of `Ef`,  
and the relative tolerance for the nonzero elements of `Af`, `Bf`, `Cf`, `Df` and `Ef`. 
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 


_Method:_ For the definition of the structure matrix, see [1]. For the
determination of the weak structure matrix, minimal realizations
are determined for each column of `Rf(λ)` if `block = true` or for 
each element of `Rf(λ)` if `block = false` and the nonzero columns or 
elements in each column are identified  (see Corollary 7.1 of [1]).
For the determination of the strong structure matrix, minimal realizations
are determined for each column of `M(λ)*Rf(λ)` if `block = true` or for 
each element of `M(λ)*Rf(λ)` if `block = false` and  the full rank of the
corresponding system matrix is checked for all frequencies in `FDfreq`
(see Corollary 7.2 in [1]) (i.e., the lack of zeros in all frequencies).

_References:_

[1] Varga A. Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017; sec.3.4.
"""
function fditspec(sysr::Union{FDFilterIF{T},FDIModel{T}}; 
    block::Bool = false, FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
    poleshift::Bool = false, FDtol::Real = 0., FDStol::Real = 0., 
    atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, 
    rtol::Real = 0, fast::Bool = true) where T

    S = fditspec_(sysr.sys[:,sysr.faults]; block, FDfreq, poleshift, FDtol, FDStol, 
                  atol1, atol2, rtol, fast)
    ismissing(FDfreq) && (return S)
    for j = 1:size(S,2)
        for i = 1:size(S,1)
            S[i,j,1] = all(view(S,i,j,:))
        end
    end
    return S[:,:,1]
end
"""
    S = fditspec(sysr::FDIFilterIF; FDfreq = missing, poleshift = false, 
                 FDtol, FDStol, atol = 0, atol1 = atol, atol2 = atol, rtol, fast = true) 

Compute the weak or strong binary structure matrix `S` of the global transfer function matrix `Rf(λ)` of 
the transfer channel from the fault inputs to residuals of 
a fault detection and isolation filter internal form object `sysr::FDIFilterIF`.  
The filter `sysr` consists of `N` individual FDI filters `sysr.sys[i]`, for `i = 1, ..., N`, where
the fault to residual channel of the `i`-th filter `sysr.sys[i][:,sysr.faults]` 
has `qi` residual outputs and `mf` fault inputs, has the descriptor system representation
`sysr.sys[i][:,sysr.faults] := (Afi-lambda*Efi,Bfi,Cfi,Dfi)` and 
`Rfi(λ)` is the corresponding `qi x mf` transfer function matrix. 
The global transfer function matrix `Rf(λ)` is formed by row concatenation of the 
transfer function matrices of the `N` individual filters, i.e.,  `Rf(λ) := [ Rf1(λ); Rf2(λ); ...; RfN(λ)]`. 
For the evaluation of the strong structure matrix, the structure matrix of the stable 
transfer function matrix `M(λ)*Rf(λ)` is determined, with a `M(λ)` block-diagonal
`M(λ) = block-diag(M1(λ), M2(λ), ..., MN(λ))`, where `Mi(λ)` is a suitable square and invertible 
transfer function matrix (see below). 

`FDtol = tol1` specifies an absolute threshold `tol1` for the magnitudes of nonzero elements in the system matrices 
`Bf` and `Df` and is used to determine the weak structure matrix. 
Its default value is `tol1 = 0.0001*max(1, norm(Bf,1), norm(Df,1))`. 

`FDStol = tol2` specifies an absolute  threshold `tol2` for the magnitudes of nonzero elements in the system matrices 
`Af`, `Ef`, `Bf`, `Cf` and `Df` and is used to determine the strong structure matrix. 
Its default value is 
`tol2 = epsm*max(1, norm(Ef,1), norm(Af,1), norm(Bf,1), norm(Cf,Inf), norm(Df,1)))`, 
where `epsm` is the working machine precision.

If `FDfreq = missing` (default), then `S` contains the weak structure matrix of `Rf(λ)`. 
`S` is determined as a `N x mf` binary matrix, whose `(i,j)`-th element is `S[i,j] = 1`, 
if the `j`-th column of `Rfi(λ)` is nonzero, and otherwise, `S[i,j] = 0`. 

If `FDfreq = freq` specifies a vector `freq` of `nf` real frequencies 
which characterize the classes of persistent fault signals, then 
for a suitable proper and invertible `M(λ)` (see below),  
`S` contains the strong structure matrix of `M(λ)*Rf(λ)` with respect to a set of `nf` complex frequencies `Ω`, defined as follows: 
if `f` is a real frequency in `freq`, then the corresponding complex frequency in `Ω` 
is `λ := im*f`, for a continuous-time system,
or `λ := exp(im*f*abs(Ts))`, for a discrete-time system with sampling-time `Ts`. 

`S` is determined as a `N x mf` binary matrix, whose `(i,j)`-th element is `S[i,j] = 1`, 
if the `j`-th column of `Mi(λ)*Rfi(λ)` is nonzero for all frequencies in `Ω`, and otherwise, `S[i,j] = 0`. 
If `poleshift = true`, `Mi(λ)` is chosen such that `Mi(λ)*Rfi(λ)` has no poles in `Ω` and
if `poleshift = false` (default), `Mi(λ) = I` is used and an error is issued if any `Rfi(λ)` has poles in `Ω`. 
   
The keyword arguments `atol1`, `atol2`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Af`, `Bf`, `Cf`, `Df`, the absolute tolerance for the nonzero elements of `Ef`,  
and the relative tolerance for the nonzero elements of `Af`, `Bf`, `Cf`, `Df` and `Ef`. 
The keyword argument `atol` can be used to simultaneously set `atol1 = atol` and `atol2 = atol`. 


_Method:_ For the definition of the structure matrix, see [1]. For the
determination of the weak structure matrix, minimal realizations
are determined for each column of `Rfi(λ)` and the nonzero columns are identified  (see Corollary 7.1 of [1]).
For the determination of the strong structure matrix, minimal realizations
are determined for each column of `Mi(λ)*Rfi(λ)` and  the full rank of the
corresponding system matrix is checked for all frequencies in `Ω`
(see Corollary 7.2 in [1]) (i.e., the lack of zeros in all frequencies in `Ω`).

_References:_

[1] Varga A. Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017; sec.3.4.
"""
function fditspec(sysr::FDIFilterIF{T}; FDfreq::Union{AbstractVector{<:Real},Real,Missing} = missing, 
                  poleshift::Bool = false, FDtol::Real = 0., FDStol::Real = 0., 
                  atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, 
                  rtol::Real = 0, fast::Bool = true) where T
    N = length(sysr.sys)
    mf = length(sysr.faults)
    S = trues(N,mf)
    if ismissing(FDfreq)
       for i = 1:N
           S[i,:] = fditspec_(sysr.sys[i][:,sysr.faults]; block = true, FDtol, atol1, atol2, fast, 
                                         rtol = (size(sysr.sys[i].A,1)+1)*eps(float(one(real(T))))*iszero(max(atol1,atol2,rtol)))
       end
    else
       isa(FDfreq,Vector) || (FDfreq = [FDfreq]) 
       lfreq = length(FDfreq);
       for i = 1:N
           Sc  = fditspec_(sysr.sys[i][:,sysr.faults]; block = true, FDfreq, poleshift, FDtol, FDStol, atol1, atol2, fast, 
                                         rtol = (size(sysr.sys[i].A,1)+1)*eps(float(one(real(T))))*iszero(max(atol1,atol2,rtol)))
           lfreq == 1 ? (S[i,:] = Sc[1,:,1]) : ([S[i,j] = all(view(Sc,1,j,:)) for j = 1:mf])           
       end
    end
    return S
end
"""
     S = fdisspec(sysr::FDFilterIF, freq; block = false, stabilize = false, FDGainTol = 0.01, 
                     atol = 0, atol1 = atol, atol2 = atol, rtol, fast = true) 

Compute, for a given set of real frequencies `freq`,  
the strong binary structure matrix `S` of the stable transfer function matrix `M(λ)*Rf(λ)`, 
where `Rf(λ)` is the transfer function matrix of the
transfer channel from the fault inputs to residuals of the
fault detection filter internal form object `sysr::FDFilterIF` and  
`M(λ)` is a suitable proper and invertible stabilizing transfer function matrix (see below).  
For a filter `sysr` with `q` residual outputs and `mf` fault inputs, 
`Rf(λ)` is the `q x mf` transfer function matrix of the fault inputs channel with the descriptor system representation
`sysr.sys[:,sysr.faults] := (Af-lambda*Ef,Bf,Cf,Df)`. 

`freq` must contain a real frequency value or a vector of `nf` real frequencies 
which characterize the classes of persistent fault signals 
(default: `freq = 0`, i.e., characterizing constant faults).
`S` contains the strong 
structure matrix of `M(λ)*Rf(λ)` with respect to a set of `nf` complex frequencies `Ω`, defined as follows: 
if `f` is a real frequency in `freq`, then the corresponding complex frequency in `Ω` 
is `λ := im*f`, for a continuous-time system,
or `λ := exp(im*f*abs(Ts))`, for a discrete-time system with sampling-time `Ts`.  

`FDGainTol = tol` specifies an absolute  threshold `tol` for the nonzero magnitudes of 
the frequency response gains (default: `tol = 0.01`). 

For `block = false`, then, if `stabilize = true`, `M(λ)` is chosen diagonal such that `M(λ)*Rf(λ)`
has only stable poles and if `stabilize = false` (default), `M(λ) = I` is used and 
an error is issued if `Rf(λ)` has poles in `Ω`. 
`S` is determined as a `q x mf` binary matrix, whose `(i,j)`-th element is `S[i,j] = 1`, 
if the `(i,j)`-th element of `M(λ)*Rf(λ)` 
evaluated for all frequencies in `freq` is larger than or equal to `tol`, and otherwise, `S[i,j] = 0`.  


For `block = true`, then, if `stabilize = true`, `M(λ)` is chosen such that `M(λ)*Rf(λ)` 
has only stable poles and if `stabilize = false` (default), `M(λ) = I` is used and 
an error is issued if `Rf(λ)` has poles in `Ω`. 
`S` is determined as an `1 x mf` binary matrix, whose `(1,j)`-th element is `S[1,j] = 1`, 
if the `j`-th column of `M(λ)*Rf(λ)` evaluated for all frequencies in `Ω`
is larger than or equal to `tol` and otherwise, `S[1,j] = 0`. 

The keyword arguments `atol1`, `atol2`, `atol3`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Af`, `Bf`, `Cf`, `Df`, the absolute tolerance for the nonzero elements of `Ef`, 
the absolute tolerance for the nonzero elements of `Cf`,   
and the relative tolerance for the nonzero elements of `Af`, `Bf`, `Cf`, `Df` and `Ef`. 
The keyword argument `atol` can be used to simultaneously set `atol1 = atol`, `atol2 = atol` and `atol3 = atol`. 

The computation of minimal realizations of individual input-output channels relies on pencil manipulation algorithms,
which employ rank determinations based on either the use of 
rank revealing QR-decomposition with column pivoting, if `fast = true`, or the SVD-decomposition.
The rank decision based on the SVD-decomposition is generally more reliable, but the involved computational effort is higher.

_Method:_ `S` is evaluated using the definition of the strong structure matrix in [1]. 

_References:_

[1] Varga A. Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017; sec. 3.4.

"""
function fdisspec(sysr::Union{FDFilterIF{T},FDIModel{T}}, freq::Union{AbstractVector{<:Real},Real} = 0; kwargs...) where T
    S = fdisspec_(sysr.sys[:,sysr.faults], freq; kwargs...)[1] 
    length(freq) == 1 && (return S[:,:,1])
    for j = 1:size(S,2)
        for i = 1:size(S,1)
            S[i,j,1] = all(view(S,i,j,:))
        end
    end
    return S[:,:,1]
end  
"""
     S = fdisspec(sysr::FDIFilterIF, freq; stabilize = false, FDGainTol = 0.01, 
                     atol = 0, atol1 = atol, atol2 = atol, rtol, fast = true) 

Compute, for a given set of real frequencies `freq`,  
the strong binary structure matrix `S` 
of the stable transfer function matrix `M(λ)*Rf(λ)`, where `Rf(λ)` is the global 
transfer function matrix of the transfer channel from the fault inputs to residuals of the
fault detection and isolation filter internal form object `sysr::FDIFilterIF` and  
`M(λ)` is a suitable block-diagonal proper and invertible stabilizing transfer function matrix (see below).  
The filter `sysr` consists of `N` individual FDI filters `sysr.sys[i]`, for `i = 1, ..., N`, where
the fault to residual channel of the `i`-th filter `sysr.sys[i][:,sysr.faults]` 
has `qi` residual outputs and `mf` fault inputs, has the descriptor system representation
`sysr.sys[i][:,sysr.faults] := (Afi-lambda*Efi,Bfi,Cfi,Dfi)` and 
`Rfi(λ)` is the corresponding `qi x mf` transfer function matrix. 
The global transfer function matrix `Rf(λ)` is formed by row concatenation of the 
transfer function matrices of the `N` individual filters, i.e.,  `Rf(λ) := [ Rf1(λ); Rf2(λ); ...; RfN(λ)]`. 
`M(λ) = block-diag(M1(λ), M2(λ), ..., MN(λ))`, where `Mi(λ)` is square and invertible 
and chosen such that `Mi(λ)Rfi(λ)` is stable (see below). 

`freq` must contain a real frequency value or a vector of `nf` real frequencies 
which characterize the classes of persistent fault signals 
(default: `freq = 0`, i.e., characterizing constant faults).
`S` contains the strong 
structure matrix of `M(λ)*Rf(λ)` with respect to a set of `nf` complex frequencies `Ω`, defined as follows: 
if `f` is a real frequency in `freq`, then the corresponding complex frequency in `Ω` 
is `λ := im*f`, for a continuous-time system,
or `λ := exp(im*f*abs(Ts))`, for a discrete-time system with sampling-time `Ts`.  

`FDGainTol = tol` specifies an absolute  threshold `tol` for the nonzero magnitudes of 
the frequency response gains (default: `tol = 0.01`). 

If `stabilize = true`, `Mi(λ)` is chosen such that `Mi(λ)*Rfi(λ)` has only stable poles and
if `stabilize = false` (default), `Mi(λ) = I` is used and an error is issued 
if any `Rfi(λ)` has poles in `Ω`. 

`S` is determined as a `N x mf` 
binary matrix, whose `(i,j)`-th element is `S[i,j] = 1`, if the norm of the 
`j`-th column of `Mi(λ)*Rfi(λ)` evaluated for all frequencies in `Ω` 
is larger than or equal to `tol`, and otherwise, `S[i,j] = 0`. 

The keyword arguments `atol1`, `atol2`, `atol3`, and `rtol`, specify, respectively, the absolute tolerance for the 
nonzero elements of matrices `Afi`, `Bfi`, `Cfi`, `Dfi`, the absolute tolerance for the nonzero elements of `Efi`, 
the absolute tolerance for the nonzero elements of `Cfi`,   
and the relative tolerance for the nonzero elements of `Afi`, `Bfi`, `Cfi`, `Dfi` and `Efi`. 
The keyword argument `atol` can be used to simultaneously set `atol1 = atol`, `atol2 = atol` and `atol3 = atol`. 

The computation of minimal realizations of individual input-output channels relies on pencil manipulation algorithms,
which employ rank determinations based on either the use of 
rank revealing QR-decomposition with column pivoting, if `fast = true`, or the SVD-decomposition.
The rank decision based on the SVD-decomposition is generally more reliable, but the involved computational effort is higher.

_Method:_ `S` is evaluated using the definition of the strong structure matrix in [1]. 

_References:_

[1] Varga A. Solving Fault Diagnosis Problems - Linear Synthesis Techniques. Springer Verlag, 2017; sec. 3.4.

"""
function fdisspec(sysr::FDIFilterIF{T}, freq::Union{AbstractVector{<:Real},Real} = 0; stabilize::Bool = false, 
                  FDGainTol::Real = 0.01, atol::Real = zero(float(real(T))), atol1::Real = atol, atol2::Real = atol, 
                  rtol::Real = 0, fast::Bool = true) where T
    isa(freq,Vector) || (freq = [freq]) 
    lfreq = length(freq);
    N = length(sysr.sys)
    mf = length(sysr.faults)
    S = trues(N,mf)
    for i = 1:N
        Sc = fdisspec_(sysr.sys[i][:,sysr.faults], freq; block = true, stabilize, FDGainTol, 
                       atol1, atol2, fast, rtol = (size(sysr.sys[i].A,1)+1)*eps(float(one(real(T))))*iszero(max(atol1,atol2,rtol)))[1]
        lfreq == 1 ? (S[i,:] = Sc[1,:,1]) : ([S[i,j] = all(view(Sc,1,j,:)) for j = 1:mf])           
    end
    return S
end
"""
     fdiscond(sysr::FDFilterIF,freq) -> (scond, β, γ)

Compute for the stable transfer function matrix `Rf(λ)` of the
transfer channel from the fault inputs to residuals of the
fault detection filter internal form object `sysr::FDFilterIF` the quantities: 
`β` - the H∞- index of `Rf(λ)`, `γ` - the maximum of the columns norms of `Rf(λ)` and 
the fault detection sensitivity condition `scond` evaluated as `scond := β/γ`. 
If `freq` is a vector of real frequency values, then `β` and `γ`
are evaluated over the frequencies contained in `freq`. 
"""
function fdiscond(sysr::FDFilterIF{T}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
    fdiscond_(sysr.sys[:,sysr.faults], freq)
end
"""
     fdiscond(sysr::FDFilterIF,SFDI,freq) -> (scond, β, γ)

Compute the detection and isolation sensitivity condition
`scond` (and related quatities `β` and `γ`) for the binary structure vector `SFDI` associated to
the stable transfer function matrix `Rf(λ)` of the
transfer channel from the fault inputs to residuals of the
fault detection filter internal form object `sysr::FDFilterIF`. 
If  `Rff(λ)` is the transfer function matrix formed of those `j`-th columns of `Rf(λ)` 
for which `SFDI[j] = 1`, then:   
`β` - the H∞- index of `Rff(λ)`, `γ` - the maximum of the columns norms of `Rff(λ)` and 
the fault detection sensitivity condition `scond` evaluated as `scond := β/γ`. 
If `freq` is a vector of real frequency values, then `β` and `γ`
are evaluated over the frequencies contained in `freq`. 
"""
function fdiscond(sysr::FDFilterIF{T}, SFDI::Union{BitVector,Vector{Bool}}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
    fdiscond_(sysr.sys[:,sysr.faults[SFDI]], freq)
end
"""
     fdiscond(sysr::FDFilterIF, SFDI, freq) -> (scond, β, γ)

Compute the detection and isolation sensitivity condition vector
`scond` (and related quatities `β` and `γ`) for the `q × mf` binary structure matrix `SFDI` associated to
the stable transfer function matrix `Rf(λ)` of the
transfer channel from the fault inputs to residuals of the
fault detection filter internal form object `sysr::FDFilterIF`. 
The filter `sysr` consists of `N` individual FDI filters `sysr.sys[i]`, for `i = 1, ..., N`, where
`sysr.sys[i][:,sysr.faults]` is the fault to residual channel of the `i`-th filter and 
`Rfi(λ)` is the corresponding transfer function matrix. 
The global transfer function matrix `Rf(λ)` is formed by row concatenation of the 
transfer function matrices of the `N` individual filters, i.e.,  `Rf(λ) := [ Rf1(λ); Rf2(λ); ...; RfN(λ)]`.
The i-th element of the vectors `scond`, `β` and `γ` contain the quantities: 
`β[i]` - the H∞- index of the nonzero columns of `Rffi(λ)`, `γ` - the maximum of the nonzero columns norms of `Rfi(λ)` and 
the correspomding fault detection sensitivity condition `scond[i]` evaluated as `scond[i] := β[i]/γ[i]`,
where `Rffi(λ)` is formed of those `j`-th  columns of the `i`-th row of `Rf(λ)` for which `S[i,j] = 1` and 
`Rfi(λ)` is the `i`-th row of `Rf(λ)`. 
It is assumed that for each `j` such that `SFDI[i,j] = 1`, the `j`-th column of `Rfi(λ)` is nonzero  and 
for each `j` such that `SFDI[i,j] = 0`, the `j`-th column of `Rfi(λ)` is zero. 
If `freq` is a vector of real frequency values, then `β[i]` and `γ[i]`
are evaluated over the frequencies contained in `freq`. 
"""
function fdiscond(sysr::FDFilterIF{T}, SFDI::Union{BitMatrix,Array{Bool,2}}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
    q = size(sysr.sys,1)   
    inpf = sysr.faults
    mf = length(inpf)
    (q, mf ) == size(SFDI) || error("missmatch between the number of fault inputs and column dimension of SFDI")
 
    mf == 0 && (return zeros(T,q,0), zeros(T,q,0), zeros(T,q,0))
    scond = similar(Array{T,1},q)
    β = similar(Array{T,1},q)
    γ = similar(Array{T,1},q)
    for i = 1:q
        scond[i], β[i], γ[i] = fdiscond_(sysr.sys[i,inpf[view(SFDI,i,:)]], freq)
    end
    return scond, β, γ
end
"""
     fdiscond(sysr::FDIFilterIF, SFDI, freq) -> (scond, β, γ)

Compute the detection and isolation sensitivity condition
`scond` (and related quatities `β` and `γ`) for the `N × mf` binary structure matrix `SFDI` associated to
the stable global transfer function matrix `Rf(λ)` of the
transfer channel from the fault inputs to residuals of the
fault detection and isolation filter internal form object `sysr::FDIFilterIF`. 
The filter `sysr` consists of `N` individual FDI filters `sysr.sys[i]`, for `i = 1, ..., N`, where
`sysr.sys[i][:,sysr.faults]` is the fault to residual channel of the `i`-th filter and 
`Rfi(λ)` is the corresponding transfer function matrix. 
The global transfer function matrix `Rf(λ)` is formed by row concatenation of the 
transfer function matrices of the `N` individual filters, i.e.,  `Rf(λ) := [ Rf1(λ); Rf2(λ); ...; RfN(λ)]`.
It is assumed that for each `j` such that `SFDI[i,j] = 1`, the `j`-th column of `Rfi(λ)` is nonzero  and 
for each `j` such that `SFDI[i,j] = 0`, the `j`-th column of `Rfi(λ)` is zero. 
The i-th element of the vectors `scond`, `β` and `γ` contain the quantities: 
`β[i]` - the H∞- index of the nonzero columns of `Rfi(λ)`, `γ` - the maximum of the nonzero columns norms of `Rf(λ)` and 
the fault detection sensitivity condition `scond` evaluated as `scond[i] := β[i]/γ[i]`. 
If `freq` is a vector of real frequency values, then `β[i]` and `γ[i]`
are evaluated over the frequencies contained in `freq`. 
"""
function fdiscond(sysr::FDIFilterIF{T}, SFDI::Union{BitMatrix,BitVector,Array{Bool,2},Array{Bool,1}}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
    N = length(sysr.sys)
    isa(SFDI,Union{BitVector,Array{Bool,1}}) ? (nb = 1; mf = length(SFDI); SFDI = reshape(SFDI,nb,mf) ) :
                                               (nb = size(SFDI,1); mf = size(SFDI,2))
    nb == N || error("missmatch between number of filters $N and number of specifications $nb")
    inpf = sysr.faults
    mf == length(inpf) || error("missmatch between number of faults and column dimension/length of SFDI")
    scond = similar(Array{T,1},N)
    β = similar(Array{T,1},N)
    γ = similar(Array{T,1},N)
    for i = 1:N
        scond[i], β[i], γ[i] = fdiscond_(sysr.sys[i][:,inpf[view(SFDI,i,:)]], freq)
    end
    return scond, β, γ
end
fdiscond(sysr::FDIFilterIF{T}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T = fdiscond(sysr,trues(1,length(sysr.faults)),freq)

"""
     fdif2ngap(sysr::FDFilterIF, freq) -> (gap, β, γ)

Compute the fault-to-noise gap `gap` (and the related quantities `β` and `γ`) 
for the stable fault detection filter internal form object `sysr::FDFilterIF`.
For the fault to residual channel of the filter `sysr.sys[:,sysr.faults]`   
with the corresponding transfer function matrix `Rf(λ)` and 
the noise to residual channel of the filter `sysr.sys[:,sysr.noise]`   
with the corresponding transfer function matrix `Rw(λ)`,    
`β` is the H∞- index of `Rf(λ)`, `γ` is the H∞-norm of `Rw(λ)` and 
`gap` is the fault-to-noise gap evaluated as `gap := β/γ`. 
If `freq` is a vector of real frequency values, then `β` and `γ`
are evaluated over the frequencies contained in `freq`. 
`gap = ∞` if there are no noise inputs and `gap = 0` if there are no fault inputs.
"""
function fdif2ngap(sysr::FDFilterIF{T}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
   mf = length(sysr.faults)
   mw = length(sysr.noise)

   mf+mw == 0 && (return T[], T[], T[])

   isstable(sysr.sys) || error("the system is unstable")
   
   β = mf == 0 ? 0 : fdhinfminus(sysr.sys[:,sysr.faults],freq)[1]
   γ = mw == 0 ? 0 : ghinfnorm(sysr.sys[:,sysr.noise])[1]
   return β/γ, β, γ
end
"""
     fdif2ngap(sysr::FDFilterIF, SFDI, freq) -> (gap, β, γ)

Compute the fault-to-noise gap `gap` (and the related quantities `β` and `γ`) 
for the stable fault detection filter internal form object `sysr::FDFilterIF`
and the associated binary structure vector `SFDI`.
`sysr.sys[:,sysr.faults]` is the fault to residual channel of the filter 
with the corresponding transfer function matrix `Rf(λ)` and 
`sysr.sys[:,sysr.noise]` is the noise to residual channel of the filter 
with the corresponding transfer function matrix `Rw(λ)`.   
If  `Rff(λ)` is the transfer function matrix formed of those `j`-th columns of `Rf(λ)` 
for which `SFDI[j] = 1` and `Rdf(λ)` is the transfer function matrix formed of 
those `j`-th columns of `Rf(λ)` for which `SFDI[j] = false`, then:   
`β` is the H∞- index of `Rff(λ)`, `γ` is the H∞-norm of `[Rdf(λ) Rw(λ)]` and 
`gap` is the fault-to-noise gap evaluated as `gap := β/γ`. 
If `freq` is a vector of real frequency values, then `β` and `γ`
are evaluated over the frequencies contained in `freq`. 
`gap = ∞` if `[Rdf(λ) Rw(λ)] = 0` and `gap = 0` if `Rff(λ) = 0`.
"""
function fdif2ngap(sysr::FDFilterIF{T}, SFDI::Union{BitVector,Vector{Bool}}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
   inpf = sysr.faults
   mf = length(inpf)
   mw = length(sysr.noise)

   mf+mw == 0 && (return T[], T[], T[])

   isstable(sysr.sys) || error("the system is unstable")
   
   β = mf == 0 ? 0 : fdhinfminus(sysr.sys[:,inpf[SFDI]],freq)[1]
   γ = ghinfnorm(sysr.sys[:,[inpf[.!SFDI]; sysr.noise]])[1]
   return β/γ, β, γ
end
"""
      fdif2ngap(sysr::FDFilterIF, SFDI, freq; atol = √ϵ) -> (gap, β, γ)

Compute the fault-to-noise gap 
`gap` (and the related quatities `β` and `γ`) for the `q × mf` binary structure matrix `SFDI` associated to
the stable transfer function matrices `Rf(λ)` and `Rw(λ)` of the
transfer channels from the fault inputs to residuals and noise inputs to residuals, respectively,  of the
fault detection filter internal form object `sysr::FDFilterIF`. 
The `i`-th element of the vectors `gap`, `β` and `γ` contain the quantities: 
`β[i]` - the H∞- index of `Rffi(λ)`, `γ[i]` - the H∞-norm of `[Rdfi(λ) Rwi(λ)]` and 
the fault-to-noise gap `gap` evaluated as `gap[i] := β[i]/γ[i]`, where
`Rffi(λ)` is formed of those `j`-th  columns of the `i`-th row of `Rf(λ)` for which `S[i,j] = 1` and 
`Rdfi(λ)` is formed of those `j`-th  columns of the `i`-th row of `Rf(λ)` for which `S[i,j] = 0`. 
`gap[i] = ∞` if `[Rdfi(λ) Rwi(λ)] = 0` and `gap[i] = 0` if `Rffi(λ) = 0`.
If `freq` is a vector of real frequency values, then `β[i]` and `γ[i]`
are evaluated over the frequencies contained in `freq`. 
`atol` is an absolute tolerance for the norms `Rwi(λ)`, such that norm values less than or equal to `atol` are 
considered zero (default:  `√ϵ`, where `ϵ` is the working machine precision.)
"""
function fdif2ngap(sysr::FDFilterIF{T}, SFDI::Union{BitMatrix,Matrix{Bool}}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing;
                    atol::Real = sqrt(eps(T))) where T
    q = size(sysr.sys,1)   
    inpf = sysr.faults
    mf = length(inpf)
    (q, mf ) == size(SFDI) || error("missmatch between the number of fault inputs and column dimension of SFDI")
    mw = length(sysr.noise)
 
    mf+mw == 0 && (return zeros(T,q,0), zeros(T,q,0), zeros(T,q,0))
 
    isstable(sysr.sys) || error("the system is unstable")
    β = similar(Array{T,1},q) 
    γ = similar(Array{T,1},q) 
    for i = 1:q
        β[i] = mf == 0 ? 0 : fdhinfminus(sysr.sys[i,inpf[SFDI[i,:]]],freq)[1]
        t = ghinfnorm(sysr.sys[i,[inpf[.!SFDI[i,:]]; sysr.noise]])[1]
        γ[i] = t <= atol ? 0 : t
    end
    return β ./ γ, β, γ
 end
 """
     fdif2ngap(sysr::FDIFilterIF, SFDI, freq; atol = √ϵ) -> (gap, β, γ)

Compute the fault-to-noise gap 
`gap` (and the related quatities `β` and `γ`) for the `N × mf` binary structure matrix `SFDI` associated to
the stable global transfer function matrices `Rf(λ)` and `Rw(λ)` of the
transfer channels from the fault inputs to residuals and noise inputs to residuals, respectively,  of the
fault detection and isolation filter internal form object `sysr::FDIFilterIF`. 
The filter `sysr` consists of `N` individual FDI filters `sysr.sys[i]`, for `i = 1, ..., N`, where
`sysr.sys[i][:,sysr.faults]` is the fault to residual channel of the `i`-th filter 
with the corresponding transfer function matrix `Rfi(λ)` and `sysr.sys[i][:,sysr.noise]` is the noise to residual channel of the `i`-th filter 
with the corresponding transfer function matrix `Rwi(λ)`.   
The global transfer function matrices `Rf(λ)` and `Rw(λ)` are formed by row concatenation of the 
transfer function matrices of the `N` individual filters, i.e.,  `Rf(λ) := [ Rf1(λ); Rf2(λ); ...; RfN(λ)]` and 
`Rw(λ) := [ Rw1(λ); Rw2(λ); ...; RwN(λ)]` 
Let  `Rffi(λ)` be the transfer function matrix formed of those `j`-th columns of `Rfi(λ)` 
for which `SFDI[i,j] = 1` and let `Rdfi(λ)` be the transfer function matrix formed of those `j`-th columns of `Rfi(λ)` 
for which `SFDI[i,j] = 0`. 
The `i`-th element of the vectors `gap`, `β` and `γ` contain the quantities: 
`β[i]` - the H∞- index of `Rffi(λ)`, `γ[i]` - the H∞-norm of `[Rdfi(λ) Rwi(λ)]` and 
the fault-to-noise gap `gap` evaluated as `gap[i] := β[i]/γ[i]`. 
`gap[i] = ∞` if `[Rdfi(λ) Rwi(λ)] = 0` and `gap[i] = 0` if `Rffi(λ) = 0`.
If `freq` is a vector of real frequency values, then `β[i]` and `γ[i]`
are evaluated over the frequencies contained in `freq`. 
`atol` is an absolute tolerance for the norms `Rwi(λ)`, such that norm values less than or equal to `atol` are 
considered zero (default:  `√ϵ`, where `ϵ` is the working machine precision.)
"""
function fdif2ngap(sysr::FDIFilterIF{T}, SFDI::Union{BitMatrix,BitVector,Array{Bool,2},Array{Bool,1}}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing; 
                   atol::Real = sqrt(eps(T))) where T
    N = length(sysr.sys)
    isa(SFDI,Union{BitVector,Array{Bool,1}}) ? (nb = 1; mf = length(SFDI); SFDI = reshape(SFDI,nb,mf) ) :
                                               (nb = size(SFDI,1); mf = size(SFDI,2))
    nb == N || error("missmatch between the number of filters $N and the number of specifications $nb")
    inpf = sysr.faults
    mf == length(inpf) || error("missmatch between the number of faults and the column dimension/length of SFDI")
    mw = length(sysr.noise)
    β = mf > 0 ? similar(Array{T,1},N) : zeros(T,N)
    γ = mf+mw > 0 ? similar(Array{T,1},N) : zeros(T,N)
    if mf > 0
        for i = 1:N
            β[i] = fdhinfminus(sysr.sys[i][:,inpf[view(SFDI,i,:)]],freq)[1]
        end
    end
    if mf+mw > 0
        for i = 1:N
            t = ghinfnorm(sysr.sys[i][:,[inpf[.!SFDI[i,:]]; sysr.noise]])[1]
            γ[i] = t <= atol ? 0 : t
        end
    end
    
    return β./γ, β, γ
end
"""
     γ = fdimmperf(sysr::FDFilterIF[, nrmflag]) 

Compute the model-matching performance `γ` of the fault detection filter internal form object `sysr::FDFilterIF`. 
If `Rw(λ)` is the transfer function matrix of the transfer channel from the noise inputs to residuals 
`sysr.sys[:,sysr.noise]`, then
`γ` is the  H∞-norm of `Rw(λ)`, if `nrmflag = Inf` (default) and the  H2-norm of `Rw(λ)`, if `nrmflag = 2`.
The value of `γ` is infinite for an unstable filter or if `nrmflag = 2` and the transfer function matrix
`Rw(λ)` of a continuous-time system is not strictly proper.
"""    
function fdimmperf(sysr::FDFilterIF{T}, nrmflag::Real = Inf) where T

    (nrmflag == Inf || nrmflag == 2) || error("only H∞- and H2-norms supported")
 
    length(sysr.noise) == 0 && (return T(0))
 
    return isinf(nrmflag) ? ghinfnorm(sysr.sys[:,sysr.noise])[1] : 
                            gh2norm(sysr.sys[:,sysr.noise])
end
"""
     γ = fdimmperf(sysr::FDFilterIF, SFDI[, nrmflag]) 

Compute the model-matching performance `γ` of the fault detection filter internal form object `sysr::FDFilterIF`
for a given binary structure vector `SFDI`. If `Rf(λ)` is the transfer function matrix of the 
transfer channel from the fault inputs to residuals `sysr.sys[:,sysr.faults]` and 
`Rw(λ)` is the transfer function matrix of the transfer channel from the noise inputs to residuals 
`sysr.sys[:,sysr.noise]`, then
`γ` is the  H∞-norm of `[Rdf(λ) Rw(λ)]`, if `nrmflag = Inf` (default) and the  H2-norm of `[Rdf(λ) Rw(λ)]`, 
if `nrmflag = 2`, where `Rdf(λ)` is the transfer function matrix formed by those `j`-th columns of `Rf(λ)` for which 
`SFDI[j] = 0`.  
The value of `γ` is infinite for an unstable filter or if `nrmflag = 2` and the transfer function matrix
`[Rdf(λ) Rw(λ)]` of a continuous-time system is not strictly proper.
"""    
function fdimmperf(sysr::FDFilterIF{T}, SFDI::Union{BitVector,Vector{Bool}}, nrmflag::Real = Inf) where T

    (nrmflag == Inf || nrmflag == 2) || error("only H∞- and H2-norms supported")
    inpf = sysr.faults
    length(SFDI) == length(inpf) || error("missmatch between the number of faults and the column dimension/length of SFDI")
 
    return isinf(nrmflag) ? ghinfnorm(sysr.sys[:,[inpf[.!SFDI]; sysr.noise]])[1] : 
                            gh2norm(sysr.sys[:,[inpf[.!SFDI]; sysr.noise]])
end
"""
     γ = fdimmperf(sysr::FDFilterIF, SFDI[, nrmflag]) 

Compute the model-matching performance `γ` of the fault detection filter internal form object `sysr::FDFilterIF`
for a given binary structure matrix `SFDI`. If `Rf(λ)` is the transfer function matrix of the 
transfer channel from the fault inputs to residuals `sysr.sys[:,sysr.faults]` and 
`Rw(λ)` is the transfer function matrix of the transfer channel from the noise inputs to residuals 
`sysr.sys[:,sysr.noise]`, then
`γ` is the  H∞-norm of `[Rdf(λ) Rw(λ)]`, if `nrmflag = Inf` (default) and the  H2-norm of `[Rdf(λ) Rw(λ)]`, 
if `nrmflag = 2`, where `Rdf(λ) = .!SFDI .* Rf(λ)` (i.e., the element-wise product of `.!SFDI` and `Rf(λ)`.
The value of `γ` is infinite for an unstable filter or if `nrmflag = 2` and the transfer function matrix
`[Rdf(λ) Rw(λ)]` of a continuous-time system is not strictly proper.
"""    
function fdimmperf(sysr::FDFilterIF{T}, SFDI::Union{BitMatrix,Matrix{Bool}}, nrmflag::Real = Inf) where T

    (nrmflag == Inf || nrmflag == 2) || error("only H∞- and H2-norms supported")

    return isinf(nrmflag) ? ghinfnorm([dssubsel(sysr.sys[:,sysr.faults], .!SFDI) sysr.sys[:,sysr.noise]])[1] : 
                            gh2norm([dssubsel(sysr.sys[:,sysr.faults], .!SFDI) sysr.sys[:,sysr.noise]])
end

"""
     γ = fdimmperf(sysr::FDIFilterIF[, nrmflag]) 

Compute the model-matching performance `γ` of the stable fault detection and isolation filter internal form object `sysr::FDIFilterIF`. 
The filter `sysr` consists of `N` individual FDI filters `sysr.sys[i]`, for `i = 1, ..., N`, where
`sysr.sys[i][:,sysr.noise]` is the noise to residual channel of the `i`-th filter 
with the corresponding transfer function matrix `Rwi(λ)`. Then, 
`γ` is an `N`-dimensional vector whose `i`-th component is the  H∞-norm of `Rwi(λ)`, if `nrmflag = Inf` (default) 
and the  H2-norm of `Rwi(λ)`, if `nrmflag = 2`.
The `i`-th component of `γ` is infinite for an unstable filter or if `nrmflag = 2` and the transfer function matrix
`Rwi(λ)` of a continuous-time system is not strictly proper.
"""    
function fdimmperf(sysr::FDIFilterIF{T}, nrmflag::Real = Inf) where T
    N = length(sysr.sys)

    (nrmflag == Inf || nrmflag == 2) || error("only H∞- and H2-norms supported")
 
    length(sysr.noise) == 0 && (return zeros(T,N))
 

    γ = similar(Array{T,1},N) 
    for i = 1:N
        γ[i] = isinf(nrmflag) ? ghinfnorm(sysr.sys[i][:,sysr.noise])[1] : 
                                gh2norm(sysr.sys[i][:,sysr.noise])
    end
    return γ
end
"""
     γ = fdimmperf(sysr::FDIFilterIF, SFDI[, nrmflag]) 

Compute the model-matching performance `γ` of the stable fault detection and isolation 
filter internal form object `sysr::FDIFilterIF` and the associated binary structure matrix `SFDI`. 
The filter `sysr` consists of `N` individual FDI filters `sysr.sys[i]`, for `i = 1, ..., N`, where
`sysr.sys[i][:,sysr.faults]` is the fault to residual channel of the `i`-th filter 
with the corresponding transfer function matrix `Rfi(λ)` and
`sysr.sys[i][:,sysr.noise]` is the noise to residual channel of the `i`-th filter 
with the corresponding transfer function matrix `Rwi(λ)`. Then, 
`γ` is an `N`-dimensional vector whose `i`-th component is the  H∞-norm of `[Rfdi(λ) Rwi(λ)]`, if `nrmflag = Inf` (default) 
or the  H2-norm of `[Rfdi(λ) Rwi(λ)]`, if `nrmflag = 2`, where `Rfdi(λ)` is the transfer function matrix whose `j`-th column is 
zero if `SFDI[i,j] = 1` and is equal to the `j`-th column of `Rfi(λ)` if `SFDI[i,j] = 0`. 
The `i`-th component of `γ` is infinite for an unstable filter or if `nrmflag = 2` and the transfer function matrix
`[Rfdi(λ) Rwi(λ)]` of a continuous-time system is not strictly proper.
"""    
function fdimmperf(sysr::FDIFilterIF{T}, SFDI::Union{BitMatrix,BitVector,Array{Bool,2},Array{Bool,1}}, nrmflag::Real = Inf) where T
    N = length(sysr.sys)

    (nrmflag == Inf || nrmflag == 2) || error("only H∞- and H2-norms supported")
    isa(SFDI,Union{BitVector,Array{Bool,1}}) ? (nb = 1; mf = length(SFDI); SFDI = reshape(SFDI,nb,mf) ) :
                                               (nb = size(SFDI,1); mf = size(SFDI,2))
    nb == N || error("missmatch between the number of filters $N and the number of specifications $nb")
    inpf = sysr.faults
    mf == length(inpf) || error("missmatch between the number of faults and the column dimension/length of SFDI")

    γ = similar(Array{T,1},N) 
    for i = 1:N
        inddfw = [inpf[.!SFDI[i,:]]; sysr.noise]
        γ[i] = isinf(nrmflag) ? ghinfnorm(sysr.sys[i][:,inddfw])[1] : 
                                gh2norm(sysr.sys[i][:,inddfw])
    end
    return γ
end
"""
     γ = fdimmperf(sysr::FDFilterIF, sysref::Union{FDFilterIF,FDIModel}[, nrmflag]) 

Compute the model-matching performance `γ` of the fault detection filter internal form object `sysr::FDFilterIF`
with respect to the fault detection reference filter internal form  `sysref::FDFilterIF`. 
If `R(λ)` is the transfer function matrix of the fault detection filter internal form   
`sysr.sys` and `Mr(λ)` is the transfer function matrix of the fault detection reference filter internal form   
`sysref.sys`, then
`γ` is the  H∞-norm of `R(λ)-Mr(λ)`, if `nrmflag = Inf` (default) or the  H2-norm of `R(λ)-Mr(λ)`, if `nrmflag = 2`.
The value of `γ` is infinite for an unstable difference `R(λ)-Mr(λ)` or if `nrmflag = 2` and the transfer function matrix
`R(λ)-Mr(λ)` of a continuous-time system is not strictly proper. In general, `R(λ)` and `Mr(λ)` are partitioned as
`R(λ) = [ Ru(λ) Rd(λ) Rf(λ) Rw(λ) Ra(λ) ]` and `Mr(λ) = [ Mru(λ) Mrd(λ) Mrf(λ) Mrw(λ) Mra(λ) ]` in accordance with 
the partitioning of the inputs in control inputs, disturbance inputs, fault inputs, noise inputs and auxiliary inputs.
Void components of `Mr(λ)` corresponding to non-void components in `R(λ)` are assumed to be zero. 
"""    
function fdimmperf(sysr::FDFilterIF{T}, sysref::Union{FDFilterIF,FDIModel}, nrmflag::Real = Inf; atolinf::Real = 0) where T

    (nrmflag == Inf || nrmflag == 2) || error("only H∞- and H2-norms supported")
    DescriptorSystems.promote_Ts(sysr.sys.Ts, sysref.sys.Ts)  # check sampling time

    inpu = sysr.controls; mu = length(inpu)  
    inpd = sysr.disturbances; md = length(inpd) 
    inpf = sysr.faults; mf = length(inpf)  
    inpw = sysr.noise;  mw = length(inpw) 
    inpaux = sysr.aux;  maux = length(inpaux)  
    mru = length(sysref.controls); 
    mru == 0 || mru == mu ||  error("Incompatible control input groups in sysr and sysref") 
    mrd = length(sysref.disturbances); 
    mrd == 0 || mrd == md ||  error("Incompatible disturbance input groups in sysr and sysref") 
    mrf = length(sysref.faults); 
    mrf == 0 || mrf == mf ||  error("Incompatible fault input groups in sysr and sysref") 
    mrw = length(sysref.noise); 
    mrw == 0 || mrw == mw ||  error("Incompatible noise input groups in sysr and sysref") 
    mra = length(sysref.aux); 
    mra > 0 || mra == maux ||  error("Incompatible auxiliary input groups in sysr and sysref") 
    
    m = mu+md+mf+mw+maux;       # total number of inputs
 
    rinp = zeros(0,m);
    mru == 0 || (rinp = [rinp; eye(mu,m)])
    mrd == 0 || (rinp = [rinp; zeros(md,mu) eye(md,m-mu)])
    mrf == 0 || (rinp = [rinp; zeros(mf,mu+md) eye(mf,m-mu-md)])
    mrw == 0 || (rinp = [rinp; zeros(mw,mu+md+mf) eye(mw,m-mu-md-mf)])
    mra == 0 || (rinp = [rinp; zeros(maux,m-maux) eye(maux)])
    return isinf(nrmflag) ? ghinfnorm(sysr.sys-sysref.sys*rinp)[1] : 
                            gh2norm(sysr.sys-sysref.sys*rinp; atolinf)
end
"""
     γ = fdimmperf(sysr::FDIFilterIF, sysref::FDIFilterIF[, nrmflag]) 

Compute the model-matching performance `γ` of the fault detection and isolation filter internal form object `sysr::FDIFilterIF`
with respect to the fault detection and isolation reference filter internal form  `sysref::FDIFilterIF`. 
If `Ri(λ)` is the transfer function matrix of the `i`-th fault detection and isolation filter internal form   
`sysr.sys[i]` and `Mri(λ)` is the transfer function matrix of the `i`-th fault detection and isolation reference filter internal form   
`sysref.sys[i]`, then `γ` is a vector whose `i`-th component
`γ[i]` is the  H∞-norm of `Ri(λ)-Mri(λ)`, if `nrmflag = Inf` (default) or the  H2-norm of `Ri(λ)-Mri(λ)`, if `nrmflag = 2`.
The value of `γ[i]` is infinite for an unstable difference `Ri(λ)-Mri(λ)` or if `nrmflag = 2` and the transfer function matrix
`Ri(λ)-Mri(λ)` of a continuous-time system is not strictly proper. In general, `Ri(λ)` and `Mri(λ)` are partitioned as
`Ri(λ) = [ Rui(λ) Rdi(λ) Rfi(λ) Rwi(λ) Rai(λ) ]` and `Mri(λ) = [ Mrui(λ) Mrdi(λ) Mrfi(λ) Mrwi(λ) Mrai(λ) ]` in accordance with 
the partitioning of the inputs in control inputs, disturbance inputs, fault inputs, noise inputs and auxiliary inputs.
Void components of `Mri(λ)` corresponding to non-void components in `Ri(λ)` are assumed to be zero. 
"""    
function fdimmperf(sysr::FDIFilterIF{T1}, sysref::FDIFilterIF{T2}, nrmflag::Real = Inf) where {T1,T2}
    N = length(sysr.sys)
    N == length(sysref.sys) || error("sysr and sysref must have the same number of component filters")
    DescriptorSystems.promote_Ts(sysr.sys[1].Ts, sysref.sys[1].Ts)  # check sampling time

    (nrmflag == Inf || nrmflag == 2) || error("only H∞- and H2-norms supported")
    inpu = sysr.controls; mu = length(inpu)  
    inpd = sysr.disturbances; md = length(inpd) 
    inpf = sysr.faults; mf = length(inpf)  
    inpw = sysr.noise;  mw = length(inpw) 
    inpaux = sysr.aux;  maux = length(inpaux)  
    mru = length(sysref.controls); 
    mru == 0 || mru == mu ||  error("Incompatible control input groups in sysr and sysref") 
    mrd = length(sysref.disturbances); 
    mrd == 0 || mrd == md ||  error("Incompatible disturbance input groups in sysr and sysref") 
    mrf = length(sysref.faults); 
    mrf == 0 || mrf == mf ||  error("Incompatible fault input groups in sysr and sysref") 
    mrw = length(sysref.noise); 
    mrw == 0 || mrw == mw ||  error("Incompatible noise input groups in sysr and sysref") 
    mra = length(sysref.aux); 
    mra > 0 || mra == maux ||  error("Incompatible auxiliary input groups in sysr and sysref") 
    
    m = mu+md+mf+mw+maux;       # total number of inputs
    rinp = zeros(0,m);
    mru == 0 || (rinp = [rinp; eye(mu,m)])
    mrd == 0 || (rinp = [rinp; zeros(md,mu) eye(md,m-mu)])
    mrf == 0 || (rinp = [rinp; zeros(mf,mu+md) eye(mf,m-mu-md)])
    mrw == 0 || (rinp = [rinp; zeros(mw,mu+md+mf) eye(mw,m-mu-md-mf)])
    mra == 0 || (rinp = [rinp; zeros(maux,m-maux) eye(maux)])

    γ = similar(Array{T1,1},N) 
    for i = 1:N
        γ[i] = isinf(nrmflag) ? ghinfnorm(sysr.sys[i]-sysref.sys[i]*rinp)[1] : 
                                gh2norm(sysr.sys[i]-sysref.sys[i]*rinp)
    end
    return γ
end

