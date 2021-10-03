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
     fdscond(sysr::FDFilterIF,freq) -> (scond, β, γ)

Compute for the stable transfer function matrix `Rf(λ)` of the
transfer channel from the fault inputs to residuals of the
fault detection filter internal form object `sysr::FDFilterIF` the quantities: 
`β` - the H∞- index of `Rf(λ)`, `γ` - the maximum of the columns norms of `Rf(λ)` and 
the fault detection sensitivity condition `scond` evaluated as `scond := β/γ`. 
If `freq` is a vector of real frequency values, then `β` and `γ`
are evaluated over the frequencies contained in `freq`. 
"""
function fdscond(sysr::FDFilterIF{T}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
    fdscond_(sysr.sys[:,sysr.faults], freq)
end
"""
     fdscond(sysr::FDIFilterIF, SFDI, freq) -> (scond, β, γ)

Compute the detection and isolation sensitivity condition
`scond` (and related quatities `β` and `γ`) for the `N × mf` structure matrix `SFDI` associated to
the stable global transfer function matrix `Rf(λ)` of the
transfer channel from the fault inputs to residuals of the
fault detection and isolation filter internal form object `sysr::FDIFilterIF`. 
The filter `sysr` consists of `N` individual FDI filters `sysr.sys[i]`, for `i = 1, ..., N`, where
`sysr.sys[i][:,sysr.faults]` is the fault to residual channel of the `i`-th filter and 
`Rfi(λ)` is the corresponding transfer function matrix. 
The global transfer function matrix `Rf(λ)` is formed by row concatenation of the 
transfer function matrices of the `N` individual filters, i.e.,  `Rf(λ) := [ Rf1(λ); Rf2(λ); ...; RfN(λ)]`.
It is assumed that for each `j` such that `SFDI[i,j] = true`, the `j`-th column of `Rfi(λ)` is nonzero  and 
for each `j` such that `SFDI[i,j] = false`, the `j`-th column of `Rfi(λ)` is zero. 
The i-th element of the vectors `scond`, `β` and `γ` contain the quantities: 
`β[i]` - the H∞- index of the nonzero columns of `Rfi(λ)`, `γ` - the maximum of the nonzero columns norms of `Rf(λ)` and 
the fault detection sensitivity condition `scond` evaluated as `scond[i] := β[i]/γ[i]`. 
If `freq` is a vector of real frequency values, then `β[i]` and `γ[i]`
are evaluated over the frequencies contained in `freq`. 
"""
function fdscond(sysr::FDIFilterIF{T}, SFDI::Union{BitMatrix,BitVector,Array{Bool,2},Array{Bool,1}}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T
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
        scond[i], β[i], γ[i] = fdscond_(sysr.sys[i][:,inpf[view(SFDI,i,:)]], freq)
    end
    return scond, β, γ
end
fdscond(sysr::FDIFilterIF{T}, freq::Union{AbstractVector{<:Real},Real,Missing} = missing) where T = fdscond(sysr,trues(1,length(sysr.faults)),freq)

"""
     fdif2ngap(sysr::FDFilterIF, freq) -> (gap, β, γ)

Compute for the stable transfer function matrices `Rf(λ)` and `Rw(λ)` of the
transfer channels from the fault inputs and noise inputs to residuals, respectively, of the
fault detection filter internal form object `sysr::FDFilterIF` the quantities:      
`β` - the H∞- index of `Rf(λ)`, `γ` - the H∞-norm of `Rw(λ)` and 
`gap` - the fault-to-noise gap evaluated as `gap := β/γ`. 
If `freq` is a vector of real frequency values, then `β` and `γ`
are evaluated over the frequencies contained in `freq`. 
`gap = ∞` if there are no noise inputs and `gap = 0` if there are no fault inputs.

If `R(λ)` is the input-output form of the fault detection filter internal form `sysr.sys`,
then

       r = Rf(λ)*f + Rw(λ)*w + Rv(λ)*v ,                           

with the Laplace- or Z-transformed residual outputs `r`, fault inputs `f`, 
noise inputs `w`, and auxiliary inputs `v`, and with `Rf(λ)`, `Rw(λ)` and `Rv(λ)` the  
corresponding transfer function matrices. The inputs `f` and `w` of `sysr.sys` 
correspond to the input groups named 'sysr.faults' and 'sysr.noise', respectively.
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
