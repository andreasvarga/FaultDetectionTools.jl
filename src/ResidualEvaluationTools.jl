"""
    tstep!(fdsys::FDSystem, y::Vector, u::Vector; init = false)

Perform a time step on the FDD system `fdsys` using the output measurements `y` and control inputs `u`. 
If `init = true`, a reinitialization of filter state vector, evaluation vector and time value is performed first.  

The following quantities are updated:
- `fdsys.r`  - the residual filter output vector `r`
- `fdsys.x`  - the residual filter state vector `x` 
- `fdsys.θ`  - the residual evaluation signal vector `θ`
- `fdsys.re` - the residual evaluation filter state `re`
- `fdsys.t`  - time value `t`

Using the threshold values `τ` in `fdsys.τ` and the updated values of the residual evaluation signal vector `θ`,
the fault detection status vector `fdsys.isig` is updated as follows:
- for fault detection: `fdsys.isig[1] = 1` if `θ[1] ≥ τ[1]`; otherwise `fdsys.isig[1] = 0`;
- for fault isolation: the `i`-th component is set as `fdsys.isig[i] = 1` if `θ[i] ≥ τ[i]`; otherwise `fdsys.isig[i] = 0`.

Based on the updated values of the fault detection status vector `fdsys.isig`, 
the decison vector `fdsys.indfault` is set as follows: 
 - for fault detection: `fdsys.indfault = [1]` if `fdsys.isig[1] = 1`; otherwise `fdsys.indfault = Int[]`;
 - for weak fault isolation: `fdsys.indfault = [k]` if `fdsys.isig` matches the `k`-th column of `fdsys.SFDI`;
 - for strong fault isolation: `fdsys.indfault` contains the indices of nonzero elements of `fdsys.isig`. 
"""
function tstep!(fds::FDSystem, y::Vector, u::Vector; init = false)
    if init
       fds.x .= fds.x0
       fds.re .= 0 .* fds.re
       fds.t = fds.t0
    end
    # residual generation block
    fds.r = fds.C*fds.x + fds.Dy*y + fds.Du*u
    fds.x = fds.A*fds.x + fds.By*y + fds.Bu*u
    # evaluation block
    Ne = size(fds.SFDI,1)
    if Ne == 1
       rnorm = norm(fds.r)
       fds.θ[1] = fds.α[1]*rnorm
       if fds.β[1] != 0
          fds.θ[1] += fds.β[1]*fds.re[1] 
          fds.re[1] = fds.γ[1]*fds.re[1] + rnorm
       end
       # decision variable
       fds.isig[1] = Int(fds.θ[1] >= fds.τ[1])
       #isempty(fds.indfault) && (fds.indfault = fds.isig[1] > 0 ? [1] : Int[])
       fds.indfault = fds.isig[1] > 0 ? [1] : Int[]
    else
       rnorm = abs.(fds.r)
       fds.θ .= fds.α.*rnorm
       for i = 1:Ne
           if fds.β[i] != 0
              fds.θ[i] += fds.β[i]*fds.re[i] 
              fds.re[i] = fds.γ[i]*fds.re[i] + rnorm[i]
           end
           # decision variables
           fds.isig[i] = Int(fds.θ[i] >= fds.τ[i])
       end
       if fds.strongfdi 
          tind = findall(!iszero,fds.isig)
          length(fds.indfault) < length(tind) && (fds.indfault = tind)
       else
          ind = findfirst(!iszero,([fds.SFDI[:,i] == fds.isig for i in 1:Ne] .== 1))
          isempty(fds.indfault) && (fds.indfault = isnothing(ind) ? Int[] : [ind])
       end
    end
    fds.t += fds.Ts
    return fds
end

"""
    tstep!(fdisys::FDISystem, y::Vector, u::Vector; init = false)

Perform a time step on the FDD system `fdisys` using the output measurements `y` and control inputs `u`. 
If `init = true`, a reinitialization of the component filter state vectors, evaluation vectors and time values is performed first.  

If `Ne` is the number of component filters, the following quantities are updated for each of the `i = 1, ..., Ne` component filters:
- `fdisys.fdsys[i].r`  - the `i`-th residual filter output vector `r`
- `fdisys.fdsys[i].x`  - the `i`-th residual filter state vector `x` 
- `fdisys.fdsys[i].θ`  - the `i`-th residual evaluation signal vector `θ`
- `fdisys.fdsys[i].re` - the `i`-th residual evaluation filter state `re`
- `fdisys.fdsys[i].t`  - `i`-th filter time value `t`

Using the threshold value `τ[i] = fdisys.fdsys[i].τ[1]` and the updated value `θ[i] = fdisys.fdsys[i].θ[1]` of the residual evaluation signal,
the fault detection status vector `fdisys.isig` is updated as follows:
the `i`-th component is set as `fdisys.isig[i] = 1` if `θ[i] ≥ τ[i]`; otherwise `fdisys.isig[i] = 0`.

Based on the updated values of the fault detection status vector `fdisys.isig`, 
the decison vector `fdisys.indfault` is set as follows: 
 - for weak fault isolation: `fdisys.indfault = [k]` if `fdisys.isig` matches the `k`-th column of `fdisys.SFDI`;
 - for strong fault isolation: `fdisys.indfault` contains the indices of nonzero elements of `fdisys.isig`. 
"""
function tstep!(fds::FDISystem, y, u; init = false)
    Ne = length(fds.fdsys)
    if init
       for i = 1:Ne
           fds.fdsys[i].x .= fds.fdsys[i].x0
           fds.fdsys[i].re .= 0 .* fds.fdsys[i].re
           fds.fdsys[i].t = fds.fdsys[i].t0
       end
    end
    fds.t += fds.Ts
    for i in 1:Ne
        tstep!(fds.fdsys[i], y, u)
        ind = fds.fdsys[i].indfault
        fds.isig[i] = isempty(ind) ? 0 : 1
        fds.θ[i] = fds.fdsys[i].θ[1]
    end
    if fds.strongfdi 
       fds.indfault = findall(!iszero,fds.isig)
    else
       ind = findfirst(!iszero,([fds.SFDI[:,i] == fds.isig for i in 1:size(fds.SFDI,2)] .== 1))
       fds.indfault = isnothing(ind) ? Int[] : [ind]
    end
    return fds
end
"""
    tstep!(mfddsys::MDSystem, y::Vector, u::Vector; init = false)

Perform a time step on the FDD system `mfddsys` using the output measurements `y` and control inputs `u`. 
If `init = true`, a reinitialization of the component filter state vectors, evaluation vectors and time values is performed first.  

If `N` is the number of component filters, the following quantities are updated for each of the `i = 1, ..., N` component filters:
- `mfddsys.mdsys[i].r`  - the `i`-th residual filter output vector `r`
- `mfddsys.mdsys[i].x`  - the `i`-th residual filter state vector `x` 
- `mfddsys.mdsys[i].θ`  - the `i`-th residual evaluation signal vector `θ`
- `mfddsys.mdsys[i].re` - the `i`-th residual evaluation filter state `re`
- `mfddsys.mdsys[i].t`  - `i`-th filter time value `t`

Using the threshold value `τ[i] = mfddsys.mdsys[i].τ[1]` and the updated value `θ[i] = mfddsys.mdsys[i].θ[1]` of the residual evaluation signal,
the model detection status vector `mfddsys.isig` is updated as follows:
the `i`-th component is set as `mfddsys.isig[i] = 1` if `θ[i] ≥ τ[i]`; otherwise `mfddsys.isig[i] = 0`.

Based on the updated values of the model detection status vector `mfddsys.isig`, 
the decison variable `mfddsys.indmodel` is set to the index of the currently 
detected model and is set to zero if no signature match occurred. `mfddsys.indminim` provides the index of best matched model, 
corresponding to the least component of the evaluation vector.
"""
function tstep!(mds::MDSystem, y, u; init = false)
    Ne = length(mds.mdsys)
    if init
       for i = 1:Ne
           mds.mdsys[i].x .= mds.mdsys[i].x0
           mds.mdsys[i].re .= 0 .* mds.mdsys[i].re
           mds.mdsys[i].t = mds.mdsys[i].t0
       end
    end
    mds.t += mds.Ts
    indminim = 0
    θmin = Inf
    for i in 1:Ne
        tstep!(mds.mdsys[i], y, u)
        ind = mds.mdsys[i].indfault
        mds.isig[i] = isempty(ind) ? 0 : 1
        temp = mds.mdsys[i].θ[1]
        if θmin > temp
           indminim = i
           θmin = temp
        end
        mds.θ[i] = temp
    end
    mds.indmodel = count(iszero,mds.isig) == 1 ? findfirst(iszero,mds.isig) : 0
    mds.indminim = indminim
    return mds
end





