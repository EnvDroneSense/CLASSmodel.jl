# Output, time stepping and the run loop.

# Names of the output variables, in the order of the columns.
const OUTPUT_NAMES = (
    :t, :h, :theta, :thetav, :dtheta, :dthetav, :wtheta, :wthetav, :wthetae, :wthetave, :q,
    :dq, :wq, :wqe, :wqM, :qsat, :e, :esat, :CO2, :dCO2, :wCO2, :wCO2A, :wCO2R, :wCO2e,
    :wCO2M, :u, :du, :uw, :v, :dv, :vw, :T2m, :q2m, :u2m, :v2m, :e2m, :esat2m, :thetasurf,
    :thetavsurf, :qsurf, :ustar, :z0m, :z0h, :Cm, :Cs, :L, :Rib, :Swin, :Swout, :Lwin,
    :Lwout, :Q, :ra, :rs, :H, :LE, :LEliq, :LEveg, :LEsoil, :LEpot, :LEref, :G, :zlcl,
    :RH_h, :ac, :M, :dz
)

"""
    ModelOutput{T}

Time series of one run: a `NamedTuple` with one vector per variable (a table in the sense
of Tables.jl). `t` is the time of day in hours. `wCO2*` are in mg CO2 m-2 s-1.
"""
const ModelOutput{T} = NamedTuple{OUTPUT_NAMES, NTuple{length(OUTPUT_NAMES), Vector{T}}}

function empty_output(::Type{T}, n::Integer) where {T}
    ModelOutput{T}(ntuple(_ -> zeros(T, n), length(OUTPUT_NAMES)))
end

"Write the current state into row `s.t + 1` of `out`."
function store!(out::ModelOutput, s::ModelState)
    i = s.t + 1
    x = s.input
    ml, w, rad, sl, ls, cu = s.mixed_layer,
    s.wind, s.radiation, s.surface_layer, s.land_surface, s.cumulus
    f = mg_to_ppm(s.constants)

    out.t[i] = s.t * x.dt / 3600 + x.tstart
    out.h[i] = ml.h

    out.theta[i] = ml.theta
    out.thetav[i] = ml.thetav
    out.dtheta[i] = ml.dtheta
    out.dthetav[i] = ml.dthetav
    out.wtheta[i] = ml.wtheta
    out.wthetav[i] = ml.wthetav
    out.wthetae[i] = ml.wthetae
    out.wthetave[i] = ml.wthetave

    out.q[i] = ml.q
    out.dq[i] = ml.dq
    out.wq[i] = ml.wq
    out.wqe[i] = ml.wqe
    out.wqM[i] = cu.wqM
    out.qsat[i] = ml.qsat
    out.e[i] = ml.e
    out.esat[i] = ml.esat

    out.CO2[i] = ml.CO2
    out.dCO2[i] = ml.dCO2
    out.wCO2[i] = ml.wCO2 / f
    out.wCO2e[i] = ml.wCO2e / f
    out.wCO2R[i] = ml.wCO2R / f
    out.wCO2A[i] = ml.wCO2A / f
    out.wCO2M[i] = cu.wCO2M / f

    out.u[i] = w.u
    out.du[i] = w.du
    out.uw[i] = sl.uw
    out.v[i] = w.v
    out.dv[i] = w.dv
    out.vw[i] = sl.vw

    out.T2m[i] = sl.T2m
    out.q2m[i] = sl.q2m
    out.u2m[i] = sl.u2m
    out.v2m[i] = sl.v2m
    out.e2m[i] = sl.e2m
    out.esat2m[i] = sl.esat2m

    out.thetasurf[i] = sl.thetasurf
    out.thetavsurf[i] = sl.thetavsurf
    out.qsurf[i] = sl.qsurf
    out.ustar[i] = sl.ustar
    out.z0m[i] = x.z0m
    out.z0h[i] = x.z0h
    out.Cm[i] = sl.Cm
    out.Cs[i] = sl.Cs
    out.L[i] = sl.L
    out.Rib[i] = sl.Rib

    out.Swin[i] = rad.Swin
    out.Swout[i] = rad.Swout
    out.Lwin[i] = rad.Lwin
    out.Lwout[i] = rad.Lwout
    out.Q[i] = rad.Q

    out.ra[i] = sl.ra
    out.rs[i] = ls.rs
    out.H[i] = ls.H
    out.LE[i] = ls.LE
    out.LEliq[i] = ls.LEliq
    out.LEveg[i] = ls.LEveg
    out.LEsoil[i] = ls.LEsoil
    out.LEpot[i] = ls.LEpot
    out.LEref[i] = ls.LEref
    out.G[i] = ls.G

    out.zlcl[i] = ml.lcl
    out.RH_h[i] = ml.RH_h

    out.ac[i] = cu.ac
    out.M[i] = cu.M
    out.dz[i] = cu.dz_h
    return out
end

"""
    compute!(s::ModelState)

Run the components for the current time: the diagnostics, then radiation, surface layer,
land surface, cumulus and mixed layer, as far as they are switched on.
"""
function compute!(s::ModelState)
    x = s.input
    statistics!(s)
    x.sw_rad && radiation!(s)
    x.sw_sl && surface_layer!(s)
    x.sw_ls && land_surface!(s)
    x.sw_cu && cumulus!(s)
    x.sw_ml && mixed_layer!(s)
    return s
end

"Advance the land surface and the mixed layer by one time step from their tendencies."
function integrate!(s::ModelState)
    x = s.input
    x.sw_ls && integrate_land_surface!(s)
    x.sw_ml && integrate_mixed_layer!(s)
    s.t += 1
    return s
end

"""
    step!(s::ModelState)

Advance the model by one time step: compute the components, then integrate the
land surface and the mixed layer. The diagnostics of the step are left in `s`.
"""
function step!(s::ModelState)
    compute!(s)
    integrate!(s)
    return s
end

"""
    run_model(input; constants = PhysicalConstants{T}()) -> ModelOutput

Run the model for `floor(runtime / dt)` steps and return the time series. Each row holds
the state before the time step is integrated.
"""
function run_model(x::ModelInput{T}; kwargs...) where {T}
    s = init_model(x; kwargs...)
    out = empty_output(T, floor(Int, x.runtime / x.dt))
    for _ in eachindex(out.t)
        compute!(s)
        store!(out, s)
        integrate!(s)
    end
    return out
end
