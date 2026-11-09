# Output, time stepping and the run loop.

"""
    ModelOutput{T}(n)

Time series of one run, one vector of length `n` per variable. `t` is the time of day
in hours. `wCO2*` are in mg CO2 m-2 s-1. `wCO2M`, `z0m` and `z0h` are not written
and stay zero.
"""
struct ModelOutput{T<:AbstractFloat}
    t::Vector{T}
    h::Vector{T}
    theta::Vector{T}
    thetav::Vector{T}
    dtheta::Vector{T}
    dthetav::Vector{T}
    wtheta::Vector{T}
    wthetav::Vector{T}
    wthetae::Vector{T}
    wthetave::Vector{T}
    q::Vector{T}
    dq::Vector{T}
    wq::Vector{T}
    wqe::Vector{T}
    wqM::Vector{T}
    qsat::Vector{T}
    e::Vector{T}
    esat::Vector{T}
    CO2::Vector{T}
    dCO2::Vector{T}
    wCO2::Vector{T}
    wCO2A::Vector{T}
    wCO2R::Vector{T}
    wCO2e::Vector{T}
    wCO2M::Vector{T}
    u::Vector{T}
    du::Vector{T}
    uw::Vector{T}
    v::Vector{T}
    dv::Vector{T}
    vw::Vector{T}
    T2m::Vector{T}
    q2m::Vector{T}
    u2m::Vector{T}
    v2m::Vector{T}
    e2m::Vector{T}
    esat2m::Vector{T}
    thetasurf::Vector{T}
    thetavsurf::Vector{T}
    qsurf::Vector{T}
    ustar::Vector{T}
    z0m::Vector{T}
    z0h::Vector{T}
    Cm::Vector{T}
    Cs::Vector{T}
    L::Vector{T}
    Rib::Vector{T}
    Swin::Vector{T}
    Swout::Vector{T}
    Lwin::Vector{T}
    Lwout::Vector{T}
    Q::Vector{T}
    ra::Vector{T}
    rs::Vector{T}
    H::Vector{T}
    LE::Vector{T}
    LEliq::Vector{T}
    LEveg::Vector{T}
    LEsoil::Vector{T}
    LEpot::Vector{T}
    LEref::Vector{T}
    G::Vector{T}
    zlcl::Vector{T}
    RH_h::Vector{T}
    ac::Vector{T}
    M::Vector{T}
    dz::Vector{T}
end

ModelOutput{T}(n::Integer) where {T} =
    ModelOutput{T}(ntuple(_ -> zeros(T, n), fieldcount(ModelOutput{T}))...)

"Write the current state into row `s.t + 1` of `out`."
function store!(out::ModelOutput, s::ModelState)
    i = s.t + 1
    x, c = s.input, s.constants
    ml, w, rad, sl, ls, cu = s.mixed_layer, s.wind, s.radiation, s.surface_layer, s.land_surface, s.cumulus
    fac = (c.rho * c.mco2) / c.mair     # ppm m s-1 to mg CO2 m-2 s-1

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
    out.wCO2[i] = ml.wCO2 * fac
    out.wCO2e[i] = ml.wCO2e * fac
    out.wCO2R[i] = ml.wCO2R * fac
    out.wCO2A[i] = ml.wCO2A * fac

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
    step!(s::ModelState, out::ModelOutput)

Advance the model by one time step: run the components, store the output, then
integrate the land surface and the mixed layer.
"""
function step!(s::ModelState, out::ModelOutput)
    x = s.input
    statistics!(s)
    x.sw_rad && radiation!(s)
    x.sw_sl && surface_layer!(s)
    x.sw_ls && land_surface!(s)
    x.sw_cu && cumulus!(s)
    x.sw_ml && mixed_layer!(s)
    store!(out, s)
    x.sw_ls && integrate_land_surface!(s)
    x.sw_ml && integrate_mixed_layer!(s)
    s.t += 1
    return s
end

"""
    run_model(input; constants = PhysicalConstants{T}()) -> ModelOutput

Run the model for `floor(runtime / dt)` steps and return the time series.
"""
function run_model(x::ModelInput{T}; kwargs...) where {T}
    s = init_model(x; kwargs...)
    out = ModelOutput{T}(floor(Int, x.runtime / x.dt))
    for _ in eachindex(out.t)
        step!(s, out)
    end
    return out
end
