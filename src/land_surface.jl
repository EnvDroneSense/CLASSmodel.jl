# Land surface: surface energy balance with a skin temperature and a force-restore soil.

"""
    jarvis_stewart!(s::ModelState)

Surface resistance `rs` from the Jarvis-Stewart scheme: the minimum resistance divided by
the leaf area index, times factors for radiation, soil water, vapor pressure deficit and
temperature.
"""
function jarvis_stewart!(s::ModelState)
    x, ml = s.input, s.mixed_layer

    f1 = x.sw_rad ?
         1 / min(
        1, (0.004 * s.radiation.Swin + 0.05) /
           (0.81 * (0.004 * s.radiation.Swin + 1))) : one(ml.theta)
    f2 = x.w2 > x.wwilt ? (x.wfc - x.wwilt) / (x.w2 - x.wwilt) : 1e8
    f2 = max(f2, 1)   # no reduction when w2 is above field capacity
    f3 = 1 / exp(-x.gD * (ml.esat - ml.e) / 100)
    f4 = 1 / (1 - 0.0016 * (298 - ml.theta)^2)

    s.land_surface.rs = x.rsmin / x.LAI * f1 * f2 * f3 * f4
    return nothing
end

"""
    land_surface!(s::ModelState)

Aerodynamic and surface resistances, the skin temperature `Ts` from the linearised energy
balance, the surface fluxes (written to the mixed layer as `wtheta` and `wq`) and the
tendencies of the soil temperature, soil moisture and wet-leaf water layer.
"""
function land_surface!(s::ModelState)
    x, c, ml, wind, sl, ls = s.input, s.constants, s.mixed_layer, s.wind, s.surface_layer,
    s.land_surface
    rho, cp, Lv = c.rho, c.cp, c.Lv

    ueff = sqrt(wind.u^2 + wind.v^2 + ml.wstar^2)
    sl.ra = x.sw_sl ? (sl.Cs * ueff)^-1 : ueff / max(1e-3, sl.ustar)^2

    # thermodynamic variables at the mixed-layer temperature
    ml.esat = esat(ml.theta)
    ml.qsat = qsat(ml.theta, x.Ps)
    dqsatdT = ε * desat(ml.theta) / x.Ps
    ml.e = ml.q * x.Ps / ε

    x.ls_type == :js ? jarvis_stewart!(s) : ags!(s)

    # soil resistance from the top-layer moisture
    f2 = ls.wg > x.wwilt ? (x.wfc - x.wwilt) / (ls.wg - x.wwilt) : 1e8
    rssoil = x.rssoilmin * f2

    cliq = min(1, ls.Wl / (x.LAI * x.Wmax))   # wet fraction of the vegetation
    ra = sl.ra

    # latent heat per unit humidity deficit of the dry vegetation, the soil and the wet vegetation
    kveg = x.cveg * (1 - cliq) * rho * Lv / (ra + ls.rs)
    ksoil = (1 - x.cveg) * rho * Lv / (ra + rssoil)
    kliq = x.cveg * cliq * rho * Lv / ra
    kLE = kveg + ksoil + kliq

    # skin temperature from the linearised energy balance
    ls.Ts = (s.radiation.Q + rho * cp / ra * ml.theta +
             kLE * (dqsatdT * ml.theta - ml.qsat + ml.q) +
             x.Lambda * ls.Tsoil) / (rho * cp / ra + kLE * dqsatdT + x.Lambda)

    dTs = ls.Ts - ml.theta
    dq = dqsatdT * dTs + ml.qsat - ml.q
    ls.LEveg = kveg * dq
    ls.LEliq = kliq * dq
    ls.LEsoil = ksoil * dq
    ls.Wltend = -ls.LEliq / (c.rhow * Lv)

    ls.LE = ls.LEsoil + ls.LEveg + ls.LEliq
    ls.H = rho * cp / ra * dTs
    ls.G = x.Lambda * (ls.Ts - ls.Tsoil)
    avail = dqsatdT * (s.radiation.Q - ls.G) + rho * cp / ra * (ml.qsat - ml.q)
    ls.LEpot = avail / (dqsatdT + cp / Lv)
    ls.LEref = avail / (dqsatdT + cp / Lv * (1 + x.rsmin / x.LAI / ra))

    # force-restore soil
    CG = x.CGsat * (x.wsat / x.w2)^(x.b / (2 * log(10)))
    ls.Tsoiltend = CG * ls.G - 2π / DAY * (ls.Tsoil - x.T2)

    d1 = 0.1
    C1 = x.C1sat * (x.wsat / ls.wg)^(x.b / 2 + 1)
    C2 = x.C2ref * (x.w2 / (x.wsat - x.w2))
    wgeq = x.w2 - x.wsat * x.a * ((x.w2 / x.wsat)^x.p * (1 - (x.w2 / x.wsat)^(8 * x.p)))
    ls.wgtend = -C1 / (c.rhow * d1) * ls.LEsoil / Lv - C2 / DAY * (ls.wg - wgeq)

    # kinematic surface fluxes for the mixed layer
    ml.wtheta = ls.H / (rho * cp)
    ml.wq = ls.LE / (rho * Lv)
    return nothing
end

"Forward Euler step of the soil temperature, soil moisture and wet-leaf water layer."
function integrate_land_surface!(s::ModelState)
    ls, dt = s.land_surface, s.input.dt
    ls.Tsoil += dt * ls.Tsoiltend
    ls.wg += dt * ls.wgtend
    ls.Wl += dt * ls.Wltend
    return nothing
end
