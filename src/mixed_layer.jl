# Mixed-layer model: diagnostics, tendencies and time integration.

"""
    statistics!(s::ModelState)

Virtual temperatures, values at the top of the mixed layer and the lifting
condensation level (LCL). The LCL is found by a fixed-point iteration: it restarts
from `h` at the first step and otherwise continues from the previous value.
"""
function statistics!(s::ModelState)
    x, c, ml = s.input, s.constants, s.mixed_layer

    ml.thetav = ml.theta + VIRT * ml.theta * ml.q
    ml.wthetav = ml.wtheta + VIRT * ml.theta * ml.wq
    ml.dthetav = (ml.theta + ml.dtheta) * (1 + VIRT * (ml.q + ml.dq)) - ml.theta * (1 + VIRT * ml.q)

    ml.P_h = x.Ps - c.rho * c.g * ml.h
    ml.T_h = ml.theta - c.g / c.cp * ml.h
    ml.RH_h = ml.q / qsat(ml.T_h, ml.P_h)

    if s.t == 0
        lcl = ml.h
        RHlcl = 0.5
    else
        lcl = ml.lcl
        RHlcl = 0.9998
    end
    itmax = 30
    it = 0
    while (RHlcl <= 0.9999 || RHlcl >= 1.0001) && it < itmax
        lcl += (1 - RHlcl) * 1000
        p_lcl = x.Ps - c.rho * c.g * lcl
        T_lcl = ml.theta - c.g / c.cp * lcl
        RHlcl = ml.q / qsat(T_lcl, p_lcl)
        it += 1
    end
    ml.lcl = lcl
    it == itmax && @warn "LCL calculation not converged" RHlcl lcl
    return s
end

# Surface momentum fluxes from the friction velocity `ustar` and the wind direction, for
# runs without the surface layer.
function prescribed_stress!(s::ModelState)
    w, sl = s.wind, s.surface_layer
    sl.uw = -sign(w.u) * sqrt(sl.ustar^4 / (w.v^2 / w.u^2 + 1))
    sl.vw = -sign(w.v) * sqrt(sl.ustar^4 / (w.u^2 / w.v^2 + 1))
    return nothing
end

"""
    mixed_layer!(s::ModelState)

Entrainment, fluxes and tendencies of the mixed layer. Without the surface layer the
momentum fluxes follow from `ustar` and the wind direction.
"""
function mixed_layer!(s::ModelState)
    x, c, ml, w, sl, cu = s.input, s.constants, s.mixed_layer, s.wind, s.surface_layer, s.cumulus

    x.sw_sl || prescribed_stress!(s)

    ml.ws = -x.divU * ml.h
    # vertical velocity that keeps the free troposphere fixed under subsidence
    ws_ft = x.sw_fixft ? ml.ws : zero(ml.ws)

    ml.wf = x.dFz / (c.rho * c.cp * ml.dtheta)

    if ml.wthetav > 0
        ml.wstar = ((c.g * ml.h * ml.wthetav) / ml.thetav)^(1 / 3)
    else
        ml.wstar = 1e-6
    end

    ml.wthetave = -x.beta * ml.wthetav

    if x.sw_shearwe
        ml.we = (-ml.wthetave + 5 * sl.ustar^3 * ml.thetav / (c.g * ml.h)) / ml.dthetav
    else
        ml.we = -ml.wthetave / ml.dthetav
    end
    # no shrinking of the boundary layer
    ml.we = max(ml.we, 0)

    ml.wthetae = -ml.we * ml.dtheta
    ml.wqe = -ml.we * ml.dq
    ml.wCO2e = -ml.we * ml.dCO2

    # growth of the layer relative to the free troposphere, and of the layer itself
    wrel = ml.we + ml.wf - cu.M
    ml.htend = wrel + ml.ws

    ml.thetatend = (ml.wtheta - ml.wthetae) / ml.h + x.advtheta
    ml.qtend = (ml.wq - ml.wqe - cu.wqM) / ml.h + x.advq
    ml.CO2tend = (ml.wCO2 - ml.wCO2e - cu.wCO2M) / ml.h + x.advCO2

    ml.dthetatend = x.gammatheta * (wrel + ws_ft) - ml.thetatend
    ml.dqtend = x.gammaq * (wrel + ws_ft) - ml.qtend
    ml.dCO2tend = x.gammaCO2 * (wrel + ws_ft) - ml.CO2tend

    # the wind above the layer is geostrophic: u + du = ug
    if x.sw_wind
        w.utend = -x.fc * w.dv + (sl.uw + ml.we * w.du) / ml.h + x.advu
        w.vtend = x.fc * w.du + (sl.vw + ml.we * w.dv) / ml.h + x.advv

        w.dutend = x.gammau * wrel - w.utend
        w.dvtend = x.gammav * wrel - w.vtend
    end

    if cu.ac > 0 || ml.lcl - ml.h < 300
        cu.dztend = ((ml.lcl - ml.h) - cu.dz_h) / 7200
    else
        cu.dztend = zero(cu.dztend)
    end
    return s
end

"Advance the mixed-layer variables by one forward Euler step."
function integrate_mixed_layer!(s::ModelState)
    x, ml, w, cu = s.input, s.mixed_layer, s.wind, s.cumulus
    dt = x.dt

    ml.h += dt * ml.htend
    ml.theta += dt * ml.thetatend
    ml.dtheta += dt * ml.dthetatend
    ml.q += dt * ml.qtend
    ml.dq += dt * ml.dqtend
    ml.CO2 += dt * ml.CO2tend
    ml.dCO2 += dt * ml.dCO2tend
    cu.dz_h += dt * cu.dztend
    # the transition layer is at least 50 m thick
    cu.dz_h = max(cu.dz_h, 50)

    if x.sw_wind
        w.u += dt * w.utend
        w.du += dt * w.dutend
        w.v += dt * w.vtend
        w.dv += dt * w.dvtend
    end
    return s
end
