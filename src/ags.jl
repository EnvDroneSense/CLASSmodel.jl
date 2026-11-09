# A-Gs: photosynthesis-based canopy resistance and surface CO2 flux.

"""
    E1(x)

Exponential integral E1 for `x > 0`, from its power series (99 terms). Accurate for
small `x`; the alternating series loses precision for large `x`.
"""
function E1(x)
    s = zero(x)
    term = one(x)   # x^k / k! with the sign (-1)^k, built up term by term
    for k in 1:99
        term *= -x / k
        s += term / k
    end
    return -0.57721566490153286060 - log(x) - s
end

"""
    ags!(s::ModelState)

Surface resistance `rs` for moisture from the CO2 assimilation of the canopy (A-Gs), and
the surface CO2 fluxes: assimilation `wCO2A`, soil respiration `wCO2R` and their sum
`wCO2`, in ppm m s-1. Uses the surface resistance `ra`, the incoming shortwave radiation
and the surface temperature of the previous step.
"""
function ags!(s::ModelState{T}) where {T}
    x, c, ml, sl, ls = s.input, s.constants, s.mixed_layer, s.surface_layer, s.land_surface
    p = x.c3c4 == :c3 ? AGsParameters{T}(C3()) : AGsParameters{T}(C4())
    rho = c.rho

    # CO2 compensation concentration and mesophyll conductance
    CO2comp = p.CO2comp298 * rho * p.Q10CO2^(0.1 * (sl.thetasurf - 298))
    gm = p.gm298 * p.Q10gm^(0.1 * (sl.thetasurf - 298)) /
         ((1 + exp(0.3 * (p.T1gm - sl.thetasurf))) * (1 + exp(0.3 * (sl.thetasurf - p.T2gm))))
    gm = gm / 1000   # mm s-1 to m s-1

    # CO2 concentration inside the leaf
    fmin0 = p.gmin / c.nuco2q - 1 / 9 * gm
    fmin = -fmin0 + (fmin0^2 + 4 * p.gmin / c.nuco2q * gm)^0.5 / (2 * gm)
    Ds = (esat(ls.Ts) - ml.e) / 1000   # vapor pressure deficit [kPa]
    D0 = (p.f0 - fmin) / p.ad
    cfrac = p.f0 * (1 - (Ds / D0)) + fmin * (Ds / D0)
    co2abs = ml.CO2 * (c.mco2 / c.mair) * rho   # ppm to mg CO2 m-3
    ci = cfrac * (co2abs - CO2comp) + CO2comp

    # maximal gross assimilation, with the effect of soil moisture stress
    Ammax = p.Ammax298 * p.Q10Am^(0.1 * (sl.thetasurf - 298)) /
            ((1 + exp(0.3 * (p.T1Am - sl.thetasurf))) * (1 + exp(0.3 * (sl.thetasurf - p.T2Am))))
    betaw = max(1e-3, min(1, (x.w2 - x.wwilt) / (x.wfc - x.wwilt)))
    if x.c_beta == 0
        fstr = betaw
    else
        # stress function of Combe et al. (2016)
        P = x.c_beta < 0.25 ? 6.4 * x.c_beta : x.c_beta < 0.5 ? 7.6 * x.c_beta - 0.3 : 2^(3.66 * x.c_beta + 0.34) - 1
        fstr = (1 - exp(-P * betaw)) / (1 - exp(-P))
    end

    # leaf assimilation and light use efficiency
    Am = Ammax * (1 - exp(-(gm * (ci - CO2comp) / Ammax)))
    Rdark = (1 / 9) * Am
    PAR = 0.5 * max(1e-1, s.radiation.Swin * x.cveg)
    alphac = p.alpha0 * (co2abs - CO2comp) / (co2abs + 2 * CO2comp)

    # upscaling from leaf to canopy: net CO2 uptake and canopy conductance
    y = alphac * p.Kx * PAR / (Am + Rdark)
    An = (Am + Rdark) * (1 - 1 / (p.Kx * x.LAI) * (E1(y * exp(-p.Kx * x.LAI)) - E1(y)))
    a1 = 1 / (1 - p.f0)
    Dstar = D0 / (a1 * (p.f0 - fmin))
    gcco2 = x.LAI * (p.gmin / c.nuco2q + a1 * fstr * An / ((co2abs - CO2comp) * (1 + Ds / Dstar)))

    # surface resistance for moisture and CO2, and the net CO2 flux into the plant
    ls.rs = 1 / (1.6 * gcco2)
    rsCO2 = 1 / gcco2
    An = -(co2abs - ci) / (sl.ra + rsCO2)

    # soil respiration
    fw = c.Cw * c.wmax / (ls.wg + c.wmin)
    Resp = c.R10 * (1 - fw) * exp(c.E0 / (283.15 * 8.314) * (1 - 283.15 / ls.Tsoil))

    ml.wCO2A = An * (c.mair / (rho * c.mco2))
    ml.wCO2R = Resp * (c.mair / (rho * c.mco2))
    ml.wCO2 = ml.wCO2A + ml.wCO2R
    return nothing
end
