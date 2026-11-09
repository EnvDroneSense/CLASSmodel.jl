# Surface layer: Monin-Obukhov similarity between the surface and the lowest tenth of the mixed layer.

"Integrated stability function for momentum at `zeta = z / L`."
function psim(zeta)
    if zeta <= 0
        x = (1 - 16 * zeta)^0.25
        return 3.14159265 / 2 - 2 * atan(x) + log((1 + x)^2 * (1 + x^2) / 8)
    else
        return -2 / 3 * (zeta - 5 / 0.35) * exp(-0.35 * zeta) - zeta - (10 / 3) / 0.35
    end
end

"Integrated stability function for heat and scalars at `zeta = z / L`."
function psih(zeta)
    if zeta <= 0
        x = (1 - 16 * zeta)^0.25
        return 2 * log((1 + x * x) / 2)
    else
        return -2 / 3 * (zeta - 5 / 0.35) * exp(-0.35 * zeta) - (1 + (2 / 3) * zeta)^1.5 - (10 / 3) / 0.35 + 1
    end
end

# Bulk Richardson number as a function of the Obukhov length `L`.
ribfun(L, zsl, z0m, z0h) =
    zsl / L * (log(zsl / z0h) - psih(zsl / L) + psih(z0h / L)) / (log(zsl / z0m) - psim(zsl / L) + psim(z0m / L))^2

"""
    ribtol(Rib, zsl, z0m, z0h)

Obukhov length [m] for the bulk Richardson number `Rib` at height `zsl`, found by Newton
iteration (central difference for the derivative) from `L = ±1`. The iteration stops when
`L` changes by less than 1 mm, when `|L|` exceeds 1e15 (neutral), or after 100 iterations
(with a warning).
"""
function ribtol(Rib, zsl, z0m, z0h)
    L, L0 = Rib > 0 ? (one(Rib), 2 * one(Rib)) : (-one(Rib), -2 * one(Rib))
    itmax = 100
    it = 0
    while abs(L - L0) > 0.001
        if it == itmax
            @warn "Obukhov length not converged" Rib L
            break
        end
        it += 1
        L0 = L
        fx = Rib - ribfun(L, zsl, z0m, z0h)
        Lstart = L - 0.001 * L
        Lend = L + 0.001 * L
        fxdif = (-ribfun(Lstart, zsl, z0m, z0h) + ribfun(Lend, zsl, z0m, z0h)) / (Lstart - Lend)
        L = L - fx / fxdif
        abs(L) > 1e15 && break
    end
    return L
end

"""
    surface_layer!(s::ModelState)

Surface values, bulk Richardson number, Obukhov length, drag coefficients, friction
velocity, momentum fluxes and the 2 m diagnostics. The surface values use the drag
coefficient and the surface resistance of the previous call.
"""
function surface_layer!(s::ModelState)
    x, c, ml, wind, sl, ls = s.input, s.constants, s.mixed_layer, s.wind, s.surface_layer, s.land_surface

    ueff = max(0.01, sqrt(wind.u^2 + wind.v^2 + ml.wstar^2))
    sl.thetasurf = ml.theta + ml.wtheta / (sl.Cs * ueff)
    qsatsurf = qsat(sl.thetasurf, x.Ps)
    cq = (1 + sl.Cs * ueff * ls.rs)^-1
    sl.qsurf = (1 - cq) * ml.q + cq * qsatsurf
    sl.thetavsurf = sl.thetasurf * (1 + VIRT * sl.qsurf)

    zsl = 0.1 * ml.h
    sl.Rib = min(c.g / ml.thetav * zsl * (ml.thetav - sl.thetavsurf) / ueff^2, 0.2)
    sl.L = ribtol(sl.Rib, zsl, x.z0m, x.z0h)

    fm = log(zsl / x.z0m) - psim(zsl / sl.L) + psim(x.z0m / sl.L)
    fh = log(zsl / x.z0h) - psih(zsl / sl.L) + psih(x.z0h / sl.L)
    sl.Cm = c.k^2 / fm^2
    sl.Cs = c.k^2 / fm / fh

    sl.ustar = sqrt(sl.Cm) * ueff
    sl.uw = -sl.Cm * ueff * wind.u
    sl.vw = -sl.Cm * ueff * wind.v

    # 2 m diagnostics
    fh2 = log(2 / x.z0h) - psih(2 / sl.L) + psih(x.z0h / sl.L)
    fm2 = log(2 / x.z0m) - psim(2 / sl.L) + psim(x.z0m / sl.L)
    sl.T2m = sl.thetasurf - ml.wtheta / sl.ustar / c.k * fh2
    sl.q2m = sl.qsurf - ml.wq / sl.ustar / c.k * fh2
    sl.u2m = -sl.uw / sl.ustar / c.k * fm2
    sl.v2m = -sl.vw / sl.ustar / c.k * fm2
    sl.esat2m = esat(sl.T2m)
    sl.e2m = sl.q2m * x.Ps / ε
    return nothing
end
