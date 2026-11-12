# Cumulus: mass-flux closure for shallow cumulus at the top of the mixed layer.

# Square root that gives NaN, not an error, for a negative argument.
nan_sqrt(x) = x < 0 ? oftype(x, NaN) : sqrt(x)

"""
    cumulus!(s::ModelState)

Cloud core fraction `ac`, mass flux `M` and the moisture, CO2, scalar and species fluxes `wqM`,
`wCO2M`, `wSM`, `species.wCM` carried by the clouds, from the humidity, CO2 and scalar variance at the top of the mixed layer.
Needs the entrainment fluxes of a previous `mixed_layer!` call. With `sw_curad` it also sets
the cloud cover of the radiation to `2 ac`.
"""
function cumulus!(s::ModelState)
    ml, cu = s.mixed_layer, s.cumulus

    # humidity and CO2 variance at the top of the mixed layer
    if ml.wthetav > 0
        q2_h = -(ml.wqe + cu.wqM) * ml.dq * ml.h / (cu.dz_h * ml.wstar)
        CO22_h = -(ml.wCO2e + cu.wCO2M) * ml.dCO2 * ml.h / (cu.dz_h * ml.wstar)
        S2_h = -(ml.wSe + cu.wSM) * ml.dS * ml.h / (cu.dz_h * ml.wstar)
    else
        q2_h = zero(ml.q)
        CO22_h = zero(ml.q)
        S2_h = zero(ml.q)
    end

    # cloud core fraction (zero when not positive or undefined), mass flux and moisture flux
    ac = 0.5 + 0.36 * atan(1.55 * ((ml.q - qsat(ml.T_h, ml.P_h)) / nan_sqrt(q2_h)))
    cu.ac = ac > 0 ? ac : zero(ac)
    cu.M = cu.ac * ml.wstar
    s.input.sw_curad && (s.radiation.cc = 2 * cu.ac)
    cu.wqM = cu.M * nan_sqrt(q2_h)

    # CO2 mass flux only if the jump at the top is negative
    cu.wCO2M = ml.dCO2 < 0 ? cu.M * nan_sqrt(CO22_h) : zero(ml.q)
    cu.wSM = ml.dS < 0 ? cu.M * nan_sqrt(S2_h) : zero(ml.q)

    # the species, as the scalar
    sp = s.species
    for i in 1:NSPECIES
        C2_h = ml.wthetav > 0 ?
               -(sp.wCe[i] + sp.wCM[i]) * sp.dC[i] * ml.h / (cu.dz_h * ml.wstar) : zero(ml.q)
        sp.wCM[i] = sp.dC[i] < 0 ? cu.M * nan_sqrt(C2_h) : zero(ml.q)
    end
    return nothing
end
