# Cumulus: mass-flux closure for shallow cumulus at the top of the mixed layer.

# Square root that gives NaN, not an error, for a negative argument.
nan_sqrt(x) = x < 0 ? oftype(x, NaN) : sqrt(x)

"""
    cumulus!(s::ModelState)

Cloud core fraction `ac`, mass flux `M` and the moisture and CO2 fluxes `wqM`, `wCO2M`
carried by the clouds, from the humidity and CO2 variance at the top of the mixed layer.
Needs the entrainment fluxes of a previous `mixed_layer!` call.
"""
function cumulus!(s::ModelState)
    ml, cu = s.mixed_layer, s.cumulus

    # humidity and CO2 variance at the top of the mixed layer
    if ml.wthetav > 0
        q2_h = -(ml.wqe + cu.wqM) * ml.dq * ml.h / (cu.dz_h * ml.wstar)
        CO22_h = -(ml.wCO2e + cu.wCO2M) * ml.dCO2 * ml.h / (cu.dz_h * ml.wstar)
    else
        q2_h = zero(ml.q)
        CO22_h = zero(ml.q)
    end

    # cloud core fraction (zero when not positive or undefined), mass flux and moisture flux
    ac = 0.5 + 0.36 * atan(1.55 * ((ml.q - qsat(ml.T_h, ml.P_h)) / nan_sqrt(q2_h)))
    cu.ac = ac > 0 ? ac : zero(ac)
    cu.M = cu.ac * ml.wstar
    cu.wqM = cu.M * nan_sqrt(q2_h)

    # CO2 mass flux only if the jump at the top is negative
    cu.wCO2M = ml.dCO2 < 0 ? cu.M * nan_sqrt(CO22_h) : zero(ml.q)
    return nothing
end
