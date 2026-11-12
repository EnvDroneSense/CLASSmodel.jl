# Cumulus: mass-flux closure for shallow cumulus at the top of the mixed layer.

# Square root that gives NaN, not an error, for a negative argument.
nan_sqrt(x) = x < 0 ? oftype(x, NaN) : sqrt(x)

# Variance at h of a scalar with jump `d`, entrainment flux `fe` and cloud flux `fM`.
function top_variance(ml::MixedLayerState, cu::CumulusState, fe, fM, d)
    return ml.wthetav > 0 ? -(fe + fM) * d * ml.h / (cu.dz_h * ml.wstar) : zero(ml.q)
end

# Cloud flux of a scalar that is not moisture: only for a negative jump at h.
function cloud_flux(ml::MixedLayerState, cu::CumulusState, fe, fM, d)
    return d < 0 ? cu.M * nan_sqrt(top_variance(ml, cu, fe, fM, d)) : zero(ml.q)
end

"""
    cumulus!(s::ModelState)

Cloud core fraction `ac`, mass flux `M` and the cloud fluxes of moisture, CO2 and the
species, from their variance at h. Needs the entrainment fluxes of a previous
`mixed_layer!` call. With `sw_curad` it sets the cloud cover of the radiation to `2 ac`.
"""
function cumulus!(s::ModelState)
    ml, cu, sp = s.mixed_layer, s.cumulus, s.species

    # cloud core fraction (zero when not positive or undefined), mass flux and moisture flux
    sigma_q = nan_sqrt(top_variance(ml, cu, ml.wqe, cu.wqM, ml.dq))
    ac = 0.5 + 0.36 * atan(1.55 * ((ml.q - qsat(ml.T_h, ml.P_h)) / sigma_q))
    cu.ac = ac > 0 ? ac : zero(ac)
    cu.M = cu.ac * ml.wstar
    s.input.sw_curad && (s.radiation.cc = 2 * cu.ac)
    cu.wqM = cu.M * sigma_q

    cu.wCO2M = cloud_flux(ml, cu, ml.wCO2e, cu.wCO2M, ml.dCO2)
    for i in 1:NSPECIES
        sp.wCM[i] = cloud_flux(ml, cu, sp.wCe[i], sp.wCM[i], sp.dC[i])
    end
    return nothing
end
