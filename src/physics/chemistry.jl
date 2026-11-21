# Gas-phase chemistry of the 22 species: rate constants and the solver. Mixing ratios in ppb.

const MW_AIR = 28.97   # molar mass of dry air [g mol-1]
const MW_H2O = 18.0    # molar mass of water [g mol-1]
const I_H2O, I_NO, I_NO2, I_O3 = species_index.((:H2O, :NO, :NO2, :O3))
# Rate constant of the reaction `name`, or NaN when the scheme does not have it.
function rate_of(K, scheme::ChemistryScheme, name)
    r = findfirst(rx -> rx.name == name, scheme.reactions)
    return r === nothing ? oftype(zero(eltype(K)), NaN) : K[r]
end

# Specific humidity [kg kg-1] as a mixing ratio of H2O [ppb].
h2o_ppb(q) = q * 1e9 * MW_AIR / MW_H2O

"""
    rate_constants!(K, scheme, p, T, q, coszen, cloud)

Fill `K` with the rate constants [s-1, ppb-1 s-1] of the reactions of `scheme` at pressure
`p` [Pa], temperature `T` [K] and specific humidity `q`. A photolysis rate is multiplied by
`cloud`. The forms are in [`Reaction`](@ref).
"""
function rate_constants!(K, scheme::ChemistryScheme, p, T, q, coszen, cloud)
    q_ppb = h2o_ppb(q)
    conv = 6.023e8 * (p / 100) / (8.314e-2 * T)   # molecules cm-3 per ppb
    for (r, rx) in enumerate(scheme.reactions)
        A, B, C, D, E, F, G = rx.k
        form = rx.form
        if form === :photolysis
            K[r] = coszen > 0 ? cloud * (A * exp(B / coszen)) : zero(coszen)
        elseif form === :constant
            K[r] = A * conv
        elseif form === :arrhenius
            K[r] = A * exp(B / T) * conv
        elseif form === :falloff || form === :decomposition
            k1 = A * (T / 300)^B * exp(C / T) * conv * 1e9
            k2 = D * (T / 300)^E * exp(F / T)
            K[r] = k1 * k2 / (k1 + k2) * G * (form === :falloff ? conv : 1)
        elseif form === :ho2
            k1 = A * exp(B / T) * conv
            k2 = C * exp(D / T) * conv * conv * 1e9
            k3 = E * exp(F / T) * conv * q_ppb
            K[r] = (k1 + k2) * (1 + k3)
        else   # :power
            K[r] = A * (T / B)^C * exp(D / T) * conv^2
        end
    end
    return K
end

# Advance the mixing ratios `c` with H2O from `q` by `dt` into `y`: four Gauss-Seidel passes
# of y = (y0 + dt P) / (1 + dt L), in the order of the species.
function solve_chemistry!(y, y0, c, q, scheme::ChemistryScheme, K, dt)
    copyto!(y0, c)
    copyto!(y, c)
    y0[I_H2O] = y[I_H2O] = h2o_ppb(q)
    for _ in 1:4
        for n in scheme.solved
            P = L = zero(eltype(y))
            for t in scheme.terms[n]
                v = t.coef * K[t.reaction] * y[t.i1]^t.e1 * y[t.i2]^t.e2
                if t.production
                    P += v
                else
                    L += v
                end
            end
            y[n] = max(zero(eltype(y)), (y0[n] + dt * P) / (1 + dt * L))
        end
    end
    return y
end

"""
    chemistry!(s::ModelState, dt)

Advance the species by `dt` in the mixed layer and in the free troposphere (`C + dC`), and
set `dC` from the two. H2O follows `q`. Also store `phi`, the photostationary state of NO,
NO2 and O3, and `k_r05`, the NO2 photolysis rate. See `ModelInput` for the switches.
"""
function chemistry!(s::ModelState, dt)
    x, c, ml, ch, sp = s.input, s.constants, s.mixed_layer, s.chemistry, s.species
    sinlea = solar_elevation(x, utc_seconds(x, s.t * x.dt))
    if x.sw_photo_constant && sinlea >= 0
        sinlea = solar_elevation(x, x.tod_ref * 3600)
    end
    if x.sw_chem_constant
        p_cbl, p_ft, T_cbl, T_ft = x.P_ref, x.P_ref, x.Tcbl_ref, x.Tft_ref
        q_cbl, q_ft = x.qcbl_ref, x.qft_ref
    else
        p_cbl, p_ft = x.Ps, x.Ps - c.rho * c.g * ml.h
        T_cbl, T_ft = ml.theta, ml.theta - c.g / c.cp * ml.h
        q_cbl, q_ft = ml.q, ml.q + ml.dq
    end
    scheme = ch.scheme
    rate_constants!(ch.Kcbl, scheme, p_cbl, T_cbl, q_cbl, sinlea, 1 - 0.5 * s.radiation.cc)
    rate_constants!(ch.Kft, scheme, p_ft, T_ft, q_ft, sinlea, one(sinlea))

    y = ch.y
    @. ch.fsc = sp.C + sp.dC   # before C changes
    solve_chemistry!(y, ch.y0, sp.C, q_cbl, scheme, ch.Kcbl, dt)
    copyto!(sp.C, y)
    K05 = rate_of(ch.Kcbl, scheme, :r05)   # NO2 photolysis
    K21 = rate_of(ch.Kcbl, scheme, :r21)   # NO + O3
    ch.phi = isnan(K05 + K21) ? K05 + K21 :
             K05 > 1e-10 && y[I_NO2] > 1e-5 ? K21 * y[I_NO] * y[I_O3] / (K05 * y[I_NO2]) : zero(K05)
    ch.k_r05 = K05

    solve_chemistry!(y, ch.y0, ch.fsc, q_ft, scheme, ch.Kft, dt)
    @. sp.dC = y - sp.C
    return s
end
