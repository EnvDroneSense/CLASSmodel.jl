# Title, unit, symbol and slider range of the model inputs, for forms and labels.
#
# Units are the model units of `ModelInput`. `min` and `max` bound a slider; `validate`
# stays the only check of a run. The default of `ModelInput()` is the typical value.

num(title, unit, symbol, min, max) = (; title, unit, symbol, min, max)
flag(title, symbol) = (; title, unit = "-", symbol, min = nothing, max = nothing)
choice(title, symbol, options) = (; title, unit = "-", symbol, options)

const INPUT_INFO = (
    runtime = num("Run time", "s", "t_run", 600, 172800),
    dt = num("Time step", "s", "Δt", 1, 600),
    sw_ml = flag("Mixed-layer model", "sw_ml"),
    sw_shearwe = flag("Shear growth of the mixed layer", "sw_shearwe"),
    sw_fixft = flag("Fixed free troposphere", "sw_fixft"),
    sw_wind = flag("Prognostic wind", "sw_wind"),
    sw_sl = flag("Surface layer", "sw_sl"),
    sw_rad = flag("Radiation", "sw_rad"),
    sw_ls = flag("Land surface", "sw_ls"),
    sw_sea = flag("Sea instead of land", "sw_sea"),
    sw_cu = flag("Cumulus parameterization", "sw_cu"),
    sw_curad = flag("Cumulus cloud cover in radiation", "sw_curad"),
    sw_chem = flag("Chemical reactions", "sw_chem"),
    sinperiod = num("Half period of the flux cycle", "s", "T_sin", 3600, 86400),
    sw_wtheta = flag("Diurnal cycle of the heat flux", "sw_wtheta"),
    sw_wq = flag("Diurnal cycle of the moisture flux", "sw_wq"),
    h = num("Initial boundary-layer height", "m", "h", 50, 3000),
    Ps = num("Surface pressure", "Pa", "Pₛ", 90000, 105000),
    divU = num("Large-scale wind divergence", "s-1", "div U", 0, 1e-5),
    fc = num("Coriolis parameter", "s-1", "f", 0, 2e-4),
    beta = num("Entrainment ratio for virtual heat", "-", "β", 0, 1),
    theta = num("Initial potential temperature", "K", "θ", 270, 320),
    dtheta = num("Initial potential temperature jump", "K", "Δθ", 0, 10),
    gammatheta = num("Potential temperature lapse rate", "K m-1", "γθ", 0, 0.02),
    advtheta = num("Heat advection", "K s-1", "advθ", -5e-4, 5e-4),
    wtheta = num("Surface kinematic heat flux", "K m s-1", "w'θ'ₛ", -0.1, 0.5),
    q = num("Initial specific humidity", "kg kg-1", "q", 0, 0.03),
    dq = num("Initial specific humidity jump", "kg kg-1", "Δq", -0.01, 0.01),
    gammaq = num("Specific humidity lapse rate", "kg kg-1 m-1", "γq", -1e-5, 0),
    advq = num("Moisture advection", "kg kg-1 s-1", "advq", -1e-6, 1e-6),
    wq = num("Surface kinematic moisture flux", "kg kg-1 m s-1", "w'q'ₛ", 0, 5e-4),
    CO2 = num("Initial CO2", "ppm", "CO2", 300, 600),
    dCO2 = num("Initial CO2 jump", "ppm", "ΔCO2", -200, 100),
    gammaCO2 = num("CO2 lapse rate", "ppm m-1", "γCO2", -0.1, 0.1),
    advCO2 = num("CO2 advection", "ppm s-1", "advCO2", -0.01, 0.01),
    wCO2 = num("Surface CO2 flux", "mg CO2 m-2 s-1", "w'CO2'ₛ", -0.5, 0.5),
    S = num("Initial scalar", "S", "S", -100, 100),
    dS = num("Initial scalar jump", "S", "ΔS", -100, 100),
    gammaS = num("Scalar lapse rate", "S m-1", "γS", -0.1, 0.1),
    advS = num("Scalar advection", "S s-1", "advS", -0.01, 0.01),
    wS = num("Surface scalar flux", "S m s-1", "w'S'ₛ", -0.5, 0.5),
    C = num("Initial species mixing ratio", "ppb", "C", 0, 1e9),
    dC = num("Initial species jump", "ppb", "ΔC", -100, 100),
    gammaC = num("Species lapse rate", "ppb m-1", "γC", -0.1, 0.1),
    advC = num("Species advection", "ppb s-1", "advC", -0.01, 0.01),
    wC = num("Surface species flux", "ppb m s-1", "w'C'ₛ", -10, 10),
    sw_wC = choice("Species flux mode", "sw_wC", (0, 1, 2)),
    sw_chem_constant = flag("Chemistry at the reference state", "sw_chem_constant"),
    sw_photo_constant = flag("Photolysis at the reference time", "sw_photo_constant"),
    sw_reactions = choice("Reaction switches", "sw_reactions", (false, true)),
    tod_ref = num("Reference time of the day", "h UTC", "tod_ref", 0, 24),
    Tcbl_ref = num("Reference temperature of the mixed layer", "K", "Tcbl_ref", 250, 320),
    Tft_ref = num("Reference temperature of the free troposphere", "K", "Tft_ref", 250, 320),
    qcbl_ref = num("Reference humidity of the mixed layer", "kg kg-1", "qcbl_ref", 0, 0.03),
    qft_ref = num("Reference humidity of the free troposphere", "kg kg-1", "qft_ref", 0, 0.03),
    P_ref = num("Reference pressure", "Pa", "P_ref", 50000, 105000),
    stocoef = num("OH yield of RO2 + HO2", "-", "stocoef", 0, 10),
    u = num("Initial u wind", "m s-1", "u", -20, 20),
    du = num("Initial u wind jump", "m s-1", "Δu", -20, 20),
    gammau = num("u wind lapse rate", "s-1", "γu", -0.01, 0.01),
    advu = num("u wind advection", "m s-2", "advu", -1e-3, 1e-3),
    v = num("Initial v wind", "m s-1", "v", -20, 20),
    dv = num("Initial v wind jump", "m s-1", "Δv", -20, 20),
    gammav = num("v wind lapse rate", "s-1", "γv", -0.01, 0.01),
    advv = num("v wind advection", "m s-2", "advv", -1e-3, 1e-3),
    ustar = num("Friction velocity", "m s-1", "u*", 0.01, 2),
    z0m = num("Roughness length for momentum", "m", "z0m", 1e-4, 5),
    z0h = num("Roughness length for scalars", "m", "z0h", 1e-5, 5),
    lat = num("Latitude", "deg", "φ", -90, 90),
    lon = num("Longitude", "deg", "λ", -180, 180),
    doy = num("Day of the year", "-", "doy", 1, 366),
    tstart = num("Start time of the day", "h UTC", "tstart", 0, 24),
    cc = num("Cloud cover fraction", "-", "cc", 0, 1),
    Q = num("Net radiation", "W m-2", "Q", 0, 1000),
    dFz = num("Cloud top radiative divergence", "W m-2", "ΔFz", 0, 200),
    ls_type = choice("Land surface type", "ls_type", (:js, :ags)),
    wg = num("Water content top soil layer", "m3 m-3", "wg", 0.01, 0.6),
    w2 = num("Water content deeper soil layer", "m3 m-3", "w2", 0.01, 0.6),
    Tsoil = num("Temperature top soil layer", "K", "Tsoil", 250, 320),
    T2 = num("Temperature deeper soil layer", "K", "T2", 250, 320),
    a = num("Retention curve parameter a", "-", "a", 0.1, 0.6),
    b = num("Retention curve parameter b", "-", "b", 2, 13),
    p = num("Retention curve parameter p", "-", "p", 1, 12),
    CGsat = num("Saturated soil conductivity for heat", "K m2 J-1", "CGsat", 1e-6, 1e-5),
    wsat = num("Saturated water content", "m3 m-3", "wsat", 0.3, 0.6),
    wfc = num("Water content at field capacity", "m3 m-3", "wfc", 0.1, 0.5),
    wwilt = num("Water content at wilting point", "m3 m-3", "wwilt", 0.05, 0.3),
    C1sat = num("Force-restore coefficient C1 at saturation", "-", "C1sat", 0.05, 0.5),
    C2ref = num("Force-restore coefficient C2 at reference", "-", "C2ref", 0, 4),
    c_beta = num("Curvature of the plant water-stress factor", "-", "cβ", 0, 1),
    LAI = num("Leaf area index", "-", "LAI", 0, 8),
    gD = num("Correction factor of transpiration for VPD", "-", "gD", 0, 0.05),
    rsmin = num("Minimum resistance of transpiration", "s m-1", "rsmin", 20, 500),
    rssoilmin = num("Minimum resistance of soil evaporation", "s m-1", "rssoilmin", 0, 500),
    alpha = num("Surface albedo", "-", "α", 0, 0.6),
    cveg = num("Vegetation fraction", "-", "cveg", 0, 1),
    Ts = num("Initial surface temperature", "K", "Ts", 250, 330),
    Wmax = num("Thickness of water layer on wet vegetation", "m", "Wmax", 0, 1e-3),
    Wl = num("Water layer depth on wet vegetation", "m", "Wl", 0, 1e-3),
    Lambda = num("Thermal diffusivity of the skin layer", "W m-2 K-1", "Λ", 1, 20),
    c3c4 = choice("Plant type", "c3c4", (:c3, :c4)),
    dz_h = num("Transition layer thickness", "m", "Δz", 10, 500)
)

"""
    input_info(name::Symbol) -> NamedTuple

Title, unit and symbol of the `ModelInput` field `name`, with the range `min` and `max` for a
slider (not a validity rule; see [`validate`](@ref)). The unit is the model unit, `"-"` for a
dimensionless field. A species field (`C`, `dC`, ...) holds one value for each species, and `sw_reactions` one
switch for each reaction; the entry gives the range of one value. A switch has `min = max = nothing`; a `Symbol` field has `options`
instead of a range (`sw_wC`: the modes 0, 1 and 2). The typical value is the default of `ModelInput()`.
"""
function input_info(name::Symbol)
    haskey(INPUT_INFO, name) || throw(ArgumentError("unknown input field :$name"))
    return INPUT_INFO[name]
end
