# Title, unit and symbol of the output variables, for axis labels and legends.
#
# Units follow the comments of `ModelInput` and the units section of the documentation
# (`wCO2*` in mg CO2 m-2 s-1, `CO2` and `dCO2` in ppm). Dimensionless variables have "-".

info(title, unit, symbol) = (; title, unit, symbol)

const BASE_INFO = (
    t = info("Time", "h", "t"),
    h = info("Boundary-layer height", "m", "h"),
    theta = info("Mixed-layer potential temperature", "K", "θ"),
    thetav = info("Mixed-layer virtual potential temperature", "K", "θv"),
    dtheta = info("Potential temperature jump", "K", "Δθ"),
    dthetav = info("Virtual potential temperature jump", "K", "Δθv"),
    wtheta = info("Surface kinematic heat flux", "K m s-1", "w'θ'ₛ"),
    wthetav = info("Surface kinematic virtual heat flux", "K m s-1", "w'θv'ₛ"),
    wthetae = info("Entrainment kinematic heat flux", "K m s-1", "w'θ'ₑ"),
    wthetave = info("Entrainment kinematic virtual heat flux", "K m s-1", "w'θv'ₑ"),
    q = info("Mixed-layer specific humidity", "kg kg-1", "q"),
    dq = info("Specific humidity jump", "kg kg-1", "Δq"),
    wq = info("Surface kinematic moisture flux", "kg kg-1 m s-1", "w'q'ₛ"),
    wqe = info("Entrainment kinematic moisture flux", "kg kg-1 m s-1", "w'q'ₑ"),
    wqM = info("Cumulus kinematic moisture flux", "kg kg-1 m s-1", "w'q'ₘ"),
    qsat = info("Saturated specific humidity", "kg kg-1", "qsat"),
    e = info("Vapour pressure", "Pa", "e"),
    esat = info("Saturated vapour pressure", "Pa", "esat"),
    CO2 = info("Mixed-layer CO2", "ppm", "CO2"),
    dCO2 = info("CO2 jump", "ppm", "ΔCO2"),
    wCO2 = info("Surface CO2 flux", "mg CO2 m-2 s-1", "w'CO2'ₛ"),
    wCO2A = info("CO2 flux of assimilation", "mg CO2 m-2 s-1", "w'CO2'ₐ"),
    wCO2R = info("CO2 flux of respiration", "mg CO2 m-2 s-1", "w'CO2'ᵣ"),
    wCO2e = info("Entrainment CO2 flux", "mg CO2 m-2 s-1", "w'CO2'ₑ"),
    wCO2M = info("Cumulus CO2 flux", "mg CO2 m-2 s-1", "w'CO2'ₘ"),
    S = info("Scalar", "S", "S"),
    dS = info("Scalar jump", "S", "ΔS"),
    wS = info("Surface scalar flux", "S m s-1", "w'S'ₛ"),
    wSe = info("Entrainment scalar flux", "S m s-1", "w'S'ₑ"),
    wSM = info("Cumulus scalar flux", "S m s-1", "w'S'ₘ"),
    u = info("Mixed-layer u wind", "m s-1", "u"),
    du = info("u wind jump", "m s-1", "Δu"),
    uw = info("Surface momentum flux in u", "m2 s-2", "u'w'ₛ"),
    v = info("Mixed-layer v wind", "m s-1", "v"),
    dv = info("v wind jump", "m s-1", "Δv"),
    vw = info("Surface momentum flux in v", "m2 s-2", "v'w'ₛ"),
    T2m = info("Temperature at 2 m", "K", "T2m"),
    q2m = info("Specific humidity at 2 m", "kg kg-1", "q2m"),
    u2m = info("u wind at 2 m", "m s-1", "u2m"),
    v2m = info("v wind at 2 m", "m s-1", "v2m"),
    e2m = info("Vapour pressure at 2 m", "Pa", "e2m"),
    esat2m = info("Saturated vapour pressure at 2 m", "Pa", "esat2m"),
    thetasurf = info("Surface potential temperature", "K", "θsurf"),
    thetavsurf = info("Surface virtual potential temperature", "K", "θv,surf"),
    qsurf = info("Surface specific humidity", "kg kg-1", "qsurf"),
    ustar = info("Friction velocity", "m s-1", "u*"),
    z0m = info("Roughness length for momentum", "m", "z0m"),
    z0h = info("Roughness length for scalars", "m", "z0h"),
    Cm = info("Drag coefficient for momentum", "-", "Cm"),
    Cs = info("Drag coefficient for scalars", "-", "Cs"),
    L = info("Obukhov length", "m", "L"),
    Rib = info("Bulk Richardson number", "-", "Rib"),
    Swin = info("Incoming shortwave radiation", "W m-2", "SW↓"),
    Swout = info("Outgoing shortwave radiation", "W m-2", "SW↑"),
    Lwin = info("Incoming longwave radiation", "W m-2", "LW↓"),
    Lwout = info("Outgoing longwave radiation", "W m-2", "LW↑"),
    Q = info("Net radiation", "W m-2", "Q"),
    cc = info("Cloud cover", "-", "cc"),
    ra = info("Aerodynamic resistance", "s m-1", "ra"),
    rs = info("Surface resistance", "s m-1", "rs"),
    H = info("Sensible heat flux", "W m-2", "H"),
    LE = info("Latent heat flux", "W m-2", "LE"),
    LEliq = info("Latent heat flux from wet vegetation", "W m-2", "LEliq"),
    LEveg = info("Latent heat flux from transpiration", "W m-2", "LEveg"),
    LEsoil = info("Latent heat flux from soil", "W m-2", "LEsoil"),
    LEpot = info("Potential latent heat flux", "W m-2", "LEpot"),
    LEref = info("Reference latent heat flux", "W m-2", "LEref"),
    G = info("Ground heat flux", "W m-2", "G"),
    zlcl = info("Lifting condensation level", "m", "zlcl"),
    RH_h = info("Relative humidity at h", "-", "RH(h)"),
    ac = info("Cloud core fraction", "-", "ac"),
    M = info("Cumulus mass flux", "m s-1", "M"),
    dz = info("Transition layer thickness", "m", "Δz"),
    phi = info("Photostationary state of NO, NO2 and O3", "-", "φ"),
    k_r05 = info("NO2 photolysis rate", "s-1", "k_r05")
)

const VARIABLE_INFO = merge(BASE_INFO,
    NamedTuple(n => info("$n mixing ratio", "ppb", String(n)) for n in SPECIES))

"""
    variable_info(name::Symbol) -> NamedTuple

Title, unit and symbol of the output variable `name`. A dimensionless variable has unit `"-"`.
"""
function variable_info(name::Symbol)
    haskey(VARIABLE_INFO, name) || throw(ArgumentError("unknown output variable :$name"))
    return VARIABLE_INFO[name]
end
