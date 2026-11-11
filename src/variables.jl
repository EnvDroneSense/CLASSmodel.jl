# Title, unit and symbol of the output variables, for axis labels and legends.
#
# Units follow the comments of `ModelInput` and the units section of the documentation
# (`wCO2*` in mg CO2 m-2 s-1, `CO2` and `dCO2` in ppm). Dimensionless variables have "-".

const VARIABLE_INFO = (
    t = (title = "Time", unit = "h", symbol = "t"),
    h = (title = "Boundary-layer height", unit = "m", symbol = "h"),
    theta = (title = "Mixed-layer potential temperature", unit = "K", symbol = "θ"),
    thetav = (title = "Mixed-layer virtual potential temperature", unit = "K", symbol = "θv"),
    dtheta = (title = "Potential temperature jump", unit = "K", symbol = "Δθ"),
    dthetav = (title = "Virtual potential temperature jump", unit = "K", symbol = "Δθv"),
    wtheta = (title = "Surface kinematic heat flux", unit = "K m s-1", symbol = "w'θ'ₛ"),
    wthetav = (title = "Surface kinematic virtual heat flux", unit = "K m s-1", symbol = "w'θv'ₛ"),
    wthetae = (title = "Entrainment kinematic heat flux", unit = "K m s-1", symbol = "w'θ'ₑ"),
    wthetave = (title = "Entrainment kinematic virtual heat flux", unit = "K m s-1", symbol = "w'θv'ₑ"),
    q = (title = "Mixed-layer specific humidity", unit = "kg kg-1", symbol = "q"),
    dq = (title = "Specific humidity jump", unit = "kg kg-1", symbol = "Δq"),
    wq = (title = "Surface kinematic moisture flux", unit = "kg kg-1 m s-1", symbol = "w'q'ₛ"),
    wqe = (title = "Entrainment kinematic moisture flux", unit = "kg kg-1 m s-1", symbol = "w'q'ₑ"),
    wqM = (title = "Cumulus kinematic moisture flux", unit = "kg kg-1 m s-1", symbol = "w'q'ₘ"),
    qsat = (title = "Saturated specific humidity", unit = "kg kg-1", symbol = "qsat"),
    e = (title = "Vapour pressure", unit = "Pa", symbol = "e"),
    esat = (title = "Saturated vapour pressure", unit = "Pa", symbol = "esat"),
    CO2 = (title = "Mixed-layer CO2", unit = "ppm", symbol = "CO2"),
    dCO2 = (title = "CO2 jump", unit = "ppm", symbol = "ΔCO2"),
    wCO2 = (title = "Surface CO2 flux", unit = "mg CO2 m-2 s-1", symbol = "w'CO2'ₛ"),
    wCO2A = (title = "CO2 flux of assimilation", unit = "mg CO2 m-2 s-1", symbol = "w'CO2'ₐ"),
    wCO2R = (title = "CO2 flux of respiration", unit = "mg CO2 m-2 s-1", symbol = "w'CO2'ᵣ"),
    wCO2e = (title = "Entrainment CO2 flux", unit = "mg CO2 m-2 s-1", symbol = "w'CO2'ₑ"),
    wCO2M = (title = "Cumulus CO2 flux", unit = "mg CO2 m-2 s-1", symbol = "w'CO2'ₘ"),
    u = (title = "Mixed-layer u wind", unit = "m s-1", symbol = "u"),
    du = (title = "u wind jump", unit = "m s-1", symbol = "Δu"),
    uw = (title = "Surface momentum flux in u", unit = "m2 s-2", symbol = "u'w'ₛ"),
    v = (title = "Mixed-layer v wind", unit = "m s-1", symbol = "v"),
    dv = (title = "v wind jump", unit = "m s-1", symbol = "Δv"),
    vw = (title = "Surface momentum flux in v", unit = "m2 s-2", symbol = "v'w'ₛ"),
    T2m = (title = "Temperature at 2 m", unit = "K", symbol = "T2m"),
    q2m = (title = "Specific humidity at 2 m", unit = "kg kg-1", symbol = "q2m"),
    u2m = (title = "u wind at 2 m", unit = "m s-1", symbol = "u2m"),
    v2m = (title = "v wind at 2 m", unit = "m s-1", symbol = "v2m"),
    e2m = (title = "Vapour pressure at 2 m", unit = "Pa", symbol = "e2m"),
    esat2m = (title = "Saturated vapour pressure at 2 m", unit = "Pa", symbol = "esat2m"),
    thetasurf = (title = "Surface potential temperature", unit = "K", symbol = "θsurf"),
    thetavsurf = (title = "Surface virtual potential temperature", unit = "K", symbol = "θv,surf"),
    qsurf = (title = "Surface specific humidity", unit = "kg kg-1", symbol = "qsurf"),
    ustar = (title = "Friction velocity", unit = "m s-1", symbol = "u*"),
    z0m = (title = "Roughness length for momentum", unit = "m", symbol = "z0m"),
    z0h = (title = "Roughness length for scalars", unit = "m", symbol = "z0h"),
    Cm = (title = "Drag coefficient for momentum", unit = "-", symbol = "Cm"),
    Cs = (title = "Drag coefficient for scalars", unit = "-", symbol = "Cs"),
    L = (title = "Obukhov length", unit = "m", symbol = "L"),
    Rib = (title = "Bulk Richardson number", unit = "-", symbol = "Rib"),
    Swin = (title = "Incoming shortwave radiation", unit = "W m-2", symbol = "SW↓"),
    Swout = (title = "Outgoing shortwave radiation", unit = "W m-2", symbol = "SW↑"),
    Lwin = (title = "Incoming longwave radiation", unit = "W m-2", symbol = "LW↓"),
    Lwout = (title = "Outgoing longwave radiation", unit = "W m-2", symbol = "LW↑"),
    Q = (title = "Net radiation", unit = "W m-2", symbol = "Q"),
    ra = (title = "Aerodynamic resistance", unit = "s m-1", symbol = "ra"),
    rs = (title = "Surface resistance", unit = "s m-1", symbol = "rs"),
    H = (title = "Sensible heat flux", unit = "W m-2", symbol = "H"),
    LE = (title = "Latent heat flux", unit = "W m-2", symbol = "LE"),
    LEliq = (title = "Latent heat flux from wet vegetation", unit = "W m-2", symbol = "LEliq"),
    LEveg = (title = "Latent heat flux from transpiration", unit = "W m-2", symbol = "LEveg"),
    LEsoil = (title = "Latent heat flux from soil", unit = "W m-2", symbol = "LEsoil"),
    LEpot = (title = "Potential latent heat flux", unit = "W m-2", symbol = "LEpot"),
    LEref = (title = "Reference latent heat flux", unit = "W m-2", symbol = "LEref"),
    G = (title = "Ground heat flux", unit = "W m-2", symbol = "G"),
    zlcl = (title = "Lifting condensation level", unit = "m", symbol = "zlcl"),
    RH_h = (title = "Relative humidity at h", unit = "-", symbol = "RH(h)"),
    ac = (title = "Cloud core fraction", unit = "-", symbol = "ac"),
    M = (title = "Cumulus mass flux", unit = "m s-1", symbol = "M"),
    dz = (title = "Transition layer thickness", unit = "m", symbol = "Δz")
)

"""
    variable_info(name::Symbol) -> NamedTuple

Title, unit and symbol of the output variable `name`; the unit of a dimensionless variable is `"-"`.
"""
function variable_info(name::Symbol)
    haskey(VARIABLE_INFO, name) || throw(ArgumentError("unknown output variable :$name"))
    return VARIABLE_INFO[name]
end
