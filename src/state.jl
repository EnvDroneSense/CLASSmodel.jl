# Model state: the input, the constants and one state struct per component. Each
# struct holds the variables that are stored in the output; the components add the
# rest when they are implemented.

"Saturated vapor pressure [Pa] at temperature `T` [K]."
esat(T) = 0.611e3 * exp(17.2694 * (T - 273.16) / (T - 35.86))

"Saturated specific humidity [kg kg-1] at temperature `T` [K] and pressure `p` [Pa]."
qsat(T, p) = 0.622 * esat(T) / p

"State of the mixed layer: bulk values, surface and entrainment fluxes, top diagnostics."
Base.@kwdef mutable struct MixedLayerState{T<:AbstractFloat}
    h::T = NaN              # ABL height [m]
    ws::T = NaN             # large-scale vertical velocity [m s-1]
    wf::T = NaN             # growth due to radiative divergence [m s-1]
    we::T = NaN             # entrainment velocity [m s-1]
    wstar::T = NaN          # convective velocity scale [m s-1]

    theta::T = NaN          # potential temperature [K]
    dtheta::T = NaN         # potential temperature jump at h [K]
    wtheta::T = NaN         # surface kinematic heat flux [K m s-1]
    wthetae::T = NaN        # entrainment kinematic heat flux [K m s-1]
    thetav::T = NaN         # virtual potential temperature [K]
    dthetav::T = NaN        # virtual potential temperature jump at h [K]
    wthetav::T = NaN        # surface kinematic virtual heat flux [K m s-1]
    wthetave::T = NaN       # entrainment kinematic virtual heat flux [K m s-1]

    q::T = NaN              # specific humidity [kg kg-1]
    dq::T = NaN             # specific humidity jump at h [kg kg-1]
    wq::T = NaN             # surface kinematic moisture flux [kg kg-1 m s-1]
    wqe::T = NaN            # entrainment moisture flux [kg kg-1 m s-1]
    qsat::T = NaN           # saturated specific humidity [kg kg-1]
    esat::T = NaN           # saturated vapor pressure [Pa]
    e::T = NaN              # vapor pressure [Pa]

    CO2::T = NaN            # CO2 [ppm]
    dCO2::T = NaN           # CO2 jump at h [ppm]
    wCO2::T = NaN           # surface CO2 flux [ppm m s-1]
    wCO2A::T = 0            # surface assimilation CO2 flux [ppm m s-1]
    wCO2R::T = 0            # surface respiration CO2 flux [ppm m s-1]
    wCO2e::T = NaN          # entrainment CO2 flux [ppm m s-1]

    P_h::T = NaN            # pressure at the top [Pa]
    T_h::T = NaN            # absolute temperature at the top [K]
    RH_h::T = NaN           # relative humidity at the top [-]
    lcl::T = NaN            # lifting condensation level [m]

    htend::T = NaN          # tendency of h [m s-1]
    thetatend::T = NaN      # tendency of theta [K s-1]
    dthetatend::T = NaN     # tendency of dtheta [K s-1]
    qtend::T = NaN          # tendency of q [kg kg-1 s-1]
    dqtend::T = NaN         # tendency of dq [kg kg-1 s-1]
    CO2tend::T = NaN        # tendency of CO2 [ppm s-1]
    dCO2tend::T = NaN       # tendency of dCO2 [ppm s-1]
end

"State of the mixed-layer wind: components and jumps at h."
Base.@kwdef mutable struct WindState{T<:AbstractFloat}
    u::T = NaN              # u-wind [m s-1]
    du::T = NaN             # u-wind jump at h [m s-1]
    v::T = NaN              # v-wind [m s-1]
    dv::T = NaN             # v-wind jump at h [m s-1]
    utend::T = NaN          # tendency of u [m s-2]
    dutend::T = NaN         # tendency of du [m s-2]
    vtend::T = NaN          # tendency of v [m s-2]
    dvtend::T = NaN         # tendency of dv [m s-2]
end

"State of the radiation scheme."
Base.@kwdef mutable struct RadiationState{T<:AbstractFloat}
    Swin::T = NaN           # incoming shortwave radiation [W m-2]
    Swout::T = NaN          # outgoing shortwave radiation [W m-2]
    Lwin::T = NaN           # incoming longwave radiation [W m-2]
    Lwout::T = NaN          # outgoing longwave radiation [W m-2]
    Q::T = NaN              # net radiation [W m-2]
end

"State of the surface layer: friction, momentum fluxes, drag, stability, surface and 2 m values."
Base.@kwdef mutable struct SurfaceLayerState{T<:AbstractFloat}
    ustar::T = NaN          # friction velocity [m s-1]
    uw::T = NaN             # surface momentum flux in u-direction [m2 s-2]
    vw::T = NaN             # surface momentum flux in v-direction [m2 s-2]
    Cm::T = 1e12            # drag coefficient for momentum [-]
    Cs::T = 1e12            # drag coefficient for scalars [-]
    L::T = NaN              # Obukhov length [m]
    Rib::T = NaN            # bulk Richardson number [-]
    ra::T = NaN             # aerodynamic resistance [s m-1]

    thetasurf::T = NaN      # surface potential temperature [K]
    thetavsurf::T = NaN     # surface virtual potential temperature [K]
    qsurf::T = NaN          # surface specific humidity [kg kg-1]

    T2m::T = NaN            # 2 m temperature [K]
    q2m::T = NaN            # 2 m specific humidity [kg kg-1]
    e2m::T = NaN            # 2 m vapor pressure [Pa]
    esat2m::T = NaN         # 2 m saturated vapor pressure [Pa]
    u2m::T = NaN            # 2 m u-wind [m s-1]
    v2m::T = NaN            # 2 m v-wind [m s-1]
end

"State of the land surface: resistance and surface fluxes."
Base.@kwdef mutable struct LandSurfaceState{T<:AbstractFloat}
    Ts::T = NaN             # surface temperature [K]
    rs::T = 1e6             # resistance of transpiration [s m-1]
    H::T = NaN              # sensible heat flux [W m-2]
    LE::T = NaN             # evapotranspiration [W m-2]
    LEliq::T = NaN          # open water evaporation [W m-2]
    LEveg::T = NaN          # transpiration [W m-2]
    LEsoil::T = NaN         # soil evaporation [W m-2]
    LEpot::T = NaN          # potential evaporation [W m-2]
    LEref::T = NaN          # reference evaporation at rs = rsmin / LAI [W m-2]
    G::T = NaN              # ground heat flux [W m-2]
end

"State of the cumulus scheme: cloud core and transition layer."
Base.@kwdef mutable struct CumulusState{T<:AbstractFloat}
    ac::T = 0               # cloud core fraction [-]
    M::T = 0                # cloud core mass flux [m s-1]
    wqM::T = 0              # cloud core moisture flux [kg kg-1 m s-1]
    wCO2M::T = 0            # cloud core CO2 flux [ppm m s-1]
    dz_h::T = NaN           # transition layer thickness [m]
    dztend::T = NaN         # tendency of dz_h [m s-1]
end

"""
    ModelState{T}

State of a running model: the (fixed) `input` and `constants`, the step counter `t`
(number of completed steps) and the state of each component.
"""
mutable struct ModelState{T<:AbstractFloat}
    input::ModelInput{T}
    constants::PhysicalConstants{T}
    mixed_layer::MixedLayerState{T}
    wind::WindState{T}
    radiation::RadiationState{T}
    surface_layer::SurfaceLayerState{T}
    land_surface::LandSurfaceState{T}
    cumulus::CumulusState{T}
    t::Int
end

# Component functions. They do nothing yet; each is replaced when its component is added.
surface_layer!(s::ModelState) = nothing
land_surface!(s::ModelState) = nothing
cumulus!(s::ModelState) = nothing
integrate_land_surface!(s::ModelState) = nothing

"""
    init_model(input; constants = PhysicalConstants{T}()) -> ModelState

Build the initial state from `input` and compute the initial diagnostics.
"""
function init_model(x::ModelInput{T}; constants::PhysicalConstants{T}=PhysicalConstants{T}()) where {T}
    fac = constants.mair / (constants.rho * constants.mco2)   # mg CO2 m-2 s-1 to ppm m s-1
    s = ModelState{T}(
        x, constants,
        MixedLayerState{T}(; h=x.h, theta=x.theta, dtheta=x.dtheta, wtheta=x.wtheta,
            q=x.q, dq=x.dq, wq=x.wq, CO2=x.CO2, dCO2=x.dCO2, wCO2=x.wCO2 * fac),
        WindState{T}(; u=x.u, du=x.du, v=x.v, dv=x.dv),
        RadiationState{T}(; Q=x.Q),
        SurfaceLayerState{T}(; ustar=x.ustar, thetasurf=x.theta),
        LandSurfaceState{T}(; Ts=x.Ts),
        CumulusState{T}(; dz_h=x.dz_h),
        0,
    )

    statistics!(s)
    x.sw_rad && radiation!(s)
    if x.sw_sl
        for _ in 1:10
            surface_layer!(s)
        end
    end
    x.sw_ls && land_surface!(s)
    if x.sw_cu
        mixed_layer!(s)
        cumulus!(s)
    end
    x.sw_ml && mixed_layer!(s)
    return s
end
