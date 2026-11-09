# Model input, physical constants and A-Gs plant parameters.

"""
    PlantType

Plant type of the A-Gs scheme: [`C3`](@ref) or [`C4`](@ref).
"""
abstract type PlantType end

"C3 plant type of the A-Gs scheme."
struct C3 <: PlantType end

"C4 plant type of the A-Gs scheme."
struct C4 <: PlantType end

# Element type for a constructor called without `{T}`: the promoted type of the numbers
# given, at least Float64 for integers. The keyword constructors made by `@kwdef` end up
# here, so `ModelInput(; dt=60, h=200.5)` is a `ModelInput{Float64}`.
_float_type(args) = float(promote_type((typeof(a) for a in args if a isa Real && !(a isa Bool))...))

"""
    PhysicalConstants{T}(; kwargs...)

Physical constants of the model, with the usual values as defaults.
"""
Base.@kwdef struct PhysicalConstants{T<:AbstractFloat}
    Lv::T = 2.5e6           # heat of vaporization [J kg-1]
    cp::T = 1005            # specific heat of dry air [J kg-1 K-1]
    rho::T = 1.2            # density of air [kg m-3]
    k::T = 0.4              # Von Karman constant [-]
    g::T = 9.81             # gravity acceleration [m s-2]
    Rd::T = 287             # gas constant for dry air [J kg-1 K-1]
    Rv::T = 461.5           # gas constant for moist air [J kg-1 K-1]
    bolz::T = 5.67e-8       # Stefan-Boltzmann constant [W m-2 K-4]
    rhow::T = 1000          # density of water [kg m-3]
    S0::T = 1368            # solar constant [W m-2]
    mco2::T = 44            # molecular weight CO2 [g mol-1]
    mair::T = 28.9          # molecular weight air [g mol-1]
    nuco2q::T = 1.6         # ratio of molecular viscosity of water to carbon dioxide [-]
    Cw::T = 0.0016          # constant of the water stress correction of soil respiration [-]
    wmax::T = 0.55          # upper reference value of soil water [-]
    wmin::T = 0.005         # lower reference value of soil water [-]
    R10::T = 0.23           # respiration at 10 C [mg CO2 m-2 s-1]
    E0::T = 53.3e3          # activation energy [J kmol-1]
end
PhysicalConstants(args...) = PhysicalConstants{_float_type(args)}(args...)

"""
    AGsParameters{T}(plant)

Parameters of the A-Gs scheme for one plant type, `plant` being [`C3`](@ref)`()`
or [`C4`](@ref)`()`.
"""
Base.@kwdef struct AGsParameters{T<:AbstractFloat}
    CO2comp298::T   # CO2 compensation concentration [mg m-3]
    Q10CO2::T       # parameter of the CO2 compensation concentration [-]
    gm298::T        # mesophyll conductance at 298 K [mm s-1]
    Ammax298::T     # maximal primary productivity [mg m-2 s-1]
    Q10gm::T        # parameter of the mesophyll conductance [-]
    T1gm::T         # reference temperature of the mesophyll conductance [K]
    T2gm::T         # reference temperature of the mesophyll conductance [K]
    Q10Am::T        # parameter of the maximal primary productivity [-]
    T1Am::T         # reference temperature of the maximal primary productivity [K]
    T2Am::T         # reference temperature of the maximal primary productivity [K]
    f0::T           # maximum value of Cfrac [-]
    ad::T           # regression coefficient of Cfrac [kPa-1]
    alpha0::T       # initial low light conditions [mg J-1]
    Kx::T           # extinction coefficient of PAR [-]
    gmin::T         # cuticular (minimum) conductance [mm s-1]
end

AGsParameters(plant::PlantType) = AGsParameters{Float64}(plant)
AGsParameters{T}(::C3) where {T} = AGsParameters{T}(;
    CO2comp298=68.5, Q10CO2=1.5, gm298=7.0, Ammax298=2.2, Q10gm=2.0, T1gm=278,
    T2gm=301, Q10Am=2.0, T1Am=281, T2Am=311, f0=0.89, ad=0.07, alpha0=0.017,
    Kx=0.7, gmin=0.25e-3)
AGsParameters{T}(::C4) where {T} = AGsParameters{T}(;
    CO2comp298=4.3, Q10CO2=1.5, gm298=17.5, Ammax298=1.7, Q10gm=2.0, T1gm=286,
    T2gm=309, Q10Am=2.0, T1Am=286, T2Am=311, f0=0.85, ad=0.15, alpha0=0.014,
    Kx=0.7, gmin=0.25e-3)

"""
    ModelInput{T}(; kwargs...)
    ModelInput(; kwargs...)    # T from the numbers given: Float64 unless Float32 only

All input of a run. The defaults describe a dry mixed layer with fixed surface
fluxes. The input is checked when it is constructed (see [`validate`](@ref)).

| Group | Fields |
|:------|:-------|
| time | `runtime`, `dt` |
| switches | `sw_ml`, `sw_shearwe`, `sw_fixft`, `sw_wind`, `sw_sl`, `sw_rad`, `sw_ls`, `sw_cu` |
| mixed layer | `h`, `Ps`, `divU`, `fc`, `beta` |
| heat, moisture, CO2 | `theta`, `dtheta`, `gammatheta`, `advtheta`, `wtheta`, `q`, `dq`, `gammaq`, `advq`, `wq`, `CO2`, `dCO2`, `gammaCO2`, `advCO2`, `wCO2` |
| wind | `u`, `du`, `gammau`, `advu`, `v`, `dv`, `gammav`, `advv` |
| surface layer | `ustar`, `z0m`, `z0h` |
| radiation | `lat`, `lon`, `doy`, `tstart`, `cc`, `Q`, `dFz` |
| land surface | `ls_type`, `c3c4`, `c_beta`, soil (`wg`, `w2`, `Tsoil`, `T2`, `a`, `b`, `p`, `CGsat`, `wsat`, `wfc`, `wwilt`, `C1sat`, `C2ref`), vegetation (`LAI`, `gD`, `rsmin`, `rssoilmin`, `alpha`, `cveg`, `Wmax`, `Wl`), `Ts`, `Lambda` |
| cumulus | `dz_h` |

Notes:
- `ls_type` is `:js` (Jarvis-Stewart) or `:ags` (A-Gs); `c3c4` is `:c3` or `:c4`.
- `wCO2` is in mg CO2 m-2 s-1 (converted with `mair / (rho * mco2)`).
- `c_beta = 0` is a linear water stress response.
"""
Base.@kwdef struct ModelInput{T<:AbstractFloat}
    # time
    runtime::T = 12 * 3600      # duration of the run [s]
    dt::T = 60                  # time step [s]

    # switches
    sw_ml::Bool = true          # mixed-layer model
    sw_shearwe::Bool = false    # shear growth of the mixed layer
    sw_fixft::Bool = false      # fixed free-troposphere
    sw_wind::Bool = false       # prognostic wind
    sw_sl::Bool = false         # surface layer
    sw_rad::Bool = false        # radiation
    sw_ls::Bool = false         # land surface
    sw_cu::Bool = false         # cumulus parameterization

    # mixed layer
    h::T = 200                  # initial ABL height [m]
    Ps::T = 101300              # surface pressure [Pa]
    divU::T = 0                 # horizontal large-scale divergence of wind [s-1]
    fc::T = 1e-4                # Coriolis parameter [s-1]
    beta::T = 0.2               # entrainment ratio for virtual heat [-]

    # heat
    theta::T = 288              # initial mixed-layer potential temperature [K]
    dtheta::T = 1               # initial temperature jump at h [K]
    gammatheta::T = 0.006       # free atmosphere potential temperature lapse rate [K m-1]
    advtheta::T = 0             # advection of heat [K s-1]
    wtheta::T = 0.1             # surface kinematic heat flux [K m s-1]

    # moisture
    q::T = 0.008                # initial mixed-layer specific humidity [kg kg-1]
    dq::T = -0.001              # initial specific humidity jump at h [kg kg-1]
    gammaq::T = 0               # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
    advq::T = 0                 # advection of moisture [kg kg-1 s-1]
    wq::T = 0.1e-3              # surface kinematic moisture flux [kg kg-1 m s-1]

    # CO2
    CO2::T = 422                # initial mixed-layer CO2 [ppm]
    dCO2::T = -44               # initial CO2 jump at h [ppm]
    gammaCO2::T = 0             # free atmosphere CO2 lapse rate [ppm m-1]
    advCO2::T = 0               # advection of CO2 [ppm s-1]
    wCO2::T = 0                 # surface kinematic CO2 flux [mg CO2 m-2 s-1, see notes]

    # wind
    u::T = 6                    # initial mixed-layer u-wind speed [m s-1]
    du::T = 4                   # initial u-wind jump at h [m s-1]
    gammau::T = 0               # free atmosphere u-wind lapse rate [s-1]
    advu::T = 0                 # advection of u-wind [m s-2]
    v::T = -4                   # initial mixed-layer v-wind speed [m s-1]
    dv::T = 4                   # initial v-wind jump at h [m s-1]
    gammav::T = 0               # free atmosphere v-wind lapse rate [s-1]
    advv::T = 0                 # advection of v-wind [m s-2]

    # surface layer
    ustar::T = 0.3              # surface friction velocity [m s-1]
    z0m::T = 0.02               # roughness length for momentum [m]
    z0h::T = 0.002              # roughness length for scalars [m]

    # radiation
    lat::T = 51.97              # latitude [deg]
    lon::T = -4.93              # longitude [deg]
    doy::T = 268                # day of the year [-]
    tstart::T = 6.8             # start time of the day [h UTC]
    cc::T = 0                   # cloud cover fraction [-]
    Q::T = 400                  # net radiation [W m-2]
    dFz::T = 0                  # cloud top radiative divergence [W m-2]

    # land surface
    ls_type::Symbol = :js       # :js (Jarvis-Stewart) or :ags (A-Gs)
    wg::T = 0.21                # volumetric water content top soil layer [m3 m-3]
    w2::T = 0.21                # volumetric water content deeper soil layer [m3 m-3]
    Tsoil::T = 285              # temperature top soil layer [K]
    T2::T = 286                 # temperature deeper soil layer [K]
    a::T = 0.219                # Clapp and Hornberger retention curve parameter a [-]
    b::T = 4.90                 # Clapp and Hornberger retention curve parameter b [-]
    p::T = 4                    # Clapp and Hornberger retention curve parameter p [-]
    CGsat::T = 3.56e-6          # saturated soil conductivity for heat [K m2 J-1]
    wsat::T = 0.472             # saturated volumetric water content [-]
    wfc::T = 0.323              # volumetric water content at field capacity [-]
    wwilt::T = 0.171            # volumetric water content at wilting point [-]
    C1sat::T = 0.132            # force-restore coefficient C1 at saturation [-]
    C2ref::T = 1.8              # force-restore coefficient C2 at reference [-]
    c_beta::T = 0               # curvature of the plant water-stress factor (0..1) [-]
    LAI::T = 2                  # leaf area index [-]
    gD::T = 0                   # correction factor of transpiration for VPD [-]
    rsmin::T = 110              # minimum resistance of transpiration [s m-1]
    rssoilmin::T = 50           # minimum resistance of soil evaporation [s m-1]
    alpha::T = 0.25             # surface albedo [-]
    cveg::T = 0.85              # vegetation fraction [-]
    Ts::T = 290                 # initial surface temperature [K]
    Wmax::T = 0.0002            # thickness of water layer on wet vegetation [m]
    Wl::T = 0.0                 # equivalent water layer depth on wet vegetation [m]
    Lambda::T = 5.9             # thermal diffusivity of the skin layer [W m-2 K-1]

    # A-Gs
    c3c4::Symbol = :c3          # plant type, :c3 or :c4

    # cumulus
    dz_h::T = 150               # transition layer thickness [m]

    function ModelInput{T}(args...) where {T<:AbstractFloat}
        input = new{T}(args...)
        validate(input)
        return input
    end
end
ModelInput(args...) = ModelInput{_float_type(args)}(args...)

"""
    validate(input::ModelInput)

Throw an `ArgumentError` if `input` cannot be run. Called by the
[`ModelInput`](@ref) constructor. Checks:

- positive time step, run time, `h`, `Ps`, `z0m`, `z0h`
- `0 <= cc <= 1`, `0 <= cveg <= 1`, `0 <= c_beta <= 1`
- `ls_type` and `c3c4` have a known value
- with `sw_ls`: `LAI > 0` and `wwilt < wfc <= wsat`
- A-Gs (`sw_ls` with `ls_type = :ags`) needs `sw_rad`, because it uses the
  incoming shortwave radiation
"""
function validate(x::ModelInput)
    positive(name, v) = v > 0 || throw(ArgumentError("$name must be > 0, got $v"))
    fraction(name, v) = 0 <= v <= 1 || throw(ArgumentError("$name must be in [0, 1], got $v"))

    positive(:runtime, x.runtime)
    positive(:dt, x.dt)
    positive(:h, x.h)
    positive(:Ps, x.Ps)
    positive(:z0m, x.z0m)
    positive(:z0h, x.z0h)
    x.runtime >= x.dt || throw(ArgumentError("runtime must be at least dt"))
    fraction(:cc, x.cc)
    fraction(:cveg, x.cveg)
    fraction(:c_beta, x.c_beta)
    x.ls_type in (:js, :ags) || throw(ArgumentError("ls_type must be :js or :ags, got :$(x.ls_type)"))
    x.c3c4 in (:c3, :c4) || throw(ArgumentError("c3c4 must be :c3 or :c4, got :$(x.c3c4)"))
    if x.sw_ls
        positive(:LAI, x.LAI)
        x.wwilt < x.wfc <= x.wsat || throw(ArgumentError("need wwilt < wfc <= wsat"))
        x.ls_type === :ags && !x.sw_rad &&
            throw(ArgumentError("ls_type = :ags needs sw_rad = true"))
    end
    return x
end
