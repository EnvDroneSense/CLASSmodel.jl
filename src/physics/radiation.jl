# Radiation: solar geometry, bulk transmissivity and net radiation.

"""
    solar_elevation(x::ModelInput, time)

Sine of the solar elevation at `time` [s UTC since midnight of the start day], from the day of
the year, the latitude and the longitude. It is negative at night.
"""
function solar_elevation(x::ModelInput, time)
    sda = 0.409 * cos(2π * (x.doy - 173) / 365)
    hour_angle = 2π * time / DAY + 2π * x.lon / 360
    return sin(2π * x.lat / 360) * sin(sda) -
           cos(2π * x.lat / 360) * cos(sda) * cos(hour_angle)
end

"""
    radiation!(s::ModelState)

Shortwave and longwave radiation at the surface and the net radiation `Q`. The solar
elevation follows from the day of the year, the latitude, the longitude and the time of
day; the day of the year is not advanced during a run. The outgoing longwave radiation
uses the surface temperature of the land surface. The cloud cover is `rad.cc`: the input
`cc`, or `2 ac` of the previous step with `sw_curad`.
"""
function radiation!(s::ModelState)
    x, c, ml, rad = s.input, s.constants, s.mixed_layer, s.radiation

    sinlea = max(solar_elevation(x, s.t * x.dt + x.tstart * 3600), 0.0001)

    Ta = ml.theta * ((x.Ps - 0.1 * ml.h * c.rho * c.g) / x.Ps)^(c.Rd / c.cp)   # at 10 % of h
    Tr = (0.6 + 0.2 * sinlea) * (1 - 0.4 * rad.cc)

    rad.Swin = c.S0 * Tr * sinlea
    rad.Swout = x.alpha * c.S0 * Tr * sinlea
    rad.Lwin = 0.8 * c.bolz * Ta^4
    rad.Lwout = c.bolz * s.land_surface.Ts^4
    rad.Q = rad.Swin - rad.Swout + rad.Lwin - rad.Lwout
    return nothing
end
