# Thermodynamic constants and helpers shared by the components.

const ε = 0.622     # ratio of the molar masses of water vapor and dry air [-]
const VIRT = 0.61   # virtual temperature coefficient, Rv / Rd - 1 [-]
const DAY = 86400   # length of a day [s]

"Saturated vapor pressure [Pa] at temperature `T` [K]."
esat(T) = 0.611e3 * exp(17.2694 * (T - 273.16) / (T - 35.86))

"Derivative of [`esat`](@ref) with respect to the temperature [Pa K-1]."
desat(T) = esat(T) * 17.2694 * (273.16 - 35.86) / (T - 35.86)^2

"Saturated specific humidity [kg kg-1] at temperature `T` [K] and pressure `p` [Pa]."
qsat(T, p) = ε * esat(T) / p

# Factor from a CO2 flux in mg CO2 m-2 s-1 to ppm m s-1.
mg_to_ppm(c) = c.mair / (c.rho * c.mco2)
