# Soil and surface types of the CLASS model (modelgui, landsoil.cpp) as keywords of `ModelInput`.

const SOILS = (
    sandy_loam = (wsat = 0.472, wfc = 0.323, wwilt = 0.171, a = 0.219, b = 4.9, p = 4.0,
        CGsat = 3.56e-6, C1sat = 0.132, C2ref = 1.8),
    sand = (wsat = 0.403, wfc = 0.244, wwilt = 0.059, a = 0.387, b = 4.05, p = 4.0,
        CGsat = 3.222e-6, C1sat = 0.082, C2ref = 3.9),
    clay = (wsat = 0.614, wfc = 0.541, wwilt = 0.335, a = 0.083, b = 11.4, p = 12.0,
        CGsat = 3.6e-6, C1sat = 0.342, C2ref = 0.3)
)

# The trees have no plant type: they have no A-Gs parameters in the model GUI.
const SURFACES = (
    short_grass = (LAI = 2.0, gD = 0.0, rsmin = 110.0, alpha = 0.25, cveg = 0.85, Lambda = 5.9,
        z0m = 0.02, z0h = 0.002, c3c4 = :c3),
    maize = (LAI = 3.0, gD = 0.0, rsmin = 180.0, alpha = 0.2, cveg = 0.9, Lambda = 10.0,
        z0m = 0.15, z0h = 0.015, c3c4 = :c4),
    broadleaf_trees = (LAI = 5.0, gD = 0.03, rsmin = 200.0, alpha = 0.25, cveg = 0.9,
        Lambda = 20.0, z0m = 2.0, z0h = 2.0),
    needleleaf_trees = (LAI = 5.0, gD = 0.03, rsmin = 500.0, alpha = 0.25, cveg = 0.9,
        Lambda = 20.0, z0m = 2.0, z0h = 2.0)
)

"""
    soil(name) -> NamedTuple

Soil parameters of `ModelInput` for `name`: `:sandy_loam` (the default of `ModelInput`), `:sand`
or `:clay`. Use it as `ModelInput(; soil(:sand)...)`. It does not change the soil moisture
`wg` and `w2`.
"""
soil(name::Symbol) = preset(SOILS, name, "soil")

"""
    surface(name) -> NamedTuple

Vegetation and roughness parameters of `ModelInput` for `name`: `:short_grass` (the default of
`ModelInput`), `:maize`, `:broadleaf_trees` or `:needleleaf_trees`. Grass
and maize also set `c3c4`. Use it as `ModelInput(; surface(:maize)...)`.
"""
surface(name::Symbol) = preset(SURFACES, name, "surface")

function preset(table, name, what)
    haskey(table, name) ||
        throw(ArgumentError("unknown $what :$name, choose one of $(keys(table))"))
    return table[name]
end
