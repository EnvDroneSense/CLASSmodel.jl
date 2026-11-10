```@meta
CurrentModule = CLASSmodel
```

# CLASSmodel

CLASS, the Chemistry Land-surface Atmosphere Soil Slab model, is a mixed-layer model of the
convective boundary layer coupled to a land surface. This package is a Julia version of the
Python model [modelpy](https://github.com/classmodel/modelpy).

## Quick start

```@example quickstart
using CLASSmodel

x = ModelInput(; runtime = 12 * 3600, dt = 60)   # the default case: a dry, growing mixed layer
out = run_model(x)
(h = out.h[end], theta = out.theta[end], q = out.q[end])
```

The output is a `NamedTuple` of vectors, so it is also a table in the sense of
[Tables.jl](https://tables.juliadata.org): `DataFrame(out)` and `CSV.write("run.csv", out)` work
once those packages are loaded. `out.t` is the time in hours since midnight of the start day.

## Switches

Every component is off by default and is switched on in the input.

| Switch | Component |
|:-------|:----------|
| `sw_ml` | mixed-layer growth, entrainment and tendencies (default on) |
| `sw_shearwe` | shear contribution to entrainment |
| `sw_fixft` | free troposphere fixed under subsidence |
| `sw_wind` | prognostic wind |
| `sw_sl` | surface layer: `ustar`, drag coefficients, 2 m diagnostics |
| `sw_rad` | radiation: solar angle, net radiation |
| `sw_ls` | land surface, with `ls_type = :js` (Jarvis-Stewart) or `:ags` (A-Gs, needs `sw_rad`) |
| `sw_cu` | shallow cumulus: mass flux and cloud core fraction |

```@example quickstart
x = ModelInput(; sw_rad = true, sw_sl = true, sw_ls = true, ls_type = :ags, c3c4 = :c4)
out = run_model(x)
(H = maximum(out.H), LE = maximum(out.LE), wCO2A = minimum(out.wCO2A))
```

## Step by step

```@example quickstart
s = init_model(ModelInput(; runtime = 600))
step!(s)
s.mixed_layer.h
```

## Units

`wCO2` and the fluxes `wCO2*` in the output are in mg CO2 m-2 s-1; the mixed-layer
concentrations `CO2` and `dCO2` are in ppm.

## Reference

```@index
```

```@autodocs
Modules = [CLASSmodel]
Private = false
```
