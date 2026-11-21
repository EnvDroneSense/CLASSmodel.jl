# Switches

```@setup quickstart
using CLASSmodel
```

## Switches

Only `sw_ml` is on by default. Switch the other components on in the input.

| Switch | Component |
|:-------|:----------|
| `sw_ml` | mixed-layer growth, entrainment and tendencies (default on) |
| `sw_shearwe` | shear contribution to entrainment |
| `sw_fixft` | free troposphere fixed under subsidence |
| `sw_wind` | prognostic wind |
| `sw_sl` | surface layer: `ustar`, drag coefficients, 2 m diagnostics |
| `sw_rad` | radiation: solar angle, net radiation |
| `sw_ls` | land surface, with `ls_type = :js` (Jarvis-Stewart) or `:ags` (A-Gs, needs `sw_rad`) |
| `sw_sea` | with `sw_ls`: open water at the fixed `Ts` instead of land. `G` and `LEref` are `NaN` |
| `sw_cu` | shallow cumulus: mass flux and cloud core fraction |
| `sw_curad` | with `sw_cu` and `sw_rad`: the cloud cover of the radiation is `2 ac` of the previous step |
| `sw_chem` | gas-phase chemistry of the 22 species (see [Species and chemistry](chemistry.md)) |

```@example quickstart
x = ModelInput(; sw_rad = true, sw_sl = true, sw_ls = true, ls_type = :ags, c3c4 = :c4)
out = run_model(x)
(H = maximum(out.H), LE = maximum(out.LE), wCO2A = minimum(out.wCO2A))
```

## Surface fluxes that change in time

Without `sw_ls`, `wtheta` and `wq` are prescribed. `wtheta`, `wq`, `wCO2` and each entry of
`wC` take a number (constant flux) or a [`FluxShape`](@ref). `t` is the time in seconds from
the start of the run. `period` and `start` are in seconds; the defaults are 43200 and 0.

| Shape | Flux |
|:------|:-----|
| `Sine(a; period, start)` | `a * sin(pi (t - start) / period)` in the window, zero outside |
| `Cosine(a; period, start)` | `a / 2 * (1 - cos(2 pi (t - start) / period))` in the window, zero outside |
| `Block(a; period, start)` | `a` in the window, zero outside |
| `Table(values; interval)` | linear interpolation of `values`, one every `interval` seconds |
| `Sun(a)` | `a` times the solar elevation divided by its noon value, zero at night |
| `Deposition(vd)` | `-vd` times the mixed-layer value, with `vd` in m s-1 |

With `sw_ls`, a shape for `wtheta` or `wq` is an error: the land surface gives these fluxes.

```@example quickstart
x = ModelInput(; wtheta = Sine(0.1; period = 12 * 3600), wq = Sine(0.1e-3; period = 12 * 3600))
out = run_model(x)
maximum(out.wtheta)
```
