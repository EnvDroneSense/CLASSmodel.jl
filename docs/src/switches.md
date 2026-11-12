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
| `sw_wtheta`, `sw_wq` | diurnal cycle of the prescribed fluxes (see below) |
| `sw_chem` | gas-phase chemistry of the 22 species (see [Species and chemistry](chemistry.md)) |

```@example quickstart
x = ModelInput(; sw_rad = true, sw_sl = true, sw_ls = true, ls_type = :ags, c3c4 = :c4)
out = run_model(x)
(H = maximum(out.H), LE = maximum(out.LE), wCO2A = minimum(out.wCO2A))
```

## Diurnal cycle of the surface fluxes

Without `sw_ls`, `wtheta` and `wq` are constant. With `sw_wtheta` or `sw_wq`, the flux is the
input value times `sin(pi t / sinperiod)`, with `t` in seconds from the start of the run.
`sinperiod` is the half period of the sine [s]. The default is 43200 s.

```@example quickstart
x = ModelInput(; sw_wtheta = true, sw_wq = true, sinperiod = 12 * 3600)
out = run_model(x)
maximum(out.wtheta)
```
