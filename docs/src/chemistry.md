# Species and chemistry

```@setup quickstart
using CLASSmodel
```

## Chemical species

The model carries 22 species. Their names are in `CLASSmodel.SPECIES`. The input tuples `C`, `dC`,
`gammaC`, `advC` and `wC` hold one value for each species, in that order. The mixing
ratios are in ppb and the surface fluxes are in ppb m s-1. The output has one column for each
species. The species are transported like CO2.

Each entry of `wC` is a number (constant flux) or a [`FluxShape`](@ref), for example
`Sun(0.1)` or `Deposition(0.01)`. See [Surface fluxes that change in time](@ref).

`Inert` takes part in no reaction. It is the passive scalar. Without `sw_chem`, all species are
passive.

```@example quickstart
x = ModelInput(; runtime = 3600, wtheta = Sine(0.1; period = 12 * 3600))
out = run_model(x)
(Inert = out.Inert[end], CH4 = out.CH4[end])
```

## List of species

The position in the table is the position in the input tuples `C`, `dC`, `gammaC`, `advC`
and `wC`.

```@eval
using CLASSmodel, Markdown
x = ModelInput()
rows = map(enumerate(CLASSmodel.SPECIES)) do (i, name)
    "| $i | `$name` | $(x.C[i]) | $(x.wC[i]) |"
end
Markdown.parse(join(["| Position | Species | Default `C` [ppb] | Default `wC` [ppb m s-1] |",
    "|:--|:--|:--|:--|", rows...], "\n"))
```

## Chemistry

`sw_chem` runs the gas-phase reactions with an implicit Gauss-Seidel solver. The scheme has 28
reactions. `sw_reactions` holds one switch for each of them. Set a switch to `false` to turn
that reaction off.

| Input | Meaning |
|:------|:--------|
| `sw_chem_constant` | the rates use the reference state, not the mixed layer (default on) |
| `sw_photo_constant` | the photolysis uses the sun at the hour `tod_ref` [h UTC] (default on) |
| `P_ref`, `Tcbl_ref`, `Tft_ref`, `qcbl_ref`, `qft_ref` | reference pressure, temperature and humidity of the mixed layer and the free troposphere |
| `stocoef` | OH yield of the reaction RO2 + HO2 |

The output has two extra columns: `phi`, the photostationary state of NO, NO2 and O3, and
`k_r05`, the NO2 photolysis rate [s-1].

```@example quickstart
x = ModelInput(; runtime = 3600, sw_chem = true)
out = run_model(x)
(NO2 = out.NO2[end], O3 = out.O3[end], phi = out.phi[end])
```

## List of reactions

The position in the table is the position in `sw_reactions`. A photolysis rate depends on the
solar elevation and is zero at night. The OH yield of r19 is the input `stocoef`.

```@eval
using CLASSmodel, Markdown
coef(c) = c == 0 ? "stocoef " : c == 1 ? "" : "$(isinteger(c) ? Int(c) : c) "   # r19 OH is stocoef
side(terms) = join(("$(coef(c))$s" for (s, c) in terms), " + ")
rows = map(enumerate(CLASSmodel.REACTIONS)) do (i, r)
    "| $i | $(r.name) | $(side(r.inputs)) -> $(side(r.outputs)) | $(r.photolysis ? "yes" : "no") |"
end
Markdown.parse(join(["| Position | Name | Reaction | Photolysis |", "|:--|:--|:--|:--|", rows...], "\n"))
```
