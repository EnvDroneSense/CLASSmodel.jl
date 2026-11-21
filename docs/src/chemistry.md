# Species and chemistry

```@setup quickstart
using CLASSmodel
```

## Chemical species

The input `species` lists the species; the default is the 22 species of `CLASSmodel.SPECIES`. The
vectors `C`, `dC`, `gammaC`, `advC` and `wC` hold one value for each species, in that order.
With another list, give all five vectors and the `reactions` that use those species. The mixing
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

The position in the table is the position in the input vectors `C`, `dC`, `gammaC`, `advC`
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

`sw_chem` runs the gas-phase reactions with an implicit Gauss-Seidel solver. The input
`reactions` is a vector of [`Reaction`](@ref). The default is `CLASSmodel.REACTIONS`, the 28
reactions of the C++ model. To turn a reaction off, leave it out of the vector. To add one,
append a `Reaction`. An empty vector gives transport only.

```julia
reactions = filter(r -> r.name != :r19, CLASSmodel.REACTIONS)
x = ModelInput(; sw_chem = true, reactions)
```

| Input | Meaning |
|:------|:--------|
| `sw_chem_constant` | the rates use the reference state, not the mixed layer (default on) |
| `sw_photo_constant` | the photolysis uses the sun at the hour `tod_ref` [h UTC] (default on) |
| `P_ref`, `Tcbl_ref`, `Tft_ref`, `qcbl_ref`, `qft_ref` | reference pressure, temperature and humidity of the mixed layer and the free troposphere |

The output has two extra columns: `phi`, the photostationary state of NO, NO2 and O3, and
`k_r05`, the NO2 photolysis rate [s-1]. Both find the reactions `:r05` and `:r21` by name. They
are `NaN` when the scheme does not have them.

```@example quickstart
x = ModelInput(; runtime = 3600, sw_chem = true)
out = run_model(x)
(NO2 = out.NO2[end], O3 = out.O3[end], phi = out.phi[end])
```

## List of reactions

The position in the table is the position in `CLASSmodel.REACTIONS`. The `form` of a reaction selects its rate expression, see [`Reaction`](@ref). The OH yield of r19 is its coefficient of OH; it is
0 in the default table. To change it, make a new `Reaction` for r19 with another coefficient.

```@eval
using CLASSmodel, Markdown
coef(c) = c == 1 ? "" : "$(isinteger(c) ? Int(c) : c) "
side(terms) = join(("$(coef(c))$s" for (s, c) in terms), " + ")
rows = map(enumerate(CLASSmodel.REACTIONS)) do (i, r)
    "| $i | $(r.name) | $(side(r.inputs)) -> $(side(r.outputs)) | $(r.form) |"
end
Markdown.parse(join(["| Position | Name | Reaction | Form |", "|:--|:--|:--|:--|", rows...], "\n"))
```
