# Outputs

```@setup quickstart
using CLASSmodel
```

## Output table

The output is a `NamedTuple` of vectors, so it is also a table in the sense of
[Tables.jl](https://tables.juliadata.org): `DataFrame(out)` and `CSV.write("run.csv", out)` work
once those packages are loaded. `out.t` is the time in hours since midnight of the start day.

## Units

`wCO2` and the fluxes `wCO2*` in the output are in mg CO2 m-2 s-1. The mixed-layer
concentrations `CO2` and `dCO2` are in ppm. The 22 species are in ppb.

## Output information

`variable_info(name)` returns the title, unit and symbol of an output variable.

```@example quickstart
variable_info(:k_r05)
```

## Output variables

The output has one column for each variable below, then one column for each species in ppb.

```@eval
using CLASSmodel, Markdown
names = filter(n -> !(n in CLASSmodel.SPECIES), CLASSmodel.OUTPUT_NAMES)
rows = map(names) do name
    info = variable_info(name)
    "| `$name` | $(info.title) | $(info.unit) |"
end
Markdown.parse(join(["| Column | Title | Unit |", "|:--|:--|:--|", rows...], "\n"))
```

## Mixed-layer profile

`mixed_layer_profile(input, out, i)` gives the corner points `(z, theta, q)` of the profiles at
output row `i`. The profile is constant up to `h`, has a jump at `h`, then follows the lapse
rates of the input up to `ztop` (default `2h`).

```@example quickstart
x = ModelInput()
out = run_model(x)
mixed_layer_profile(x, out, length(out.t))
```
