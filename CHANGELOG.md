# Changelog

Notable changes to CLASSmodel.jl, in the spirit of
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Each version uses
these sections, in this order, and leaves out the empty ones:

| Section | Holds |
|:--|:--|
| Breaking changes | What a caller has to change: renamed or removed names, new signatures, new errors |
| Changes to model results | Every change to the numbers a run returns |
| Features | New things a user can do |
| Changes | Other changes a user can notice |
| Fixes | Wrong behaviour corrected, results unchanged |
| Documentation | The README, the docs site and the docstrings |
| Tools and tests | CI, the test suite and maintainer scripts |

An entry outside "Changes to model results" leaves results unchanged.

## [Unreleased]

### Breaking changes

- The inputs `sw_reactions` and `stocoef` are removed. The input `reactions` is a vector of `Reaction`
  (default `CLASSmodel.REACTIONS`). Leave a reaction out to switch it off; the OH yield of r19 is the
  coefficient of OH in its `Reaction`. `phi` and `k_r05` are `NaN` when the scheme has no `:r05` or `:r21`.
- The inputs `sinperiod`, `sw_wtheta`, `sw_wq` and `sw_wC` are removed. Give `wtheta`, `wq`,
  `wCO2` or an entry of `wC` a `FluxShape` instead: `sw_wtheta = true` with `sinperiod = P`
  becomes `wtheta = Sine(wtheta; period = P)`, and `sw_wC` 2 becomes `Sun(wC)`.

### Features

- The output columns `qsat`, `e` and `esat` have values in runs without `sw_ls`. Before, they were `NaN`.
- `Reaction` is exported: a user can give `reactions` a scheme of their own. Its `form` is a symbol such as
  `:photolysis` or `:arrhenius`. `ModelInput` checks the reactions, finite flux shapes and `Segments`, and
  `gammatheta > 0` in a run without an inversion.
- Encroachment: a run with `dtheta = 0`, `beta = 0` and `dq = 0` has no inversion and the layer grows
  into the stable air above it, `h^2 = h0^2 + 2 wtheta t / gammatheta`. Before, `dtheta = 0` was an error.
  It is still an error with `sw_cu`.
- `Segments(z, gamma)`: a lapse rate that changes with height for `gammatheta`, `gammaq`,
  `gammaCO2`, `gammau` and `gammav`. `mixed_layer_profile` follows the segments.
- Surface flux shapes `Sine`, `Cosine`, `Block`, `Table`, `Sun` and `Deposition` for `wtheta`,
  `wq`, `wCO2` and each entry of `wC`. `Sine` is zero outside its window.

## [0.2.0] - 2026-11-18

### Breaking changes

- `run_model` with `ls_type = :ags` throws a `DomainError` with a clear message when the vapour
  pressure deficit is above `f0 / ad`. Before, it threw a bare `DomainError` from `log`.

### Changes to model results

- A-Gs: the minimum leaf conductance `fmin` now uses the root of its quadratic, as in MXLCH and DALES.
  Results with `ls_type = :ags` move by 3e-11 or less (relative to the column maximum).

### Features

- Diurnal cycle of the prescribed surface fluxes: `sw_wtheta`, `sw_wq` and `sinperiod` make `wtheta`
  and `wq` follow a half sine (only without `sw_ls`).
- 22 chemical species: inputs `C`, `dC`, `gammaC`, `advC`, `wC` and `sw_wC` (tuples with one
  value for each species of `CLASSmodel.SPECIES`, in ppb) and 22 output columns, one for each species.
  The species are transported like CO2; `Inert` takes part in no reaction and is the passive
  scalar. `sw_wC` is 0 (constant), 1 (sine) or 2 (sun).
- Gas-phase chemistry of the C++ model: `sw_chem` runs its 28 reactions (`sw_reactions` switches each one)
  with an implicit Gauss-Seidel solver. `sw_chem_constant` and `sw_photo_constant` fix the rates at a
  reference state (`P_ref`, `Tcbl_ref`, `Tft_ref`, `qcbl_ref`, `qft_ref`, `tod_ref`); `stocoef` is the OH yield of
  RO2 + HO2. New output columns `phi` (photostationary state) and `k_r05` (NO2 photolysis rate).
  Compared with the C++ model, the powers use double precision (the C++ model uses single precision
  in a few terms, a relative difference of 1e-8) and the photolysis of the free troposphere is zero at night
  (the C++ model keeps its last daytime value). With `sw_chem_constant` and `sw_photo_constant` the
  difference after sunset is up to 10 % in some species.
- `sw_curad`: with `sw_cu` and `sw_rad`, the cloud cover of the radiation is `2 ac` of the previous step.
  New output column `cc`.
- `sw_sea`: with `sw_ls`, the land surface becomes open water at the fixed `Ts` (no soil, no
  land-input checks). `G` and `LEref` are `NaN`.
- `run_model(input; progress, progress_every)`: a callback `progress(i, n, out)` reports the
  progress and stops the run when it returns `false`.
- `input_info(name)`: title, unit, symbol and slider range of each `ModelInput` field.
- `variable_info(name)`: title, unit and symbol of each output variable.
- `mixed_layer_profile(input, out, i)`: theta and q profiles of the mixed layer at output row `i`.
- `run_model` and `init_model` take a `plant` keyword: an `AGsParameters` that replaces the
  preset of `c3c4`.

### Tools and tests

- CompatHelper installs Julia before it runs. The daily run failed with "julia: command not found".
- `[compat]` entries for the test and docs environments (Aqua, JET, Documenter, CLASSmodel). The CI
  actions `checkout`, `julia-actions/cache` and `codecov-action` move to v7, v3 and v7.

## [0.1.0] - 2026-11-10

First release. A port of modelpy, the Python version of CLASS, that reproduces its
results to round-off.

### Features

- `ModelInput`, `PhysicalConstants` and `AGsParameters`, checked by `validate`;
  the numeric type can be any `T <: Real`.
- `run_model(input)` returns a `ModelOutput`, a table with one column per variable.
  `init_model` and `step!` run the model one step at a time.
- Mixed-layer model, radiation, surface layer, land surface, the A-Gs
  photosynthesis model (C3 and C4) and the cumulus scheme, each with its switch.

### Documentation

- README with install and usage, and a documentation site.

### Tools and tests

- Aqua and JET tests, and a CI check of the formatting (JuliaFormatter).
