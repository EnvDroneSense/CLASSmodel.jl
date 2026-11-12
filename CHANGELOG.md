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

- `run_model` with `ls_type = :ags` throws a `DomainError` with a clear message when the vapour
  pressure deficit is above `f0 / ad`. Before, it threw a bare `DomainError` from `log`.

### Changes to model results

- A-Gs: the minimum leaf conductance `fmin` now uses the root of its quadratic, as in MXLCH and DALES.
  Results with `ls_type = :ags` move by 3e-11 or less (relative to the column maximum).

### Features

- Diurnal cycle of the prescribed surface fluxes: `sw_wtheta`, `sw_wq` and `sinperiod` make `wtheta`
  and `wq` follow a half sine (only without `sw_ls`).
- `run_model(input; progress, progress_every)`: a callback `progress(i, n, out)` reports the
  progress and stops the run when it returns `false`.
- `input_info(name)`: title, unit, symbol and slider range of each `ModelInput` field.
- `variable_info(name)`: title, unit and symbol of each output variable.
- `mixed_layer_profile(input, out, i)`: theta and q profiles of the mixed layer at output row `i`.

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
