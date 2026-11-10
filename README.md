# CLASSmodel

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://EnvDroneSense.github.io/CLASSmodel.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://EnvDroneSense.github.io/CLASSmodel.jl/dev/)
[![Build Status](https://github.com/EnvDroneSense/CLASSmodel.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/EnvDroneSense/CLASSmodel.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/EnvDroneSense/CLASSmodel.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/EnvDroneSense/CLASSmodel.jl)

A Julia implementation of CLASS, the Chemistry Land-surface Atmosphere Soil Slab model: a mixed-layer model of the convective boundary layer coupled to a land surface. It is a port of [modelpy](https://github.com/classmodel/modelpy), the Python version of CLASS, and reproduces its results to round-off.

## Installation

```julia
using Pkg
Pkg.add("CLASSmodel")
```

## Usage

```julia
using CLASSmodel

x = ModelInput(; runtime = 12 * 3600, dt = 60, sw_rad = true, sw_sl = true, sw_ls = true)
out = run_model(x)       # a NamedTuple of vectors, also a table (Tables.jl)
out.h[end]               # boundary layer height at the end of the run [m]
```

`init_model` and `step!` advance the model one time step at a time. See the [documentation](https://EnvDroneSense.github.io/CLASSmodel.jl/stable/) for the switches, the output variables and their units.

## Licence

GPL-3.0-or-later, the licence of CLASS and modelpy. The copyright of the CLASS authors (Meteorology and Air Quality section, Wageningen University and Research centre; Jordi Vila-Guerau de Arellano, Chiel van Heerwaarden, Bart van Stratum, Kees van den Dries) is kept. See [LICENSE](LICENSE).
