```@meta
CurrentModule = CLASSmodel
```

# CLASSmodel

CLASS, the Chemistry Land-surface Atmosphere Soil Slab model, is a mixed-layer model of the
convective boundary layer coupled to a land surface. This package is a Julia version of the
Python model [modelpy](https://github.com/classmodel/modelpy), with added chemistry.

## Quick start

```@example quickstart
using CLASSmodel

x = ModelInput(; runtime = 12 * 3600, dt = 60)   # the default case: a dry, growing mixed layer
out = run_model(x)
(h = out.h[end], theta = out.theta[end], q = out.q[end])
```

## Contents

| Page | Content |
|:-----|:--------|
| [Switches](switches.md) | the components and their switches, the diurnal cycle of the fluxes |
| [Species and chemistry](chemistry.md) | the 22 species, their surface fluxes, the reactions |
| [Inputs](inputs.md) | all fields of `ModelInput` with unit and default |
| [Outputs](outputs.md) | all output columns with unit, the mixed-layer profile |
| [Running the model](running.md) | a run step by step, progress and early stop |
| [API](api.md) | the exported names and their docstrings |
