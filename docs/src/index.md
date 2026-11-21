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
