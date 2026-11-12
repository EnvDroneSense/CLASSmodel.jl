# CLASSmodel.jl, a Julia implementation of the CLASS mixed-layer model
#
# Copyright (c) 2010-2015 Meteorology and Air Quality section, Wageningen University and Research centre
# Copyright (c) 2011-2015 Jordi Vila-Guerau de Arellano, Chiel van Heerwaarden, Bart van Stratum, Kees van den Dries
# Copyright (c) 2026 Thomas Struys
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. It is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.

"""
    CLASSmodel

CLASS, the Chemistry Land-surface Atmosphere Soil Slab model: a mixed-layer model of the
convective boundary layer coupled to a land surface.

```julia
using CLASSmodel
out = run_model(ModelInput(; runtime = 12 * 3600))   # a table with one column per variable
out.h                                                # boundary-layer height [m]
```

For a run step by step use [`init_model`](@ref) and [`step!`](@ref).
"""
module CLASSmodel

export ModelInput, PhysicalConstants, AGsParameters, PlantType, C3, C4, validate
export ModelOutput, run_model, mixed_layer_profile, variable_info, input_info
export ModelState, init_model, step!

include("io/input.jl")
include("io/input_info.jl")
include("core/state.jl")
include("physics/mixed_layer.jl")
include("physics/radiation.jl")
include("physics/surface_layer.jl")
include("physics/ags.jl")
include("physics/land_surface.jl")
include("physics/cumulus.jl")
include("io/output.jl")
include("io/output_info.jl")
include("core/model.jl")

end
