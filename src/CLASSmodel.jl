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

module CLASSmodel

export ModelInput, PhysicalConstants, AGsParameters, PlantType, C3, C4, validate
export ModelOutput, run_model

include("input.jl")
include("state.jl")
include("mixed_layer.jl")
include("radiation.jl")
include("surface_layer.jl")
include("model.jl")

end
