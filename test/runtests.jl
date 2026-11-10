using CLASSmodel
using Test

@testset "CLASSmodel.jl" begin
    include("test_input.jl")
    include("test_model.jl")
    include("test_mixed_layer.jl")
    include("test_radiation.jl")
    include("test_surface_layer.jl")
    include("test_land_surface.jl")
    include("test_ags.jl")
    include("test_cumulus.jl")
    include("test_switches.jl")
    include("test_quality.jl")
end
