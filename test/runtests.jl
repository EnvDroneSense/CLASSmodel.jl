using CLASSmodel
using Test

@testset "CLASSmodel.jl" begin
    include("test_input.jl")
    include("test_model.jl")
    include("test_mixed_layer.jl")
    include("test_radiation.jl")
    include("compare_modelpy.jl")   # temporary, see the file
end
