using CLASSmodel
using Test

@testset "CLASSmodel.jl" begin
    include("test_input.jl")
    include("test_model.jl")
    include("compare_modelpy.jl")   # temporary, see the file
end
