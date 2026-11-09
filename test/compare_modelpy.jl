# Comparison with the modelpy reference output. Temporary: remove this file and its
# include in runtests.jl when the port is validated. The reference files are made by
# tools/make_reference.py and are not tracked.
using TOML

const REFERENCE = joinpath(@__DIR__, "reference", "cases.toml")

@testset "modelpy reference" begin
    if isfile(REFERENCE)
        cases = TOML.parsefile(REFERENCE)["cases"]
        value(v) = v isa String ? Symbol(v) : v
        @testset "input: $name" for (name, case) in cases
            x = ModelInput(; (Symbol(k) => value(v) for (k, v) in case)...)
            @test all(getfield(x, Symbol(k)) == value(v) for (k, v) in case)
        end
        # the defaults are the "ml" case
        @test all(getfield(ModelInput(), Symbol(k)) == value(v) for (k, v) in cases["ml"])
    else
        @info "test/reference not found (run tools/make_reference.py); skipping the modelpy comparison"
    end
end
