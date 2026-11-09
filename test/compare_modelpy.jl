# Comparison with the modelpy reference output. Temporary: remove this file and its
# include in runtests.jl when the port is validated. The reference files are made by
# tools/make_reference.py and are not tracked.
using TOML

const REFERENCE = joinpath(@__DIR__, "reference", "cases.toml")
const IMPLEMENTED = ("ml", "ml_wind_shear_fixft", "rad")   # cases compared in full; grows with each step

@testset "modelpy reference" begin
    if isfile(REFERENCE)
        cases = TOML.parsefile(REFERENCE)["cases"]
        value(v) = v isa String ? Symbol(v) : v
        @testset "input and time axis: $name" for (name, case) in cases
            x = ModelInput(; (Symbol(k) => value(v) for (k, v) in case)...)
            @test all(getfield(x, Symbol(k)) == value(v) for (k, v) in case)

            lines = filter(l -> !startswith(l, "#"), readlines(joinpath(@__DIR__, "reference", "$name.csv")))
            names = Symbol.(split(lines[1], ','))
            ref = [parse.(Float64, split(l, ',')) for l in lines[2:end]]
            out = run_model(x)
            @test out.t == getindex.(ref, 1)

            # all columns, for the cases whose physics is implemented
            if name in IMPLEMENTED
                for (j, n) in enumerate(names)
                    @test isapprox(getfield(out, n), getindex.(ref, j); rtol=1e-12, nans=true)
                end
            end
        end
        # the defaults are the "ml" case
        @test all(getfield(ModelInput(), Symbol(k)) == value(v) for (k, v) in cases["ml"])
    else
        @info "test/reference not found (run tools/make_reference.py); skipping the modelpy comparison"
    end
end
