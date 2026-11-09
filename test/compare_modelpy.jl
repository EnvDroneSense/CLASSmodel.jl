# Comparison with the modelpy reference output. Temporary: remove this file and its
# include in runtests.jl when the port is validated. The reference files are made by
# tools/make_reference.py and are not tracked.
using TOML

const REFERENCE = joinpath(@__DIR__, "reference", "cases.toml")
const IMPLEMENTED = ("ml", "ml_wind_shear_fixft", "rad", "sl", "ls_js", "ls_ags_c3", "ls_ags_c4")   # cases compared in full; grows with each step

# Every element within rtol of the reference, or within `atol_frac` of the largest value of
# the column (for values that pass through zero, where a difference of two large numbers
# loses relative precision); equal NaN patterns. (isapprox on vectors would compare norms.)
function elementwise_close(a, b, rtol, atol_frac)
    atol = atol_frac * maximum(abs, filter(!isnan, b); init=0.0)
    return all(zip(a, b)) do (x, y)
        isnan(x) || isnan(y) ? isnan(x) && isnan(y) : abs(x - y) <= max(rtol * max(abs(x), abs(y)), atol)
    end
end

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
                    # L: the Newton iteration stops at |dL| < 1 mm, so round-off can add or skip
                    # one last step (about 1e-8 m near neutral, where |L| is large)
                    rtol, atol_frac = n == :L ? (1e-10, 0.0) : (1e-12, 1e-13)
                    @test elementwise_close(getfield(out, n), getindex.(ref, j), rtol, atol_frac)
                end
            end
        end
        # the defaults are the "ml" case
        @test all(getfield(ModelInput(), Symbol(k)) == value(v) for (k, v) in cases["ml"])
    else
        @info "test/reference not found (run tools/make_reference.py); skipping the modelpy comparison"
    end
end
