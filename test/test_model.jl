using CLASSmodel: esat, desat, qsat, init_model, step!

@testset "model" begin
    @testset "saturation" begin
        @test esat(273.16) ≈ 611 rtol = 1e-12
        @test esat(293.15) ≈ 2339 rtol = 0.01     # textbook value at 20 C
        @test qsat(293.15, 101300) ≈ 0.622 * esat(293.15) / 101300
        @test desat(285.0) ≈ (esat(285.0 + 1e-4) - esat(285.0 - 1e-4)) / 2e-4 rtol = 1e-8
    end

    @testset "run loop" begin
        x = ModelInput(; runtime=1000, dt=60, tstart=6, wCO2=0.2, sw_ml=false)
        s = init_model(x)
        @test s.t == 0
        @test s.mixed_layer.wCO2 ≈ 0.2 * 28.9 / (1.2 * 44)

        out = run_model(x)
        @test length(out.t) == 16                                  # floor(1000 / 60)
        @test out.t ≈ 6 .+ (0:15) .* 60 / 3600
        @test all(==(x.h), out.h) && all(==(x.u), out.u)           # mixed layer off: no change
        @test out.wCO2 ≈ fill(0.2, 16)                             # unit conversion round trip
        @test all(isnan, out.Swin)                                 # not computed
        @test all(==(x.z0m), out.z0m) && all(==(x.z0h), out.z0h) && all(==(0), out.wCO2M)
        @test run_model(ModelInput{Float32}()) isa ModelOutput{Float32}

        # stepping by hand gives the run's state: the row of step n holds the state before it
        y = ModelInput(; runtime=600)
        s = init_model(y)
        for _ in 1:3
            step!(s)
        end
        @test s.t == 3 && s.mixed_layer.h ≈ run_model(y).h[4]
    end
end
