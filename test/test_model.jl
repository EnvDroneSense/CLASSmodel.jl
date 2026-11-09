using CLASSmodel: esat, qsat, init_model

@testset "model" begin
    @testset "saturation" begin
        @test esat(273.16) ≈ 611 rtol = 1e-12
        @test esat(293.15) ≈ 2339 rtol = 0.01     # textbook value at 20 C
        @test qsat(293.15, 101300) ≈ 0.622 * esat(293.15) / 101300
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
        @test all(isnan, out.Swin) && all(==(0), out.z0m)          # not computed, not written
        @test run_model(ModelInput{Float32}()) isa ModelOutput{Float32}
    end
end
