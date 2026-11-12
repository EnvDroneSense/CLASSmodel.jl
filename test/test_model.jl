using CLASSmodel: esat, desat, qsat, init_model, step!

@testset "model" begin
    @testset "saturation" begin
        @test esat(273.16)≈611 rtol=1e-12
        @test esat(293.15)≈2339 rtol=0.01     # textbook value at 20 C
        @test qsat(293.15, 101300) ≈ 0.622 * esat(293.15) / 101300
        @test desat(285.0)≈(esat(285.0 + 1e-4) - esat(285.0 - 1e-4)) / 2e-4 rtol=1e-8
    end

    @testset "variable_info" begin
        @test keys(CLASSmodel.VARIABLE_INFO) == CLASSmodel.OUTPUT_NAMES    # every output, no extra
        @test variable_info(:wCO2).unit == "mg CO2 m-2 s-1"                # as in the docs
        @test_throws ArgumentError variable_info(:nonsense)
    end

    @testset "mixed_layer_profile" begin
        x = ModelInput(; runtime = 3600, gammatheta = 0.005, gammaq = -2e-6, dq = -0.002)
        out = run_model(x)
        p = mixed_layer_profile(x, out, 1; ztop = 1000)
        @test p.z == [0, 200, 200, 1000]
        @test p.theta ≈ [288, 288, 289, 289 + 0.005 * 800]
        @test p.q ≈ [0.008, 0.008, 0.006, 0.006 - 2e-6 * 800]

        i = length(out.h)                                    # a later row uses that row's state
        p = mixed_layer_profile(x, out, i)
        @test p.theta[3] - p.theta[2] ≈ out.dtheta[i] && p.z[2] == out.h[i]
        @test_throws ArgumentError mixed_layer_profile(x, out, 1; ztop = 100)
    end

    @testset "progress callback" begin
        x = ModelInput(; runtime = 3600, dt = 60)
        full = run_model(x)
        calls = Int[]
        out = run_model(x; progress_every = 7,
            progress = (i, n, o) -> (push!(calls, i); i < 10))
        @test calls == [7, 14] && isequal(out, map(v -> v[1:14], full))
        calls = Int[]
        run_model(x; progress_every = 25, progress = (i, n, o) -> push!(calls, i))
        @test calls == [25, 50, 60]
    end

    @testset "run loop" begin
        x = ModelInput(; runtime = 1000, dt = 60, tstart = 6, wCO2 = 0.2, sw_ml = false)
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
        y = ModelInput(; runtime = 600)
        s = init_model(y)
        for _ in 1:3
            step!(s)
        end
        @test s.t == 3 && s.mixed_layer.h ≈ run_model(y).h[4]
    end

    @testset "every output is written" begin
        # a column the run loop forgets to write would stay at zero
        out = run_model(ModelInput(;
            runtime = 3600, sw_rad = true, sw_sl = true, sw_ls = true,
            sw_cu = true, sw_chem = true, sw_wind = true, ls_type = :ags, wCO2 = 0.1, Wl = 1e-4,
            dC = ntuple(_ -> 1.0, 22), wC = ntuple(_ -> 1.0, 22)))
        for name in keys(out)
            name in (:wqM, :wCO2M, :M, :ac, :cc) && continue    # zero without clouds
            @test any(!=(0), out[name])
        end
    end
end
