# One test per model switch that the other test files do not cover.

@testset "switches" begin
    run(; kw...) = run_model(ModelInput(; runtime = 3600, kw...))

    @testset "sw_wind" begin
        # the wind above the layer is geostrophic: u + du and v + dv stay constant
        x = (; sw_wind = true, ustar = 0)
        out = run(; x...)
        @test all(≈(10), out.u .+ out.du) && all(≈(0; atol = 1e-12), out.v .+ out.dv)
        @test out.v[end] > out.v[1]                    # the Coriolis force turns the wind
        @test run(; x..., advu = 1e-5).u[end] > out.u[end]
        # switched off: the wind does not change
        @test all(==(6), run().u)
    end

    @testset "sw_shearwe" begin
        # shear adds to the entrainment; without friction velocity it does not
        @test run(; sw_shearwe = true).h[end] > run().h[end]
        @test run(; sw_shearwe = true, ustar = 0).h == run(; ustar = 0).h
    end

    @testset "sw_fixft" begin
        # with subsidence the fixed free troposphere changes the jump tendency by gamma * ws
        x = (; divU = 1e-5, gammatheta = 0.006)
        out, ref = run(; x..., sw_fixft = true), run(; x...)
        @test out.dtheta[2] - ref.dtheta[2] ≈ x.gammatheta * (-x.divU * 200) * 60
        @test run(; sw_fixft = true).dtheta == run().dtheta   # no subsidence, no difference
    end

    @testset "A-Gs C4 and water stress" begin
        ags(; kw...) = run(;
            sw_rad = true, sw_sl = true, sw_ls = true, ls_type = :ags, kw...)
        c4 = ags(c3c4 = :c4)
        @test all(c4.wCO2A .< 0) && c4.wCO2A != ags().wCO2A
        # curvature of the water stress: the stronger the curvature, the lower the resistance
        rs(c_beta) = ags(; w2 = 0.2, c_beta).rs[1]
        @test rs(0) > rs(0.5) > rs(1)
    end

    @testset "cumulus CO2 flux" begin
        # clouds carry CO2 out of the layer only when the jump at its top is negative
        out = run_model(ModelInput(; runtime = 6 * 3600, sw_cu = true))
        @test any(>(0), out.ac) && any(>(0), out.wCO2M) && all(>=(0), out.wCO2M)
        out = run_model(ModelInput(; runtime = 6 * 3600, sw_cu = true, dCO2 = 10))
        @test any(>(0), out.ac) && all(==(0), out.wCO2M)
    end

    @testset "lifting condensation level" begin
        # a layer that stays far below saturation does not converge, and says so
        x = ModelInput(; runtime = 60, h = 3000, q = 0.0001, theta = 300)
        @test_logs (:warn, "LCL calculation not converged") match_mode = :any run_model(x)
    end

    @testset "every output is written" begin
        # a column the run loop forgets to write would stay at zero
        out = run_model(ModelInput(;
            runtime = 3600, sw_rad = true, sw_sl = true, sw_ls = true,
            sw_cu = true, sw_wind = true, ls_type = :ags, wCO2 = 0.1, Wl = 1e-4))
        for name in keys(out)
            name in (:wqM, :wCO2M, :M, :ac) && continue    # zero without clouds
            @test any(!=(0), out[name])
        end
    end
end
