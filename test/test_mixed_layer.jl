using CLASSmodel: qsat

@testset "mixed layer" begin
    # without entrainment (beta = 0) and moisture flux, the heat flux warms the layer linearly
    x = ModelInput(; runtime = 3600, dt = 60, beta = 0, wq = 0, wtheta = 0.1, h = 500)
    out = run_model(x)
    @test all(==(500), out.h)
    @test out.theta ≈ x.theta .+ (0:59) .* x.dt .* x.wtheta ./ x.h
    @test all(==(x.q), out.q)

    # no fluxes and no entrainment: nothing changes
    out = run_model(ModelInput(; runtime = 3600, wtheta = 0, wq = 0, wCO2 = 0))
    @test all(==(200), out.h) && all(==(288), out.theta) && all(==(422), out.CO2)

    # subsidence alone makes h decay like exp(-divU t)
    x = ModelInput(; runtime = 3600, dt = 1, wtheta = 0, wq = 0, wCO2 = 0, divU = 1e-4)
    out = run_model(x)
    @test out.h[end]≈x.h * exp(-x.divU * (3600 - 1)) rtol=1e-4

    # without surface layer the stress has the size ustar^2 and the direction of the wind, also for u = 0
    out = run_model(ModelInput(; runtime = 60, ustar = 0.3, u = 3, v = -4))
    @test hypot(out.uw[1], out.vw[1]) ≈ 0.3^2 && out.uw[1] ≈ -0.09 * 3 / 5 &&
          out.vw[1] ≈ 0.09 * 4 / 5
    out = run_model(ModelInput(; runtime = 60, ustar = 0.3, u = 0, v = -4))
    @test out.uw[1] == 0 && out.vw[1] ≈ 0.09

    # the lifting condensation level is where the mixed-layer air is saturated
    out = run_model(ModelInput(; runtime = 600))
    x = ModelInput()
    T_lcl = out.theta[1] - 9.81 / 1005 * out.zlcl[1]
    @test out.q[1] / qsat(T_lcl, 101300 - 1.2 * 9.81 * out.zlcl[1])≈1 atol=2e-4
    # a layer that stays far below saturation does not converge, and says so
    x = ModelInput(; runtime = 60, h = 3000, q = 0.0001, theta = 300)
    @test_logs (:warn, "LCL calculation not converged") match_mode=:any run_model(x)
end

@testset "mixed-layer switches" begin
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

    @testset "lapse rate in segments" begin
        # one segment, or equal segments, is the constant lapse rate, bit for bit
        @test run(; gammatheta = Segments([1000.0], [0.006, 0.006])).h == run().h
        # the layer grows faster once it is above 250 m, where the air is less stable
        out = run(; gammatheta = Segments([250.0], [0.006, 0.002]))
        i = findfirst(>(250), out.h)
        @test out.h[1:i] == run().h[1:i] && out.h[end] > run().h[end]
        @test_throws ArgumentError Segments([2.0, 1.0], [1.0, 2.0, 3.0])
        p = mixed_layer_profile(ModelInput(; gammatheta = Segments([300.0], [0.01, 0.0])), run(), 1; ztop = 400)
        @test p.theta[end] ≈ 288 + 1 + 0.01 * 100
    end

    @testset "encroachment" begin
        # dry air without an inversion: h^2 = h0^2 + 2 wtheta t / gammatheta
        x = ModelInput(; runtime = 3600, dt = 1, beta = 0, dtheta = 0, dq = 0, wq = 0)
        out = run_model(x)
        @test out.h[end] ≈ sqrt(200^2 + 2 * 0.1 * 3599 / 0.006) rtol = 1e-3
        @test maximum(abs, out.dtheta) < 1e-9
        @test_throws ArgumentError ModelInput(; beta = 0, dtheta = 0, dq = 0, sw_cu = true)
    end

    @testset "diurnal cycle of the fluxes" begin
        # half sine, taken at the end of each step: the amplitude at period / 2, zero after
        x = (; runtime = 2500, wtheta = Sine(0.1; period = 1200), wq = Sine(1e-4; period = 1200))
        out = run(; x...)
        @test out.wtheta[1] ≈ 0.1 * sin(pi * 60 / 1200)
        @test out.wtheta[10] ≈ 0.1 && out.wq[10] ≈ 1e-4        # end of the step at 600 s
        @test out.wtheta[31] == 0
        @test run(; x..., wq = 1e-4).wq == fill(1e-4, 41)      # each flux has its own shape
        # the land surface sets the fluxes itself
        ls = (; sw_sl = true, sw_ls = true)
        @test_throws ArgumentError run(; ls..., x...)
    end
end
