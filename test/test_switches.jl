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

    @testset "diurnal cycle of the fluxes" begin
        # half sine, taken at the end of each step: the input value at sinperiod / 2, negative after
        x = (; runtime = 2500, sinperiod = 1200, sw_wtheta = true, sw_wq = true)
        out = run(; x...)
        @test out.wtheta[1] ≈ 0.1 * sin(pi * 60 / 1200)
        @test out.wtheta[10] ≈ 0.1 && out.wq[10] ≈ 1e-4        # end of the step at 600 s
        @test out.wtheta[31] < 0
        @test run(; x..., sw_wq = false).wq == fill(1e-4, 41)  # each switch acts alone
        # the land surface sets the fluxes itself
        ls = (; sw_sl = true, sw_ls = true)
        @test isequal(run(; ls..., x...).wtheta, run(; ls..., runtime = x.runtime).wtheta)
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
        @test_logs (:warn, "LCL calculation not converged") match_mode=:any run_model(x)
    end

    @testset "cumulus cloud cover in radiation" begin
        x = (; runtime = 6 * 3600, sw_rad = true, sw_cu = true, cc = 0.1, q = 0.010, dq = -0.002)
        off = run_model(ModelInput(; x...))
        on = run_model(ModelInput(; x..., sw_curad = true))
        @test all(==(0.1), off.cc)
        @test on.cc == 2 .* on.ac    # radiation reads the value of the step before
        @test any(>(0), on.ac) && on.Swin != off.Swin
    end

    @testset "sea" begin
        x = (; runtime = 3600, sw_rad = true, sw_sl = true, sw_ls = true)
        out = run_model(ModelInput(; x..., sw_sea = true, w2 = 0.6))   # w2 > wsat is not checked
        @test all(==(0), out.rs) && all(==(0), out.LEsoil) && all(==(0), out.LEveg)
        @test out.LE == out.LEliq == out.LEpot && all(isnan, out.G) && all(isnan, out.LEref)
        @test out.H[1] > 0 && out.LE[1] > 0
        @test_throws ArgumentError ModelInput(; x..., w2 = 0.6)
    end

    @testset "passive scalar" begin
        # the surface flux fills the layer: dS/dt = wS / h at the first step
        out = run(; wS = 0.01, dS = 0, wtheta = 0)
        @test out.S[2] - out.S[1] ≈ 60 * 0.01 / out.h[1]
        # a jump with entrainment mixes into the layer, like CO2
        out = run(; S = 1, dS = 2, dCO2 = 2, CO2 = 1)
        @test out.S ≈ out.CO2 && out.dS ≈ out.dCO2
    end

    @testset "chemical species" begin
        # passive for now: the flux of Inert (default 1 ppb m s-1) fills the layer
        out = run(; wtheta = 0)
        @test out.Inert[2] - out.Inert[1] ≈ 60 * 1 / out.h[1]
        @test out.O3[end] == 30 && out.N2[end] == 0.8e9    # no flux, no jump, no change
        # the sine mode of a species, and a species with a jump behaves like the scalar S
        modes = ntuple(i -> i == 1 ? 1 : 0, 22)
        sine = run(; sw_wC = modes, wC = ntuple(_ -> 0.1, 22))
        @test sine.Inert[2] - sine.Inert[1] ≈ 60 * 0.1 * sin(pi * 60 / 43200) / sine.h[1]
        jump = run(; S = 30, dS = 2, dC = ntuple(i -> i == 2 ? 2.0 : 0.0, 22))
        @test jump.O3 ≈ jump.S
    end

    @testset "chemistry" begin
        chem(; kw...) = run(; sw_chem = true, wtheta = 0, kw...)
        # no reactions: only the transport and H2O, which is set from q
        none = ntuple(_ -> false, 28)
        a, b = chem(; sw_reactions = none), run(; wtheta = 0)
        @test all(getproperty(a, n) == getproperty(b, n) for n in CLASSmodel.SPECIES if n != :H2O)
        @test a.H2O[end] ≈ 0.010 * 1e9 * 28.97 / 18    # qcbl_ref
        @test all(==(0), a.phi) && all(==(0), a.k_r05)
        # NO2 + sun -> NO + O3 alone, without transport: the implicit step is exact
        only = ntuple(i -> i == 5, 28)
        c = chem(; sw_reactions = only, sw_ml = false)
        K = c.k_r05[2]
        @test K > 0 && c.NO2[2] ≈ c.NO2[1] / (1 + 60K) && c.NO[2] + c.NO2[2] ≈ c.NO[1] + c.NO2[1]
        # photolysis follows the sun when it is not fixed: zero at night, and positive by day
        d = run_model(ModelInput(; runtime = 86400, sw_chem = true, sw_photo_constant = false,
            tstart = 0, lat = 0, lon = 0))
        @test d.k_r05[1] == 0 && maximum(d.k_r05) > 0 && count(==(0), d.k_r05) > 400
        # RO2 + HO2 gives OH only with stocoef
        iso = (; wC = ntuple(i -> i == 10 ? 0.5 : 0.0, 22), C = ntuple(i -> i in (10, 12, 13) ? 1.0 : 0.0, 22))
        @test chem(; iso..., stocoef = 2).OH[end] > chem(; iso..., stocoef = 0).OH[end]
    end

    @testset "every output is written" begin
        # a column the run loop forgets to write would stay at zero
        out = run_model(ModelInput(;
            runtime = 3600, sw_rad = true, sw_sl = true, sw_ls = true,
            sw_cu = true, sw_chem = true, sw_wind = true, ls_type = :ags, wCO2 = 0.1, Wl = 1e-4,
            S = 1, dS = 1, wS = 0.01,
            dC = ntuple(_ -> 1.0, 22), wC = ntuple(_ -> 1.0, 22)))
        for name in keys(out)
            name in (:wqM, :wCO2M, :wSM, :M, :ac, :cc) && continue    # zero without clouds
            @test any(!=(0), out[name])
        end
    end
end
