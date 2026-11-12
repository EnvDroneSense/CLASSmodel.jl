@testset "species and chemistry" begin
    run(; kw...) = run_model(ModelInput(; runtime = 3600, kw...))
    one_species(i, v) = ntuple(j -> j == i ? v : zero(v), 22)

    @testset "transport" begin
        # the flux of Inert (default 1 ppb m s-1) fills the layer
        out = run(; wtheta = 0)
        @test out.Inert[2] - out.Inert[1] ≈ 60 * 1 / out.h[1]
        @test out.O3[end] == 30 && out.N2[end] == 0.8e9    # no flux, no jump, no change
        # a species with a jump mixes into the layer like CO2
        jump = run(; CO2 = 30, dCO2 = 2, dC = one_species(2, 2.0))
        @test jump.O3 ≈ jump.CO2
        # flux mode 1 follows the sine, mode 2 the sun: zero at night, peak at noon
        sine = run(; sw_wC = one_species(1, 1), wC = ntuple(_ -> 0.1, 22))
        @test sine.Inert[2] - sine.Inert[1] ≈ 60 * 0.1 * sin(pi * 60 / 43200) / sine.h[1]
        sun = run(; sw_wC = one_species(1, 2), wC = ntuple(_ -> 0.1, 22), tstart = 0,
            runtime = 86400)
        noon = argmax(diff(sun.Inert))
        @test sun.Inert[2] == sun.Inert[1] && 0.45 < noon / 1440 < 0.55
    end

    @testset "reactions" begin
        chem(; kw...) = run(; sw_chem = true, wtheta = 0, kw...)
        # no reactions: only the transport and H2O, which is set from q
        a, b = chem(; sw_reactions = ntuple(_ -> false, 28)), run(; wtheta = 0)
        @test all(getproperty(a, n) == getproperty(b, n)
        for n in CLASSmodel.SPECIES if n != :H2O)
        @test a.H2O[end] ≈ 0.010 * 1e9 * 28.97 / 18    # qcbl_ref
        @test all(==(0), a.phi) && all(==(0), a.k_r05)
        # NO2 + sun -> NO + O3 alone, without transport: the implicit step is exact
        c = chem(; sw_reactions = ntuple(==(5), 28), sw_ml = false)
        K = c.k_r05[2]
        @test K > 0 && c.NO2[2] ≈ c.NO2[1] / (1 + 60K) &&
              c.NO[2] + c.NO2[2] ≈ c.NO[1] + c.NO2[1]
        # photolysis follows the sun when it is not fixed: zero at night, and positive by day
        d = run(; runtime = 86400, sw_chem = true, sw_photo_constant = false, tstart = 0,
            lat = 0, lon = 0)
        @test d.k_r05[1] == 0 && maximum(d.k_r05) > 0 && count(==(0), d.k_r05) > 400
        # RO2 + HO2 gives OH only with stocoef
        iso = (;
            wC = one_species(10, 0.5), C = ntuple(i -> i in (10, 12, 13) ? 1.0 : 0.0, 22))
        @test chem(; iso..., stocoef = 2).OH[end] > chem(; iso..., stocoef = 0).OH[end]
    end
end
