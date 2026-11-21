@testset "species and chemistry" begin
    run(; kw...) = run_model(ModelInput(; runtime = 3600, kw...))
    one_species(i, v) = map(j -> j == i ? v : zero(v), 1:22)

    @testset "transport" begin
        # the flux of Inert (default 1 ppb m s-1) fills the layer
        out = run(; wtheta = 0)
        @test out.Inert[2] - out.Inert[1] ≈ 60 * 1 / out.h[1]
        @test out.O3[end] == 30 && out.N2[end] == 0.8e9    # no flux, no jump, no change
        # a species with a jump mixes into the layer like CO2
        jump = run(; CO2 = 30, dCO2 = 2, dC = one_species(2, 2.0))
        @test jump.O3 ≈ jump.CO2
        # a Sine flux follows the sine, a Sun flux the sun: zero at night, peak at noon
        only(shape) = map(j -> j == 1 ? shape : 0.0, 1:22)
        sine = run(; wC = only(Sine(0.1)))
        @test sine.Inert[2] - sine.Inert[1] ≈ 60 * 0.1 * sin(pi * 60 / 43200) / sine.h[1]
        sun = run(; wC = only(Sun(0.1)), tstart = 0, runtime = 86400)
        noon = argmax(diff(sun.Inert))
        @test sun.Inert[2] == sun.Inert[1] && 0.45 < noon / 1440 < 0.55
    end

    @testset "reactions" begin
        chem(; kw...) = run(; sw_chem = true, wtheta = 0, kw...)
        # no reactions: only the transport and H2O, which is set from q
        a, b = chem(; reactions = Reaction[]), run(; wtheta = 0)
        @test all(getproperty(a, n) == getproperty(b, n)
        for n in CLASSmodel.SPECIES if n != :H2O)
        @test a.H2O[end] ≈ 0.010 * 1e9 * 28.97 / 18    # qcbl_ref
        @test all(isnan, a.phi) && all(isnan, a.k_r05)   # no r05 and r21 in the scheme
        # NO2 + sun -> NO + O3 alone, without transport: the implicit step is exact
        c = chem(; reactions = CLASSmodel.REACTIONS[5:5], sw_ml = false)
        K = c.k_r05[2]
        @test K > 0 && c.NO2[2] ≈ c.NO2[1] / (1 + 60K) &&
              c.NO[2] + c.NO2[2] ≈ c.NO[1] + c.NO2[1]
        # photolysis follows the sun when it is not fixed: zero at night, and positive by day
        d = run(; runtime = 86400, sw_chem = true, sw_photo_constant = false, tstart = 0,
            lat = 0, lon = 0)
        @test d.k_r05[1] == 0 && maximum(d.k_r05) > 0 && count(==(0), d.k_r05) > 400
        # RO2 + HO2 gives OH only with an OH yield in r19
        iso = (;
            wC = one_species(10, 0.5), C = map(i -> i in (10, 12, 13) ? 1.0 : 0.0, 1:22))
        r = CLASSmodel.REACTIONS[19]
        yield = copy(CLASSmodel.REACTIONS)
        yield[19] = Reaction(r.name, r.form, r.inputs, [:OH => 2.0, :Product => 1.0], r.k)
        @test chem(; iso..., reactions = yield).OH[end] > chem(; iso...).OH[end]
        # a scheme of one new reaction: O3 decays to a product at the rate A
        decay = [Reaction(:d1, :constant, [:O3 => 1], [:Product => 1.0], (1e-13, 0, 0, 0, 0, 0, 0))]
        e = chem(; reactions = decay, sw_ml = false)
        @test 0 < e.O3[end] < e.O3[1] && all(isnan, e.k_r05)
        @test_throws ArgumentError chem(; reactions = [Reaction(:bad, :constant,
            [:NO => 1, :O3 => 1, :OH => 1], [:Product => 1.0], r.k)])
    end
end

@testset "scheme with two species" begin
    rx = [Reaction(:a, :constant, [:A => 1], [:B => 1.0], (1e-12, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0))]
    x = ModelInput(; species = [:A, :B], reactions = rx, sw_chem = true, runtime = 600,
        C = [10.0, 0.0], dC = [0.0, 0.0], gammaC = [0.0, 0.0], advC = [0.0, 0.0], wC = [0.0, 0.0])
    out = run_model(x)
    @test keys(out)[end-1:end] == (:A, :B) && out.B[end] > 0
    @test isnan(out.phi[end]) && out.A[end] + out.B[end] ≈ 10
    @test_throws ArgumentError ModelInput(; species = [:A, :B], reactions = rx, sw_chem = true)   # default C has 22 values
    s = CLASSmodel.init_model(x)
    @inferred CLASSmodel.run_loop(s, CLASSmodel.empty_output(Float64, 3, x.species), 3, nothing, 100)
end
