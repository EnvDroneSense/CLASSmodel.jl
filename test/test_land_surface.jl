@testset "land surface" begin
    x = (; runtime = 3600, sw_rad = true, sw_sl = true, sw_ls = true)
    out = run_model(ModelInput(; x...))

    # the skin temperature closes the energy balance
    @test out.Q≈out.H + out.LE + out.G rtol=1e-9
    @test out.LE ≈ out.LEveg + out.LEsoil + out.LEliq
    @test all(==(0), out.LEliq)                      # dry canopy: no open-water evaporation
    # the fluxes are the ones the mixed layer sees
    @test out.wtheta ≈ out.H ./ (1.2 * 1005)
    @test out.wq ≈ out.LE ./ (1.2 * 2.5e6)

    # more leaves give a lower surface resistance; a soil at wilting point blocks transpiration
    rs(; kw...) = run_model(ModelInput(; x..., runtime = 600, kw...)).rs[1]
    @test rs(LAI = 4) < rs(LAI = 2)
    @test rs(w2 = 0.1) > 1e6

    # sw_sea: open water at the fixed Ts, without the soil and its input checks
    out = run_model(ModelInput(; x..., sw_sea = true, w2 = 0.6))   # w2 > wsat
    @test all(==(0), out.rs) && all(==(0), out.LEsoil) && all(==(0), out.LEveg)
    @test out.LE == out.LEliq == out.LEpot && all(isnan, out.G) && all(isnan, out.LEref)
    @test out.H[1] > 0 && out.LE[1] > 0
    @test_throws ArgumentError ModelInput(; x..., w2 = 0.6)
end
