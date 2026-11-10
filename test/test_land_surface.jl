@testset "land surface" begin
    x = ModelInput(; runtime = 3600, sw_rad = true, sw_sl = true, sw_ls = true)
    out = run_model(x)

    # the skin temperature closes the energy balance
    @test out.Q≈out.H + out.LE + out.G rtol=1e-9
    @test out.LE ≈ out.LEveg + out.LEsoil + out.LEliq
    @test all(==(0), out.LEliq)                      # dry canopy: no open-water evaporation
    # the fluxes are the ones the mixed layer sees
    @test out.wtheta ≈ out.H ./ (1.2 * 1005)
    @test out.wq ≈ out.LE ./ (1.2 * 2.5e6)

    # more leaves give a lower surface resistance; a soil at wilting point blocks transpiration
    rs(; kw...) = run_model(ModelInput(;
    runtime = 600, sw_rad = true, sw_sl = true, sw_ls = true, kw...)).rs[1]
    @test rs(LAI = 4) < rs(LAI = 2)
    @test rs(w2 = 0.1) > 1e6
end
