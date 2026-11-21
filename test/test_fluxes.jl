using CLASSmodel: flux

@testset "flux shapes" begin
    x = ModelInput()
    at(f, t) = flux(f, x, t, 1.0)
    # every shape is zero outside its window and has its amplitude inside
    @test at(Sine(2; period = 100, start = 50), 100) ≈ 2 && at(Sine(2; period = 100), 101) == 0
    @test at(Cosine(2; period = 100), 50) ≈ 2 && at(Cosine(2; period = 100), 0) == 0
    @test at(Block(2; period = 100, start = 50), 60) == 2 && at(Block(2; start = 50), 49) == 0
    # a table is linear between its values and constant after the last one
    @test at(Table([0, 2, 4]), 1800) ≈ 1 && at(Table([0, 2, 4]), 5400) ≈ 3
    @test at(Table([0, 2, 4]), 1e5) == 4
    @test at(Sun(2.0), 5.5 * 3600) > at(Sun(2.0), 2 * 3600) > 0 == at(Sun(2.0), 17 * 3600)
    @test at(0.3, 10) == 0.3 && flux(Deposition(0.01), x, 0, 5.0) == -0.05

    # a constant table is the constant flux, bit for bit
    run(; kw...) = run_model(ModelInput(; runtime = 7200, kw...))
    @test run(; wtheta = Table([0.1, 0.1, 0.1])).h == run().h
    # dry deposition in a layer of fixed height: C decays as exp(-vd t / h)
    dep = map(j -> j == 2 ? Deposition(0.01) : 0.0, 1:22)
    out = run_model(ModelInput(; runtime = 3600, dt = 1, wtheta = 0, wq = 0, wC = dep))
    @test out.O3[end] ≈ 30 * exp(-0.01 * 3599 / 200) rtol = 1e-4
    @test all(==(200), out.h)
    # the same for CO2, whose flux is in ppm m s-1 inside the layer
    x0 = ModelInput()
    co2 = run_model(ModelInput(; runtime = 3600, dt = 1, wtheta = 0, wq = 0, wCO2 = Deposition(0.01)))
    @test co2.CO2[end] ≈ x0.CO2 * exp(-0.01 * 3599 / 200) rtol = 1e-4
end
