using CLASSmodel: qsat

@testset "mixed layer" begin
    # without entrainment (beta = 0) and moisture flux, the heat flux warms the layer linearly
    x = ModelInput(; runtime=3600, dt=60, beta=0, wq=0, wtheta=0.1, h=500)
    out = run_model(x)
    @test all(==(500), out.h)
    @test out.theta ≈ x.theta .+ (0:59) .* x.dt .* x.wtheta ./ x.h
    @test all(==(x.q), out.q)

    # no fluxes and no entrainment: nothing changes
    out = run_model(ModelInput(; runtime=3600, wtheta=0, wq=0, wCO2=0))
    @test all(==(200), out.h) && all(==(288), out.theta) && all(==(422), out.CO2)

    # subsidence alone makes h decay like exp(-divU t)
    x = ModelInput(; runtime=3600, dt=1, wtheta=0, wq=0, wCO2=0, divU=1e-4)
    out = run_model(x)
    @test out.h[end] ≈ x.h * exp(-x.divU * (3600 - 1)) rtol = 1e-4

    # the lifting condensation level is where the mixed-layer air is saturated
    out = run_model(ModelInput(; runtime=600))
    x = ModelInput()
    T_lcl = out.theta[1] - 9.81 / 1005 * out.zlcl[1]
    @test out.q[1] / qsat(T_lcl, 101300 - 1.2 * 9.81 * out.zlcl[1]) ≈ 1 atol = 2e-4
end
