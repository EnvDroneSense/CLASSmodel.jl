@testset "cumulus" begin
    moist = ModelInput(; runtime = 6 * 3600, sw_rad = true, sw_sl = true,
        sw_ls = true, sw_cu = true, q = 0.01, dq = -0.002, wq = 1e-4)
    out = run_model(moist)

    # clouds form in the moist run
    @test maximum(out.ac) > 0 && maximum(out.M) > 0
    @test all(out.M .>= 0) && all(out.wqM .>= 0)
    @test all(out.dz .>= 50)

    # a dry mixed layer stays below saturation at its top: no clouds
    dry = run_model(ModelInput(; runtime = 3600, sw_rad = true, sw_sl = true, sw_ls = true,
        sw_cu = true, q = 0.002, dq = -0.001, wq = 0))
    @test all(==(0), dry.ac) && all(==(0), dry.M)

    # the cloud mass flux slows the growth of the mixed layer
    off = run_model(ModelInput(;
        runtime = 6 * 3600, sw_rad = true, sw_sl = true, sw_ls = true,
        sw_cu = false, q = 0.01, dq = -0.002, wq = 1e-4))
    @test out.h[end] < off.h[end]

    # clouds carry CO2 out of the layer only when the jump at its top is negative
    out = run_model(ModelInput(; runtime = 6 * 3600, sw_cu = true))
    @test any(>(0), out.ac) && any(>(0), out.wCO2M) && all(>=(0), out.wCO2M)
    out = run_model(ModelInput(; runtime = 6 * 3600, sw_cu = true, dCO2 = 10))
    @test any(>(0), out.ac) && all(==(0), out.wCO2M)

    # sw_curad: the radiation uses the cloud cover 2 ac of the step before
    x = (;
        runtime = 6 * 3600, sw_rad = true, sw_cu = true, cc = 0.1, q = 0.010, dq = -0.002)
    off = run_model(ModelInput(; x...))
    on = run_model(ModelInput(; x..., sw_curad = true))
    @test all(==(0.1), off.cc)
    @test on.cc == 2 .* on.ac
    @test any(>(0), on.ac) && on.Swin != off.Swin
end
