@testset "cumulus" begin
    moist = ModelInput(; runtime = 6 * 3600, sw_rad = true, sw_sl = true,
        sw_ls = true, sw_cu = true, q = 0.01, dq = -0.002, wq = 1e-4)
    out = run_model(moist)

    # the closure bounds the cloud core fraction by 0.5 + 0.18 pi (slightly above 1), clouds form in the moist run
    @test all(0 .<= out.ac .<= 0.5 + 0.18π)
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
end
