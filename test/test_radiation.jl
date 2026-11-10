@testset "radiation" begin
    # lon = 0, doy = 81 (near equinox), equator: the sun is overhead at 12 h UTC
    x = ModelInput(; sw_ml = false, sw_rad = true, runtime = 24 * 3600,
        dt = 1800, tstart = 0, lat = 0, lon = 0,
        doy = 81, cc = 0, alpha = 0.25)
    out = run_model(x)
    @test argmax(out.Swin) == 25                      # t = 12 h
    @test out.Swin[25] ≈ 1368 * 0.8 atol = 15         # sinlea ≈ 1, Tr = 0.8
    @test out.Swin[1] < 1                             # night: floor on the solar elevation
    @test out.Swout ≈ x.alpha .* out.Swin
    @test out.Q ≈ out.Swin - out.Swout + out.Lwin - out.Lwout

    # clouds reduce the shortwave radiation by the factor 1 - 0.4 cc
    cloudy = run_model(ModelInput(; sw_ml = false, sw_rad = true, runtime = 3600,
        doy = 81, cc = 0.5, tstart = 12, lon = 0, lat = 0))
    clear = run_model(ModelInput(; sw_ml = false, sw_rad = true, runtime = 3600,
        doy = 81, cc = 0, tstart = 12, lon = 0, lat = 0))
    @test cloudy.Swin ≈ clear.Swin .* 0.8
end
