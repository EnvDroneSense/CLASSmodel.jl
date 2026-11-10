using CLASSmodel: E1

@testset "A-Gs" begin
    # exponential integral at known values
    @test E1(1.0) ≈ 0.21938393439552029
    @test E1(0.1) ≈ 1.8229239584193906
    @test E1(2.0) < E1(1.0)
    # series and continued fraction agree, also for large arguments and for Float32
    @test E1(nextfloat(1.0))≈E1(1.0) rtol=1e-14
    @test E1(5.0)≈1.1482955912753257e-3 rtol=1e-14
    @test E1(10.0)≈4.156968929685325e-6 rtol=1e-14
    @test E1(40.0)≈1.036773261451657e-19 rtol=1e-14
    @test E1(3.0f0) isa Float32
    @test E1(3.0f0)≈E1(3.0) rtol=1e-6

    ags(; kw...) = run_model(ModelInput(;
        runtime = 600, sw_rad = true, sw_sl = true, sw_ls = true, ls_type = :ags, kw...))

    # by day the canopy takes up CO2 and the soil releases it
    out = ags()
    @test all(out.wCO2A .< 0) && all(out.wCO2R .> 0)
    @test out.wCO2 ≈ out.wCO2A + out.wCO2R
    @test out.wCO2A ≈ ags(c3c4 = :c3).wCO2A && out.wCO2A != ags(c3c4 = :c4).wCO2A

    # drier soil, stronger water stress, higher surface resistance
    @test ags(w2 = 0.2).rs[1] > ags(w2 = 0.3).rs[1]
    @test ags(w2 = 0.2, c_beta = 0.5).rs[1] != ags(w2 = 0.2).rs[1]
end
