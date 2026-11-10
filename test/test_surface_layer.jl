using CLASSmodel: psim, psih, ribtol, ribfun

@testset "surface layer" begin
    # no stability correction at the surface, continuous across neutral
    @test abs(psim(0.0)) < 1e-7 && abs(psih(0.0)) < 1e-7
    @test abs(psim(-1e-8) - psim(1e-8)) < 1e-6 && abs(psih(-1e-8) - psih(1e-8)) < 1e-6
    # stable air mixes less (psi < 0), unstable air more (psi > 0)
    @test psim(1.0) < 0 && psih(1.0) < 0
    @test psim(-1.0) > 0 && psih(-1.0) > 0

    # ribtol inverts the bulk Richardson number, with the sign of L following Rib
    zsl, z0m, z0h = 20.0, 0.02, 0.002
    for Rib in (-1.0, -0.1, 0.05, 0.15)
        L = ribtol(Rib, zsl, z0m, z0h)
        @test sign(L) == sign(Rib)
        @test ribfun(L, zsl, z0m, z0h) ≈ Rib atol = 1e-4
    end
    @test abs(ribtol(0.0, zsl, z0m, z0h)) > 1e3     # neutral: very large |L|

    # surface layer in a run: fluxes follow from the drag coefficient, warm surface is unstable
    out = run_model(ModelInput(; runtime = 600, sw_sl = true, wtheta = 0.1, ustar = 0.1))
    @test out.L[1] < 0 && out.Rib[1] < 0
    @test out.ustar[1] ≈ sqrt(out.Cm[1]) * hypot(out.u[1], out.v[1], 0.0) atol = 0.5 *
                                                                                 out.ustar[1]
    @test out.uw[1] ≈ -out.Cm[1] * hypot(out.u[1], out.v[1]) * out.u[1] rtol = 0.2
end
