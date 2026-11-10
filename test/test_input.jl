@testset "input" begin
    @testset "element type" begin
        @test ModelInput() isa ModelInput{Float64}
        @test ModelInput(; h = 250).h === 250.0
        y = ModelInput{Float32}(; h = 300)
        @test y isa ModelInput{Float32} && y.h === 300.0f0
        @test PhysicalConstants{Float32}().g === 9.81f0
        @test AGsParameters(C3()) isa AGsParameters{Float64}
        @test AGsParameters{Float32}(C4()) isa AGsParameters{Float32}
    end

    @testset "validation" begin
        invalid = (
            (; dt = 0), (; runtime = 30, dt = 60), (; h = -1), (; Ps = 0), (; z0m = 0), (;
                z0h = -0.1),
            (; cc = 1.5), (; cveg = -0.1), (; c_beta = 2),
            (; ls_type = :other), (; c3c4 = :c5),
            (; sw_ls = true, LAI = 0), (; sw_ls = true, wfc = 0.1, wwilt = 0.2),
            (; sw_ls = true, ls_type = :ags),   # A-Gs needs radiation
            (; dtheta = 0), (; sw_ml = false, sw_cu = true, dtheta = 0),   # division by dtheta
            (; sw_ls = true, Wmax = 0), (; sw_ls = true, wg = 0), (; sw_ls = true, w2 = 0),
            (; sw_ls = true, w2 = 0.472), (; sw_cu = true, dz_h = 0),
            (; theta = NaN), (; fc = Inf), (; dq = NaN)
        )
        for kw in invalid
            @test_throws ArgumentError ModelInput(; kw...)
        end
        @test ModelInput(; sw_ls = true, ls_type = :ags, sw_rad = true) isa ModelInput
        @test ModelInput(; LAI = 0) isa ModelInput   # land surface checks need sw_ls
        @test ModelInput(; dtheta = 0, sw_ml = false) isa ModelInput
        @test ModelInput(; Wmax = 0, dz_h = 0) isa ModelInput
    end
end
