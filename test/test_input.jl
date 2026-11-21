@testset "input" begin
    @testset "element type" begin
        y = ModelInput{Float32}(; h = 300)
        @test y isa ModelInput{Float32} && y.h === 300.0f0
        @test ModelInput{Float32}(; wtheta = Sine(0.1)).wtheta isa Sine{Float64}
        @test ModelInput{Float32}(; gammatheta = Segments([1.0], [0.0, 1.0])).gammatheta isa Segments{Float32}
        @test AGsParameters(C3()) isa AGsParameters{Float64}
        @test AGsParameters{Float32}(C4()) isa AGsParameters{Float32}
    end

    @testset "input_info" begin
        x = ModelInput()
        @test keys(CLASSmodel.INPUT_INFO) == fieldnames(ModelInput)   # every field, no extra
        for name in fieldnames(ModelInput)
            info, v = input_info(name), getfield(x, name)
            if v isa Tuple    # one value for each species
                @test hasproperty(info, :options) ? all(in(info.options), v) :
                      all(e -> info.min <= e <= info.max, v)
            elseif v isa Symbol
                @test v in info.options
            elseif v isa Vector    # the reactions: a title only
                @test info.title isa String
            elseif v isa Bool
                @test info.min === nothing && info.max === nothing
            else
                @test info.min <= v <= info.max
            end
        end
        @test_throws ArgumentError input_info(:nonsense)
    end

    @testset "validation" begin
        invalid = (
            (; dt = 0), (; runtime = 30, dt = 60), (; h = -1),
            (; cc = 1.5), (; cveg = -0.1), (; c_beta = 2),
            (; ls_type = :other), (; c3c4 = :c5),
            (; sw_ls = true, LAI = 0), (; sw_ls = true, wfc = 0.1, wwilt = 0.2),
            (; sw_ls = true, ls_type = :ags),   # A-Gs needs radiation
            (; dtheta = 0), (; sw_ml = false, sw_cu = true, dtheta = 0),   # division by dtheta
            (; sw_ls = true, w2 = 0.472), (; sw_cu = true, dz_h = 0),
            (; theta = NaN), (; wtheta = Sine(NaN)),
            (; dtheta = 0, beta = 0, dq = 0, gammatheta = 0),   # no inversion needs a lapse rate
            (; sw_chem = true, reactions = [Reaction(:x, :nonsense, [:O3 => 1], [:Product => 1.0], ntuple(_ -> 1.0, 7))]),
            (; C = ntuple(_ -> -1.0, 22)), (; wC = ntuple(_ -> NaN, 22)),
            (; tod_ref = 25), (; P_ref = 0), (; qcbl_ref = -0.1)
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
