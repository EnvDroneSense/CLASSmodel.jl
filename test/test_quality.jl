using Aqua
using JET

@testset "package quality" begin
    Aqua.test_all(CLASSmodel)
    # JET follows the compiler closely, so it only runs on released Julia versions.
    if isempty(VERSION.prerelease)
        JET.test_package(CLASSmodel; target_modules = (CLASSmodel,))
    end
end
