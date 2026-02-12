using HomotopyContinuation, Test

@testset "Code quality" begin
    using Aqua
    Aqua.test_all(HomotopyContinuation)
end
