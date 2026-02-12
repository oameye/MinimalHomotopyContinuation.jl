using MinimalHomotopyContinuation, Test

@testset "Code quality" begin
    using Aqua
    Aqua.test_all(MinimalHomotopyContinuation)
end
