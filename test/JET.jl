using MinimalHomotopyContinuation, Test

@testset "Code linting" begin
    using JET
    JET.test_package(
        MinimalHomotopyContinuation;
        target_modules = (
            MinimalHomotopyContinuation,
            MinimalHomotopyContinuation.ModelKit,
            MinimalHomotopyContinuation.DoubleDouble,
        ),
    )

    rep = JET.report_package(MinimalHomotopyContinuation)
    reports = JET.get_reports(rep)
    @test isempty(reports)
end

@testset "Concretely typed" begin
    using CheckConcreteStructs
    @test all_concrete(MinimalHomotopyContinuation)

end
