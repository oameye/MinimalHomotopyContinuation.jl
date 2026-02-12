using HomotopyContinuation, Test

@testset "Code linting" begin
    using JET
    JET.test_package(
        HomotopyContinuation;
        target_modules = (
            HomotopyContinuation,
            HomotopyContinuation.ModelKit,
            HomotopyContinuation.DoubleDouble,
        ),
    )

    rep = JET.report_package(HomotopyContinuation)
    reports = JET.get_reports(rep)
    @test isempty(reports)
end

@testset "Concretely typed" begin
    using CheckConcreteStructs
    @test all_concrete(HomotopyContinuation)

end
