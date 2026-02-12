using HomotopyContinuation
using LinearAlgebra, Test, Parameters, Random
using Arblib
using HomotopyContinuation.DoubleDouble: ComplexDF64
const HC = HomotopyContinuation

@testset "Public API contract" begin
    @test Base.isexported(HC, :SystemProblem)
    @test Base.isexported(HC, :ParameterHomotopyProblem)
    @test Base.isexported(HC, :ParameterSweepProblem)
    @test Base.isexported(HC, :HomotopyProblem)
    @test Base.isexported(HC, :PolyhedralAlgorithm)
    @test Base.isexported(HC, :TotalDegreeAlgorithm)
    @test Base.isexported(HC, :PathTrackingAlgorithm)
    @test Base.isexported(HC, :CompileAll)
    @test Base.isexported(HC, :CompileMixed)
    @test Base.isexported(HC, :CompileNone)
    @test Base.isexported(HC, :init)
    @test Base.isexported(HC, :solve)
    @test Base.isexported(HC, :solve!)
    @test Base.isexported(HC, :Result)
    @test Base.isexported(HC, :ResultIterator)
    @test Base.isexported(HC, :paths_to_track)

    @test !Base.isexported(HC, :Tracker)
    @test !Base.isexported(HC, :EndgameTracker)
    @test !Base.isexported(HC, :total_degree)
    @test !Base.isexported(HC, :polyhedral)
    @test !Base.isexported(HC, :newton)
    @test !Base.isexported(HC, :MatrixWorkspace)
    @test !Base.isexported(HC, :track)
end


@testset "ExplicitImports" begin
    using ExplicitImports

    non_public_explicit_imports = (:init, :MPFR, :MPFRRoundingMode, :MPZ, :solve, :solve!)
    non_public_qualified_accesses = (
        HomotopyContinuation.ModelKit,
        HomotopyContinuation.DoubleDouble,
        :FastMath,
        :HasEltype,
        :HasLength,
        :MPFRRoundNearest,
        :ROUNDING_MODE,
        :RefValue,
        :SizeUnknown,
        :add!,
        :add_ui!,
        :broadcastable,
        :checked_add,
        :checked_mul,
        :div_fast,
        :foreachfield,
        :gcdext!,
        :inv!,
        :max_fast,
        :mpz_t,
        :mul!,
        :mul_2exp!,
        :mul_si!,
        :qrfactUnblocked!,
        :set!,
        :set_si!,
        :set_ui!,
        :sizeinbase,
        :tdiv_q!,
        :tty_width,
    )

    @test check_no_implicit_imports(HomotopyContinuation) == nothing
    @test check_all_explicit_imports_via_owners(HomotopyContinuation) == nothing
    @test check_all_explicit_imports_are_public(
        HomotopyContinuation; ignore = non_public_explicit_imports
    ) == nothing
    @test check_no_stale_explicit_imports(HomotopyContinuation) == nothing
    @test check_all_qualified_accesses_via_owners(HomotopyContinuation) == nothing
    @test check_all_qualified_accesses_are_public(
        HomotopyContinuation; ignore = non_public_qualified_accesses
    ) == nothing
    @test check_no_self_qualified_accesses(HomotopyContinuation) == nothing
end
