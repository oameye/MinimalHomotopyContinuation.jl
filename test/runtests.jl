using HomotopyContinuation
using LinearAlgebra, Test, Parameters, Random
using Arblib
using HomotopyContinuation.DoubleDouble: ComplexDF64
const HC = HomotopyContinuation


@testset "Concretely typed" begin
    using CheckConcreteStructs
    @test all_concrete(HomotopyContinuation)

end

@testset "ExplicitImports" begin
    using ExplicitImports

    non_public_explicit_imports = (
        :init,
        :MPFR,
        :MPFRRoundingMode,
        :MPZ,
        :solve,
        :solve!,
    )
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
        HomotopyContinuation;
        ignore = non_public_explicit_imports,
    ) == nothing
    @test check_no_stale_explicit_imports(HomotopyContinuation) == nothing
    @test check_all_qualified_accesses_via_owners(HomotopyContinuation) == nothing
    @test check_all_qualified_accesses_are_public(
        HomotopyContinuation;
        ignore = non_public_qualified_accesses,
    ) == nothing
    @test check_no_self_qualified_accesses(HomotopyContinuation) == nothing
end

@testset "Code quality" begin
    using Aqua
    Aqua.test_all(HomotopyContinuation)
end

@testset "Code linting" begin
    using JET
    JET.test_package(HomotopyContinuation; target_modules = (HomotopyContinuation, HomotopyContinuation.ModelKit, HomotopyContinuation.DoubleDouble))
end


include("test_systems.jl")

Random.seed!(0x8b868a97)

# include("./model_kit/symbolic_test.jl")
# include("./model_kit/operations_test.jl")
# include("./model_kit/e2e_test.jl")
# include("./model_kit/slp_test.jl")

include("systems_test.jl")
include("tracker_test.jl")
include("endgame_tracker_test.jl")
include("polyhedral_test.jl")
include("solve_test.jl")
include("result_test.jl")
include("endgame_test.jl")
