using MinimalHomotopyContinuation
using Test
using Random, LinearAlgebra

using MinimalHomotopyContinuation.DoubleDouble: ComplexDF64

const HC = MinimalHomotopyContinuation
Random.seed!(0x8b868a97)

using JET

@testset "Type stability" begin
    @testset "Utils" begin
        @test (@inferred HC.fast_abs(3.0 + 4.0im)) == 5.0
        @test (@inferred HC.fast_abs(-3.0)) == 3.0
        @test (@inferred HC.nanmin(1.0, 2.0)) == 1.0
        @test (@inferred HC.nanmax(1.0, 2.0)) == 2.0
        @test (@inferred HC.always_false(1)) === false
        @test (@inferred HC.nthroot(16.0, 4)) == 2.0
        @test (@inferred HC.nthroot(8.0, 3)) == 2.0
        @test (@inferred HC.all2(==, [1, 2], [1, 2])) === true
        @test (@inferred HC.all2(==, [1, 2], [1, 3])) === false

        stepper = HC.SegmentStepper(0.0, 2.0)
        stepper_t(s) = s.t
        stepper_t′(s) = s.t′
        stepper_Δt(s) = s.Δt
        stepper_Δs(s) = s.Δs
        @test (@inferred HC.propose_step!(stepper, 0.5)) isa HC.SegmentStepper
        @test (@inferred HC.step_success!(stepper)) isa HC.SegmentStepper
        @test (@inferred HC.dist_to_target(stepper)) isa Float64
        @test (@inferred HC.is_done(stepper)) isa Bool
        @test (@inferred stepper_t(stepper)) isa ComplexF64
        @test (@inferred stepper_t′(stepper)) isa ComplexF64
        @test (@inferred stepper_Δt(stepper)) isa ComplexF64
        @test (@inferred stepper_Δs(stepper)) isa Float64
    end

    @testset "Fixed system dispatch" begin
        @var x y
        sys = System([x^2 + y^2 - 1, x - y])
        @test (@inferred MixedSystem fixed(sys)) isa MixedSystem
        @test (@inferred InterpretedSystem fixed(sys; compile_mode = CompileNone())) isa
            InterpretedSystem
        @test (@inferred CompiledSystem fixed(sys; compile_mode = CompileAll())) isa
            CompiledSystem
    end

    @testset "Norms" begin
        inf_norm = @inferred HC.InfNorm()
        weighted = @inferred HC.WeightedNorm(inf_norm, 2)
        @test (@inferred HC.init!(weighted, [1.0, 2.0])) isa typeof(weighted)
        @test (@inferred HC.update!(weighted, [1.5, 2.5])) isa typeof(weighted)
        @test (@inferred HC.norm([1.0, 2.0], inf_norm)) isa Float64
        @test (@inferred HC.distance([1.0, 2.0], [2.0, 3.0], inf_norm)) isa Float64
        @test (@inferred HC.norm([1.0, 2.0], weighted)) isa Float64
        @test (@inferred HC.distance([1.0, 2.0], [2.0, 3.0], weighted)) isa Float64
    end

    @testset "Linear algebra helpers" begin
        workspace = @inferred HC.MatrixWorkspace(2, 2)
        @test (@inferred HC.matrix(workspace)) isa Matrix{ComplexF64}
        @test (@inferred HC.updated!(workspace)) isa HC.MatrixWorkspace
    end

    @testset "Result and solve helpers" begin
        @var x y
        prob = SystemProblem(System([x^2 + y^2 - 1, x - y]))
        alg = TotalDegreeAlgorithm(seed = UInt32(0x11223344))
        cache = @inferred HC.PathSolveCache HC.init(
            prob, alg; show_progress = false, threading = false
        )
        cache_iter = @inferred HC.PathIteratorSolveCache HC.init(
            prob, alg, ResultIterator; show_progress = false
        )
        @test cache_iter isa HC.PathIteratorSolveCache
        iter = @inferred solve(cache.solver, cache.starts, ResultIterator)
        @test (@inferred seed(iter)) isa UInt32
        @test (@inferred start_solutions(iter)) === cache.starts
        @test (@inferred solver(iter)) === cache.solver
        @test (@inferred length(iter)) isa Int
        @test (@inferred HC.bitmask(HC.is_success, iter)) isa BitVector
        @test (@inferred HC.bitmask_filter(HC.is_success, iter)) isa HC.ResultIterator

        result = @inferred HC.Result(iter)
        @test (@inferred path_results(result)) isa Vector{HC.PathResult}
        @test (@inferred nresults(result)) isa Int
        @test (@inferred nsolutions(result)) isa Int
        @test (@inferred real_solutions(result)) isa Vector{Vector{Float64}}

        opts = @inferred HC._result_filter_options()
        @test (@inferred HC._with_multiple_results(opts, true)) isa HC.ResultFilterOptions
    end

    @testset "Interpreted fast paths" begin
        @var x y
        sys = InterpretedSystem(System([x^2 + y^2 - 1, x - y]))
        u = zeros(ComplexF64, 2)
        J = zeros(ComplexF64, 2, 2)
        x0 = ComplexF64[0.2 + 0im, 0.4 + 0im]
        @test (@inferred evaluate!(u, sys, x0)) isa Vector{ComplexF64}
        @test (@inferred evaluate_and_jacobian!(u, J, sys, x0)) === nothing
        @test (@inferred jacobian!(J, sys, x0)) isa Matrix{ComplexF64}

        @var s a b
        f = [x^2 - a, x * y - a + b]
        g = subs(f, [a, b] => s .* [a, b] .+ (1 .- s) .* [a, b])
        h = Homotopy(g, [x, y], s, [a, b, a, b])
        H = InterpretedHomotopy(h)
        t = 0.5 + 0.0im
        p = [1.0, 2.0, 1.0, 2.0]
        @test (@inferred evaluate!(u, H, x0, t, p)) isa Vector{ComplexF64}
        @test (@inferred evaluate_and_jacobian!(u, J, H, x0, t, p)) === nothing
        @test (@inferred jacobian!(J, H, x0, t, p)) isa Matrix{ComplexF64}
    end

    @testset "JET report_opt regressions" begin
        no_reports(r) = isempty(JET.get_reports(r))

        @var x y
        sys = System([x^2 + y^2 - 1, x - y])
        prob = SystemProblem(sys)
        alg = TotalDegreeAlgorithm(seed = UInt32(0x11223344))
        _ = prob
        _ = alg

        isys = InterpretedSystem(sys)
        u = zeros(ComplexF64, 2)
        J = zeros(ComplexF64, 2, 2)
        x0 = ComplexF64[0.2 + 0im, 0.4 + 0im]
        @test no_reports(@report_opt evaluate!(u, isys, x0))
        @test no_reports(@report_opt evaluate_and_jacobian!(u, J, isys, x0))
        @test no_reports(@report_opt jacobian!(J, isys, x0))

        @test no_reports(@report_opt HC.all2(==, [1, 2], [1, 2]))
    end

    @testset "Newton helpers" begin
        newton_result = HC.NewtonResult(
            HC.NewtonCode.success, ComplexF64[1.0 + 0.0im], 1.0, 1.0, 1, 0.1
        )
        @test (@inferred HC.is_success(newton_result)) === true
        @test (@inferred HC.solution(newton_result)) isa Vector{ComplexF64}

        corrector_result = HC.NewtonCorrectorResult(
            HC.NEWT_CONVERGED, 1.0, 1, 1.0, 1.0, 1.0, 1.0
        )
        @test (@inferred HC.is_converged(corrector_result)) === true

        @test (@inferred HC.t_to_s_plane(1.0 + 1.0im, 2)) isa ComplexF64
    end
end
