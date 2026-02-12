const total_degree = HC.total_degree
const polyhedral = HC.polyhedral
const track = HC.track
const mixed_volume = HC.MixedSubdivisions.mixed_volume

@testset "solve" begin

    @testset "total degree (simple)" begin
        @var x y
        affine_sqr = System(
            [
                2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3,
                2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5,
            ]
        )
        @test count(is_success, track.(total_degree(affine_sqr; compile_mode = CompileNone())...)) == 2

        @var x y
        affine_ov = System(
            [
                (x^2 + y^2 + x * y - 3) * (x + 3),
                (x^2 + y^2 + x * y - 3) * (y - x + 2),
                2x + 5y - 3,
            ]
        )
        @test_throws ArgumentError total_degree(affine_ov; compile_mode = CompileNone())

        @var x y
        affine_ov_reordering = System(
            [
                (x^2 + y^2 + x * y - 3) * (x + 3),
                2x + 5y - 3,
                (x^2 + y^2 + x * y - 3) * (y^2 - x + 2),
            ]
        )
        @test_throws ArgumentError total_degree(affine_ov_reordering; compile_mode = CompileNone())

        @var x y z
        homogeneous = System(
            [
                x^2 + y^2 + z^2,
                x * z + y * z,
            ]
        )
        @test_throws ArgumentError total_degree(homogeneous; compile_mode = CompileNone())

        @var x y
        affine_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3])
        @test_throws HC.FiniteException total_degree(affine_underdetermined)
    end

    @testset "polyhedral" begin
        @var x y
        affine_sqr = System(
            [
                2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3,
                2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5,
            ]
        )
        @test count(
            is_success,
            track.(polyhedral(affine_sqr; compile_mode = CompileNone(), show_progress = false)...),
        ) == 2

        @var x y z
        homogeneous = System(
            [
                x^2 + y^2 + z^2,
                x * z + y * z,
            ]
        )
        @test_throws ArgumentError polyhedral(homogeneous; compile_mode = CompileNone())

        @var x y
        affine_ov = System(
            [
                (x^2 + y^2 + x * y - 3) * (x + 3),
                (x^2 + y^2 + x * y - 3) * (y - x + 2),
                2x + 5y - 3,
            ]
        )
        @test_throws ArgumentError polyhedral(affine_ov; compile_mode = CompileNone())

        @var x y
        affine_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3])
        @test_throws HC.FiniteException polyhedral(affine_underdetermined)
    end

    @testset "paths to track" begin
        @var x y
        f = System([2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y])

        @test paths_to_track(SystemProblem(f), TotalDegreeAlgorithm()) == 16
        @test paths_to_track(SystemProblem(f), TotalDegreeAlgorithm()) == 16
        @test paths_to_track(SystemProblem(f), PolyhedralAlgorithm()) == 8
        @test paths_to_track(SystemProblem(f), PolyhedralAlgorithm(only_non_zero = true)) == 3
        @test_throws MethodError paths_to_track(f)
        @test mixed_volume(f) == 3

        @var x y a
        g = System([2y + a * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y], parameters = [a])
        @test paths_to_track(SystemProblem(g), TotalDegreeAlgorithm()) == 16
        @test paths_to_track(SystemProblem(g), PolyhedralAlgorithm()) == 8
    end

    @testset "solve (parameter homotopy)" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        s = [1, 1]
        prob = ParameterHomotopyProblem(
            F,
            [s];
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        res = solve(prob, PathTrackingAlgorithm(); show_progress = false)
        @test nsolutions(res) == 1

        prob = ParameterHomotopyProblem(
            InterpretedSystem(F),
            [s];
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        res = solve(prob, PathTrackingAlgorithm(); threading = false, show_progress = false)
        @test nsolutions(res) == 1

        prob = ParameterHomotopyProblem(
            F,
            s;
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        res = solve(
            prob,
            PathTrackingAlgorithm(compile_mode = CompileNone());
            threading = false,
            show_progress = false,
        )
        @test nsolutions(res) == 1

        @var x a y b
        F = System([x^2 - a], [x, y], [a, b])
        s = [1, 1]
        prob = ParameterHomotopyProblem(
            F,
            [s];
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        @test_throws FiniteException(1) solve(prob, PathTrackingAlgorithm(); show_progress = false)

        @var x a y b z
        F_homogeneous = System([x^2 - a * z^2, x * y + (b - a) * z^2], [x, y, z], [a, b])
        s_homogeneous = [1, 1, 1]
        prob = ParameterHomotopyProblem(
            F_homogeneous,
            [s_homogeneous];
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        @test_throws ArgumentError solve(prob, PathTrackingAlgorithm(); show_progress = false)

        @var x y v w a b
        @test_throws MethodError System(
            [x * y - a * v * w, x^2 - b * v^2],
            parameters = [a, b],
            variable_groups = [[x, v], [y, w]],
        )
    end

    @testset "solve (Homotopy)" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        s = [1, 1]
        H = ParameterHomotopy(F, [1, 0], [2, 4])
        res = solve(HomotopyProblem(H, [s]), PathTrackingAlgorithm(); show_progress = false)
        @test nsolutions(res) == 1
    end

    @testset "solve (Vector{Expression})" begin
        @var x a y b
        F = [x^2 - a, x * y - a + b]
        s = [1, 1]
        prob = ParameterHomotopyProblem(
            F,
            [s];
            parameters = [a, b],
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        res = solve(prob, PathTrackingAlgorithm(); show_progress = false)
        @test nsolutions(res) == 1
    end

    @testset "solve (DynamicPolynomials)" begin
        @polyvar x y
        f₁ = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
        f₂ = x^2 + 2x * y^2 - 2 * y^2 - 1 / 2
        result = solve(SystemProblem([f₁, f₂]), PolyhedralAlgorithm(); show_progress = false)
        @test nsolutions(result) == 18

        @polyvar x a y b
        F = [x^2 - a, x * y - a + b]
        s = [1, 1]
        prob = ParameterHomotopyProblem(
            F,
            [s];
            variables = [x, y],
            parameters = [a, b],
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        res = solve(prob, PathTrackingAlgorithm(compile_mode = CompileNone()); show_progress = false)
        @test nsolutions(res) == 1

        prob2 = ParameterHomotopyProblem(
            F,
            [s];
            variable_ordering = [y, x],
            parameters = [a, b],
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        res2 = solve(
            prob2,
            PathTrackingAlgorithm(compile_mode = CompileNone());
            show_progress = false,
        )
        s = solutions(res)[1]
        s2 = solutions(res2)[1]
        @test s ≈ [s2[2], s2[1]]
    end

    @testset "change parameters" begin
        @var x a y b
        F = System([x^2 - a, x * y - a + b]; parameters = [a, b])
        s = [1.0, 1.0 + 0im]
        prob = ParameterHomotopyProblem(
            F,
            [s];
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        cache = HC.init(
            prob,
            PathTrackingAlgorithm(seed = 0x12345678);
            show_progress = false,
        )
        S = cache.solver
        start_parameters!(S, [1, 0])
        target_parameters!(S, [2, 4])
        @test is_success(track(S, s))
    end

    @testset "solve (threading)" begin
        prob = SystemProblem(cyclic(7))
        alg = PolyhedralAlgorithm()
        res = solve(prob, alg; threading = true, show_progress = false)
        @test nsolutions(res) == 924
    end

    @testset "stop early callback" begin
        @var x
        first_result = nothing
        prob = SystemProblem([(x - 3) * (x + 6) * (x + 2)])
        alg = TotalDegreeAlgorithm()
        results = solve(
            prob,
            alg;
            stop_early_cb = r -> begin
                first_result = r
                true
            end,
            threading = false,
            show_progress = false,
        )
        @test length(results) == 1
        @test first(results) === first_result

        result = let k = 0
            solve(
                prob,
                alg;
                stop_early_cb = r -> (k += 1) == 2,
                show_progress = false,
                threading = false,
            )
        end
        @test length(result) == 2

        @var x y z
        first_result = nothing
        prob = SystemProblem(
            [
                (x - 3) * (x + 6) * (x + 2) * (x - 2) * (x + 2.5),
                (y + 2) * (y - 2) * (y + 3) * (y + 5) * (y - 1),
                (z + 2) * (z - 2) * (z + 3) * (z + 5) * (z - 2.1),
            ],
        )
        results = solve(
            prob,
            TotalDegreeAlgorithm();
            stop_early_cb = r -> begin
                first_result = r
                true
            end,
            threading = true,
            show_progress = false,
        )
        @test length(results) < 125
    end

    @testset "Many parameters solver" begin
        @var x y
        f = x^2 + y^2 - 1

        @var a b c
        l = a * x + b * y + c
        F = [f, l]

        p₀ = randn(ComplexF64, 3)
        S₀ = solutions(
            solve(
                SystemProblem(subs(F, [a, b, c] => p₀)),
                PolyhedralAlgorithm();
                show_progress = false,
            ),
        )
        params = [rand(3) for i in 1:100]

        prob = ParameterSweepProblem(
            F,
            S₀;
            start_parameters = p₀,
            targets = params,
            parameters = [a, b, c],
        )
        result1 = solve(prob, PathTrackingAlgorithm(); threading = true, show_progress = false)
        @test eltype(result1) <: Result

        result1 = solve(prob, PathTrackingAlgorithm(); show_progress = false, threading = false)
        @test eltype(result1) <: Result

        prob2 = ParameterSweepProblem(
            F,
            S₀;
            start_parameters = p₀,
            targets = params,
            parameters = [a, b, c],
            reducer = MapReducer((r, p) -> real_solutions(r)),
        )
        result2 = solve(prob2, PathTrackingAlgorithm(); threading = true, show_progress = false)
        @test typeof(result2) == Vector{Vector{Vector{Float64}}}
        @test !isempty(result2)

        prob3 = ParameterSweepProblem(
            F,
            S₀;
            start_parameters = p₀,
            targets = params,
            parameters = [a, b, c],
            reducer = FlatMapReducer((r, p) -> real_solutions(r)),
        )
        result3 = solve(prob3, PathTrackingAlgorithm(); threading = false, show_progress = false)
        @test typeof(result3) == Vector{Vector{Float64}}
        @test !isempty(result3)

        prob4 = ParameterSweepProblem(
            F,
            S₀;
            start_parameters = p₀,
            targets = [rand(3) for _ in 1:100],
            parameters = [a, b, c],
            reducer = MapReducer((r, p) -> (real_solutions(r), p)),
        )
        result4 = solve(prob4, PathTrackingAlgorithm(); show_progress = false)
        @test typeof(result4) == Vector{Tuple{Vector{Vector{Float64}}, Vector{Float64}}}

        @testset "Many parameters threaded" begin
            @var u1, v1, ω, α, γ, λ, ω0

            eqs = [
                -u1 * ω^2 +
                    u1 * ω0^2 +
                    (3 / 4) * u1^3 * α +
                    (3 / 4) * u1 * v1^2 * α +
                    (-1 / 2) * u1 * λ * ω0^2 +
                    v1 * γ * ω,
                -v1 * ω^2 + v1 * ω0^2 + (3 / 4) * v1^3 * α - u1 * γ * ω +
                    (3 / 4) * u1^2 * v1 * α +
                    (1 / 2) * v1 * λ * ω0^2,
            ]

            F = System(eqs, parameters = [ω, α, γ, λ, ω0], variables = [u1, v1])

            input_array = [
                [0.9, 1.0, 0.01, 0.01, 1.1],
                [0.9105263157894737, 1.0, 0.01, 0.01, 1.1],
                [0.9210526315789473, 1.0, 0.01, 0.01, 1.1],
                [0.9315789473684211, 1.0, 0.01, 0.01, 1.1],
                [0.9421052631578948, 1.0, 0.01, 0.01, 1.1],
                [0.9526315789473684, 1.0, 0.01, 0.01, 1.1],
            ]

            generic_parameters = randn(ComplexF64, 5)

            R0 = solve(
                SystemProblem(F; target_parameters = generic_parameters),
                PolyhedralAlgorithm();
                threading = true,
                show_progress = false,
            )
            R1 = solve(
                ParameterSweepProblem(
                    F,
                    solutions(R0);
                    start_parameters = generic_parameters,
                    targets = input_array,
                ),
                PathTrackingAlgorithm();
                threading = true,
                show_progress = false,
            )

            @test length(R1) == 6
        end
    end

    @testset "CommonSolve interface contract" begin
        @var x y
        sysprob = SystemProblem(System([x^2 + y^2 - 1, x - y]))
        alg = TotalDegreeAlgorithm(seed = 0x11112222)
        cache = HC.init(sysprob, alg; show_progress = false, threading = false)
        res1 = HC.solve!(cache)
        res2 = solve(sysprob, alg; show_progress = false, threading = false)
        @test nsolutions(res1) == nsolutions(res2)

        cache2 = HC.init(sysprob, alg; show_progress = false, threading = false)
        rA = HC.solve!(cache2)
        rB = HC.solve!(cache2)
        @test nsolutions(rA) == nsolutions(rB)

        @var u v a
        param_prob = ParameterHomotopyProblem(
            System([u^2 - a, u + v], [u, v], [a]),
            [[1.0, -1.0]];
            start_parameters = [1.0],
            target_parameters = [2.0],
        )
        c = HC.init(
            param_prob,
            PathTrackingAlgorithm();
            show_progress = false,
            threading = false,
        )
        r = HC.solve!(c)
        @test r isa Result
    end

    @testset "Inference" begin
        @var x y
        base_prob = SystemProblem(System([x^2 + y^2 - 1, x - y]))
        base_alg = TotalDegreeAlgorithm(seed = UInt32(0x11223344))

        cache = @inferred HC.PathSolveCache HC.init(
            base_prob,
            base_alg;
            show_progress = false,
            threading = false,
        )
        @test cache isa HC.PathSolveCache

        res = @inferred solve(base_prob, base_alg; threading = false, show_progress = false)
        @test res isa Result

        npaths = @inferred paths_to_track(base_prob, base_alg)
        @test npaths isa Int

        @var a b
        Fp = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        starts = [[1.0, 1.0], [-1.0, -1.0]]
        targets = [[2.0, 4.0], [3.0, 5.0]]
        sweep_alg = PathTrackingAlgorithm(seed = UInt32(0x55667788))

        sweep_id = ParameterSweepProblem(
            Fp,
            starts;
            start_parameters = [1.0, 0.0],
            targets = targets,
            reducer = IdentityReducer(),
        )
        r_id = @inferred solve(sweep_id, sweep_alg; threading = false, show_progress = false)
        @test r_id isa Vector{<:Result}

        sweep_map = ParameterSweepProblem(
            Fp,
            starts;
            start_parameters = [1.0, 0.0],
            targets = targets,
            reducer = MapReducer((r, p) -> nsolutions(r)),
        )
        r_map = @inferred solve(sweep_map, sweep_alg; threading = false, show_progress = false)
        @test r_map isa Vector{Int}

        sweep_flat = ParameterSweepProblem(
            Fp,
            starts;
            start_parameters = [1.0, 0.0],
            targets = targets,
            reducer = FlatMapReducer((r, p) -> real_solutions(r)),
        )
        r_flat =
            @inferred solve(sweep_flat, sweep_alg; threading = false, show_progress = false)
        @test r_flat isa Vector{Vector{Float64}}
    end

    @testset "Removed symbol-based APIs" begin
        @var x y
        f = System([x^2 + y^2 - 1, x - y])

        @test_throws MethodError PolyhedralAlgorithm(compile_mode = :none)
        @test_throws MethodError TotalDegreeAlgorithm(compile_mode = :none)
        @test_throws MethodError PathTrackingAlgorithm(compile_mode = :none)

        @test_throws TypeError HC.TrackerOptions(parameter_preset = :fast)

        @test_throws MethodError paths_to_track(SystemProblem(f), PolyhedralAlgorithm(); start_system = :total_degree)
        @test_throws MethodError paths_to_track(f; start_system = :polyhedral)

        starts = [[1.0, 1.0]]
        @test_throws MethodError ParameterSweepProblem(
            f,
            starts;
            start_parameters = [1.0, 0.0],
            targets = [[2.0, 4.0]],
            transform_result = identity,
        )
        @test_throws MethodError ParameterSweepProblem(
            f,
            starts;
            start_parameters = [1.0, 0.0],
            targets = [[2.0, 4.0]],
            flatten = true,
        )
    end

end
