@testset "solve" begin

    @testset "total degree (simple)" begin
        @var x y
        affine_sqr = System(
            [
                2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3,
                2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5,
            ]
        )
        @test count(is_success, track.(total_degree(affine_sqr; compile = false)...)) == 2

        @var x y
        affine_ov = System(
            [
                (x^2 + y^2 + x * y - 3) * (x + 3),
                (x^2 + y^2 + x * y - 3) * (y - x + 2),
                2x + 5y - 3,
            ]
        )
        @test_throws ArgumentError total_degree(affine_ov; compile = false)

        @var x y
        affine_ov_reordering = System(
            [
                (x^2 + y^2 + x * y - 3) * (x + 3),
                2x + 5y - 3,
                (x^2 + y^2 + x * y - 3) * (y^2 - x + 2),
            ]
        )
        @test_throws ArgumentError total_degree(affine_ov_reordering; compile = false)

        @var x y z
        homogeneous = System(
            [
                x^2 + y^2 + z^2,
                x * z + y * z,
            ]
        )
        @test_throws ArgumentError total_degree(homogeneous; compile = false)

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
            track.(polyhedral(affine_sqr; compile = false, show_progress = false)...),
        ) == 2

        @var x y z
        homogeneous = System(
            [
                x^2 + y^2 + z^2,
                x * z + y * z,
            ]
        )
        @test_throws ArgumentError polyhedral(homogeneous; compile = false)

        @var x y
        affine_ov = System(
            [
                (x^2 + y^2 + x * y - 3) * (x + 3),
                (x^2 + y^2 + x * y - 3) * (y - x + 2),
                2x + 5y - 3,
            ]
        )
        @test_throws ArgumentError polyhedral(affine_ov; compile = false)

        @var x y
        affine_underdetermined = System([2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3])
        @test_throws HC.FiniteException polyhedral(affine_underdetermined)
    end

    @testset "paths to track" begin
        @var x y
        f = System([2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y])

        @test paths_to_track(f; start_system = :total_degree) == 16
        @test paths_to_track(f; start_system = :total_degree) == 16
        @test paths_to_track(f; start_system = :polyhedral) == 8
        @test paths_to_track(f; start_system = :polyhedral, only_non_zero = true) == 3
        @test paths_to_track(f) == 8
        @test mixed_volume(f) == 3

        @var x y a
        g = System([2y + a * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y], parameters = [a])
        @test paths_to_track(g; start_system = :total_degree) == 16
        @test paths_to_track(g; start_system = :polyhedral) == 8
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
        res = solve(prob; show_progress = false)
        @test nsolutions(res) == 1

        prob = ParameterHomotopyProblem(
            InterpretedSystem(F),
            [s];
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        res = solve(prob; threading = false, show_progress = false)
        @test nsolutions(res) == 1

        prob = ParameterHomotopyProblem(
            F,
            s;
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        res = solve(prob; threading = false, compile = false, show_progress = false)
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
        @test_throws FiniteException(1) solve(prob; show_progress = false)

        @var x a y b z
        F_homogeneous = System([x^2 - a * z^2, x * y + (b - a) * z^2], [x, y, z], [a, b])
        s_homogeneous = [1, 1, 1]
        prob = ParameterHomotopyProblem(
            F_homogeneous,
            [s_homogeneous];
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        @test_throws ArgumentError solve(prob; show_progress = false)

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
        res = solve(HomotopyProblem(H, [s]); show_progress = false)
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
        res = solve(prob; show_progress = false)
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
        res = solve(prob; compile = false, show_progress = false)
        @test nsolutions(res) == 1

        prob2 = ParameterHomotopyProblem(
            F,
            [s];
            variable_ordering = [y, x],
            parameters = [a, b],
            start_parameters = [1, 0],
            target_parameters = [2, 4],
        )
        res2 = solve(prob2; compile = false, show_progress = false)
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
        cache = HC.init(prob; seed = 0x12345678, show_progress = false)
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
        result1 = solve(prob; threading = true, show_progress = false)
        @test eltype(result1) <: Tuple{<:Result, Vector{Float64}}

        result1 = solve(prob; show_progress = false, threading = false)
        @test eltype(result1) <: Tuple{<:Result, Vector{Float64}}

        prob2 = ParameterSweepProblem(
            F,
            S₀;
            start_parameters = p₀,
            targets = params,
            parameters = [a, b, c],
            transform_result = (r, p) -> real_solutions(r),
        )
        result2 = solve(prob2; threading = true, show_progress = false)
        @test typeof(result2) == Vector{Vector{Vector{Float64}}}
        @test !isempty(result2)

        prob3 = ParameterSweepProblem(
            F,
            S₀;
            start_parameters = p₀,
            targets = params,
            parameters = [a, b, c],
            transform_result = (r, p) -> real_solutions(r),
            flatten = true,
        )
        result3 = solve(prob3; threading = false, show_progress = false)
        @test typeof(result3) == Vector{Vector{Float64}}
        @test !isempty(result3)

        prob4 = ParameterSweepProblem(
            F,
            S₀;
            start_parameters = p₀,
            targets = 1:100,
            parameters = [a, b, c],
            transform_result = (r, p) -> (real_solutions(r), p),
            transform_parameters = _ -> rand(3),
        )
        result4 = solve(prob4; show_progress = false)
        @test typeof(result4) == Vector{Tuple{Vector{Vector{Float64}}, Int64}}

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
                );
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
        c = HC.init(param_prob; show_progress = false, threading = false)
        r = HC.solve!(c)
        @test r isa Result
    end

    @testset "Breaking API" begin
        @var x
        F = System([x^2 - 1])
        @test_throws MethodError solve(F)
        @test_throws MethodError solve(F, [[1.0]])
        @test_throws MethodError solve(SystemProblem(F), TotalDegreeAlgorithm(); iterator_only = true)
    end
end
