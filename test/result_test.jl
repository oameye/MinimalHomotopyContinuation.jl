@testset "Result" begin
    @testset "Basic functionality of Result" begin
        d = 2
        @var x y a[1:6]
        F = System(
            [
                (a[1] * x^d + a[2] * y) * (a[3] * x + a[4] * y) + 1,
                (a[1] * x^d + a[2] * y) * (a[5] * x + a[6] * y) + 1,
            ];
            parameters = a,
        )
        prob = SystemProblem(
            F; target_parameters = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32]
        )
        res = solve(prob, PolyhedralAlgorithm(); show_progress = false)

        @test startswith(sprint(show, res), "Result with 3 solutions")
        @test seed(res) isa UInt32

        seeded_res = solve(prob, PolyhedralAlgorithm(seed = seed(res)), show_progress = false)
        @test seed(seeded_res) == seed(res)

        @test length(path_results(res)) == ntracked(res) == 7
        @test length(results(res)) == nresults(res) == 3
        @test length(solutions(res)) == 3
        @test length(findall(is_success, res)) == 3
        @test real_solutions(res) isa Vector{Vector{Float64}}
        @test length(real_solutions(res)) == nreal(res) == 1
        @test length(real_solutions(res; atol = 1.0e-16)) == 1
        @test length(real_solutions(res; atol = 0.0)) == 0
        @test length(real_solutions(res; atol = 0.0, rtol = 1.0e-16)) == 1
        @test length(real_solutions(res; atol = Inf, rtol = Inf)) == 3
        @test length(real_solutions(res; atol = 1.0, rtol = 1.0e-16)) == 1
        @test length(nonsingular(res)) == nnonsingular(res) == 3
        @test isempty(singular(res))
        @test nsingular(res) == 0
        @test nexcess_solutions(res) == 0
        @test !isempty(sprint(show, statistics(res)))

        @var x y
        g = System([(29 / 16) * x^3 - 2 * x * y, x^2 - y])
        res = solve(SystemProblem(g), TotalDegreeAlgorithm(); show_progress = false)
        @test startswith(sprint(show, res), "Result with 1 solution")
        @test seed(res) isa UInt32
        @test !isempty(sprint(show, statistics(res)))
        @test !isempty(sprint(show, res))
    end
end

@testset "ResultIterator" begin

    @testset "Result tests" begin
        d = 2
        @var x y a[1:6]
        F = System(
            [
                (a[1] * x^d + a[2] * y) * (a[3] * x + a[4] * y) + 1,
                (a[1] * x^d + a[2] * y) * (a[5] * x + a[6] * y) + 1,
            ];
            parameters = a,
        )
        param = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32]
        res = solve(
            SystemProblem(F; target_parameters = param),
            PolyhedralAlgorithm(),
            ResultIterator;
            show_progress = false,
        )

        @test startswith(sprint(show, res), "ResultIterator")
        @test seed(res) isa UInt32

        @test length(path_results(res)) == ntracked(res) == 7
        @test nresults(res) == 3
        @test nsolutions(res) == 3
        real_sols = collect(real_solutions(res))
        @test length(real_sols) == nreal(res) == 1
        @test nnonsingular(res) == 3
        @test isempty(singular(res))
        @test nsingular(res) == 0
        @test nexcess_solutions(res) == 0

        B = BitVector([1, 0, 0, 0, 0, 0, 0])
        res_B = solve(
            SystemProblem(F; target_parameters = param),
            PolyhedralAlgorithm(),
            ResultIterator;
            bitmask = B,
            show_progress = false,
        )
        @test length(res_B) == 1

        @var x y
        g = System([(29 / 16) * x^3 - 2 * x * y, x^2 - y])
        res = solve(SystemProblem(g), TotalDegreeAlgorithm(), ResultIterator)
        @test startswith(sprint(show, res), "ResultIterator")
        @test seed(res) isa UInt32
        @test !isempty(sprint(show, res))

        @var x y p
        f₁ = y - x^2 + p
        f₂ = y - x^3 - p
        F = System([f₁, f₂]; variables = [x; y], parameters = [p])
        R = solve(
            ParameterHomotopyProblem(
                F, [[1, 1], [-1, 1]]; start_parameters = [0], target_parameters = [-1]
            ),
            PathTrackingAlgorithm(),
            ResultIterator,
        )
        RR = solve(
            ParameterHomotopyProblem(F, R; start_parameters = [-1], target_parameters = [-2]),
            PathTrackingAlgorithm(),
            ResultIterator,
        )
        @test length(collect(RR)) == 2
    end

    @testset "Basic functionality of ResultIterator" begin
        @var x y
        f₁ = y - x^2
        f₂ = y - x^3
        F = System([f₁, f₂])
        tsi_polyhedral = solve(
            SystemProblem(F), PolyhedralAlgorithm(), ResultIterator; show_progress = false
        )
        tsi_total_degree = solve(SystemProblem(F), TotalDegreeAlgorithm(), ResultIterator)

        @test isa(tsi_polyhedral, ResultIterator)
        @test isa(tsi_total_degree, ResultIterator)

        @test nsolutions(tsi_polyhedral; only_nonsingular = false, multiple_results = true) == 3
        @test nsolutions(tsi_total_degree; multiple_results = true) == 1
        @test length(tsi_polyhedral) == 3
        @test length(tsi_total_degree) == 6

        BM = bitmask_filter(isfinite, tsi_total_degree)
        @test length(BM) == sum(bitmask(isfinite, tsi_total_degree))
        @test length(BM) > 0

        t = trace(BM)
        @test norm([1.0 + 0.0im, 1.0 + 0.0im] - t) < 1.0e-12
    end

    @testset "Manual start solutions" begin
        @var x y p
        f₁ = y - x^2 + p
        f₂ = y - x^3 - p
        F = System([f₁, f₂]; variables = [x; y], parameters = [p])
        R = solve(
            ParameterHomotopyProblem(
                F, [1, 1]; start_parameters = [0], target_parameters = [-1]
            ),
            PathTrackingAlgorithm(),
            ResultIterator,
        )

        @test isa(R, ResultIterator)
        @test nsolutions(R) == 1

        R2 = solve(
            ParameterHomotopyProblem(
                F, [[1, 1], [1, 1]]; start_parameters = [0], target_parameters = [-1]
            ),
            PathTrackingAlgorithm(),
            ResultIterator,
        )

        @test isa(R2, ResultIterator)
        @test nsolutions(R2) == 2
    end

    @testset "Many parameters" begin
        @var x y
        f = x^2 + y^2 - 1

        @var a b c
        l = a * x + b * y + c
        F = System([f, l]; parameters = [a, b, c])

        p₀ = randn(ComplexF64, 3)
        S₀ = solutions(
            solve(
                SystemProblem(System(subs(F.expressions, [a, b, c] => p₀))),
                PolyhedralAlgorithm();
                show_progress = false,
            ),
        )
        params = [rand(3) for i in 1:100]

        result1 = solve(
            ParameterSweepProblem(F, S₀; start_parameters = p₀, targets = params),
            PathTrackingAlgorithm(),
            ResultIterator,
        )
        r1 = collect.(result1)
        @test !isempty(r1)

        result1 = solve(
            ParameterSweepProblem(F, S₀; start_parameters = p₀, targets = params),
            PathTrackingAlgorithm(),
            ResultIterator,
        )
        r1 = collect.(result1)
        @test !isempty(r1)

        result2 = solve(
            ParameterSweepProblem(
                F,
                S₀;
                start_parameters = p₀,
                targets = params,
                reducer = MapReducer((r, p) -> real_solutions(r)),
            ),
            PathTrackingAlgorithm(),
            ResultIterator,
        )
        r2 = collect.(result2)
        @test !isempty(r2)

        result3 = solve(
            ParameterSweepProblem(
                F,
                S₀;
                start_parameters = p₀,
                targets = params,
                reducer = FlatMapReducer((r, p) -> real_solutions(r)),
            ),
            PathTrackingAlgorithm(),
            ResultIterator,
        )
        r3 = collect.(result3)
        @test !isempty(r3)

        result4 = solve(
            ParameterSweepProblem(
                F,
                S₀;
                start_parameters = p₀,
                targets = [rand(3) for _ in 1:100],
                reducer = MapReducer((r, p) -> (real_solutions(r), p)),
            ),
            PathTrackingAlgorithm(),
            ResultIterator,
        )
        r4 = collect.(result4)
        @test !isempty(r4)
    end

end
