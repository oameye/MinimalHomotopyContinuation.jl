using MinimalHomotopyContinuation
using Test
using Random, LinearAlgebra

using MinimalHomotopyContinuation.DoubleDouble: ComplexDF64

const HC = MinimalHomotopyContinuation
Random.seed!(0x8b868a97)

const track = HC.track
const EndgameOptions = HC.EndgameOptions
const TotalDegreeAlgorithm = HC.TotalDegreeAlgorithm

@testset "Endgame" begin
    if !isdefined(@__MODULE__, :TEST_SYSTEM_COLLECTION)
        include("test_systems.jl")
    end

    @testset "Cyclic 7" begin
        res = solve(SystemProblem(cyclic(7)), TotalDegreeAlgorithm(); show_progress = false)
        @test nsolutions(res) == 924
    end

    @testset "Hyperbolic - 6,6" begin
        @var x z
        y = 1
        F = [
            0.75 * x^4 + 1.5 * x^2 * y^2 - 2.5 * x^2 * z^2 + 0.75 * y^4 - 2.5 * y^2 * z^2 +
                0.75 * z^4
            10 * x^2 * z + 10 * y^2 * z - 6 * z^3
        ]
        res = solve(SystemProblem(System(F)), TotalDegreeAlgorithm(); show_progress = false)
        @test count(r -> r.winding_number == 3, path_results(res)) == 12
        @test nresults(res) == 2
        @test nsingular(res) == 2
        @test all(r -> multiplicity(r) == 6, results(res))
    end

    @testset "Singular 1" begin
        @var x y
        z = 1
        F = [
            x^2 + 2 * y^2 + 2 * im * y * z,
            (18 + 3 * im) * x * y + 7 * im * y^2 - (3 - 18 * im) * x * z - 14 * y * z -
                7 * im * z^2,
        ]
        result = solve(
            SystemProblem(System(F)), TotalDegreeAlgorithm(); show_progress = false
        )
        @test nsingular(result) == 1
        @test nresults(result; only_nonsingular = true) == 1
    end

    @testset "Wilkinson $d" for d in [12]
        @var x
        f = System([expand(prod(x - i for i in 1:d))])
        cache = HC.init(
            SystemProblem(f),
            TotalDegreeAlgorithm(endgame_options = EndgameOptions(only_nonsingular = true));
            show_progress = false,
            threading = false,
        )
        res = track.(cache.solver.trackers[1], cache.starts)
        @test all(is_success, res)
        @test round.(Int, real.(sort(first.(solution.(res)); by = abs))) == 1:d
        @test maximum(abs.(imag.(first.(solution.(res))))) < 1.0e-4
    end

    @testset "(x-10)^$d" for d in [2, 6]
        @var x
        f = System([(x - 10)^d])
        cache = HC.init(
            SystemProblem(f), TotalDegreeAlgorithm(); show_progress = false, threading = false
        )
        res = track.(cache.solver.trackers[1], cache.starts)
        @test count(r -> r.winding_number == d, res) == d
    end

    @testset "Beyond Polyhedral Homotopy Example" begin
        @var x y
        f = [2.3 * x^2 + 1.2 * y^2 + 3x - 2y + 3, 2.3 * x^2 + 1.2 * y^2 + 5x + 2y - 5]
        cache = HC.init(
            SystemProblem(System(f)),
            TotalDegreeAlgorithm();
            show_progress = false,
            threading = false,
        )
        res = track.(cache.solver.trackers[1], cache.starts)
        @test count(is_success, res) == 2
        @test count(is_at_infinity, res) == 2
    end

    @testset "Winding Number Family d=$d" for d in 2:2:6
        @var x y
        a = [0.257, -0.139, -1.73, -0.199, 1.79, -1.32]
        f1 = (a[1] * x^d + a[2] * y) * (a[3] * x + a[4] * y) + 1
        f2 = (a[1] * x^d + a[2] * y) * (a[5] * x + a[6] * y) + 1
        cache = HC.init(
            SystemProblem(System([f1, f2])),
            TotalDegreeAlgorithm();
            show_progress = false,
            threading = false,
        )
        res = track.(cache.solver.trackers[1], cache.starts)
        @test count(is_success, res) == d + 1
    end

    @testset "Mohab" begin
        @var x y z
        F = [
            -9091098778555951517 * x^3 * y^4 * z^2 +
                5958442613080401626 * y^2 * z^7 +
                17596733865548170996 * x^2 * z^6 - 17979170986378486474 * x * y * z^6 -
                2382961149475678300 * x^4 * y^3 - 15412758154771986214 * x * y^3 * z^3 + 133,
            -10798198881812549632 * x^6 * y^3 * z - 11318272225454111450 * x * y^9 -
                14291416869306766841 * y^9 * z - 5851790090514210599 * y^2 * z^8 +
                15067068695242799727 * x^2 * y^3 * z^4 +
                7716112995720175148 * x^3 * y * z^3 +
                171,
            13005416239846485183 * x^7 * y^3 + 4144861898662531651 * x^5 * z^4 -
                8026818640767362673 * x^6 - 6882178109031199747 * x^2 * y^4 +
                7240929562177127812 * x^2 * y^3 * z +
                5384944853425480296 * x * y * z^4 +
                88,
        ]
        @time res = solve(
            SystemProblem(System(F, [x, z, y])), TotalDegreeAlgorithm(); show_progress = false
        )
        @test nnonsingular(res) == 693

        @time res = solve(
            SystemProblem(System(F, [x, z, y])),
            TotalDegreeAlgorithm(gamma = -0.9132549847010242 + 0.4073884300256109im),
            show_progress = false,
        )
        @test nnonsingular(res) == 693
    end
end
