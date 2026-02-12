const track = HC.track
const TrackerOptions = HC.TrackerOptions
const PolyhedralAlgorithm = HC.PolyhedralAlgorithm

@testset "Polyhedral" begin
    @testset "affine + torus solutions" begin
        @var x y
        f = System([2y + 3 * y^2 - x * y^3, x + 4 * x^2 - 2 * x^3 * y])

        cache = HC.init(
            SystemProblem(f),
            PolyhedralAlgorithm(only_torus = false);
            show_progress = false,
            threading = false,
        )
        tracker, starts = cache.solver.trackers[1], cache.starts
        @test length(collect(starts)) == 8
        @test count(is_success, track.(tracker, starts)) == 6

        cache = HC.init(
            SystemProblem(f),
            PolyhedralAlgorithm(only_torus = true);
            show_progress = false,
            threading = false,
        )
        tracker, starts = cache.solver.trackers[1], cache.starts
        @test length(collect(starts)) == 3
        @test count(is_success, track.(tracker, starts)) == 3

        cache = HC.init(
            SystemProblem(f),
            PolyhedralAlgorithm(
                only_torus = true, tracker_options = TrackerOptions(automatic_differentiation = 3)
            );
            show_progress = false,
            threading = false,
        )
        tracker, starts = cache.solver.trackers[1], cache.starts
        @test length(collect(starts)) == 3
        @test count(is_success, track.(tracker, starts)) == 3
    end

    @testset "only torus" begin
        @var x₁ x₂ s

        HC_I = [
            x₁^3 * s^15 + x₁ * x₂ * s + x₂^3 + s^12,
            x₁^2 * s^9 + x₁ * x₂^2 + x₂ * s^3,
            x₁^2 * x₂ * s^5 + x₁ * s^8 + x₂^2,
        ]
        F = System(HC_I)
        starts = HC.init(
            SystemProblem(F), PolyhedralAlgorithm(only_torus = false); show_progress = false
        ).starts
        @test length(collect(starts)) ==
            paths_to_track(SystemProblem(F), PolyhedralAlgorithm(only_torus = false)) ==
            92
        starts = HC.init(
            SystemProblem(F), PolyhedralAlgorithm(only_torus = true); show_progress = false
        ).starts
        @test length(collect(starts)) ==
            paths_to_track(SystemProblem(F), PolyhedralAlgorithm(only_torus = true)) ==
            54
    end

    @testset "cyclic" begin
        cache = HC.init(
            SystemProblem(cyclic(5)),
            PolyhedralAlgorithm();
            show_progress = false,
            threading = false,
        )
        tracker, starts = cache.solver.trackers[1], cache.starts
        res = track.(tracker, starts)
        @test count(is_success, res) == 70

        cache = HC.init(
            SystemProblem(cyclic(7)),
            PolyhedralAlgorithm();
            show_progress = false,
            threading = false,
        )
        tracker, starts = cache.solver.trackers[1], cache.starts
        res = track.(tracker, starts)
        @test count(is_success, res) == 924
    end
end
