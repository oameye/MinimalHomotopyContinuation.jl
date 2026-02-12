using BenchmarkTools
using HomotopyContinuation
using Statistics

function _median_ns(b::BenchmarkTools.Benchmark)
    tr = run(b)
    return median(tr).time
end

function _fmt_ms(ns)
    return round(ns / 1.0e6; sigdigits = 5)
end

function bench_small_total_degree()
    @var x y
    f = System([x^2 + y^2 - 1, x - y])
    prob = SystemProblem(f)
    alg = TotalDegreeAlgorithm(seed = UInt32(0x11112222))
    warm = _median_ns(@benchmarkable solve($prob, $alg; show_progress = false, threading = false))
    steady = _median_ns(@benchmarkable solve($prob, $alg; show_progress = false, threading = false) samples = 50)
    return (warm = warm, steady = steady)
end

function bench_medium_polyhedral()
    @polyvar x y
    f₁ = (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y
    f₂ = x^2 + 2x * y^2 - 2 * y^2 - 1 / 2
    prob = SystemProblem([f₁, f₂])
    alg = PolyhedralAlgorithm(seed = UInt32(0x33334444))
    warm = _median_ns(@benchmarkable solve($prob, $alg; show_progress = false, threading = false))
    steady = _median_ns(@benchmarkable solve($prob, $alg; show_progress = false, threading = false) samples = 30)
    return (warm = warm, steady = steady)
end

function bench_parameter_homotopy()
    @var x y a b
    F = System([x^2 - a, x * y - a + b], [x, y], [a, b])
    prob = ParameterHomotopyProblem(
        F,
        [[1.0, 1.0], [-1.0, -1.0]];
        start_parameters = [1.0, 0.0],
        target_parameters = [2.0, 4.0],
    )
    alg = PathTrackingAlgorithm(seed = UInt32(0x55556666))
    warm = _median_ns(@benchmarkable solve($prob, $alg; show_progress = false, threading = false))
    steady = _median_ns(@benchmarkable solve($prob, $alg; show_progress = false, threading = false) samples = 50)
    return (warm = warm, steady = steady)
end

function bench_ttfx()
    @var x y
    prob = SystemProblem(System([x^2 + y^2 - 1, x - y]))
    alg = TotalDegreeAlgorithm(seed = UInt32(0x77778888))
    t = @elapsed solve(prob, alg; show_progress = false, threading = false)
    return t
end

function run_refactor_baseline()
    t_first = bench_ttfx()
    small = bench_small_total_degree()
    medium = bench_medium_polyhedral()
    param = bench_parameter_homotopy()

    println("Refactor baseline (times in ms)")
    println("first_solve: ", round(t_first * 1000; sigdigits = 5))
    println(
        "small_total_degree: warm=",
        _fmt_ms(small.warm),
        ", steady=",
        _fmt_ms(small.steady),
    )
    println(
        "medium_polyhedral: warm=",
        _fmt_ms(medium.warm),
        ", steady=",
        _fmt_ms(medium.steady),
    )
    println(
        "parameter_homotopy: warm=",
        _fmt_ms(param.warm),
        ", steady=",
        _fmt_ms(param.steady),
    )

    return (
        first_solve_seconds = t_first,
        small_total_degree = small,
        medium_polyhedral = medium,
        parameter_homotopy = param,
    )
end

run_refactor_baseline()
