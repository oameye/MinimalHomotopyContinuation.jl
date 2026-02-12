#
# The stable public API surface is defined here.
# Low-level internals (tracker/newton/start-system constructors) are intentionally not exported.
#

export AbstractCompileMode, CompileAll, CompileMixed, CompileNone, DEFAULT_COMPILE_MODE

export AbstractHCAlgorithm,
    PolyhedralAlgorithm,
    TotalDegreeAlgorithm,
    PathTrackingAlgorithm,
    AbstractSweepReducer,
    IdentityReducer,
    MapReducer,
    FlatMapReducer

export AbstractHCProblem,
    SystemProblem, ParameterHomotopyProblem, ParameterSweepProblem, HomotopyProblem

export init, solve, solve!, paths_to_track

export Result,
    ResultStatistics,
    statistics,
    seed,
    algorithm,
    path_results,
    results,
    nresults,
    solutions,
    real_solutions,
    nonsingular,
    singular,
    at_infinity,
    failed,
    nfinite,
    nsolutions,
    nsingular,
    nat_infinity,
    nexcess_solutions,
    nfailed,
    nnonsingular,
    nreal,
    ntracked

export ResultIterator, bitmask, bitmask_filter, solver, start_solutions, trace

export PathResult,
    PathResultCode,
    solution,
    accuracy,
    residual,
    steps,
    accepted_steps,
    rejected_steps,
    winding_number,
    path_number,
    start_solution,
    multiplicity,
    last_path_point,
    valuation,
    is_success,
    is_at_infinity,
    is_excess_solution,
    is_failed,
    is_finite,
    is_singular,
    is_nonsingular,
    is_real
