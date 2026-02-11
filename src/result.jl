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

export ResultIterator, bitmask, bitmask_filter, solver, start_solutions
export trace

include("result/types.jl")
include("result/iterate.jl")
include("result/filter.jl")
include("result/show.jl")
