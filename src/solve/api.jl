const PathProblems = Union{SystemProblem, ParameterHomotopyProblem, HomotopyProblem}

function _cache_from_startsolutions(
        startsolutions::Tuple,
        ::PathProblems;
        stop_early_cb,
        show_progress::Bool,
        threading::Bool,
        catch_interrupt::Bool,
    )
    solver, starts = startsolutions
    return PathSolveCache(
        solver, starts, stop_early_cb, show_progress, threading, catch_interrupt
    )
end

function _cache_from_startsolutions(
        startsolutions::Tuple,
        prob::ParameterSweepProblem;
        stop_early_cb,
        show_progress::Bool,
        threading::Bool,
        catch_interrupt::Bool,
    )
    _ = stop_early_cb
    solver, starts = startsolutions
    return SweepSolveCache(
        solver,
        starts,
        prob.targets,
        prob.reducer,
        show_progress,
        threading,
        catch_interrupt,
    )
end

function _cache_from_startsolutions(startsolutions::Tuple, prob::AbstractHCProblem; kwargs...)
    _ = startsolutions
    _ = kwargs
    throw(ArgumentError("Unsupported cache construction for problem type $(typeof(prob))."))
end

function _iter_cache_from_startsolutions(
        startsolutions::Tuple,
        ::PathProblems,
        ::AbstractHCAlgorithm;
        bitmask = nothing,
    )
    solver, starts = startsolutions
    return PathIteratorSolveCache(solver, starts, bitmask)
end

function _iter_cache_from_startsolutions(
        startsolutions::Tuple,
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm;
        bitmask = nothing,
    )
    _, starts = startsolutions
    solvers = map(p -> _sweep_solver(prob, alg, p), prob.targets)
    return SweepIteratorSolveCache(solvers, starts, bitmask)
end

function _iter_cache_from_startsolutions(
        startsolutions::Tuple,
        prob::AbstractHCProblem,
        alg::AbstractHCAlgorithm;
        bitmask = nothing,
    )
    _ = startsolutions
    _ = alg
    _ = bitmask
    throw(
        ArgumentError(
            "Unsupported iterator cache construction for problem type $(typeof(prob)).",
        ),
    )
end

function init(
        prob::AbstractHCProblem,
        alg::AbstractHCAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    startsolutions = _solver_startsolutions(prob, alg; show_progress)
    return _cache_from_startsolutions(
        startsolutions,
        prob;
        stop_early_cb,
        show_progress,
        threading,
        catch_interrupt,
    )
end

function _cache_solve_kwargs(cache::PathSolveCache)
    return (
        stop_early_cb = cache.stop_early_cb,
        show_progress = cache.show_progress,
        threading = cache.threading,
        catch_interrupt = cache.catch_interrupt,
    )
end

function _cache_solve_kwargs(cache::SweepSolveCache)
    return (
        reducer = cache.reducer,
        show_progress = cache.show_progress,
        threading = cache.threading,
        catch_interrupt = cache.catch_interrupt,
    )
end

function init(
        prob::AbstractHCProblem,
        alg::AbstractHCAlgorithm,
        ::Type{ResultIterator};
        show_progress::Bool = true,
        bitmask = nothing,
    )
    startsolutions = _solver_startsolutions(prob, alg; show_progress)
    return _iter_cache_from_startsolutions(startsolutions, prob, alg; bitmask)
end

solve(prob::AbstractHCProblem, alg::AbstractHCAlgorithm; kwargs...) = solve!(
    init(prob, alg; kwargs...)
)

solve!(cache::PathSolveCache) = solve(
    cache.solver, cache.starts; _cache_solve_kwargs(cache)...
)

solve!(cache::SweepSolveCache) = solve(
    cache.solver, cache.starts, cache.targets; _cache_solve_kwargs(cache)...
)

solve!(cache::PathIteratorSolveCache) = solve(
    cache.solver, cache.starts, ResultIterator; bitmask = cache.bitmask
)

function solve!(cache::SweepIteratorSolveCache)
    return map(cache.solvers) do solverᵢ
        solve(solverᵢ, cache.starts, ResultIterator; bitmask = cache.bitmask)
    end
end

solve(prob::AbstractHCProblem, alg::AbstractHCAlgorithm, ::Type{ResultIterator}; kwargs...) = solve!(
    init(prob, alg, ResultIterator; kwargs...)
)
