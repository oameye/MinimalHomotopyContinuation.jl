const PathProblems = Union{SystemProblem, ParameterHomotopyProblem, HomotopyProblem}

function _cache_from_startsolutions(
        startsolutions::Tuple;
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
        show_progress::Bool,
        threading::Bool,
        catch_interrupt::Bool,
    )
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

_iter_cache_from_startsolutions(startsolutions::Tuple; bitmask = nothing) = begin
    solver, starts = startsolutions
    PathIteratorSolveCache(solver, starts, bitmask)
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

function init(
        prob::PathProblems,
        alg::AbstractHCAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    startsolutions = _solver_startsolutions(prob, alg; show_progress)
    return _cache_from_startsolutions(
        startsolutions; stop_early_cb, show_progress, threading, catch_interrupt
    )
end

function init(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    _ = stop_early_cb
    startsolutions = _solver_startsolutions(prob, alg; show_progress)
    return _cache_from_startsolutions(
        startsolutions, prob; show_progress, threading, catch_interrupt
    )
end

function init(
        prob::PathProblems,
        alg::AbstractHCAlgorithm,
        ::Type{ResultIterator};
        show_progress::Bool = true,
        bitmask = nothing,
    )
    startsolutions = _solver_startsolutions(prob, alg; show_progress)
    return _iter_cache_from_startsolutions(startsolutions; bitmask)
end

function init(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm,
        ::Type{ResultIterator};
        show_progress::Bool = true,
        bitmask = nothing,
    )
    startsolutions = _solver_startsolutions(prob, alg; show_progress)
    return _iter_cache_from_startsolutions(startsolutions, prob, alg; bitmask)
end

function solve(
        prob::AbstractHCProblem,
        alg::AbstractHCAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    return solve!(init(prob, alg; stop_early_cb, show_progress, threading, catch_interrupt))
end

function solve!(cache::PathSolveCache)
    return solve(
        cache.solver,
        cache.starts;
        stop_early_cb = cache.stop_early_cb,
        show_progress = cache.show_progress,
        threading = cache.threading,
        catch_interrupt = cache.catch_interrupt,
    )
end

function solve!(cache::SweepSolveCache)
    return solve(
        cache.solver,
        cache.starts,
        cache.targets;
        reducer = cache.reducer,
        show_progress = cache.show_progress,
        threading = cache.threading,
        catch_interrupt = cache.catch_interrupt,
    )
end

solve!(cache::PathIteratorSolveCache) = solve(
    cache.solver, cache.starts, ResultIterator; bitmask = cache.bitmask
)

function solve!(cache::SweepIteratorSolveCache)
    return map(cache.solvers) do solverᵢ
        solve(solverᵢ, cache.starts, ResultIterator; bitmask = cache.bitmask)
    end
end

function solve(
        prob::PathProblems,
        alg::AbstractHCAlgorithm,
        ::Type{ResultIterator};
        show_progress::Bool = true,
        bitmask = nothing,
    )
    return solve!(init(prob, alg, ResultIterator; show_progress, bitmask))
end

function solve(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm,
        ::Type{ResultIterator};
        show_progress::Bool = true,
        bitmask = nothing,
    )
    return solve!(init(prob, alg, ResultIterator; show_progress, bitmask))
end
