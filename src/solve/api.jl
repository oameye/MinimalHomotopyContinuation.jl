function _path_cache(
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

function _path_iter_cache(startsolutions::Tuple; bitmask = nothing)
    solver, starts = startsolutions
    return PathIteratorSolveCache(solver, starts, bitmask)
end

function _sweep_cache(
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

function _sweep_iter_cache(
        prob::ParameterSweepProblem, alg::PathTrackingAlgorithm; bitmask = nothing
    )
    _, starts = _sweep_solver_startsolutions(prob, alg)
    solvers = map(p -> _sweep_solver(prob, alg, p), prob.targets)
    return SweepIteratorSolveCache(solvers, starts, bitmask)
end

function _system_startsolutions(
        prob::SystemProblem, alg::PolyhedralAlgorithm; show_progress::Bool = true
    )
    return _system_solver_startsolutions(prob, alg; show_progress)
end

function _system_startsolutions(
        prob::SystemProblem, alg::TotalDegreeAlgorithm; show_progress::Bool = true
    )
    return _system_solver_startsolutions(prob, alg)
end

_path_startsolutions(prob::ParameterHomotopyProblem, alg::PathTrackingAlgorithm) = _parameter_solver_startsolutions(
    prob, alg
)
_path_startsolutions(prob::HomotopyProblem, alg::PathTrackingAlgorithm) = _homotopy_solver_startsolutions(
    prob, alg
)

function _solve_from_init(prob, alg; kwargs...)
    return solve!(init(prob, alg; kwargs...))
end

function _solve_from_iter_init(prob, alg; kwargs...)
    return solve!(init(prob, alg, ResultIterator; kwargs...))
end

function init(
        prob::SystemProblem,
        alg::Union{PolyhedralAlgorithm, TotalDegreeAlgorithm};
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    sol = _system_startsolutions(prob, alg; show_progress)
    return _path_cache(sol; stop_early_cb, show_progress, threading, catch_interrupt)
end

function init(
        prob::SystemProblem,
        alg::Union{PolyhedralAlgorithm, TotalDegreeAlgorithm},
        ::Type{ResultIterator};
        show_progress::Bool = true,
        bitmask = nothing,
    )
    return _path_iter_cache(_system_startsolutions(prob, alg; show_progress); bitmask)
end

function init(
        prob::Union{ParameterHomotopyProblem, HomotopyProblem},
        alg::PathTrackingAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    sol = _path_startsolutions(prob, alg)
    return _path_cache(sol; stop_early_cb, show_progress, threading, catch_interrupt)
end

function init(
        prob::Union{ParameterHomotopyProblem, HomotopyProblem},
        alg::PathTrackingAlgorithm,
        ::Type{ResultIterator};
        bitmask = nothing,
    )
    return _path_iter_cache(_path_startsolutions(prob, alg); bitmask)
end

function init(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm;
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    )
    sol = _sweep_solver_startsolutions(prob, alg)
    return _sweep_cache(sol, prob; show_progress, threading, catch_interrupt)
end

function init(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm,
        ::Type{ResultIterator};
        bitmask = nothing,
    )
    return _sweep_iter_cache(prob, alg; bitmask)
end

function solve(
        prob::SystemProblem,
        alg::Union{PolyhedralAlgorithm, TotalDegreeAlgorithm};
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    return _solve_from_init(
        prob, alg; stop_early_cb, show_progress, threading, catch_interrupt
    )
end

function solve(
        prob::Union{ParameterHomotopyProblem, HomotopyProblem},
        alg::PathTrackingAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    return _solve_from_init(
        prob, alg; stop_early_cb, show_progress, threading, catch_interrupt
    )
end

function solve(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm;
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    )
    return _solve_from_init(prob, alg; show_progress, threading, catch_interrupt)
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
        prob::SystemProblem,
        alg::Union{PolyhedralAlgorithm, TotalDegreeAlgorithm},
        ::Type{ResultIterator};
        show_progress::Bool = true,
        bitmask = nothing,
    )
    return _solve_from_iter_init(prob, alg; show_progress, bitmask)
end

function solve(
        prob::Union{ParameterHomotopyProblem, ParameterSweepProblem, HomotopyProblem},
        alg::PathTrackingAlgorithm,
        ::Type{ResultIterator};
        bitmask = nothing,
    )
    return _solve_from_iter_init(prob, alg; bitmask)
end
