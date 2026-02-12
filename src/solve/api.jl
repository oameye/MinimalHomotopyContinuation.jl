function init(
        prob::SystemProblem,
        alg::PolyhedralAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    solver, starts = _system_solver_startsolutions(prob, alg; show_progress = show_progress)
    return PathSolveCache(
        solver,
        starts,
        stop_early_cb,
        show_progress,
        threading,
        catch_interrupt,
    )
end

function init(
        prob::SystemProblem,
        alg::PolyhedralAlgorithm,
        ::Type{ResultIterator};
        show_progress::Bool = true,
        bitmask = nothing,
    )
    solver, starts = _system_solver_startsolutions(prob, alg; show_progress = show_progress)
    return PathIteratorSolveCache(solver, starts, bitmask)
end

function init(
        prob::SystemProblem,
        alg::TotalDegreeAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    solver, starts = _system_solver_startsolutions(prob, alg)
    return PathSolveCache(
        solver,
        starts,
        stop_early_cb,
        show_progress,
        threading,
        catch_interrupt,
    )
end

function init(
        prob::SystemProblem,
        alg::TotalDegreeAlgorithm,
        ::Type{ResultIterator};
        show_progress::Bool = true,
        bitmask = nothing,
    )
    solver, starts = _system_solver_startsolutions(prob, alg)
    return PathIteratorSolveCache(solver, starts, bitmask)
end

function init(
        prob::ParameterHomotopyProblem,
        alg::PathTrackingAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    solver, starts = _parameter_solver_startsolutions(prob, alg)
    return PathSolveCache(
        solver,
        starts,
        stop_early_cb,
        show_progress,
        threading,
        catch_interrupt,
    )
end

function init(
        prob::ParameterHomotopyProblem,
        alg::PathTrackingAlgorithm,
        ::Type{ResultIterator};
        bitmask = nothing,
    )
    solver, starts = _parameter_solver_startsolutions(prob, alg)
    return PathIteratorSolveCache(solver, starts, bitmask)
end

function init(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm;
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    )
    solver, starts = _sweep_solver_startsolutions(prob, alg)
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

function init(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm,
        ::Type{ResultIterator};
        bitmask = nothing,
    )
    _, starts = _sweep_solver_startsolutions(prob, alg)
    solvers = map(p -> _sweep_solver(prob, alg, p), prob.targets)
    return SweepIteratorSolveCache(solvers, starts, bitmask)
end

function init(
        prob::HomotopyProblem,
        alg::PathTrackingAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    solver, starts = _homotopy_solver_startsolutions(prob, alg)
    return PathSolveCache(
        solver,
        starts,
        stop_early_cb,
        show_progress,
        threading,
        catch_interrupt,
    )
end

function init(
        prob::HomotopyProblem,
        alg::PathTrackingAlgorithm,
        ::Type{ResultIterator};
        bitmask = nothing,
    )
    solver, starts = _homotopy_solver_startsolutions(prob, alg)
    return PathIteratorSolveCache(solver, starts, bitmask)
end

function solve(
        prob::SystemProblem,
        alg::Union{PolyhedralAlgorithm, TotalDegreeAlgorithm};
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    return solve!(
        init(
            prob,
            alg;
            stop_early_cb = stop_early_cb,
            show_progress = show_progress,
            threading = threading,
            catch_interrupt = catch_interrupt,
        ),
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
    return solve!(
        init(
            prob,
            alg;
            stop_early_cb = stop_early_cb,
            show_progress = show_progress,
            threading = threading,
            catch_interrupt = catch_interrupt,
        ),
    )
end

function solve(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm;
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    )
    return solve!(
        init(
            prob,
            alg;
            show_progress = show_progress,
            threading = threading,
            catch_interrupt = catch_interrupt,
        ),
    )
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

solve!(cache::PathIteratorSolveCache) =
    solve(cache.solver, cache.starts, ResultIterator; bitmask = cache.bitmask)

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
    return solve!(
        init(prob, alg, ResultIterator; show_progress = show_progress, bitmask = bitmask),
    )
end

function solve(
        prob::ParameterHomotopyProblem,
        alg::PathTrackingAlgorithm,
        ::Type{ResultIterator};
        bitmask = nothing,
    )
    return solve!(init(prob, alg, ResultIterator; bitmask = bitmask))
end

function solve(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm,
        ::Type{ResultIterator};
        bitmask = nothing,
    )
    return solve!(init(prob, alg, ResultIterator; bitmask = bitmask))
end

function solve(
        prob::HomotopyProblem,
        alg::PathTrackingAlgorithm,
        ::Type{ResultIterator};
        bitmask = nothing,
    )
    return solve!(init(prob, alg, ResultIterator; bitmask = bitmask))
end
