export solve, Solver, paths_to_track, parameter_homotopy

struct SolveStats
    regular::Threads.Atomic{Int}
    regular_real::Threads.Atomic{Int}
    singular::Threads.Atomic{Int}
    singular_real::Threads.Atomic{Int}
end
SolveStats() = SolveStats(Threads.Atomic{Int}.((0, 0, 0, 0))...)

function init!(SS::SolveStats)
    SS.regular[] = SS.regular_real[] = SS.singular[] = SS.singular_real[] = 0
    return SS
end

function update!(stats::SolveStats, R::PathResult)
    is_success(R) || return stats

    if is_singular(R)
        Threads.atomic_add!(stats.singular_real, Int(is_real(R)))
        Threads.atomic_add!(stats.singular, 1)
    else
        Threads.atomic_add!(stats.regular_real, Int(is_real(R)))
        Threads.atomic_add!(stats.regular, 1)
    end
    return stats
end

struct Solver{T <: AbstractPathTracker} <: AbstractSolver
    trackers::Vector{T}
    seed::Union{Nothing, UInt32}
    stats::SolveStats
    start_system::Union{Nothing, Symbol}
end
Solver(
    tracker::AbstractPathTracker;
    seed::Union{Nothing, UInt32} = nothing,
    start_system = nothing,
) = Solver([tracker], seed, SolveStats(), start_system)

Base.show(io::IO, solver::Solver) = print(io, typeof(solver), "()")

# Internal helper kept for transition support.
function solver_startsolutions(
        F::AbstractVector{Expression},
        starts = nothing;
        parameters = Variable[],
        variables = setdiff(variables(F), parameters),
        variable_ordering = variables,
        kwargs...,
    )
    sys = System(
        F,
        variables = variable_ordering,
        parameters = parameters,
    )
    return solver_startsolutions(sys, starts; kwargs...)
end
function solver_startsolutions(
        F::AbstractVector{<:MP.AbstractPolynomial},
        starts = nothing;
        parameters = similar(MP.variables(F), 0),
        variables = setdiff(MP.variables(F), parameters),
        variable_ordering = variables,
        target_parameters = nothing,
        kwargs...,
    )
    if isnothing(target_parameters) && isempty(parameters)
        sys, target_parameters = ModelKit.system_with_coefficents_as_params(
            F,
            variables = variable_ordering,
        )
    else
        sys = System(
            F,
            variables = variable_ordering,
            parameters = parameters,
        )
    end
    return solver_startsolutions(sys, starts; target_parameters = target_parameters, kwargs...)
end
function solver_startsolutions(
        F::Union{System, AbstractSystem},
        starts = nothing;
        seed = rand(UInt32),
        start_system = :polyhedral,
        generic_parameters = nothing,
        p₁ = generic_parameters,
        start_parameters = p₁,
        p₀ = generic_parameters,
        target_parameters = p₀,
        compile::Union{Bool, Symbol} = COMPILE_DEFAULT[],
        kwargs...,
    )
    !isnothing(seed) && Random.seed!(seed)

    used_start_system = nothing
    if start_parameters !== nothing
        tracker = parameter_homotopy_tracker(
            F;
            start_parameters = start_parameters,
            target_parameters = target_parameters,
            compile = compile,
            kwargs...,
        )
    elseif start_system == :polyhedral
        used_start_system = :polyhedral
        tracker, starts = polyhedral(
            F;
            compile = compile,
            target_parameters = target_parameters,
            kwargs...,
        )
    elseif start_system == :total_degree
        used_start_system = :total_degree
        tracker, starts = total_degree(
            F;
            compile = compile,
            target_parameters = target_parameters,
            kwargs...,
        )
        try
            first(starts)
        catch
            throw("The number of start solutions is zero (total degree is zero).")
        end
    else
        throw(
            KeywordArgumentException(
                :start_system,
                start_system,
                "Possible values are: `:polyhedral` and `:total_degree`.",
            ),
        )
    end

    return Solver(tracker; seed = seed, start_system = used_start_system), starts
end

function solver_startsolutions(
        H::Union{Homotopy, AbstractHomotopy},
        starts = nothing;
        compile::Union{Bool, Symbol} = COMPILE_DEFAULT[],
        seed = nothing,
        kwargs...,
    )
    unsupported_kwargs(kwargs)
    !isnothing(seed) && Random.seed!(seed)
    return Solver(EndgameTracker(fixed(H; compile = compile)); seed = seed), starts
end

function solver_startsolutions(args...; kwargs...)
    throw(MethodError(solver_startsolutions, args))
end

function parameter_homotopy_tracker(
        F::Union{System, AbstractSystem};
        tracker_options = TrackerOptions(),
        endgame_options = EndgameOptions(),
        kwargs...,
    )
    H = parameter_homotopy(F; kwargs...)
    return EndgameTracker(H; tracker_options = tracker_options, options = endgame_options)
end

function parameter_homotopy(
        F::Union{System, AbstractSystem};
        generic_parameters = randn(ComplexF64, nparameters(F)),
        p₁ = generic_parameters,
        start_parameters = p₁,
        p₀ = generic_parameters,
        target_parameters = p₀,
        compile::Union{Bool, Symbol} = COMPILE_DEFAULT[],
        kwargs...,
    )
    unsupported_kwargs(kwargs)
    m, n = size(F)
    H = ParameterHomotopy(fixed(F; compile = compile), start_parameters, target_parameters)
    f = System(F)
    if is_homogeneous(f)
        throw(
            ArgumentError(
                "Homogeneous/projective systems are not supported in affine-only mode.",
            ),
        )
    end
    if m < n
        throw(FiniteException(n - m))
    elseif m > n
        throw(
            ArgumentError(
                "Only square systems are supported in this minimal build. Got $m equations, expected $n.",
            ),
        )
    end

    return H
end

struct PathSolveCache{SolverT, StartsT, StopEarlyT, BitmaskT}
    solver::SolverT
    starts::StartsT
    stop_early_cb::StopEarlyT
    show_progress::Bool
    threading::Bool
    catch_interrupt::Bool
    iterator_only::Bool
    bitmask::BitmaskT
end

struct SweepSolveCache{
        SolverT,
        StartsT,
        TargetsT,
        TransformResultT,
        TransformParamsT,
        BitmaskT,
    }
    solver::SolverT
    starts::StartsT
    targets::TargetsT
    transform_result::TransformResultT
    transform_parameters::TransformParamsT
    flatten::Bool
    show_progress::Bool
    threading::Bool
    catch_interrupt::Bool
    iterator_only::Bool
    bitmask::BitmaskT
end

function _seed!(seed)
    !isnothing(seed) && Random.seed!(seed)
    return nothing
end

function _system_solver_startsolutions(
        prob::SystemProblem,
        alg::PolyhedralAlgorithm;
        show_progress::Bool = true,
    )
    _seed!(alg.seed)
    only_non_zero = something(alg.only_non_zero, alg.only_torus)
    tracker, starts = polyhedral(
        prob.system;
        compile = alg.compile,
        target_parameters = prob.target_parameters,
        tracker_options = alg.tracker_options,
        endgame_options = alg.endgame_options,
        only_torus = alg.only_torus,
        only_non_zero = only_non_zero,
        show_progress = show_progress,
    )
    return Solver(tracker; seed = alg.seed, start_system = :polyhedral), starts
end

function _system_solver_startsolutions(prob::SystemProblem, alg::TotalDegreeAlgorithm)
    _seed!(alg.seed)
    tracker, starts = total_degree(
        prob.system;
        compile = alg.compile,
        target_parameters = prob.target_parameters,
        gamma = alg.gamma,
        tracker_options = alg.tracker_options,
        endgame_options = alg.endgame_options,
    )
    try
        first(starts)
    catch
        throw("The number of start solutions is zero (total degree is zero).")
    end
    return Solver(tracker; seed = alg.seed, start_system = :total_degree), starts
end

function _parameter_solver_startsolutions(
        prob::ParameterHomotopyProblem;
        compile::Union{Bool, Symbol} = COMPILE_DEFAULT[],
        seed::Union{Nothing, UInt32} = rand(UInt32),
        tracker_options = TrackerOptions(),
        endgame_options = EndgameOptions(),
    )
    _seed!(seed)
    tracker = parameter_homotopy_tracker(
        prob.system;
        start_parameters = prob.start_parameters,
        target_parameters = prob.target_parameters,
        compile = compile,
        tracker_options = tracker_options,
        endgame_options = endgame_options,
    )
    return Solver(tracker; seed = seed, start_system = nothing), prob.start_solutions
end

function _sweep_solver_startsolutions(
        prob::ParameterSweepProblem;
        compile::Union{Bool, Symbol} = COMPILE_DEFAULT[],
        seed::Union{Nothing, UInt32} = rand(UInt32),
        tracker_options = TrackerOptions(),
        endgame_options = EndgameOptions(),
    )
    _seed!(seed)
    isempty(prob.targets) && throw(ArgumentError("`targets` must be non-empty."))
    first_target = first(prob.targets)
    first_params = prob.transform_parameters(first_target)
    tracker = parameter_homotopy_tracker(
        prob.system;
        start_parameters = prob.start_parameters,
        target_parameters = first_params,
        compile = compile,
        tracker_options = tracker_options,
        endgame_options = endgame_options,
    )
    return Solver(tracker; seed = seed, start_system = nothing), prob.start_solutions
end

function _homotopy_solver_startsolutions(
        prob::HomotopyProblem;
        compile::Union{Bool, Symbol} = COMPILE_DEFAULT[],
        seed::Union{Nothing, UInt32} = rand(UInt32),
        tracker_options = TrackerOptions(),
        endgame_options = EndgameOptions(),
    )
    _seed!(seed)
    tracker = EndgameTracker(
        fixed(prob.homotopy; compile = compile);
        tracker_options = tracker_options,
        options = endgame_options,
    )
    return Solver(tracker; seed = seed, start_system = nothing), prob.start_solutions
end

function init(
        prob::SystemProblem,
        alg::PolyhedralAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
        iterator_only::Bool = false,
        bitmask = nothing,
    ) where {F}
    solver, starts = _system_solver_startsolutions(prob, alg; show_progress = show_progress)
    return PathSolveCache(
        solver,
        starts,
        stop_early_cb,
        show_progress,
        threading,
        catch_interrupt,
        iterator_only,
        bitmask,
    )
end

function init(
        prob::SystemProblem,
        alg::TotalDegreeAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
        iterator_only::Bool = false,
        bitmask = nothing,
    ) where {F}
    solver, starts = _system_solver_startsolutions(prob, alg)
    return PathSolveCache(
        solver,
        starts,
        stop_early_cb,
        show_progress,
        threading,
        catch_interrupt,
        iterator_only,
        bitmask,
    )
end

function init(
        prob::ParameterHomotopyProblem,
        ;
        compile::Union{Bool, Symbol} = COMPILE_DEFAULT[],
        seed::Union{Nothing, UInt32} = rand(UInt32),
        tracker_options = TrackerOptions(),
        endgame_options = EndgameOptions(),
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
        iterator_only::Bool = false,
        bitmask = nothing,
    ) where {F}
    solver, starts = _parameter_solver_startsolutions(
        prob;
        compile = compile,
        seed = seed,
        tracker_options = tracker_options,
        endgame_options = endgame_options,
    )
    return PathSolveCache(
        solver,
        starts,
        stop_early_cb,
        show_progress,
        threading,
        catch_interrupt,
        iterator_only,
        bitmask,
    )
end

function init(
        prob::ParameterSweepProblem,
        ;
        compile::Union{Bool, Symbol} = COMPILE_DEFAULT[],
        seed::Union{Nothing, UInt32} = rand(UInt32),
        tracker_options = TrackerOptions(),
        endgame_options = EndgameOptions(),
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
        iterator_only::Bool = false,
        bitmask = nothing,
    )
    solver, starts = _sweep_solver_startsolutions(
        prob;
        compile = compile,
        seed = seed,
        tracker_options = tracker_options,
        endgame_options = endgame_options,
    )
    return SweepSolveCache(
        solver,
        starts,
        prob.targets,
        prob.transform_result,
        prob.transform_parameters,
        prob.flatten,
        show_progress,
        threading,
        catch_interrupt,
        iterator_only,
        bitmask,
    )
end

function init(
        prob::HomotopyProblem,
        ;
        compile::Union{Bool, Symbol} = COMPILE_DEFAULT[],
        seed::Union{Nothing, UInt32} = rand(UInt32),
        tracker_options = TrackerOptions(),
        endgame_options = EndgameOptions(),
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
        iterator_only::Bool = false,
        bitmask = nothing,
    ) where {F}
    solver, starts = _homotopy_solver_startsolutions(
        prob;
        compile = compile,
        seed = seed,
        tracker_options = tracker_options,
        endgame_options = endgame_options,
    )
    return PathSolveCache(
        solver,
        starts,
        stop_early_cb,
        show_progress,
        threading,
        catch_interrupt,
        iterator_only,
        bitmask,
    )
end

solve(prob::SystemProblem; kwargs...) = throw(
    ArgumentError(
        "A `SystemProblem` requires an explicit algorithm. Use `solve(prob, PolyhedralAlgorithm())` or `solve(prob, TotalDegreeAlgorithm())`.",
    ),
)
function solve(
        prob::SystemProblem,
        alg::PolyhedralAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
        iterator_only::Bool = false,
        bitmask = nothing,
    ) where {F}
    return solve!(
        init(
            prob,
            alg;
            stop_early_cb = stop_early_cb,
            show_progress = show_progress,
            threading = threading,
            catch_interrupt = catch_interrupt,
            iterator_only = iterator_only,
            bitmask = bitmask,
        ),
    )
end

function solve(
        prob::SystemProblem,
        alg::TotalDegreeAlgorithm;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
        iterator_only::Bool = false,
        bitmask = nothing,
    ) where {F}
    return solve!(
        init(
            prob,
            alg;
            stop_early_cb = stop_early_cb,
            show_progress = show_progress,
            threading = threading,
            catch_interrupt = catch_interrupt,
            iterator_only = iterator_only,
            bitmask = bitmask,
        ),
    )
end
solve(prob::ParameterHomotopyProblem; kwargs...) = solve!(init(prob; kwargs...))
solve(prob::ParameterSweepProblem; kwargs...) = solve!(init(prob; kwargs...))
solve(prob::HomotopyProblem; kwargs...) = solve!(init(prob; kwargs...))

function solve!(cache::PathSolveCache)
    if cache.iterator_only
        return ResultIterator(cache.starts, cache.solver; bitmask = cache.bitmask)
    end
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
    if cache.iterator_only
        return map(cache.targets) do p
            solverᵢ = deepcopy(cache.solver)
            target_parameters!(solverᵢ, cache.transform_parameters(p))
            ResultIterator(cache.starts, solverᵢ; bitmask = cache.bitmask)
        end
    end
    return solve(
        cache.solver,
        cache.starts,
        cache.targets;
        show_progress = cache.show_progress,
        threading = cache.threading,
        catch_interrupt = cache.catch_interrupt,
        transform_result = cache.transform_result,
        transform_parameters = cache.transform_parameters,
        flatten = cache.flatten,
    )
end

solve(S::Solver, R::Result; kwargs...) =
    solve(S, solutions(R; only_nonsingular = true); kwargs...)
solve(S::Solver, s::AbstractVector{<:Number}; kwargs...) = solve(S, [s]; kwargs...)

function solve(
        S::Solver,
        starts;
        iterator_only::Bool = false,
        bitmask = nothing,
        threading::Bool = Threads.nthreads() > 1,
        kwargs...,
    )
    if iterator_only
        return ResultIterator(starts, S; bitmask = bitmask)
    else
        return solve(S, collect(starts); threading = threading, kwargs...)
    end
end

function solve(
        S::Solver,
        starts::AbstractArray;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}

    n = length(starts)
    progress = show_progress ? make_progress(n; delay = 0.3) : nothing
    init!(S.stats)
    return if threading
        threaded_solve(
            S,
            starts,
            progress,
            stop_early_cb;
            catch_interrupt = catch_interrupt,
        )
    else
        serial_solve(S, starts, progress, stop_early_cb; catch_interrupt = catch_interrupt)
    end
end

(solver::Solver)(starts; kwargs...) = solve(solver, starts; kwargs...)
track(solver::Solver, s; kwargs...) = track(solver.trackers[1], s; kwargs...)

function make_progress(n::Integer; delay::Float64 = 0.0)
    desc = "Tracking $n paths... "
    barlen = min(ProgressMeter.tty_width(desc, stdout, false), 40)
    progress =
        ProgressMeter.Progress(n; dt = 0.2, desc = desc, barlen = barlen, output = stdout)
    progress.tlast += delay
    return progress
end

function update_progress!(progress, stats, ntracked)
    t = time()
    if ntracked == progress.n || t > progress.tlast + progress.dt
        showvalues = make_showvalues(stats, ntracked)
        ProgressMeter.update!(progress, ntracked; showvalues = showvalues)
    end
    return nothing
end

@noinline function make_showvalues(stats, ntracked)
    nsols = stats.regular[] + stats.singular[]
    nreal = stats.regular_real[] + stats.singular_real[]
    return (
        ("# paths tracked", ntracked),
        ("# non-singular solutions (real)", "$(stats.regular[]) ($(stats.regular_real[]))"),
        ("# singular endpoints (real)", "$(stats.singular[]) ($(stats.singular_real[]))"),
        ("# total solutions (real)", "$nsols ($nreal)"),
    )
end

update_progress!(::Nothing, stats, ntracked) = nothing

function serial_solve(
        solver::Solver,
        starts,
        progress = nothing,
        stop_early_cb::F = always_false;
        catch_interrupt::Bool = true,
    ) where {F}
    path_results = Vector{PathResult}()
    tracker = solver.trackers[1]
    try
        for (k, s) in enumerate(starts)
            r = track(tracker, s; path_number = k)
            push!(path_results, r)
            update!(solver.stats, r)
            update_progress!(progress, solver.stats, k)
            if is_success(r) && stop_early_cb(r)
                break
            end
        end
    catch e
        (catch_interrupt && isa(e, InterruptException)) || rethrow(e)
    end

    return Result(path_results; seed = solver.seed, start_system = solver.start_system)
end

function threaded_solve(
        solver::Solver,
        S::AbstractArray,
        progress = nothing,
        stop_early_cb::F = always_false;
        catch_interrupt::Bool = true,
    ) where {F}

    N = length(S)
    path_results = Vector{PathResult}(undef, N)
    interrupted = Threads.Atomic{Bool}(false)
    finished = Threads.Atomic{Int}(0)
    next_k = Threads.Atomic{Int}(1)

    tracker = solver.trackers[1]
    ntrackers = length(solver.trackers)
    nthr = Threads.nthreads()
    resize!(solver.trackers, nthr)
    for i in (ntrackers + 1):nthr
        solver.trackers[i] = deepcopy(tracker)
    end

    progress_lock = ReentrantLock()

    try
        Base.@sync begin
            for lt in solver.trackers
                let tracker = lt
                    Threads.@spawn begin
                        while true
                            idx = Threads.atomic_add!(next_k, 1)
                            k = idx
                            if k > N || interrupted[]
                                break
                            end
                            r = track(tracker, S[k]; path_number = k)
                            path_results[k] = r
                            nfinished = Threads.atomic_add!(finished, 1) + 1
                            lock(progress_lock) do
                                update!(solver.stats, r)
                                update_progress!(progress, solver.stats, nfinished)
                            end
                            if is_success(r) && stop_early_cb(r)
                                interrupted[] = true
                            end
                        end
                    end
                end
            end
        end
    catch e
        if (
                isa(e, InterruptException) ||
                    (isa(e, TaskFailedException) && isa(e.task.exception, InterruptException))
            )
            interrupted[] = true
        end
        if !interrupted[] || !catch_interrupt
            rethrow(e)
        end
    end

    return if interrupted[]
        assigned_results = Vector{PathResult}()
        for i in eachindex(path_results)
            if isassigned(path_results, i)
                push!(assigned_results, path_results[i])
            end
        end
        Result(assigned_results; seed = solver.seed, start_system = solver.start_system)
    else
        Result(path_results; seed = solver.seed, start_system = solver.start_system)
    end

end

function start_parameters!(solver::Solver, p)
    for tracker in solver.trackers
        start_parameters!(tracker, p)
    end
    return solver
end

function target_parameters!(solver::Solver, p)
    for tracker in solver.trackers
        target_parameters!(tracker, p)
    end
    return solver
end

function parameters!(solver::Solver, p, q)
    for tracker in solver.trackers
        parameters!(tracker, p, q)
    end
    return solver
end

function paths_to_track(
        f::Union{System, AbstractSystem};
        start_system::Symbol = :polyhedral,
        kwargs...,
    )
    if start_system == :polyhedral
        f_sys = f isa System ? f : System(f)
        only_torus = get(kwargs, :only_torus, false)
        only_torus isa Bool || throw(KeywordArgumentException(:only_torus, only_torus))
        only_non_zero = get(kwargs, :only_non_zero, only_torus)
        only_non_zero isa Bool ||
            throw(KeywordArgumentException(:only_non_zero, only_non_zero))
        extra = Dict{Symbol, Any}()
        for (k, v) in pairs(kwargs)
            if k ∉ (:only_torus, :only_non_zero)
                extra[k] = v
            end
        end
        unsupported_kwargs(extra)
        return paths_to_track(
            f_sys,
            Val(:polyhedral);
            only_torus = only_torus,
            only_non_zero = only_non_zero,
        )
    elseif start_system == :total_degree
        unsupported_kwargs(kwargs)
        return paths_to_track(f, Val(:total_degree))
    else
        throw(
            KeywordArgumentException(
                :start_system,
                start_system,
                "Possible values are: `:polyhedral` and `:total_degree`.",
            ),
        )
    end
end

function solve(
        S::Solver,
        starts,
        target_parameters;
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
        transform_result = nothing,
        transform_parameters = nothing,
        flatten = nothing,
    )
    transform_result = something(transform_result, tuple)
    transform_parameters = something(transform_parameters, identity)
    flatten = something(flatten, false)

    n = length(target_parameters)

    progress = show_progress ? make_many_progress(n; delay = 0.3) : nothing

    return many_solve(
        S,
        starts,
        target_parameters,
        progress,
        transform_result,
        transform_parameters,
        Val(flatten);
        catch_interrupt = catch_interrupt,
        threading = threading,
    )
end

function make_many_progress(n::Integer; delay::Float64 = 0.0)
    desc = "Solving for $n parameters... "
    barlen = min(ProgressMeter.tty_width(desc, stdout, false), 40)
    progress =
        ProgressMeter.Progress(n; dt = 0.3, desc = desc, barlen = barlen, output = stdout)
    progress.tlast += delay
    return progress
end

function update_many_progress!(progress, results, k, paths_per_param; flatten::Bool)
    t = time()
    if k == progress.n || t > progress.tlast + progress.dt
        showvalues = make_many_showvalues(results, k, paths_per_param; flatten = flatten)
        ProgressMeter.update!(progress, k; showvalues = showvalues)
    end
    return nothing
end

@noinline function make_many_showvalues(results, k, paths_per_param; flatten::Bool)
    return if flatten
        [
            ("# parameters solved", k),
            ("# paths tracked", paths_per_param * k),
            ("# results", length(results)),
        ]
    else
        [("# parameters solved", k), ("# paths tracked", paths_per_param * k)]
    end
end

update_many_progress!(::Nothing, results, k, paths_per_param; kwargs...) = nothing

function many_solve(
        solver::Solver,
        starts,
        many_target_parameters,
        progress,
        transform_result,
        transform_parameters,
        ::Val{flatten};
        threading::Bool,
        catch_interrupt::Bool,
    ) where {flatten}


    if isa(starts, TotalDegreeStartSolutionsIterator)
        @warn "Solving for many parameters with total degree start system is not recommended. Instead, one should use a two-step approach: first solve a system with generic parameters, and track its solutions to the desired parameters. See https://www.juliahomotopycontinuation.org/guides/many-systems/."
    elseif isa(starts, PolyhedralStartSolutionsIterator)
        @error "Solving for many parameters with polyhedral start system is not implemented and also not recommended. Instead, one should use a two-step approach: solve a system with generic parameters, and track its solutions to the desired parameters. See https://www.juliahomotopycontinuation.org/guides/many-systems/."
    end
    q = first(many_target_parameters)
    target_parameters!(solver, transform_parameters(q))
    if threading
        res = threaded_solve(solver, collect(starts); catch_interrupt = false)
    else
        res = serial_solve(solver, starts; catch_interrupt = false)
    end
    if flatten
        results = transform_result(res, q)
        if !(results isa AbstractArray)
            throw(ArgumentError("Cannot flatten arguments of type `$(typeof(results))`"))
        end
    else
        results = [transform_result(res, q)]
    end
    k = 1
    m = length(starts)
    update_many_progress!(progress, results, k, m; flatten = flatten)


    try
        for q in Iterators.drop(many_target_parameters, 1)
            target_parameters!(solver, transform_parameters(q))
            if threading
                res = threaded_solve(solver, collect(starts); catch_interrupt = false)
            else
                res = serial_solve(solver, starts; catch_interrupt = false)
            end

            if flatten
                append!(results, transform_result(res, q))
            else
                push!(results, transform_result(res, q))
            end
            k += 1
            update_many_progress!(progress, results, k, m; flatten = flatten)
        end
    catch e
        if !(
                isa(e, InterruptException) ||
                    (isa(e, TaskFailedException) && isa(e.task.exception, InterruptException))
            )
            rethrow(e)
        end
    end

    return results
end
