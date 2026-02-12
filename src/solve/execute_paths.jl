_normalize_starts(R::Result) = solutions(R; only_nonsingular = true)
_normalize_starts(s::AbstractVector{<:Number}) = [s]
_normalize_starts(starts::AbstractArray) = starts
_normalize_starts(starts) = collect(starts)

solve(S::Solver, starts, ::Type{ResultIterator}; bitmask = nothing) = ResultIterator(
    starts, S; bitmask
)

function solve(
        S::Solver,
        starts;
        stop_early_cb::F = always_false,
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    ) where {F}
    starts_buffer = _normalize_starts(starts)
    n = length(starts_buffer)
    progress = show_progress ? make_progress(n; delay = 0.3) : nothing
    init!(S.stats)
    return _run_paths(S, starts_buffer, progress, stop_early_cb; threading, catch_interrupt)
end

track(solver::Solver, s; path_number::Union{Nothing, Int} = nothing) =
if isnothing(path_number)
    track(solver.trackers[1], s)
else
    track(solver.trackers[1], s; path_number)
end

function make_progress(n::Integer; delay::Float64 = 0.0)
    desc = "Tracking $n paths... "
    barlen = min(ProgressMeter.tty_width(desc, stdout, false), 40)
    progress = ProgressMeter.Progress(n; dt = 0.2, desc, barlen, output = stdout)
    progress.tlast += delay
    return progress
end

function update_progress!(progress, stats, ntracked)
    t = time()
    if ntracked == progress.n || t > progress.tlast + progress.dt
        showvalues = make_showvalues(stats, ntracked)
        ProgressMeter.update!(progress, ntracked; showvalues)
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

function _run_paths(
        solver::Solver,
        starts::AbstractArray,
        progress = nothing,
        stop_early_cb::F = always_false;
        threading::Bool,
        catch_interrupt::Bool,
    ) where {F}
    return if threading
        threaded_solve(solver, starts, progress, stop_early_cb; catch_interrupt)
    else
        serial_solve(solver, starts, progress, stop_early_cb; catch_interrupt)
    end
end

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

    return Result(path_results; seed = solver.seed, algorithm = solver.algorithm)
end

function _ensure_tracker_pool!(solver::Solver)
    tracker = solver.trackers[1]
    ntrackers = length(solver.trackers)
    nthr = Threads.nthreads()
    resize!(solver.trackers, nthr)
    for i in (ntrackers + 1):nthr
        solver.trackers[i] = deepcopy(tracker)
    end
    return solver
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

    _ensure_tracker_pool!(solver)

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
        Result(assigned_results; seed = solver.seed, algorithm = solver.algorithm)
    else
        Result(path_results; seed = solver.seed, algorithm = solver.algorithm)
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
