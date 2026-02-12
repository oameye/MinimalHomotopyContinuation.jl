function solve(
        S::Solver,
        starts,
        target_parameters;
        reducer::AbstractSweepReducer = IdentityReducer(),
        show_progress::Bool = true,
        threading::Bool = Threads.nthreads() > 1,
        catch_interrupt::Bool = true,
    )
    n = length(target_parameters)

    progress = show_progress ? make_many_progress(n; delay = 0.3) : nothing

    return many_solve(
        S, starts, target_parameters, progress, reducer; catch_interrupt, threading
    )
end

function make_many_progress(n::Integer; delay::Float64 = 0.0)
    desc = "Solving for $n parameters... "
    barlen = min(ProgressMeter.tty_width(desc, stdout, false), 40)
    progress = ProgressMeter.Progress(n; dt = 0.3, desc, barlen, output = stdout)
    progress.tlast += delay
    return progress
end

function update_many_progress!(progress, results, k, paths_per_param)
    t = time()
    if k == progress.n || t > progress.tlast + progress.dt
        showvalues = make_many_showvalues(results, k, paths_per_param)
        ProgressMeter.update!(progress, k; showvalues)
    end
    return nothing
end

@noinline function make_many_showvalues(results, k, paths_per_param)
    return [
        ("# parameters solved", k),
        ("# paths tracked", paths_per_param * k),
        ("# results", length(results)),
    ]
end

update_many_progress!(::Nothing, results, k, paths_per_param) = nothing

function _init_reduced_results(::FlatMapReducer, value)
    value isa AbstractArray ||
        throw(ArgumentError("FlatMapReducer must return an iterable collection."))
    return collect(value)
end

_init_reduced_results(::AbstractSweepReducer, value) = [value]

function _append_reduced_result!(results, ::FlatMapReducer, value)
    append!(results, value)
    return results
end

function _append_reduced_result!(results, ::AbstractSweepReducer, value)
    push!(results, value)
    return results
end

function many_solve(
        solver::Solver,
        starts,
        many_target_parameters,
        progress,
        reducer::AbstractSweepReducer;
        threading::Bool,
        catch_interrupt::Bool,
    )

    if isa(starts, TotalDegreeStartSolutionsIterator)
        @warn "Solving for many parameters with total degree start system is not recommended. Instead, one should use a two-step approach: first solve a system with generic parameters, and track its solutions to the desired parameters. See https://www.juliahomotopycontinuation.org/guides/many-systems/."
    elseif isa(starts, PolyhedralStartSolutionsIterator)
        @error "Solving for many parameters with polyhedral start system is not implemented and also not recommended. Instead, one should use a two-step approach: solve a system with generic parameters, and track its solutions to the desired parameters. See https://www.juliahomotopycontinuation.org/guides/many-systems/."
    end
    starts_buffer = starts isa AbstractArray ? starts : collect(starts)
    run_kwargs = (; threading, catch_interrupt = false)
    results = nothing
    k = 0
    m = length(starts_buffer)

    try
        for (kᵢ, q) in enumerate(many_target_parameters)
            k = kᵢ
            target_parameters!(solver, q)
            res = _run_paths(solver, starts_buffer; run_kwargs...)
            value = reducer_apply(reducer, res, q)
            if results === nothing
                results = _init_reduced_results(reducer, value)
            else
                _append_reduced_result!(results, reducer, value)
            end
            update_many_progress!(progress, results, kᵢ, m)
        end
    catch e
        interrupted = isa(e, InterruptException) ||
            (isa(e, TaskFailedException) && isa(e.task.exception, InterruptException))
        if !interrupted || !catch_interrupt
            rethrow(e)
        end
    end
    results === nothing && throw(ArgumentError("`targets` must be non-empty."))

    return results
end
