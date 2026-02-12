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
        S,
        starts,
        target_parameters,
        progress,
        reducer;
        catch_interrupt = catch_interrupt,
        threading = threading,
    )
end

function make_many_progress(n::Integer; delay::Float64 = 0.0)
    desc = "Solving for $n parameters... "
    barlen = min(ProgressMeter.tty_width(desc, stdout, false), 40)
    progress = ProgressMeter.Progress(n; dt = 0.3, desc = desc, barlen = barlen, output = stdout)
    progress.tlast += delay
    return progress
end

function update_many_progress!(progress, results, k, paths_per_param)
    t = time()
    if k == progress.n || t > progress.tlast + progress.dt
        showvalues = make_many_showvalues(results, k, paths_per_param)
        ProgressMeter.update!(progress, k; showvalues = showvalues)
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

    first_target_it = iterate(many_target_parameters)
    first_target_it === nothing && throw(ArgumentError("`targets` must be non-empty."))
    q = first(first_target_it)
    state = last(first_target_it)

    target_parameters!(solver, q)
    if threading
        res = threaded_solve(solver, starts_buffer; catch_interrupt = false)
    else
        res = serial_solve(solver, starts_buffer; catch_interrupt = false)
    end

    first_result = reducer_apply(reducer, res, q)
    results = if reducer isa FlatMapReducer
        first_result isa AbstractArray ||
            throw(ArgumentError("FlatMapReducer must return an iterable collection."))
        collect(first_result)
    else
        [first_result]
    end

    k = 1
    m = length(starts_buffer)
    update_many_progress!(progress, results, k, m)

    try
        for q in Iterators.rest(many_target_parameters, state)
            target_parameters!(solver, q)
            if threading
                res = threaded_solve(solver, starts_buffer; catch_interrupt = false)
            else
                res = serial_solve(solver, starts_buffer; catch_interrupt = false)
            end

            value = reducer_apply(reducer, res, q)
            if reducer isa FlatMapReducer
                append!(results, value)
            else
                push!(results, value)
            end
            k += 1
            update_many_progress!(progress, results, k, m)
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
