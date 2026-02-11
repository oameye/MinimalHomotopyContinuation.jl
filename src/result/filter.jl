Base.@kwdef struct ResultStatistics
    total::Int
    nonsingular::Int
    singular::Int
    singular_with_multiplicity::Int = singular
    real::Int
    real_nonsingular::Int
    real_singular::Int
    real_singular_with_multiplicity::Int = real_singular
    at_infinity::Int
    excess_solution::Int = 0
    failed::Int = 0
end
Base.show(io::IO, stats::ResultStatistics) = print_fieldnames(io, stats)

function ResultStatistics(
        result::Result;
        real_atol::Float64 = 1.0e-6,
        real_rtol::Float64 = 0.0,
    )
    failed = at_infinity = excess_solution = 0
    nonsingular = singular = real_nonsingular = real_singular = 0
    singular_with_multiplicity = real_singular_with_multiplicity = 0
    for r in result
        is_multiple_result(r, result) && continue
        if is_failed(r)
            failed += 1
        elseif is_at_infinity(r)
            at_infinity += 1
        elseif is_excess_solution(r)
            excess_solution += 1
        elseif is_singular(r)
            if is_real(r, real_atol, real_rtol)
                real_singular += 1
                real_singular_with_multiplicity += something(multiplicity(r), 1)
            end
            singular += 1
            singular_with_multiplicity += something(multiplicity(r), 1)
        else
            if is_real(r, real_atol, real_rtol)
                real_nonsingular += 1
            end
            nonsingular += 1
        end
    end
    return ResultStatistics(
        nonsingular = nonsingular,
        singular = singular,
        singular_with_multiplicity = singular_with_multiplicity,
        real_nonsingular = real_nonsingular,
        real_singular = real_singular,
        real_singular_with_multiplicity = real_singular_with_multiplicity,
        real = real_nonsingular + real_singular,
        at_infinity = at_infinity,
        excess_solution = excess_solution,
        failed = failed,
        total = result.tracked_paths,
    )
end

statistics(r; kwargs...) = ResultStatistics(r; kwargs...)

const Results = Union{Result, AbstractVector{<:PathResult}}
const AbstractResults = Union{AbstractResult, AbstractVector{<:PathResult}}

results(R::Result; kwargs...) = results(identity, R; kwargs...)
results(R::AbstractVector{<:PathResult}; kwargs...) = results(identity, R; kwargs...)
results(R::AbstractResult; kwargs...) = results(identity, R; kwargs...)

function _filtered_results(
        f::Function,
        R;
        only_real::Bool = false,
        real_atol::Float64 = 1.0e-6,
        real_rtol::Float64 = 0.0,
        only_nonsingular::Bool = false,
        only_singular::Bool = false,
        only_finite::Bool = true,
        multiple_results::Bool = false,
    )
    filter_function =
        r ->
    (!only_real || is_real(r, real_atol, real_rtol)) &&
        (!only_nonsingular || is_nonsingular(r)) &&
        (!only_singular || is_singular(r)) &&
        (!only_finite || is_finite(r)) &&
        (multiple_results || !is_multiple_result(r, R))
    return Iterators.map(f, Iterators.filter(filter_function, R))
end

function results(
        f::Function,
        R::Result;
        only_real::Bool = false,
        real_atol::Float64 = 1.0e-6,
        real_rtol::Float64 = 0.0,
        only_nonsingular::Bool = false,
        only_singular::Bool = false,
        only_finite::Bool = true,
        multiple_results::Bool = false,
    )
    return collect(
        _filtered_results(
            f,
            R;
            only_real = only_real,
            real_atol = real_atol,
            real_rtol = real_rtol,
            only_nonsingular = only_nonsingular,
            only_singular = only_singular,
            only_finite = only_finite,
            multiple_results = multiple_results,
        ),
    )
end

function results(
        f::Function,
        R::AbstractVector{<:PathResult};
        only_real::Bool = false,
        real_atol::Float64 = 1.0e-6,
        real_rtol::Float64 = 0.0,
        only_nonsingular::Bool = false,
        only_singular::Bool = false,
        only_finite::Bool = true,
        multiple_results::Bool = false,
    )
    return collect(
        _filtered_results(
            f,
            R;
            only_real = only_real,
            real_atol = real_atol,
            real_rtol = real_rtol,
            only_nonsingular = only_nonsingular,
            only_singular = only_singular,
            only_finite = only_finite,
            multiple_results = multiple_results,
        ),
    )
end

function results(
        f::Function,
        R::AbstractResult;
        only_real::Bool = false,
        real_atol::Float64 = 1.0e-6,
        real_rtol::Float64 = 0.0,
        only_nonsingular::Bool = false,
        only_singular::Bool = false,
        only_finite::Bool = true,
        multiple_results::Bool = false,
    )
    multiple_results = true
    return _filtered_results(
        f,
        R;
        only_real = only_real,
        real_atol = real_atol,
        real_rtol = real_rtol,
        only_nonsingular = only_nonsingular,
        only_singular = only_singular,
        only_finite = only_finite,
        multiple_results = multiple_results,
    )
end

function _count_results(
        R;
        only_real::Bool = false,
        real_atol::Float64 = 1.0e-6,
        real_rtol::Float64 = 0.0,
        only_nonsingular::Bool = false,
        only_singular::Bool = false,
        only_finite::Bool = true,
        multiple_results::Bool = false,
    )
    return count(R) do r
        (!only_real || is_real(r, real_atol, real_rtol)) &&
            (!only_nonsingular || is_nonsingular(r)) &&
            (!only_singular || is_singular(r)) &&
            (!only_finite || isfinite(r)) &&
            (multiple_results || !is_multiple_result(r, R))
    end
end

function nresults(
        R::Result;
        only_real::Bool = false,
        real_atol::Float64 = 1.0e-6,
        real_rtol::Float64 = 0.0,
        only_nonsingular::Bool = false,
        only_singular::Bool = false,
        only_finite::Bool = true,
        multiple_results::Bool = false,
    )
    return _count_results(
        R;
        only_real = only_real,
        real_atol = real_atol,
        real_rtol = real_rtol,
        only_nonsingular = only_nonsingular,
        only_singular = only_singular,
        only_finite = only_finite,
        multiple_results = multiple_results,
    )
end

function nresults(
        R::AbstractVector{<:PathResult};
        only_real::Bool = false,
        real_atol::Float64 = 1.0e-6,
        real_rtol::Float64 = 0.0,
        only_nonsingular::Bool = false,
        only_singular::Bool = false,
        only_finite::Bool = true,
        multiple_results::Bool = false,
    )
    return _count_results(
        R;
        only_real = only_real,
        real_atol = real_atol,
        real_rtol = real_rtol,
        only_nonsingular = only_nonsingular,
        only_singular = only_singular,
        only_finite = only_finite,
        multiple_results = multiple_results,
    )
end

function nresults(
        R::AbstractResult;
        only_real::Bool = false,
        real_atol::Float64 = 1.0e-6,
        real_rtol::Float64 = 0.0,
        only_nonsingular::Bool = false,
        only_singular::Bool = false,
        only_finite::Bool = true,
        multiple_results::Bool = false,
    )
    multiple_results = true
    return _count_results(
        R;
        only_real = only_real,
        real_atol = real_atol,
        real_rtol = real_rtol,
        only_nonsingular = only_nonsingular,
        only_singular = only_singular,
        only_finite = only_finite,
        multiple_results = multiple_results,
    )
end

solutions(result::AbstractResults; only_nonsingular = true, kwargs...) =
    results(solution, result; only_nonsingular = only_nonsingular, kwargs...)

function real_solutions(
        result::AbstractResults;
        atol::Float64 = 1.0e-6,
        rtol::Float64 = 0.0,
        kwargs...,
    )
    return results(
        real ∘ solution,
        result;
        only_real = true,
        real_atol = atol,
        real_rtol = rtol,
        kwargs...,
    )
end

nonsingular(R::AbstractResults; kwargs...) = results(R; only_nonsingular = true, kwargs...)

function singular(R::AbstractResults; kwargs...)
    return results(R; only_singular = true, kwargs...)
end

Base.real(R::Results; kwargs...) = filter(r -> is_real(r; kwargs...), path_results(R))

failed(R::Results) = filter(is_failed, path_results(R))

at_infinity(R::Results) = filter(is_at_infinity, path_results(R))

nsolutions(R::AbstractResults; only_nonsingular = true, options...) =
    nresults(R; only_nonsingular = only_nonsingular, options...)

nfinite(R::AbstractResults; kwargs...) = nresults(R; only_finite = true, kwargs...)

function nsingular(R::AbstractResult; counting_multiplicities::Bool = false, kwargs...)
    if counting_multiplicities
        return (nresults(R; only_singular = true, multiple_results = true, kwargs...))
    else
        return (nresults(R; only_singular = true, multiple_results = false, kwargs...))
    end
end

nat_infinity(R::AbstractResults) = count(is_at_infinity, R)

nexcess_solutions(R::AbstractResults) = count(is_excess_solution, R)

nfailed(R::AbstractResults) = count(is_failed, R)

nnonsingular(R::AbstractResults) = count(is_nonsingular, R)

nreal(R::Results; kwargs...) = nresults(R, only_real = true, kwargs...)
nreal(R::AbstractResult; kwargs...) =
    nresults(R, only_real = true, multiple_results = true, kwargs...)

ntracked(R::AbstractResult) = length(R)
ntracked(R::Result) = R.tracked_paths
