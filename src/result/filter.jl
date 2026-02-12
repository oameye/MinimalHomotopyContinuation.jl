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

function ResultStatistics(result::Result; real_atol::Float64 = 1.0e-6, real_rtol::Float64 = 0.0)
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
    return ResultStatistics(;
        nonsingular,
        singular,
        singular_with_multiplicity,
        real_nonsingular,
        real_singular,
        real_singular_with_multiplicity,
        real = real_nonsingular + real_singular,
        at_infinity,
        excess_solution,
        failed,
        total = result.tracked_paths,
    )
end

statistics(r; kwargs...) = ResultStatistics(r; kwargs...)

const Results = Union{Result, AbstractVector{<:PathResult}}
const AbstractResults = Union{AbstractResult, AbstractVector{<:PathResult}}

results(R::AbstractResults; kwargs...) = results(identity, R; kwargs...)

Base.@kwdef struct ResultFilterOptions
    only_real::Bool = false
    real_atol::Float64 = 1.0e-6
    real_rtol::Float64 = 0.0
    only_nonsingular::Bool = false
    only_singular::Bool = false
    only_finite::Bool = true
    multiple_results::Bool = false
end

_result_filter_options(; kwargs...) = ResultFilterOptions(; kwargs...)

function _with_multiple_results(opts::ResultFilterOptions, multiple_results::Bool)
    return ResultFilterOptions(;
        only_real = opts.only_real,
        real_atol = opts.real_atol,
        real_rtol = opts.real_rtol,
        only_nonsingular = opts.only_nonsingular,
        only_singular = opts.only_singular,
        only_finite = opts.only_finite,
        multiple_results = multiple_results,
    )
end

_effective_filter_options(::AbstractResult, opts::ResultFilterOptions) = _with_multiple_results(
    opts, true
)
_effective_filter_options(::AbstractVector{<:PathResult}, opts::ResultFilterOptions) = opts
_effective_filter_options(::Result, opts::ResultFilterOptions) = opts

function _matches_result(r, R, opts::ResultFilterOptions)
    return (!opts.only_real || is_real(r, opts.real_atol, opts.real_rtol)) &&
        (!opts.only_nonsingular || is_nonsingular(r)) &&
        (!opts.only_singular || is_singular(r)) &&
        (!opts.only_finite || is_finite(r)) &&
        (opts.multiple_results || !is_multiple_result(r, R))
end

function _filtered_results(f::Function, R, opts::ResultFilterOptions)
    filter_function = r -> _matches_result(r, R, opts)
    return Iterators.map(f, Iterators.filter(filter_function, R))
end

_results_output(::Result, iter) = collect(iter)
_results_output(::AbstractVector{<:PathResult}, iter) = collect(iter)
_results_output(::AbstractResult, iter) = iter

function _results_impl(f::Function, R; kwargs...)
    opts = _effective_filter_options(R, _result_filter_options(; kwargs...))
    iter = _filtered_results(f, R, opts)
    return _results_output(R, iter)
end

results(f::Function, R::AbstractResults; kwargs...) = _results_impl(f, R; kwargs...)

function _count_results(R; kwargs...)
    opts = _effective_filter_options(R, _result_filter_options(; kwargs...))
    return count(R) do r
        _matches_result(r, R, opts)
    end
end

nresults(R::AbstractResults; kwargs...) = _count_results(R; kwargs...)

solutions(result::AbstractResults; only_nonsingular = true, kwargs...) = results(
    solution, result; only_nonsingular, kwargs...
)

function real_solutions(
        result::AbstractResults; atol::Float64 = 1.0e-6, rtol::Float64 = 0.0, kwargs...
    )
    return results(
        real ∘ solution, result; only_real = true, real_atol = atol, real_rtol = rtol, kwargs...
    )
end

nonsingular(R::AbstractResults; kwargs...) = results(R; only_nonsingular = true, kwargs...)

function singular(R::AbstractResults; kwargs...)
    return results(R; only_singular = true, kwargs...)
end

Base.real(R::Result; kwargs...) = [r for r in path_results(R) if is_real(r; kwargs...)]
Base.real(R::AbstractVector{<:PathResult}; kwargs...) = [
    r for r in R if is_real(r; kwargs...)
]

failed(R::Result) = [r for r in path_results(R) if is_failed(r)]
failed(R::AbstractVector{<:PathResult}) = [r for r in R if is_failed(r)]

at_infinity(R::Result) = [r for r in path_results(R) if is_at_infinity(r)]
at_infinity(R::AbstractVector{<:PathResult}) = [r for r in R if is_at_infinity(r)]

nsolutions(R::AbstractResults; only_nonsingular = true, options...) = nresults(
    R; only_nonsingular, options...
)

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
nreal(R::AbstractResult; kwargs...) = nresults(
    R, only_real = true, multiple_results = true, kwargs...
)

ntracked(R::AbstractResult) = length(R)
ntracked(R::Result) = R.tracked_paths
