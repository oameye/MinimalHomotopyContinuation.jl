abstract type AbstractResult end
abstract type AbstractSolver end

######################
## MultiplicityInfo ##
######################

"""
    MultiplicityInfo

This contains information about the multiplicities of the solutions.
"""
struct MultiplicityInfo
    multiplicities::Dict{Int, Vector{Vector{Int}}}
    multiple_indicator::Set{Int}
end

function MultiplicityInfo(pathresults::Vector{<:PathResult})
    multiple_indicator = Set{Int}()
    multiplicities = compute_multiplicities(pathresults)
    for clusters in values(multiplicities), cluster in clusters
        for i in 2:length(cluster)
            push!(multiple_indicator, path_number(pathresults[cluster[i]]))
        end
    end
    for k in keys(multiplicities)
        multiplicities[k] = map(multiplicities[k]) do cluster
            map(i -> path_number(pathresults[i]), cluster)
        end
    end
    return MultiplicityInfo(multiplicities, multiple_indicator)
end

function compute_multiplicities(result::Vector{<:PathResult})
    D = Dict{Int, Vector{Vector{Int}}}()
    for m in _multiplicity_clusters(solution, result)
        if haskey(D, length(m))
            push!(D[length(m)], m)
        else
            D[length(m)] = [m]
        end
    end
    return D
end

function _multiplicity_clusters(
        f::F,
        values;
        distance_metric = EuclideanNorm(),
        atol::Float64 = 1.0e-14,
        rtol::Float64 = 1.0e-8,
    ) where {F <: Function}
    isempty(values) && return Vector{Vector{Int}}()

    mapped = map(f, values)
    clusters = Vector{Vector{Int}}()

    for i in eachindex(mapped)
        w = mapped[i]
        threshold = max(atol, rtol * _point_norm(distance_metric, w))
        assigned = false
        for cluster in clusters
            w_ref = mapped[first(cluster)]
            if _point_distance(distance_metric, w, w_ref) <= threshold
                push!(cluster, i)
                assigned = true
                break
            end
        end
        assigned || push!(clusters, [i])
    end

    filter!(c -> length(c) > 1, clusters)
    return clusters
end

_point_distance(metric::AbstractNorm, x, y) = distance(x, y, metric)
_point_distance(metric, x, y) = metric(x, y)

_point_norm(metric::AbstractNorm, x) = norm(x, metric)
_point_norm(metric, x) = metric(x, zero(x))

function assign_multiplicities!(path_results::Vector{<:PathResult}, I::MultiplicityInfo)
    for (k, clusters) in I.multiplicities, cluster in clusters, vᵢ in cluster
        path_results[vᵢ].multiplicity =
            max(k, something(path_results[vᵢ].winding_number, 1))
    end
    return path_results
end


"""
    Result

The result of [`solve`](@ref). This is a wrapper around the results of each single path
([`PathResult`](@ref)).
"""
struct Result{PathResultT <: PathResult, AlgT <: AbstractHCAlgorithm} <: AbstractResult
    path_results::Vector{PathResultT}
    tracked_paths::Int
    seed::UInt32
    algorithm::AlgT
    multiplicity_info::MultiplicityInfo
end
function Result(
        path_results::Vector{PathResultT};
        seed::UInt32 = rand(UInt32),
        tracked_paths = length(path_results),
        algorithm::AbstractHCAlgorithm = UnknownAlgorithm(),
    ) where {PathResultT <: PathResult}
    multiplicity_info = MultiplicityInfo(filter(is_singular, path_results))
    assign_multiplicities!(path_results, multiplicity_info)
    return Result{PathResultT, typeof(algorithm)}(
        path_results,
        tracked_paths,
        seed,
        algorithm,
        multiplicity_info,
    )
end

Base.size(r::Result) = (length(r),)
Base.length(r::Result) = length(r.path_results)
Base.getindex(r::Result, I) = getindex(r.path_results, I)

Base.iterate(r::Result) = iterate(r.path_results)
Base.iterate(r::Result, state) = iterate(r.path_results, state)
Base.lastindex(r::Result) = lastindex(r.path_results)
Base.eltype(::Type{Result{PathResultT, AlgT}}) where {PathResultT <: PathResult, AlgT <: AbstractHCAlgorithm} =
    PathResultT
Base.eltype(::Type{Result}) = PathResult
Base.keys(r::Result) = 1:length(r.path_results)
Base.values(r::Result) = r.path_results

seed(R::Result) = R.seed
algorithm(R::Result) = R.algorithm

path_results(R::Result) = R.path_results
path_results(R::AbstractVector{<:PathResult}) = R

multiple_indicator(::AbstractResult) =
    error("[`multiple_indicator`] Not defined for abstract results yet")
multiple_indicator(R::Result) = R.multiplicity_info.multiple_indicator
is_multiple_result(r::PathResult, R::AbstractResult) =
    path_number(r) ∈ multiple_indicator(R)
is_multiple_result(r::PathResult, R::AbstractVector{<:PathResult}) = false
