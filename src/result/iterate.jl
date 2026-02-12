"""
    ResultIterator{Iter} <: AbstractResult

A result iterator tracking paths lazily.
"""
struct ResultIterator{
        Iter, SolverT <: AbstractSolver, BitmaskT <: Union{Nothing, AbstractVector{Bool}},
    } <: AbstractResult
    starts::Iter
    S::SolverT
    bitmask::BitmaskT
end
function ResultIterator(
        starts::Iter, S::SolverT; bitmask::Union{Nothing, AbstractVector{Bool}} = nothing
    ) where {Iter, SolverT <: AbstractSolver}
    return ResultIterator{Iter, SolverT, typeof(bitmask)}(starts, S, bitmask)
end
ResultIterator(starts::AbstractVector{<:Number}, S::SolverT; bitmask::Union{Nothing, AbstractVector{Bool}} = nothing) where {SolverT <: AbstractSolver} = ResultIterator(
    [starts], S; bitmask
)

seed(ri::ResultIterator) = ri.S.seed

path_results(ri::ResultIterator) = collect(ri)

start_solutions(ri::ResultIterator) = ri.starts

solver(ri::ResultIterator) = ri.S

bitmask(ri::ResultIterator) = ri.bitmask

function Base.show(io::IO, ri::ResultIterator{Iter}) where {Iter}
    header = "ResultIterator"
    println(io, header)
    println(io, "="^(length(header)))
    println(io, "•  start solutions: $(nameof(Iter))")
    tracker = ri.S.trackers[1]
    if tracker isa EndgameTracker
        n = typeof(tracker.tracker.homotopy)
        println(io, "•  homotopy: $(nameof(n))")
    elseif tracker isa PolyhedralTracker
        println(io, "•  homotopy: Polyhedral")
    end
    println(io, "•  algorithm: $(nameof(typeof(ri.S.algorithm)))")
    return !isnothing(ri.bitmask) && println(io, "•  filtering bitmask")
end

function Base.IteratorSize(ri::ResultIterator)
    if ri.starts isa Vector
        return Base.HasLength()
    else
        ri.bitmask === nothing ? Base.IteratorSize(ri.starts) : Base.HasLength()
    end
end
Base.IteratorEltype(::ResultIterator) = Base.HasEltype()
Base.eltype(::ResultIterator) = PathResult

function _next_start_state(ri::ResultIterator, state)
    native_index, next_ss = if state === nothing
        0, iterate(ri.starts)
    else
        state[1], iterate(ri.starts, state[2])
    end

    while next_ss !== nothing
        start_value, ss_state = next_ss
        native_index += 1
        new_state = (native_index, ss_state)
        if isnothing(ri.bitmask) || ri.bitmask[native_index]
            return start_value, new_state
        end
        next_ss = iterate(ri.starts, ss_state)
    end
    return nothing
end


function Base.iterate(ri::ResultIterator, state = nothing)
    next_item = _next_start_state(ri, state)
    next_item === nothing && return nothing
    start_value, new_state = next_item
    return track(ri.S.trackers[1], start_value), new_state
end

function _iterator_length(iter)
    return if Base.IteratorSize(iter) == Base.SizeUnknown()
        count(_ -> true, iter)
    else
        length(iter)
    end
end

function Base.length(ri::ResultIterator)
    return isnothing(ri.bitmask) ? _iterator_length(ri.starts) : sum(ri.bitmask)
end

function bitmask(f::Function, ri::ResultIterator)
    return BitVector(map(f, ri))
end

function bitmask_filter(f::Function, ri::ResultIterator)
    bm = bitmask(f, ri)
    return ResultIterator(ri.starts, ri.S; bitmask = bm)
end

function trace(iter::ResultIterator)
    s = solution(first(iter))
    return mapreduce(
        x -> isfinite(x) ? solution(x) : zeros(ComplexF64, length(s)),
        +,
        iter;
        init = 0.0 .* s,
    )
end

function Result(ri::ResultIterator)
    C = collect(ri)
    i = 1
    n = length(C)
    while i <= n
        C[i].path_number = i
        i += 1
    end
    return Result(C; seed = ri.S.seed, algorithm = ri.S.algorithm)
end
