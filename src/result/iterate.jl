"""
    ResultIterator{Iter} <: AbstractResult

A result iterator tracking paths lazily.
"""
struct ResultIterator{Iter, SolverT <: AbstractSolver} <: AbstractResult
    starts::Iter
    S::SolverT
    bitmask::Union{BitVector, Nothing}
end
function ResultIterator(
        starts::Iter,
        S::SolverT;
        bitmask = nothing,
    ) where {Iter, SolverT <: AbstractSolver}
    return ResultIterator{Iter, SolverT}(starts, S, bitmask)
end
ResultIterator(starts::AbstractVector{<:Number}, S::SolverT; bitmask = nothing) where {SolverT <: AbstractSolver} =
    ResultIterator([starts], S; bitmask = bitmask)

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


function Base.iterate(ri::ResultIterator, state = nothing)
    native_state = state === nothing ? 0 : state[1]
    next_ss = state === nothing ? iterate(ri.starts) : iterate(ri.starts, state[2])
    next_ss === nothing && return nothing

    start_value, new_ss_state = next_ss
    new_state = (native_state + 1, new_ss_state)

    if ri.bitmask === nothing
        return (track(ri.S.trackers[1], start_value), new_state)
    else
        while !ri.bitmask[new_state[1]] && next_ss !== nothing
            next_ss = iterate(ri.starts, new_state[2])
            next_ss === nothing && return nothing
            start_value, new_ss_state = next_ss
            new_state = (new_state[1] + 1, new_ss_state)
        end

        return (track(ri.S.trackers[1], start_value), new_state)
    end
end


function Base.length(ri::ResultIterator)
    if Base.IteratorSize(ri) == Base.SizeUnknown()
        k = 0
        for _ in ri.starts
            k += 1
        end
        k
    elseif Base.IteratorSize(ri) == Base.HasLength()
        if ri.bitmask !== nothing
            return (sum(ri.bitmask))
        else
            return (length(ri.starts))
        end
    elseif ri.starts isa Vector
        return length(ri.starts)
    end
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
        iter,
        init = 0.0 .* s,
    )
end

function Result(ri::ResultIterator)
    C = collect(ri)
    for i in 1:length(C)
        C[i].path_number = i
    end
    return Result(C; seed = ri.S.seed, algorithm = ri.S.algorithm)
end
