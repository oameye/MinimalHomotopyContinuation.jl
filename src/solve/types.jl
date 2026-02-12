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

struct Solver{T <: AbstractPathTracker, AlgT <: AbstractHCAlgorithm} <: AbstractSolver
    trackers::Vector{T}
    seed::UInt32
    stats::SolveStats
    algorithm::AlgT
end

function Solver(
        tracker::AbstractPathTracker, algorithm::AlgT; seed::UInt32 = rand(UInt32)
    ) where {AlgT <: AbstractHCAlgorithm}
    return Solver([tracker], seed, SolveStats(), algorithm)
end

Base.show(io::IO, solver::Solver) = print(io, typeof(solver), "()")

struct PathSolveCache{SolverT, StartsT, StopEarlyT}
    solver::SolverT
    starts::StartsT
    stop_early_cb::StopEarlyT
    show_progress::Bool
    threading::Bool
    catch_interrupt::Bool
end

struct PathIteratorSolveCache{SolverT, StartsT, BitmaskT}
    solver::SolverT
    starts::StartsT
    bitmask::BitmaskT
end

struct SweepSolveCache{SolverT, StartsT, TargetsT, ReducerT}
    solver::SolverT
    starts::StartsT
    targets::TargetsT
    reducer::ReducerT
    show_progress::Bool
    threading::Bool
    catch_interrupt::Bool
end

struct SweepIteratorSolveCache{SolversT, StartsT, BitmaskT}
    solvers::SolversT
    starts::StartsT
    bitmask::BitmaskT
end
