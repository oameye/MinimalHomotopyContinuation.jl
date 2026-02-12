abstract type AbstractHCAlgorithm end

struct UnknownAlgorithm <: AbstractHCAlgorithm end

abstract type AbstractSweepReducer end

struct IdentityReducer <: AbstractSweepReducer end

struct MapReducer{F} <: AbstractSweepReducer
    f::F
end

struct FlatMapReducer{F} <: AbstractSweepReducer
    f::F
end

reducer_apply(::IdentityReducer, result, params) = result
reducer_apply(reducer::MapReducer, result, params) = reducer.f(result, params)
reducer_apply(reducer::FlatMapReducer, result, params) = reducer.f(result, params)

struct AlgorithmOptions{CompileModeT <: AbstractCompileMode}
    compile_mode::CompileModeT
    seed::UInt32
    tracker_options::TrackerOptions
    endgame_options::EndgameOptions
end

struct PolyhedralAlgorithm{CompileModeT <: AbstractCompileMode} <: AbstractHCAlgorithm
    options::AlgorithmOptions{CompileModeT}
    only_torus::Bool
    only_non_zero::Bool
end

function PolyhedralAlgorithm(
        ;
        compile_mode::CompileModeT = DEFAULT_COMPILE_MODE,
        seed::UInt32 = rand(UInt32),
        tracker_options::TrackerOptions = TrackerOptions(),
        endgame_options::EndgameOptions = EndgameOptions(),
        only_torus::Bool = false,
        only_non_zero::Bool = only_torus,
    ) where {CompileModeT <: AbstractCompileMode}
    opts = AlgorithmOptions(compile_mode, seed, tracker_options, endgame_options)
    return PolyhedralAlgorithm{CompileModeT}(opts, only_torus, only_non_zero)
end

struct TotalDegreeAlgorithm{CompileModeT <: AbstractCompileMode} <: AbstractHCAlgorithm
    options::AlgorithmOptions{CompileModeT}
    gamma::ComplexF64
end

function TotalDegreeAlgorithm(
        ;
        compile_mode::CompileModeT = DEFAULT_COMPILE_MODE,
        seed::UInt32 = rand(UInt32),
        tracker_options::TrackerOptions = TrackerOptions(),
        endgame_options::EndgameOptions = EndgameOptions(),
        gamma::ComplexF64 = cis(2π * rand()),
    ) where {CompileModeT <: AbstractCompileMode}
    opts = AlgorithmOptions(compile_mode, seed, tracker_options, endgame_options)
    return TotalDegreeAlgorithm{CompileModeT}(opts, gamma)
end

struct PathTrackingAlgorithm{CompileModeT <: AbstractCompileMode} <: AbstractHCAlgorithm
    options::AlgorithmOptions{CompileModeT}
end

function PathTrackingAlgorithm(
        ;
        compile_mode::CompileModeT = DEFAULT_COMPILE_MODE,
        seed::UInt32 = rand(UInt32),
        tracker_options::TrackerOptions = TrackerOptions(),
        endgame_options::EndgameOptions = EndgameOptions(),
    ) where {CompileModeT <: AbstractCompileMode}
    opts = AlgorithmOptions(compile_mode, seed, tracker_options, endgame_options)
    return PathTrackingAlgorithm{CompileModeT}(opts)
end

@inline function _algorithm_getproperty(alg, sym::Symbol)
    if sym === :compile_mode || sym === :seed || sym === :tracker_options ||
       sym === :endgame_options
        return getfield(getfield(alg, :options), sym)
    end
    return getfield(alg, sym)
end

Base.getproperty(alg::PolyhedralAlgorithm, sym::Symbol) = _algorithm_getproperty(alg, sym)
Base.getproperty(alg::TotalDegreeAlgorithm, sym::Symbol) = _algorithm_getproperty(alg, sym)
Base.getproperty(alg::PathTrackingAlgorithm, sym::Symbol) = _algorithm_getproperty(alg, sym)
