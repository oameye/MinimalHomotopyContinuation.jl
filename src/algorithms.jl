export AbstractHCAlgorithm,
    UnknownAlgorithm,
    PolyhedralAlgorithm,
    TotalDegreeAlgorithm,
    PathTrackingAlgorithm,
    AbstractSweepReducer,
    IdentityReducer,
    MapReducer,
    FlatMapReducer,
    reducer_apply

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

struct PolyhedralAlgorithm{CompileModeT <: AbstractCompileMode} <: AbstractHCAlgorithm
    compile_mode::CompileModeT
    seed::UInt32
    tracker_options::TrackerOptions
    endgame_options::EndgameOptions
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
    return PolyhedralAlgorithm{CompileModeT}(
        compile_mode,
        seed,
        tracker_options,
        endgame_options,
        only_torus,
        only_non_zero,
    )
end

struct TotalDegreeAlgorithm{CompileModeT <: AbstractCompileMode} <: AbstractHCAlgorithm
    compile_mode::CompileModeT
    seed::UInt32
    tracker_options::TrackerOptions
    endgame_options::EndgameOptions
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
    return TotalDegreeAlgorithm{CompileModeT}(
        compile_mode,
        seed,
        tracker_options,
        endgame_options,
        gamma,
    )
end

struct PathTrackingAlgorithm{CompileModeT <: AbstractCompileMode} <: AbstractHCAlgorithm
    compile_mode::CompileModeT
    seed::UInt32
    tracker_options::TrackerOptions
    endgame_options::EndgameOptions
end

function PathTrackingAlgorithm(
        ;
        compile_mode::CompileModeT = DEFAULT_COMPILE_MODE,
        seed::UInt32 = rand(UInt32),
        tracker_options::TrackerOptions = TrackerOptions(),
        endgame_options::EndgameOptions = EndgameOptions(),
    ) where {CompileModeT <: AbstractCompileMode}
    return PathTrackingAlgorithm{CompileModeT}(
        compile_mode,
        seed,
        tracker_options,
        endgame_options,
    )
end
