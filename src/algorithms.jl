export AbstractHCAlgorithm,
    PolyhedralAlgorithm,
    TotalDegreeAlgorithm

abstract type AbstractHCAlgorithm end

Base.@kwdef struct PolyhedralAlgorithm <: AbstractHCAlgorithm
    compile::Union{Bool, Symbol} = COMPILE_DEFAULT[]
    seed::Union{Nothing, UInt32} = rand(UInt32)
    tracker_options::TrackerOptions = TrackerOptions()
    endgame_options::EndgameOptions = EndgameOptions()
    only_torus::Bool = false
    only_non_zero::Union{Nothing, Bool} = nothing
end

Base.@kwdef struct TotalDegreeAlgorithm <: AbstractHCAlgorithm
    compile::Union{Bool, Symbol} = COMPILE_DEFAULT[]
    seed::Union{Nothing, UInt32} = rand(UInt32)
    tracker_options::TrackerOptions = TrackerOptions()
    endgame_options::EndgameOptions = EndgameOptions()
    gamma::ComplexF64 = cis(2π * rand())
end
