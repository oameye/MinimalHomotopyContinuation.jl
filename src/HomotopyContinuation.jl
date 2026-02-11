module HomotopyContinuation

export ModelKit

using Arblib: Arblib
using DynamicPolynomials: @polyvar
import ElasticArrays: ElasticArrays, ElasticArray
import LinearAlgebra
import LoopVectorization
import MixedSubdivisions: MixedSubdivisions, MixedCell, mixed_volume
import MultivariatePolynomials
const MP = MultivariatePolynomials
using Parameters: @unpack
import ProgressMeter
const PM = ProgressMeter
import Random
import Printf
using Reexport: @reexport
import StructArrays
import CommonSolve: init, solve, solve!
import Base: push!

const LA = LinearAlgebra

# To ensure that ModelKit and HomotopyContinuation define methods of the same `is_real` function:
function is_real end

include("DoubleDouble.jl")
using .DoubleDouble
Arblib.Arb(a::DoubleF64) = Arblib.Arb(BigFloat(a))
Arblib.Arf(a::DoubleF64) = Arblib.Arf(BigFloat(a))

include("ModelKit.jl")
export @polyvar

include("compile_modes.jl")

include("utils.jl")
include("norm.jl")
include("linear_algebra.jl")
include("systems.jl")
include("homotopies.jl")
include("predictor.jl")
include("newton_corrector.jl")
include("newton.jl")
include("tracker.jl")
include("valuation.jl")
include("path_result.jl")
include("endgame_tracker.jl")
include("algorithms.jl")
include("problems.jl")
include("total_degree.jl")
include("binomial_system.jl")
include("polyhedral.jl")
include("result.jl")
include("solve.jl")

end #module
