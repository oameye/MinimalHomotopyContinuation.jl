module HomotopyContinuation

# TODO:
# - use use the interface of [EnumX.jl](https://github.com/fredrikekre/EnumX.jl)

export ModelKit

using Arblib: Arblib
using DynamicPolynomials: @polyvar
import ElasticArrays: ElasticArrays, ElasticArray
import LinearAlgebra
import LoopVectorization
import MixedSubdivisions: MixedSubdivisions, MixedCell
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

include("model_kit/ModelKit.jl")
export @polyvar

include("utils.jl")
include("norm.jl")
include("linear_algebra.jl")

include("systems/compile_modes.jl")
include("systems/mixed_system.jl")
include("systems/fixed_parameter_system.jl")
include("systems/systems.jl")

include("homotopies/toric_homotopy.jl")
include("homotopies/mixed_homotopy.jl")
include("homotopies/parameter_homotopy.jl")
include("homotopies/coefficient_homotopy.jl")
include("homotopies/straight_line_homotopy.jl")
include("homotopies/fixed_parameter_homotopy.jl")
include("homotopies/homotopies.jl")

include("newton/predictor.jl")
include("newton/corrector.jl")
include("newton/newton.jl")

include("tracking/tracker.jl")
include("tracking/valuation.jl")
include("tracking/path_result.jl")
include("tracking/endgame_tracker.jl")

include("algorithms/algorithms.jl")
include("algorithms/total_degree.jl")
include("algorithms/binomial_system.jl")
include("algorithms/polyhedral.jl")

include("problems/problems_core.jl")
include("problems/problems_frontend.jl")

include("result/types.jl")
include("result/iterate.jl")
include("result/filter.jl")
include("result/show.jl")

include("solve/types.jl")
include("solve/prepare.jl")
include("solve/execute_paths.jl")
include("solve/execute_sweeps.jl")
include("solve/api.jl")

include("public_api.jl")

end #module
