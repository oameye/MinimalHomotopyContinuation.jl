@reexport module ModelKit

export @var,
    @unique_var,
    AbstractSystem,
    AbstractHomotopy,
    Expression,
    Variable,
    Interpreter,
    System,
    InterpretedSystem,
    CompiledSystem,
    Homotopy,
    InterpretedHomotopy,
    CompiledHomotopy,
    TaylorVector,
    TruncatedTaylorSeries,
    Variable,
    coefficients,
    coeffs_as_dense_poly,
    degree,
    degrees,
    differentiate,
    dense_poly,
    evaluate,
    evaluate!,
    evaluate_and_jacobian!,
    expand,
    exponents_coefficients,
    expressions,
    horner,
    interpreter,
    is_homogeneous,
    jacobian,
    jacobian!,
    parameters,
    nparameters,
    nvariables,
    monomials,
    optimize,
    parameters,
    poly_from_exponents_coefficients,
    subs,
    support_coefficients,
    rand_poly,
    taylor!,
    to_dict,
    to_number,
    variables,
    vectors

using ..DoubleDouble: ComplexDF64
import ..is_real    # required to ensure that ModelKit and HomotopyContinuation define methods of the same `is_real` function

import Arblib: Arblib, Acb, AcbRef, AcbRefVector
import SimpleGraphs
import LinearAlgebra
import MultivariatePolynomials:
    MultivariatePolynomials,
    coefficients,
    degree,
    differentiate,
    nvariables,
    monomials,
    subs,
    variables
const MP = MultivariatePolynomials

include("symengine.jl")
include("symbolic.jl")
include("operations.jl")
include("intermediate_representation.jl")
include("taylor.jl")
include("acb.jl")
include("instruction_sequence.jl")
include("instruction_interpreter.jl")
include("abstract_system_homotopy.jl")
include("interpreted_system.jl")
include("interpreted_homotopy.jl")
include("compiled_system_homotopy.jl")
end # module
