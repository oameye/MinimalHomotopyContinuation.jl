export FixedParameterHomotopy, fix_parameters

"""
    FixedParameterHomotopy(H:AbstractHomotopy, parameters)

Construct a homotopy from the given [`AbstractHomotopy`](@ref) `H` with the given `parameters`
fixed.
"""
struct FixedParameterHomotopy{S <: AbstractHomotopy, T} <: AbstractHomotopy
    homotopy::S
    parameters::Vector{T}
end
FixedParameterHomotopy(H::Homotopy, p; compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE) = FixedParameterHomotopy(
    fixed(H; compile_mode), p
)
FixedParameterHomotopy(H::AbstractHomotopy, p; compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE) = FixedParameterHomotopy(
    H, p
)
Base.size(H::FixedParameterHomotopy) = size(H.homotopy)

ModelKit.variables(H::FixedParameterHomotopy) = variables(H.homotopy)
ModelKit.parameters(H::FixedParameterHomotopy) = Variable[]

(H::FixedParameterHomotopy)(x, t, p = nothing) = H.homotopy(x, t, H.parameters)

ModelKit.evaluate!(u, H::FixedParameterHomotopy, x, t) = evaluate!(
    u, H.homotopy, x, t, H.parameters
)
ModelKit.evaluate_and_jacobian!(u, U, H::FixedParameterHomotopy, x, t) = evaluate_and_jacobian!(
    u, U, H.homotopy, x, t, H.parameters
)
ModelKit.taylor!(u, v::Val, H::FixedParameterHomotopy, tx, t) = taylor!(
    u, v, H.homotopy, tx, t, H.parameters
)

"""
    fix_parameters(H::Union{Homotopy,AbstractHomotopy}, p; compile_mode::AbstractCompileMode = CompileMixed())

Fix the parameters of the given homotopy `H`. Returns a [`FixedParameterHomotopy`](@ref).
"""
fix_parameters(H::Union{Homotopy, AbstractHomotopy}, p; compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE) = FixedParameterHomotopy(
    H, p; compile_mode
)
