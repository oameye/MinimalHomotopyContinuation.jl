export FixedParameterSystem, fix_parameters

"""
    FixedParameterSystem(F:AbstractSystem, parameters)

Construct a system from the given [`AbstractSystem`](@ref) `F` with the given `parameters`
fixed.
"""
struct FixedParameterSystem{S <: AbstractSystem, T} <: AbstractSystem
    system::S
    parameters::Vector{T}
end
FixedParameterSystem(F::AbstractSystem, p; compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE) = FixedParameterSystem(
    F, p
)
FixedParameterSystem(F::System, p; compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE) = FixedParameterSystem(
    fixed(F; compile_mode), p
)
Base.size(F::FixedParameterSystem) = size(F.system)

ModelKit.variables(F::FixedParameterSystem) = variables(F.system)
ModelKit.parameters(F::FixedParameterSystem) = Variable[]

(F::FixedParameterSystem)(x, p = nothing) = F.system(x, F.parameters)

ModelKit.evaluate!(u, F::FixedParameterSystem, x, p = nothing) = evaluate!(
    u, F.system, x, F.parameters
)
ModelKit.evaluate_and_jacobian!(u, U, F::FixedParameterSystem, x, p = nothing) = evaluate_and_jacobian!(
    u, U, F.system, x, F.parameters
)
ModelKit.taylor!(u, v::Val, F::FixedParameterSystem, tx, p = nothing) = taylor!(
    u, v, F.system, tx, F.parameters
)

# Helper for fixed-parameter wrappers.
ModelKit.jacobian!(U, F::FixedParameterSystem{<:InterpretedSystem}, x, p = nothing) = ModelKit.jacobian!(
    U, F.system, x, F.parameters
)
ModelKit.jacobian!(U, F::FixedParameterSystem{<:InterpretedSystem}, x, p, cache) = ModelKit.jacobian!(
    U, F.system, x, F.parameters
)

"""
    fix_parameters(F::Union{System,AbstractSystem}, p; compile_mode::AbstractCompileMode = CompileMixed())

Fix the parameters of the given system `F`. Returns a [`FixedParameterSystem`](@ref).
"""
fix_parameters(F::Union{System, AbstractSystem}, p; compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE) = FixedParameterSystem(
    F, p; compile_mode
)
