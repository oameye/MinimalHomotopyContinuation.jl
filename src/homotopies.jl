# internal only
include("homotopies/toric_homotopy.jl")

# public, these should be linked on the top
include("homotopies/mixed_homotopy.jl")
include("homotopies/parameter_homotopy.jl")
include("homotopies/coefficient_homotopy.jl")
include("homotopies/straight_line_homotopy.jl")
include("homotopies/fixed_parameter_homotopy.jl")

"""
    fixed(H::Homotopy; compile_mode::AbstractCompileMode = CompileMixed())

Construct either a [`CompiledHomotopy`](@ref), an [`InterpretedHomotopy`](@ref) or a
[`MixedHomotopy`](@ref) based on `compile_mode`.
"""
function fixed(
        H::Homotopy;
        compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE,
    )
    return fixed(H, compile_mode)
end

fixed(H::Homotopy, ::CompileAll) = CompiledHomotopy(H)
fixed(H::Homotopy, ::CompileNone) = InterpretedHomotopy(H)
fixed(H::Homotopy, ::CompileMixed) = MixedHomotopy(H)

fixed(H::AbstractHomotopy; compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE) = H

function set_solution!(x::AbstractVector, H::AbstractHomotopy, y::AbstractVector, t)
    return x .= y
end
get_solution(H::AbstractHomotopy, x::AbstractVector, t) = copy(x)

start_parameters!(H::AbstractHomotopy, p) = H
target_parameters!(H::AbstractHomotopy, p) = H

ModelKit.taylor!(u, v::Val, H::AbstractHomotopy, tx, t, incremental::Bool) =
    taylor!(u, v, H, tx, t)
