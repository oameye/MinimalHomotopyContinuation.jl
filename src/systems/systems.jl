export fixed

"""
    fixed(F::System; compile_mode::AbstractCompileMode = CompileMixed())

Construct either a [`CompiledSystem`](@ref), an [`InterpretedSystem`](@ref) or a
[`MixedSystem`](@ref) based on `compile_mode`.
"""
function fixed(F::System; compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE)
    return fixed(F, compile_mode)
end

fixed(F::System, ::CompileAll) = CompiledSystem(F)
fixed(F::System, ::CompileNone) = InterpretedSystem(F)
fixed(F::System, ::CompileMixed) = MixedSystem(F)

fixed(F::AbstractSystem; compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE) = F

set_solution!(x, ::AbstractSystem, y) = (x .= y; x)
