abstract type AbstractCompileMode end

struct CompileAll <: AbstractCompileMode end
struct CompileMixed <: AbstractCompileMode end
struct CompileNone <: AbstractCompileMode end

const DEFAULT_COMPILE_MODE = CompileMixed()

compile_mode_symbol(::CompileAll) = :all
compile_mode_symbol(::CompileMixed) = :mixed
compile_mode_symbol(::CompileNone) = :none

function Base.show(io::IO, mode::AbstractCompileMode)
    return print(io, Symbol(mode))
end

Base.Symbol(mode::AbstractCompileMode) = compile_mode_symbol(mode)
