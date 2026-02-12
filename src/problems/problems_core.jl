abstract type AbstractHCProblem end

struct SystemProblem{SystemT, TargetT} <: AbstractHCProblem
    system::SystemT
    target_parameters::TargetT
end

SystemProblem(F::Union{System, AbstractSystem}; target_parameters = nothing) =
    SystemProblem{typeof(F), typeof(target_parameters)}(F, target_parameters)

struct ParameterHomotopyProblem{SystemT, StartsT, PStartT, PTargetT} <: AbstractHCProblem
    system::SystemT
    start_solutions::StartsT
    start_parameters::PStartT
    target_parameters::PTargetT
end

ParameterHomotopyProblem(
    F::Union{System, AbstractSystem},
    start_solutions;
    start_parameters,
    target_parameters,
) =
    ParameterHomotopyProblem{
    typeof(F), typeof(start_solutions),
    typeof(start_parameters), typeof(target_parameters),
}(
    F,
    start_solutions,
    start_parameters,
    target_parameters,
)

struct ParameterSweepProblem{SystemT, StartsT, PStartT, TargetsT, ReducerT} <: AbstractHCProblem
    system::SystemT
    start_solutions::StartsT
    start_parameters::PStartT
    targets::TargetsT
    reducer::ReducerT
end

function ParameterSweepProblem(
        F::Union{System, AbstractSystem},
        start_solutions;
        start_parameters,
        targets,
        reducer::AbstractSweepReducer = IdentityReducer(),
    )
    return ParameterSweepProblem{
        typeof(F),
        typeof(start_solutions),
        typeof(start_parameters),
        typeof(targets),
        typeof(reducer),
    }(
        F,
        start_solutions,
        start_parameters,
        targets,
        reducer,
    )
end

struct HomotopyProblem{HomotopyT, StartsT} <: AbstractHCProblem
    homotopy::HomotopyT
    start_solutions::StartsT
end

HomotopyProblem(H::Union{Homotopy, AbstractHomotopy}, starts) =
    HomotopyProblem{typeof(H), typeof(starts)}(H, starts)
