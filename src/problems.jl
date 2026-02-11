export AbstractHCProblem,
    SystemProblem,
    ParameterHomotopyProblem,
    ParameterSweepProblem,
    HomotopyProblem

abstract type AbstractHCProblem end

struct SystemProblem{SystemT, TargetT} <: AbstractHCProblem
    system::SystemT
    target_parameters::TargetT
end

SystemProblem(F::Union{System, AbstractSystem}; target_parameters = nothing) =
    SystemProblem{typeof(F), typeof(target_parameters)}(F, target_parameters)

function SystemProblem(
        F::AbstractVector{Expression};
        parameters = Variable[],
        variables = setdiff(variables(F), parameters),
        variable_ordering = variables,
        target_parameters = nothing,
    )
    sys = System(F, variables = variable_ordering, parameters = parameters)
    return SystemProblem(sys; target_parameters = target_parameters)
end

function SystemProblem(
        F::AbstractVector{<:MP.AbstractPolynomial};
        parameters = similar(MP.variables(F), 0),
        variables = setdiff(MP.variables(F), parameters),
        variable_ordering = variables,
        target_parameters = nothing,
    )
    if isnothing(target_parameters) && isempty(parameters)
        sys, target_parameters = ModelKit.system_with_coefficents_as_params(
            F;
            variables = variable_ordering,
        )
    else
        sys = System(F, variables = variable_ordering, parameters = parameters)
    end
    return SystemProblem(sys; target_parameters = target_parameters)
end

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

function ParameterHomotopyProblem(
        F::AbstractVector{Expression},
        start_solutions;
        parameters = Variable[],
        variables = setdiff(variables(F), parameters),
        variable_ordering = variables,
        start_parameters,
        target_parameters,
    )
    sys = System(F, variables = variable_ordering, parameters = parameters)
    return ParameterHomotopyProblem(
        sys,
        start_solutions;
        start_parameters = start_parameters,
        target_parameters = target_parameters,
    )
end

function ParameterHomotopyProblem(
        F::AbstractVector{<:MP.AbstractPolynomial},
        start_solutions;
        parameters = similar(MP.variables(F), 0),
        variables = setdiff(MP.variables(F), parameters),
        variable_ordering = variables,
        start_parameters,
        target_parameters,
    )
    sys = System(F, variables = variable_ordering, parameters = parameters)
    return ParameterHomotopyProblem(
        sys,
        start_solutions;
        start_parameters = start_parameters,
        target_parameters = target_parameters,
    )
end

struct ParameterSweepProblem{
        SystemT,
        StartsT,
        PStartT,
        TargetsT,
        TransformParamsT,
        TransformResultT,
    } <: AbstractHCProblem
    system::SystemT
    start_solutions::StartsT
    start_parameters::PStartT
    targets::TargetsT
    transform_parameters::TransformParamsT
    transform_result::TransformResultT
    flatten::Bool
end

function ParameterSweepProblem(
        F::Union{System, AbstractSystem},
        start_solutions;
        start_parameters,
        targets,
        transform_parameters = identity,
        transform_result = tuple,
        flatten::Bool = false,
    )
    return ParameterSweepProblem{
        typeof(F),
        typeof(start_solutions),
        typeof(start_parameters),
        typeof(targets),
        typeof(transform_parameters),
        typeof(transform_result),
    }(
        F,
        start_solutions,
        start_parameters,
        targets,
        transform_parameters,
        transform_result,
        flatten,
    )
end

function ParameterSweepProblem(
        F::AbstractVector{Expression},
        start_solutions;
        parameters = Variable[],
        variables = setdiff(variables(F), parameters),
        variable_ordering = variables,
        start_parameters,
        targets,
        transform_parameters = identity,
        transform_result = tuple,
        flatten::Bool = false,
    )
    sys = System(F, variables = variable_ordering, parameters = parameters)
    return ParameterSweepProblem(
        sys,
        start_solutions;
        start_parameters = start_parameters,
        targets = targets,
        transform_parameters = transform_parameters,
        transform_result = transform_result,
        flatten = flatten,
    )
end

function ParameterSweepProblem(
        F::AbstractVector{<:MP.AbstractPolynomial},
        start_solutions;
        parameters = similar(MP.variables(F), 0),
        variables = setdiff(MP.variables(F), parameters),
        variable_ordering = variables,
        start_parameters,
        targets,
        transform_parameters = identity,
        transform_result = tuple,
        flatten::Bool = false,
    )
    sys = System(F, variables = variable_ordering, parameters = parameters)
    return ParameterSweepProblem(
        sys,
        start_solutions;
        start_parameters = start_parameters,
        targets = targets,
        transform_parameters = transform_parameters,
        transform_result = transform_result,
        flatten = flatten,
    )
end

struct HomotopyProblem{HomotopyT, StartsT} <: AbstractHCProblem
    homotopy::HomotopyT
    start_solutions::StartsT
end

HomotopyProblem(H::Union{Homotopy, AbstractHomotopy}, starts) =
    HomotopyProblem{typeof(H), typeof(starts)}(H, starts)
