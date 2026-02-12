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

function ParameterSweepProblem(
        F::AbstractVector{Expression},
        start_solutions;
        parameters = Variable[],
        variables = setdiff(variables(F), parameters),
        variable_ordering = variables,
        start_parameters,
        targets,
        reducer::AbstractSweepReducer = IdentityReducer(),
    )
    sys = System(F, variables = variable_ordering, parameters = parameters)
    return ParameterSweepProblem(
        sys,
        start_solutions;
        start_parameters = start_parameters,
        targets = targets,
        reducer = reducer,
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
        reducer::AbstractSweepReducer = IdentityReducer(),
    )
    sys = System(F, variables = variable_ordering, parameters = parameters)
    return ParameterSweepProblem(
        sys,
        start_solutions;
        start_parameters = start_parameters,
        targets = targets,
        reducer = reducer,
    )
end
