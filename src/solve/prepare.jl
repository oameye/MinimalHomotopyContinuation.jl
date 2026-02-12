_seed_rng(seed::UInt32) = Random.Xoshiro(seed)
_solver_from_tracker(tracker, alg::AbstractHCAlgorithm) = Solver(
    tracker, alg; seed = alg.seed
)

function _parameter_homotopy_tracker(
        F::Union{System, AbstractSystem},
        alg::PathTrackingAlgorithm;
        start_parameters,
        target_parameters = start_parameters,
    )
    _validate_affine_square_system(F)
    H = ParameterHomotopy(
        fixed(F; compile_mode = alg.compile_mode), start_parameters, target_parameters
    )
    return EndgameTracker(
        H; tracker_options = alg.tracker_options, options = alg.endgame_options
    )
end

function _first_target(targets)
    target_it = iterate(targets)
    target_it === nothing && throw(ArgumentError("`targets` must be non-empty."))
    return first(target_it)
end

function _sweep_solver(
        prob::ParameterSweepProblem, alg::PathTrackingAlgorithm, target_parameters
    )
    tracker = _parameter_homotopy_tracker(
        prob.system, alg; start_parameters = prob.start_parameters, target_parameters
    )
    return _solver_from_tracker(tracker, alg)
end

function _polyhedral_startsystem(
        prob::SystemProblem, alg::PolyhedralAlgorithm; show_progress::Bool = true
    )
    F = prob.system
    if prob.target_parameters !== nothing
        F = fix_parameters(F, prob.target_parameters; compile_mode = alg.compile_mode)
    end
    return _polyhedral_kernel(F, alg; show_progress, rng = _seed_rng(alg.seed))
end

function _total_degree_startsystem(
        prob::SystemProblem, alg::TotalDegreeAlgorithm; show_progress::Bool = true
    )
    F = prob.system
    if prob.target_parameters !== nothing
        F = fix_parameters(F, prob.target_parameters; compile_mode = alg.compile_mode)
    end
    return _total_degree_kernel(F, alg; show_progress)
end

function _solver_startsolutions(
        prob::SystemProblem, alg::PolyhedralAlgorithm; show_progress::Bool = true
    )
    tracker, starts = _polyhedral_startsystem(prob, alg; show_progress)
    return _solver_from_tracker(tracker, alg), starts
end

function _solver_startsolutions(
        prob::SystemProblem, alg::TotalDegreeAlgorithm; show_progress::Bool = true
    )
    tracker, starts = _total_degree_startsystem(prob, alg; show_progress)
    return _solver_from_tracker(tracker, alg), starts
end

function _solver_startsolutions(
        prob::ParameterHomotopyProblem, alg::PathTrackingAlgorithm; show_progress::Bool = true
    )
    _ = show_progress
    tracker = _parameter_homotopy_tracker(
        prob.system,
        alg;
        start_parameters = prob.start_parameters,
        target_parameters = prob.target_parameters,
    )
    return _solver_from_tracker(tracker, alg), prob.start_solutions
end

function _solver_startsolutions(
        prob::HomotopyProblem, alg::PathTrackingAlgorithm; show_progress::Bool = true
    )
    _ = show_progress
    tracker = EndgameTracker(
        fixed(prob.homotopy; compile_mode = alg.compile_mode);
        tracker_options = alg.tracker_options,
        options = alg.endgame_options,
    )
    return _solver_from_tracker(tracker, alg), prob.start_solutions
end

function _solver_startsolutions(
        prob::ParameterSweepProblem, alg::PathTrackingAlgorithm; show_progress::Bool = true
    )
    _ = show_progress
    return _sweep_solver(prob, alg, _first_target(prob.targets)), prob.start_solutions
end

function _solver_startsolutions(
        prob::AbstractHCProblem, alg::AbstractHCAlgorithm; show_progress::Bool = true
    )
    _ = show_progress
    throw(
        ArgumentError(
            "Unsupported problem/algorithm combination: $(typeof(prob)) with $(typeof(alg)).",
        ),
    )
end

function paths_to_track(prob::SystemProblem, alg::PolyhedralAlgorithm)
    f = prob.system isa System ? prob.system : System(prob.system)
    return _paths_to_track_polyhedral(f, alg)
end

function paths_to_track(prob::SystemProblem, alg::TotalDegreeAlgorithm)
    f = prob.system isa System ? prob.system : System(prob.system)
    return _paths_to_track_total_degree(f, alg)
end
