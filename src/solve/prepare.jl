function _seed_rng(seed::UInt32)
    return Random.Xoshiro(seed)
end

function parameter_homotopy_tracker(
        F::Union{System, AbstractSystem};
        compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE,
        tracker_options::TrackerOptions = TrackerOptions(),
        endgame_options::EndgameOptions = EndgameOptions(),
        start_parameters = randn(ComplexF64, nparameters(F)),
        target_parameters = start_parameters,
    )
    H = parameter_homotopy(
        F;
        start_parameters = start_parameters,
        target_parameters = target_parameters,
        compile_mode = compile_mode,
    )
    return EndgameTracker(H; tracker_options = tracker_options, options = endgame_options)
end

function parameter_homotopy(
        F::Union{System, AbstractSystem};
        generic_parameters = randn(ComplexF64, nparameters(F)),
        p₁ = generic_parameters,
        start_parameters = p₁,
        p₀ = generic_parameters,
        target_parameters = p₀,
        compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE,
    )
    m, n = size(F)
    H = ParameterHomotopy(fixed(F; compile_mode = compile_mode), start_parameters, target_parameters)
    f = System(F)
    if is_homogeneous(f)
        throw(
            ArgumentError(
                "Homogeneous/projective systems are not supported in affine-only mode.",
            ),
        )
    end
    if m < n
        throw(FiniteException(n - m))
    elseif m > n
        throw(
            ArgumentError(
                "Only square systems are supported in this minimal build. Got $m equations, expected $n.",
            ),
        )
    end

    return H
end

function _system_solver_startsolutions(
        prob::SystemProblem,
        alg::PolyhedralAlgorithm;
        show_progress::Bool = true,
    )
    tracker, starts = polyhedral(
        prob.system;
        compile_mode = alg.compile_mode,
        target_parameters = prob.target_parameters,
        tracker_options = alg.tracker_options,
        endgame_options = alg.endgame_options,
        only_torus = alg.only_torus,
        only_non_zero = alg.only_non_zero,
        show_progress = show_progress,
        rng = _seed_rng(alg.seed),
    )
    return Solver(tracker, alg; seed = alg.seed), starts
end

function _system_solver_startsolutions(prob::SystemProblem, alg::TotalDegreeAlgorithm)
    tracker, starts = total_degree(
        prob.system;
        compile_mode = alg.compile_mode,
        target_parameters = prob.target_parameters,
        gamma = alg.gamma,
        tracker_options = alg.tracker_options,
        endgame_options = alg.endgame_options,
    )
    try
        first(starts)
    catch
        throw("The number of start solutions is zero (total degree is zero).")
    end
    return Solver(tracker, alg; seed = alg.seed), starts
end

function _parameter_solver_startsolutions(
        prob::ParameterHomotopyProblem,
        alg::PathTrackingAlgorithm,
    )
    tracker = parameter_homotopy_tracker(
        prob.system;
        start_parameters = prob.start_parameters,
        target_parameters = prob.target_parameters,
        compile_mode = alg.compile_mode,
        tracker_options = alg.tracker_options,
        endgame_options = alg.endgame_options,
    )
    return Solver(tracker, alg; seed = alg.seed), prob.start_solutions
end

function _first_target(targets)
    target_it = iterate(targets)
    target_it === nothing && throw(ArgumentError("`targets` must be non-empty."))
    return first(target_it)
end

function _sweep_solver_startsolutions(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm,
    )
    first_target = _first_target(prob.targets)
    return _sweep_solver(prob, alg, first_target), prob.start_solutions
end

function _sweep_solver(
        prob::ParameterSweepProblem,
        alg::PathTrackingAlgorithm,
        target_parameters,
    )
    tracker = parameter_homotopy_tracker(
        prob.system;
        start_parameters = prob.start_parameters,
        target_parameters = target_parameters,
        compile_mode = alg.compile_mode,
        tracker_options = alg.tracker_options,
        endgame_options = alg.endgame_options,
    )
    return Solver(tracker, alg; seed = alg.seed)
end

function _homotopy_solver_startsolutions(
        prob::HomotopyProblem,
        alg::PathTrackingAlgorithm,
    )
    tracker = EndgameTracker(
        fixed(prob.homotopy; compile_mode = alg.compile_mode);
        tracker_options = alg.tracker_options,
        options = alg.endgame_options,
    )
    return Solver(tracker, alg; seed = alg.seed), prob.start_solutions
end

function paths_to_track(prob::SystemProblem, alg::PolyhedralAlgorithm)
    f = prob.system isa System ? prob.system : System(prob.system)
    return paths_to_track(f, alg)
end

function paths_to_track(prob::SystemProblem, alg::TotalDegreeAlgorithm)
    f = prob.system isa System ? prob.system : System(prob.system)
    return paths_to_track(f, alg)
end
