module WarntypeAudit

    using HomotopyContinuation
    using InteractiveUtils
    using Random
    using Printf
    using LinearAlgebra

    const HC = HomotopyContinuation

    struct MethodTarget
        key::String
        module_name::String
        function_name::String
        signature::String
        source_file::String
        source_line::Int
        status::Symbol
        reason::String
    end

    struct FixtureSpec
        id::String
        module_name::String
        function_name::String
        hot_path::Bool
        build::Function
    end

    struct AuditFinding
        severity::Symbol
        module_name::String
        function_name::String
        signature::String
        source_file::String
        source_line::Int
        return_type::String
        body_summary::String
        evidence::String
        recommendation::String
        covered_by_inferred_tests::Bool
    end

    struct AuditRecord
        status::Symbol
        severity::Symbol
        fixture_id::String
        module_name::String
        function_name::String
        signature::String
        source_file::String
        source_line::Int
        return_type::String
        body_summary::String
        evidence::String
        hot_path::Bool
        error_message::Union{Nothing, String}
    end

    const SEVERITY_ORDER = Dict(:none => 0, :low => 1, :medium => 2, :high => 3)

    function _repo_root()
        return normpath(joinpath(@__DIR__, ".."))
    end

    function _src_root()
        return joinpath(_repo_root(), "src")
    end

    function _normalize_path(path)
        s = String(path)
        if isempty(s) || s == "none"
            return s
        end
        return isabspath(s) ? normpath(s) : normpath(joinpath(_repo_root(), s))
    end

    function _relpath_or_self(path)
        p = _normalize_path(path)
        root = _repo_root()
        if isempty(p) || p == "none"
            return p
        end
        return startswith(p, root) ? relpath(p, root) : p
    end

    function _method_key(m::Method)
        file = _normalize_path(m.file)
        sig = sprint(show, Base.unwrap_unionall(m.sig))
        return string(m.module, "|", m.name, "|", sig, "|", file, ":", m.line)
    end

    function _is_src_method(m::Method)
        file = _normalize_path(m.file)
        src = _src_root()
        return startswith(file, src)
    end

    function _candidate_modules(root::Module)
        mods = Module[root]
        seen = Set{Module}(mods)
        i = 1
        while i <= length(mods)
            m = mods[i]
            i += 1
            for sym in names(m; all = true, imported = true)
                isdefined(m, sym) || continue
                val = getfield(m, sym)
                if val isa Module
                    modname = string(val)
                    if startswith(modname, string(root)) && !(val in seen)
                        push!(seen, val)
                        push!(mods, val)
                    end
                end
            end
        end
        return mods
    end

    function collect_method_targets()
        targets = Dict{String, MethodTarget}()
        for mod in _candidate_modules(HC)
            for sym in names(mod; all = true, imported = true)
                isdefined(mod, sym) || continue
                val = getfield(mod, sym)
                val isa Function || continue
                for m in methods(val)
                    _is_src_method(m) || continue
                    key = _method_key(m)
                    haskey(targets, key) && continue
                    fname = String(m.name)
                    unsupported = startswith(fname, "#")
                    reason =
                        unsupported ? "generated/internal method name" : "no fixture registered"
                    targets[key] = MethodTarget(
                        key,
                        string(m.module),
                        fname,
                        sprint(show, Base.unwrap_unionall(m.sig)),
                        _relpath_or_self(m.file),
                        m.line,
                        unsupported ? :unsupported : :blocked,
                        reason,
                    )
                end
            end
        end
        return sort!(
            collect(values(targets)); by = t -> (t.module_name, t.function_name, t.signature)
        )
    end

    function _extract_inferred_targets()
        files = [
            joinpath(@__DIR__, "solve_test.jl"), joinpath(@__DIR__, "type_stability_test.jl")
        ]
        out = Set{String}()
        pat = r"@inferred\s+([A-Za-z_][A-Za-z0-9_\.!]*)"
        for f in files
            isfile(f) || continue
            txt = read(f, String)
            for m in eachmatch(pat, txt)
                token = m.captures[1]
                if occursin(".", token)
                    parts = split(token, ".")
                    push!(out, last(parts))
                else
                    push!(out, token)
                end
            end
        end
        return out
    end

    function build_context()
        Random.seed!(0x22223333)

        @var x y a b t

        base_system = System([x^2 + y^2 - 1, x - y])
        base_prob = SystemProblem(base_system)
        base_alg = TotalDegreeAlgorithm(; seed = UInt32(0x11223344))
        poly_alg = PolyhedralAlgorithm(; seed = UInt32(0x11224455))
        path_cache = HC.init(base_prob, base_alg; show_progress = false, threading = false)
        poly_cache = HC.init(base_prob, poly_alg; show_progress = false, threading = false)
        path_result = solve!(path_cache)

        tracker = path_cache.solver.trackers[1]
        starts = collect(path_cache.starts)
        start1 = first(starts)

        sweep_system = System([x^2 - a, x * y - a + b], [x, y], [a, b])
        sweep_starts = [[1.0, 1.0], [-1.0, -1.0]]
        sweep_targets = [[2.0, 4.0], [3.0, 5.0]]
        sweep_alg = PathTrackingAlgorithm(; seed = UInt32(0x55667788))

        sweep_prob_map = ParameterSweepProblem(
            sweep_system,
            sweep_starts;
            start_parameters = [1.0, 0.0],
            targets = sweep_targets,
            reducer = MapReducer((r, p) -> nsolutions(r)),
        )
        sweep_cache_map = HC.init(
            sweep_prob_map, sweep_alg; show_progress = false, threading = false
        )

        sweep_prob_flat = ParameterSweepProblem(
            sweep_system,
            sweep_starts;
            start_parameters = [1.0, 0.0],
            targets = sweep_targets,
            reducer = FlatMapReducer((r, p) -> real_solutions(r)),
        )
        sweep_cache_flat = HC.init(
            sweep_prob_flat, sweep_alg; show_progress = false, threading = false
        )
        path_start = ComplexF64[1.0 + 0im, 1.0 + 0im]
        param_prob = ParameterHomotopyProblem(
            sweep_system,
            [path_start];
            start_parameters = [1.0, 0.0],
            target_parameters = [2.0, 4.0],
        )
        param_cache = HC.init(
            param_prob, sweep_alg; show_progress = false, threading = false
        )
        parameter_homotopy = ParameterHomotopy(sweep_system, [1.0, 0.0], [2.0, 4.0])
        homotopy_prob = HomotopyProblem(parameter_homotopy, [path_start])
        homotopy_cache = HC.init(
            homotopy_prob, sweep_alg; show_progress = false, threading = false
        )

        result_iter = solve(path_cache.solver, path_cache.starts, ResultIterator)

        stepper = HC.SegmentStepper(0.0, 2.0)
        weighted_norm = HC.WeightedNorm(HC.InfNorm(), 2)
        HC.init!(weighted_norm, [1.0, 2.0])

        workspace = HC.MatrixWorkspace(2, 2)
        workspace[1, 1] = 2.0 + 0.0im
        workspace[1, 2] = 1.0 + 0.0im
        workspace[2, 1] = 1.0 + 0.0im
        workspace[2, 2] = 3.0 + 0.0im

        bvec = ComplexF64[1.0 + 0im, 2.0 + 0im]
        xvec = zeros(ComplexF64, 2)

        compiled_system = CompiledSystem(base_system)
        mixed_system = MixedSystem(base_system)
        interp_system = InterpretedSystem(base_system)
        symbolic_homotopy = Homotopy([x^2 + t * y - 1, y - t * x], [x, y], t)
        interp_homotopy = InterpretedHomotopy(symbolic_homotopy)
        compiled_homotopy = CompiledHomotopy(symbolic_homotopy)
        mixed_homotopy = MixedHomotopy(symbolic_homotopy)
        eval_x = ComplexF64[0.25 + 0im, 0.75 + 0im]
        eval_r = zeros(ComplexF64, 2)
        eval_j = zeros(ComplexF64, 2, 2)
        eval_t = 0.5 + 0im

        newton_result = HC.NewtonResult(
            HC.NewtonCode.success, ComplexF64[1.0 + 0im], 1.0, 1.0, 1, 0.1
        )
        corrector_result = HC.NewtonCorrectorResult(
            HC.NEWT_CONVERGED, 1.0, 1, 1.0, 1.0, 1.0, 1.0
        )

        return (
            base_prob = base_prob,
            base_alg = base_alg,
            poly_alg = poly_alg,
            base_system = base_system,
            path_cache = path_cache,
            poly_cache = poly_cache,
            path_result = path_result,
            tracker = tracker,
            starts = starts,
            start1 = start1,
            sweep_system = sweep_system,
            sweep_alg = sweep_alg,
            sweep_prob_map = sweep_prob_map,
            sweep_prob_flat = sweep_prob_flat,
            sweep_cache_map = sweep_cache_map,
            sweep_cache_flat = sweep_cache_flat,
            param_prob = param_prob,
            param_cache = param_cache,
            homotopy_prob = homotopy_prob,
            homotopy_cache = homotopy_cache,
            result_iter = result_iter,
            stepper = stepper,
            weighted_norm = weighted_norm,
            workspace = workspace,
            bvec = bvec,
            xvec = xvec,
            compiled_system = compiled_system,
            mixed_system = mixed_system,
            interp_system = interp_system,
            symbolic_homotopy = symbolic_homotopy,
            interp_homotopy = interp_homotopy,
            compiled_homotopy = compiled_homotopy,
            mixed_homotopy = mixed_homotopy,
            eval_x = eval_x,
            eval_r = eval_r,
            eval_j = eval_j,
            eval_t = eval_t,
            newton_result = newton_result,
            corrector_result = corrector_result,
        )
    end

    function default_fixture_specs()
        return FixtureSpec[
            FixtureSpec(
                "utils.fast_abs.complex",
                "HomotopyContinuation",
                "fast_abs",
                false,
                ctx -> (f = HC.fast_abs, args = (1.0 + 2.0im,), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.fast_abs.real",
                "HomotopyContinuation",
                "fast_abs",
                false,
                ctx -> (f = HC.fast_abs, args = (-3.0,), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.nanmin",
                "HomotopyContinuation",
                "nanmin",
                false,
                ctx -> (f = HC.nanmin, args = (1.0, 2.0), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.nanmax",
                "HomotopyContinuation",
                "nanmax",
                false,
                ctx -> (f = HC.nanmax, args = (1.0, 2.0), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.always_false",
                "HomotopyContinuation",
                "always_false",
                false,
                ctx -> (f = HC.always_false, args = (1,), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.nthroot",
                "HomotopyContinuation",
                "nthroot",
                false,
                ctx -> (f = HC.nthroot, args = (8.0, 3), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.all2",
                "HomotopyContinuation",
                "all2",
                false,
                ctx -> (f = HC.all2, args = ((==), [1, 2], [1, 2]), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.unpack",
                "HomotopyContinuation",
                "unpack",
                false,
                ctx -> (f = HC.unpack, args = (nothing, 3), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.plural",
                "HomotopyContinuation",
                "plural",
                false,
                ctx -> (f = HC.plural, args = ("path", 2), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.init.stepper",
                "HomotopyContinuation",
                "init!",
                false,
                ctx -> (f = HC.init!, args = (ctx.stepper, 0.0, 2.0), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.propose_step",
                "HomotopyContinuation",
                "propose_step!",
                false,
                ctx -> (f = HC.propose_step!, args = (ctx.stepper, 0.5), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.step_success",
                "HomotopyContinuation",
                "step_success!",
                false,
                ctx -> (f = HC.step_success!, args = (ctx.stepper,), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.dist_to_target",
                "HomotopyContinuation",
                "dist_to_target",
                false,
                ctx -> (f = HC.dist_to_target, args = (ctx.stepper,), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.is_done",
                "HomotopyContinuation",
                "is_done",
                false,
                ctx -> (f = HC.is_done, args = (ctx.stepper,), kwargs = (;)),
            ),
            FixtureSpec(
                "utils.getproperty",
                "HomotopyContinuation",
                "getproperty",
                false,
                ctx -> (f = getproperty, args = (ctx.stepper, :t), kwargs = (;)),
            ),
            FixtureSpec(
                "norm.init",
                "HomotopyContinuation",
                "init!",
                false,
                ctx -> (f = HC.init!, args = (ctx.weighted_norm, [1.0, 2.0]), kwargs = (;)),
            ),
            FixtureSpec(
                "norm.update",
                "HomotopyContinuation",
                "update!",
                false,
                ctx -> (f = HC.update!, args = (ctx.weighted_norm, [1.5, 2.5]), kwargs = (;)),
            ),
            FixtureSpec(
                "norm.norm.inf",
                "HomotopyContinuation",
                "norm",
                false,
                ctx -> (f = HC.norm, args = ([1.0, 2.0], HC.InfNorm()), kwargs = (;)),
            ),
            FixtureSpec(
                "norm.distance.inf",
                "HomotopyContinuation",
                "distance",
                false,
                ctx -> (f = HC.distance, args = ([1.0, 2.0], [2.0, 3.0], HC.InfNorm()), kwargs = (;)),
            ),
            FixtureSpec(
                "norm.norm.weighted",
                "HomotopyContinuation",
                "norm",
                false,
                ctx -> (f = HC.norm, args = ([1.0, 2.0], ctx.weighted_norm), kwargs = (;)),
            ),
            FixtureSpec(
                "norm.distance.weighted",
                "HomotopyContinuation",
                "distance",
                false,
                ctx -> (
                    f = HC.distance,
                    args = ([1.0, 2.0], [2.0, 3.0], ctx.weighted_norm),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "la.workspace",
                "HomotopyContinuation",
                "MatrixWorkspace",
                true,
                ctx -> (f = HC.MatrixWorkspace, args = (2, 2), kwargs = (;)),
            ),
            FixtureSpec(
                "la.updated",
                "HomotopyContinuation",
                "updated!",
                true,
                ctx -> (f = HC.updated!, args = (ctx.workspace,), kwargs = (;)),
            ),
            FixtureSpec(
                "la.matrix",
                "HomotopyContinuation",
                "matrix",
                true,
                ctx -> (f = HC.matrix, args = (ctx.workspace,), kwargs = (;)),
            ),
            FixtureSpec(
                "la.ldiv",
                "HomotopyContinuation",
                "ldiv!",
                true,
                ctx -> (
                    f = LinearAlgebra.ldiv!,
                    args = (ctx.xvec, ctx.workspace, ctx.bvec),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "solve.init",
                "HomotopyContinuation",
                "init",
                true,
                ctx -> (
                    f = HC.init,
                    args = (ctx.base_prob, ctx.base_alg),
                    kwargs = (; show_progress = false, threading = false),
                ),
            ),
            FixtureSpec(
                "solve.init.poly",
                "HomotopyContinuation",
                "init",
                true,
                ctx -> (
                    f = HC.init,
                    args = (ctx.base_prob, ctx.poly_alg),
                    kwargs = (; show_progress = false, threading = false),
                ),
            ),
            FixtureSpec(
                "solve.init.parameter_homotopy",
                "HomotopyContinuation",
                "init",
                true,
                ctx -> (
                    f = HC.init,
                    args = (ctx.param_prob, ctx.sweep_alg),
                    kwargs = (; show_progress = false, threading = false),
                ),
            ),
            FixtureSpec(
                "solve.init.homotopy",
                "HomotopyContinuation",
                "init",
                true,
                ctx -> (
                    f = HC.init,
                    args = (ctx.homotopy_prob, ctx.sweep_alg),
                    kwargs = (; show_progress = false, threading = false),
                ),
            ),
            FixtureSpec(
                "solve.init.sweep",
                "HomotopyContinuation",
                "init",
                true,
                ctx -> (
                    f = HC.init,
                    args = (ctx.sweep_prob_map, ctx.sweep_alg),
                    kwargs = (; show_progress = false, threading = false),
                ),
            ),
            FixtureSpec(
                "solve.init.iter.total_degree",
                "HomotopyContinuation",
                "init",
                true,
                ctx -> (
                    f = HC.init,
                    args = (ctx.base_prob, ctx.base_alg, ResultIterator),
                    kwargs = (; show_progress = false),
                ),
            ),
            FixtureSpec(
                "solve.init.iter.poly",
                "HomotopyContinuation",
                "init",
                true,
                ctx -> (
                    f = HC.init,
                    args = (ctx.base_prob, ctx.poly_alg, ResultIterator),
                    kwargs = (; show_progress = false),
                ),
            ),
            FixtureSpec(
                "solve.init.iter.parameter_homotopy",
                "HomotopyContinuation",
                "init",
                true,
                ctx -> (
                    f = HC.init,
                    args = (ctx.param_prob, ctx.sweep_alg, ResultIterator),
                    kwargs = (; show_progress = false),
                ),
            ),
            FixtureSpec(
                "solve.init.iter.homotopy",
                "HomotopyContinuation",
                "init",
                true,
                ctx -> (
                    f = HC.init,
                    args = (ctx.homotopy_prob, ctx.sweep_alg, ResultIterator),
                    kwargs = (; show_progress = false),
                ),
            ),
            FixtureSpec(
                "solve.init.iter.sweep",
                "HomotopyContinuation",
                "init",
                true,
                ctx -> (
                    f = HC.init,
                    args = (ctx.sweep_prob_map, ctx.sweep_alg, ResultIterator),
                    kwargs = (; show_progress = false),
                ),
            ),
            FixtureSpec(
                "solve.start_parameters!",
                "HomotopyContinuation",
                "start_parameters!",
                false,
                ctx -> (
                    f = HC.start_parameters!,
                    args = (ctx.param_cache.solver, [1.0, 0.0]),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "solve.target_parameters!",
                "HomotopyContinuation",
                "target_parameters!",
                false,
                ctx -> (
                    f = HC.target_parameters!,
                    args = (ctx.param_cache.solver, [2.0, 4.0]),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "solve.solve!",
                "CommonSolve",
                "solve!",
                true,
                ctx -> (f = (solve!), args = (ctx.path_cache,), kwargs = (;)),
            ),
            FixtureSpec(
                "solve.solve",
                "CommonSolve",
                "solve",
                true,
                ctx -> (
                    f = solve,
                    args = (ctx.base_prob, ctx.base_alg),
                    kwargs = (; show_progress = false, threading = false),
                ),
            ),
            FixtureSpec(
                "solve.paths_to_track",
                "HomotopyContinuation",
                "paths_to_track",
                true,
                ctx -> (f = HC.paths_to_track, args = (ctx.base_prob, ctx.base_alg), kwargs = (;)),
            ),
            FixtureSpec(
                "solve.cache_kwargs",
                "HomotopyContinuation",
                "_cache_solve_kwargs",
                false,
                ctx -> (f = HC._cache_solve_kwargs, args = (ctx.path_cache,), kwargs = (;)),
            ),
            FixtureSpec(
                "solve.iter",
                "CommonSolve",
                "solve",
                true,
                ctx -> (
                    f = solve,
                    args = (ctx.path_cache.solver, ctx.path_cache.starts, ResultIterator),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "solve.sweep.solve!",
                "CommonSolve",
                "solve!",
                true,
                ctx -> (f = (solve!), args = (ctx.sweep_cache_map,), kwargs = (;)),
            ),
            FixtureSpec(
                "solve.sweep.solveflat",
                "CommonSolve",
                "solve!",
                true,
                ctx -> (f = (solve!), args = (ctx.sweep_cache_flat,), kwargs = (;)),
            ),
            FixtureSpec(
                "tracking.track",
                "HomotopyContinuation",
                "track",
                true,
                ctx -> (f = HC.track, args = (ctx.tracker, ctx.start1), kwargs = (;)),
            ),
            FixtureSpec(
                "result.seed",
                "HomotopyContinuation",
                "seed",
                false,
                ctx -> (f = HC.seed, args = (ctx.path_result,), kwargs = (;)),
            ),
            FixtureSpec(
                "result.algorithm",
                "HomotopyContinuation",
                "algorithm",
                false,
                ctx -> (f = HC.algorithm, args = (ctx.path_result,), kwargs = (;)),
            ),
            FixtureSpec(
                "result.path_results",
                "HomotopyContinuation",
                "path_results",
                false,
                ctx -> (f = HC.path_results, args = (ctx.path_result,), kwargs = (;)),
            ),
            FixtureSpec(
                "result.nresults",
                "HomotopyContinuation",
                "nresults",
                false,
                ctx -> (f = HC.nresults, args = (ctx.path_result,), kwargs = (;)),
            ),
            FixtureSpec(
                "result.nsolutions",
                "HomotopyContinuation",
                "nsolutions",
                false,
                ctx -> (f = HC.nsolutions, args = (ctx.path_result,), kwargs = (;)),
            ),
            FixtureSpec(
                "result.real_solutions",
                "HomotopyContinuation",
                "real_solutions",
                false,
                ctx -> (f = HC.real_solutions, args = (ctx.path_result,), kwargs = (;)),
            ),
            FixtureSpec(
                "result.nreal",
                "HomotopyContinuation",
                "nreal",
                false,
                ctx -> (f = HC.nreal, args = (ctx.path_result,), kwargs = (;)),
            ),
            FixtureSpec(
                "result.ntracked",
                "HomotopyContinuation",
                "ntracked",
                false,
                ctx -> (f = HC.ntracked, args = (ctx.path_result,), kwargs = (;)),
            ),
            FixtureSpec(
                "result.length_iter",
                "Base",
                "length",
                false,
                ctx -> (f = length, args = (ctx.result_iter,), kwargs = (;)),
            ),
            FixtureSpec(
                "result.bitmask",
                "HomotopyContinuation",
                "bitmask",
                false,
                ctx -> (f = HC.bitmask, args = (HC.is_success, ctx.result_iter), kwargs = (;)),
            ),
            FixtureSpec(
                "result.bitmask_filter",
                "HomotopyContinuation",
                "bitmask_filter",
                false,
                ctx -> (f = HC.bitmask_filter, args = (HC.is_success, ctx.result_iter), kwargs = (;)),
            ),
            FixtureSpec(
                "result.trace",
                "HomotopyContinuation",
                "trace",
                false,
                ctx -> (f = HC.trace, args = (ctx.result_iter,), kwargs = (;)),
            ),
            FixtureSpec(
                "result.Result",
                "HomotopyContinuation",
                "Result",
                false,
                ctx -> (f = HC.Result, args = (ctx.result_iter,), kwargs = (;)),
            ),
            FixtureSpec(
                "result.result_filter_options",
                "HomotopyContinuation",
                "_result_filter_options",
                false,
                ctx -> (f = HC._result_filter_options, args = (), kwargs = (;)),
            ),
            FixtureSpec(
                "newton.is_success",
                "HomotopyContinuation",
                "is_success",
                false,
                ctx -> (f = HC.is_success, args = (ctx.newton_result,), kwargs = (;)),
            ),
            FixtureSpec(
                "newton.solution",
                "HomotopyContinuation",
                "solution",
                false,
                ctx -> (f = HC.solution, args = (ctx.newton_result,), kwargs = (;)),
            ),
            FixtureSpec(
                "corrector.is_converged",
                "HomotopyContinuation",
                "is_converged",
                false,
                ctx -> (f = HC.is_converged, args = (ctx.corrector_result,), kwargs = (;)),
            ),
            FixtureSpec(
                "predictor.t_to_s",
                "HomotopyContinuation",
                "t_to_s_plane",
                false,
                ctx -> (f = HC.t_to_s_plane, args = (1.0 + 1.0im, 2), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.evaluate!",
                "HomotopyContinuation.ModelKit",
                "evaluate!",
                true,
                ctx -> (
                    f = (evaluate!),
                    args = (ctx.eval_r, ctx.interp_system, ctx.eval_x),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate_and_jacobian!",
                "HomotopyContinuation.ModelKit",
                "evaluate_and_jacobian!",
                true,
                ctx -> (
                    f = (evaluate_and_jacobian!),
                    args = (ctx.eval_r, ctx.eval_j, ctx.interp_system, ctx.eval_x),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.jacobian!",
                "HomotopyContinuation.ModelKit",
                "jacobian!",
                true,
                ctx -> (
                    f = (jacobian!),
                    args = (ctx.eval_j, ctx.interp_system, ctx.eval_x),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate!.compiled_system",
                "HomotopyContinuation.ModelKit",
                "evaluate!",
                true,
                ctx -> (
                    f = (evaluate!),
                    args = (ctx.eval_r, ctx.compiled_system, ctx.eval_x),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate_and_jacobian!.compiled_system",
                "HomotopyContinuation.ModelKit",
                "evaluate_and_jacobian!",
                true,
                ctx -> (
                    f = (evaluate_and_jacobian!),
                    args = (ctx.eval_r, ctx.eval_j, ctx.compiled_system, ctx.eval_x),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.jacobian!.compiled_system",
                "HomotopyContinuation.ModelKit",
                "jacobian!",
                true,
                ctx -> (
                    f = (jacobian!),
                    args = (ctx.eval_j, ctx.compiled_system, ctx.eval_x),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate!.mixed_system",
                "HomotopyContinuation.ModelKit",
                "evaluate!",
                true,
                ctx -> (
                    f = (evaluate!),
                    args = (ctx.eval_r, ctx.mixed_system, ctx.eval_x),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate_and_jacobian!.mixed_system",
                "HomotopyContinuation.ModelKit",
                "evaluate_and_jacobian!",
                true,
                ctx -> (
                    f = (evaluate_and_jacobian!),
                    args = (ctx.eval_r, ctx.eval_j, ctx.mixed_system, ctx.eval_x),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate!.interpreted_homotopy",
                "HomotopyContinuation.ModelKit",
                "evaluate!",
                true,
                ctx -> (
                    f = (evaluate!),
                    args = (ctx.eval_r, ctx.interp_homotopy, ctx.eval_x, ctx.eval_t),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate_and_jacobian!.interpreted_homotopy",
                "HomotopyContinuation.ModelKit",
                "evaluate_and_jacobian!",
                true,
                ctx -> (
                    f = (evaluate_and_jacobian!),
                    args = (
                        ctx.eval_r,
                        ctx.eval_j,
                        ctx.interp_homotopy,
                        ctx.eval_x,
                        ctx.eval_t,
                    ),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.jacobian!.interpreted_homotopy",
                "HomotopyContinuation.ModelKit",
                "jacobian!",
                true,
                ctx -> (
                    f = (jacobian!),
                    args = (ctx.eval_j, ctx.interp_homotopy, ctx.eval_x, ctx.eval_t),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate!.compiled_homotopy",
                "HomotopyContinuation.ModelKit",
                "evaluate!",
                true,
                ctx -> (
                    f = (evaluate!),
                    args = (ctx.eval_r, ctx.compiled_homotopy, ctx.eval_x, ctx.eval_t),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate_and_jacobian!.compiled_homotopy",
                "HomotopyContinuation.ModelKit",
                "evaluate_and_jacobian!",
                true,
                ctx -> (
                    f = (evaluate_and_jacobian!),
                    args = (
                        ctx.eval_r,
                        ctx.eval_j,
                        ctx.compiled_homotopy,
                        ctx.eval_x,
                        ctx.eval_t,
                    ),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.jacobian!.compiled_homotopy",
                "HomotopyContinuation.ModelKit",
                "jacobian!",
                true,
                ctx -> (
                    f = (jacobian!),
                    args = (ctx.eval_j, ctx.compiled_homotopy, ctx.eval_x, ctx.eval_t),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate!.mixed_homotopy",
                "HomotopyContinuation.ModelKit",
                "evaluate!",
                true,
                ctx -> (
                    f = (evaluate!),
                    args = (ctx.eval_r, ctx.mixed_homotopy, ctx.eval_x, ctx.eval_t),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate_and_jacobian!.mixed_homotopy",
                "HomotopyContinuation.ModelKit",
                "evaluate_and_jacobian!",
                true,
                ctx -> (
                    f = (evaluate_and_jacobian!),
                    args = (
                        ctx.eval_r,
                        ctx.eval_j,
                        ctx.mixed_homotopy,
                        ctx.eval_x,
                        ctx.eval_t,
                    ),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate.system",
                "HomotopyContinuation.ModelKit",
                "evaluate",
                false,
                ctx -> (
                    f = HC.ModelKit.evaluate,
                    args = (ctx.base_system, ctx.eval_x),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.evaluate.homotopy",
                "HomotopyContinuation.ModelKit",
                "evaluate",
                false,
                ctx -> (
                    f = HC.ModelKit.evaluate,
                    args = (ctx.symbolic_homotopy, ctx.eval_x, ctx.eval_t),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "modelkit.variables.system",
                "HomotopyContinuation.ModelKit",
                "variables",
                false,
                ctx -> (f = variables, args = (ctx.base_system,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.variables.interpreted_system",
                "HomotopyContinuation.ModelKit",
                "variables",
                false,
                ctx -> (f = variables, args = (ctx.interp_system,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.variables.compiled_system",
                "HomotopyContinuation.ModelKit",
                "variables",
                false,
                ctx -> (f = variables, args = (ctx.compiled_system,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.variables.homotopy",
                "HomotopyContinuation.ModelKit",
                "variables",
                false,
                ctx -> (f = variables, args = (ctx.symbolic_homotopy,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.variables.interpreted_homotopy",
                "HomotopyContinuation.ModelKit",
                "variables",
                false,
                ctx -> (f = variables, args = (ctx.interp_homotopy,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.nvariables.system",
                "HomotopyContinuation.ModelKit",
                "nvariables",
                false,
                ctx -> (f = nvariables, args = (ctx.base_system,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.nvariables.interpreted_system",
                "HomotopyContinuation.ModelKit",
                "nvariables",
                false,
                ctx -> (f = nvariables, args = (ctx.interp_system,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.nvariables.compiled_system",
                "HomotopyContinuation.ModelKit",
                "nvariables",
                false,
                ctx -> (f = nvariables, args = (ctx.compiled_system,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.nvariables.homotopy",
                "HomotopyContinuation.ModelKit",
                "nvariables",
                false,
                ctx -> (f = nvariables, args = (ctx.symbolic_homotopy,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.nvariables.interpreted_homotopy",
                "HomotopyContinuation.ModelKit",
                "nvariables",
                false,
                ctx -> (f = nvariables, args = (ctx.interp_homotopy,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.nvariables.compiled_homotopy",
                "HomotopyContinuation.ModelKit",
                "nvariables",
                false,
                ctx -> (f = nvariables, args = (ctx.compiled_homotopy,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.nparameters.system",
                "HomotopyContinuation.ModelKit",
                "nparameters",
                false,
                ctx -> (f = nparameters, args = (ctx.sweep_system,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.nparameters.homotopy",
                "HomotopyContinuation.ModelKit",
                "nparameters",
                false,
                ctx -> (f = nparameters, args = (ctx.symbolic_homotopy,), kwargs = (;)),
            ),
            FixtureSpec(
                "modelkit.is_homogeneous",
                "HomotopyContinuation.ModelKit",
                "is_homogeneous",
                false,
                ctx -> (f = is_homogeneous, args = (ctx.base_system,), kwargs = (;)),
            ),
            FixtureSpec(
                "systems.fixed",
                "HomotopyContinuation",
                "fixed",
                false,
                ctx -> (f = fixed, args = (ctx.base_system,), kwargs = (;)),
            ),
            FixtureSpec(
                "systems.fixed.compile_all",
                "HomotopyContinuation",
                "fixed",
                false,
                ctx -> (f = fixed, args = (ctx.base_system, CompileAll()), kwargs = (;)),
            ),
            FixtureSpec(
                "systems.fixed.compile_none",
                "HomotopyContinuation",
                "fixed",
                false,
                ctx -> (f = fixed, args = (ctx.base_system, CompileNone()), kwargs = (;)),
            ),
            FixtureSpec(
                "systems.fixed.compile_mixed",
                "HomotopyContinuation",
                "fixed",
                false,
                ctx -> (f = fixed, args = (ctx.base_system, CompileMixed()), kwargs = (;)),
            ),
            FixtureSpec(
                "systems.fixed.abstract",
                "HomotopyContinuation",
                "fixed",
                false,
                ctx -> (
                    f = fixed,
                    args = (ctx.interp_system,),
                    kwargs = (; compile_mode = CompileAll()),
                ),
            ),
            FixtureSpec(
                "homotopies.fixed",
                "HomotopyContinuation",
                "fixed",
                false,
                ctx -> (f = fixed, args = (ctx.symbolic_homotopy,), kwargs = (;)),
            ),
            FixtureSpec(
                "homotopies.fixed.compile_all",
                "HomotopyContinuation",
                "fixed",
                false,
                ctx -> (
                    f = fixed,
                    args = (ctx.symbolic_homotopy, CompileAll()),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "homotopies.fixed.compile_none",
                "HomotopyContinuation",
                "fixed",
                false,
                ctx -> (
                    f = fixed,
                    args = (ctx.symbolic_homotopy, CompileNone()),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "homotopies.fixed.compile_mixed",
                "HomotopyContinuation",
                "fixed",
                false,
                ctx -> (
                    f = fixed,
                    args = (ctx.symbolic_homotopy, CompileMixed()),
                    kwargs = (;),
                ),
            ),
            FixtureSpec(
                "homotopies.fixed.abstract",
                "HomotopyContinuation",
                "fixed",
                false,
                ctx -> (
                    f = fixed,
                    args = (ctx.interp_homotopy,),
                    kwargs = (; compile_mode = CompileAll()),
                ),
            ),
        ]
    end

    function _extract_warntype(f, argtypes)
        io = IOContext(IOBuffer(), :color => false)
        code_warntype(io, f, argtypes; debuginfo = :none)
        raw = String(take!(io.io))
        return replace(raw, r"\\e\\[[0-9;]*m" => "")
    end

    function _match_any(line::AbstractString)
        return occursin(r"::Any(?:\\W|$)", line)
    end

    function _body_lines(wt::String)
        lines = split(wt, '\n')
        idx = findfirst(l -> startswith(l, "Body::"), lines)
        idx === nothing && return String[]
        return lines[idx:end]
    end

    function _body_summary(wt::String)
        lines = _body_lines(wt)
        has_any = any(_match_any, lines)
        has_union = any(l -> occursin("UNION{", l), lines)
        if has_any
            return "body contains Any"
        elseif has_union
            return "body contains union"
        else
            return "body concrete"
        end
    end

    function _is_small_concrete_union(T)
        T isa Union || return false
        members = Base.uniontypes(T)
        return !isempty(members) && length(members) <= 4 && all(isconcretetype, members)
    end

    function _is_broad_return_type(T)
        if T === Any || T isa UnionAll || T isa TypeVar
            return true
        elseif T isa DataType
            return !isconcretetype(T)
        elseif T isa Union
            return !_is_small_concrete_union(T)
        else
            return false
        end
    end

    function _is_dynamic_segment_stepper_getproperty(method::Method, rettype)
        if String(method.name) != "getproperty"
            return false
        end
        sig = string(Base.unwrap_unionall(method.sig))
        return occursin("SegmentStepper", sig) && _is_small_concrete_union(rettype)
    end

    function classify_record(method::Method, rettype, warntype_text::String, hot_path::Bool)
        if _is_dynamic_segment_stepper_getproperty(method, rettype)
            return :none, "suppressed intentional dynamic getproperty"
        end
        if _is_small_concrete_union(rettype)
            return :none, "suppressed finite concrete return union"
        end

        if _is_broad_return_type(rettype)
            return :high, "escaping broad/abstract return type"
        end

        body_lines = _body_lines(warntype_text)
        has_any = any(_match_any, body_lines)
        has_union = any(l -> occursin("UNION{", l), body_lines)

        if has_any && hot_path
            return :medium, "dynamic Any in hot-path body"
        elseif has_any || has_union
            return :low, "local non-concrete body types"
        else
            return :none, "concrete body"
        end
    end

    function _evidence_snippet(wt::String)
        lines = split(wt, '\n')
        selected = String[]
        for l in lines
            if startswith(l, "Body::")
                push!(selected, l)
            elseif _match_any(l) || occursin("UNION{", l)
                push!(selected, l)
            end
            length(selected) >= 4 && break
        end
        return isempty(selected) ? "Body concrete in sampled call." : join(selected, "\\n")
    end

    function _recommendation(severity::Symbol, function_name::String)
        if severity == :high
            return "Stabilize return type for `$(function_name)` by separating mixed-type branches or adding concretely-typed wrappers."
        elseif severity == :medium
            return "Reduce dynamic dispatch in `$(function_name)` hot path by concretizing intermediates and callable arguments."
        else
            return "Check whether local unions in `$(function_name)` can be narrowed with concrete temporaries."
        end
    end

    function _method_signature(method::Method)
        return sprint(show, Base.unwrap_unionall(method.sig))
    end

    _argument_type(arg) = arg isa Type ? Type{arg} : typeof(arg)

    function _record_for_error(spec::FixtureSpec, err)
        return AuditRecord(
            :blocked,
            :none,
            spec.id,
            spec.module_name,
            spec.function_name,
            "",
            "",
            0,
            "",
            "call failed",
            "",
            spec.hot_path,
            sprint(showerror, err),
        )
    end

    function audit_specs(ctx, specs::Vector{FixtureSpec})
        records = AuditRecord[]
        for spec in specs
            payload = try
                spec.build(ctx)
            catch err
                push!(records, _record_for_error(spec, err))
                continue
            end

            f = payload.f
            args = payload.args
            kwargs = payload.kwargs

            argtypes = Tuple{map(_argument_type, args)...}
            m = try
                which(f, argtypes)
            catch
                nothing
            end

            if m === nothing
                push!(
                    records,
                    AuditRecord(
                        :blocked,
                        :none,
                        spec.id,
                        spec.module_name,
                        spec.function_name,
                        "",
                        "",
                        0,
                        "",
                        "unable to resolve method",
                        "",
                        spec.hot_path,
                        "which failed",
                    ),
                )
                continue
            end

            try
                # Compile and run once with deterministic fixture.
                Base.invokelatest(f, args...; kwargs...)

                typed = code_typed(f, argtypes; optimize = true)
                isempty(typed) && error("empty code_typed output")
                ci, rettype = typed[1]
                _ = ci

                wt = _extract_warntype(f, argtypes)

                sev, _reason = classify_record(m, rettype, wt, spec.hot_path)

                status = sev == :none ? :auditable : :finding
                push!(
                    records,
                    AuditRecord(
                        status,
                        sev,
                        spec.id,
                        string(m.module),
                        String(m.name),
                        _method_signature(m),
                        _relpath_or_self(m.file),
                        m.line,
                        sprint(show, rettype),
                        _body_summary(wt),
                        _evidence_snippet(wt),
                        spec.hot_path,
                        nothing,
                    ),
                )
            catch err
                push!(
                    records,
                    AuditRecord(
                        :blocked,
                        :none,
                        spec.id,
                        string(m.module),
                        String(m.name),
                        _method_signature(m),
                        _relpath_or_self(m.file),
                        m.line,
                        "",
                        "call failed",
                        "",
                        spec.hot_path,
                        sprint(showerror, err),
                    ),
                )
            end
        end
        return records
    end

    function _target_status_map(targets::Vector{MethodTarget}, records::Vector{AuditRecord})
        status = Dict(t.key => t for t in targets)

        for r in records
            if isempty(r.signature) || isempty(r.source_file)
                continue
            end
            key = string(
                r.module_name,
                "|",
                r.function_name,
                "|",
                r.signature,
                "|",
                _normalize_path(r.source_file),
                ":",
                r.source_line,
            )
            haskey(status, key) || continue
            t = status[key]
            if r.status == :auditable || r.status == :finding
                status[key] = MethodTarget(
                    t.key,
                    t.module_name,
                    t.function_name,
                    t.signature,
                    t.source_file,
                    t.source_line,
                    :auditable,
                    "covered by fixture: $(r.fixture_id)",
                )
            end
        end

        return sort!(
            collect(values(status)); by = t -> (t.module_name, t.function_name, t.signature)
        )
    end

    function build_findings(records::Vector{AuditRecord}, inferred_targets::Set{String})
        by_method = Dict{Tuple{String, String, String, String, Int}, AuditFinding}()
        for r in records
            r.status == :finding || continue
            key = (r.module_name, r.function_name, r.signature, r.source_file, r.source_line)
            candidate = AuditFinding(
                r.severity,
                r.module_name,
                r.function_name,
                r.signature,
                r.source_file,
                r.source_line,
                r.return_type,
                r.body_summary,
                r.evidence,
                _recommendation(r.severity, r.function_name),
                r.function_name in inferred_targets,
            )

            if !haskey(by_method, key) ||
                    SEVERITY_ORDER[candidate.severity] > SEVERITY_ORDER[by_method[key].severity]
                by_method[key] = candidate
            end
        end

        findings = collect(values(by_method))
        order = Dict(:high => 1, :medium => 2, :low => 3)
        sort!(
            findings; by = f -> (order[f.severity], f.module_name, f.function_name, f.signature)
        )
        return findings
    end

    function _severity_counts(findings)
        d = Dict(:high => 0, :medium => 0, :low => 0)
        for f in findings
            d[f.severity] += 1
        end
        return d
    end

    function _module_counts(findings)
        out = Dict{String, Int}()
        for f in findings
            out[f.module_name] = get(out, f.module_name, 0) + 1
        end
        return sort!(collect(out); by = x -> (-x[2], x[1]))
    end

    function _status_counts(targets)
        d = Dict(:auditable => 0, :blocked => 0, :unsupported => 0)
        for t in targets
            d[t.status] = get(d, t.status, 0) + 1
        end
        return d
    end

    function _blocked_reasons(targets)
        out = Dict{String, Int}()
        for t in targets
            t.status == :blocked || t.status == :unsupported || continue
            out[t.reason] = get(out, t.reason, 0) + 1
        end
        return sort!(collect(out); by = x -> (-x[2], x[1]))
    end

    function _audit_digest(targets, findings)
        io = IOBuffer()
        println(io, "targets:", length(targets))
        for t in targets
            println(io, t.key, "|", t.status)
        end
        println(io, "findings:", length(findings))
        for f in findings
            println(io, f.severity, "|", f.module_name, "|", f.function_name, "|", f.signature)
        end
        return String(take!(io))
    end

    function write_report(path::String, targets, records, findings)
        mkpath(dirname(path))
        sev = _severity_counts(findings)
        modcounts = _module_counts(findings)
        scounts = _status_counts(targets)
        reasons = _blocked_reasons(targets)

        io = IOBuffer()
        println(io, "# Exhaustive `@code_warntype` Audit Report")
        println(io)
        println(io, "Deterministic entry command:")
        println(io)
        println(io, "```bash")
        println(io, "julia --project=test test/warntype_audit.jl")
        println(io, "```")
        println(io)
        println(io, "## Executive Summary")
        println(io)
        println(io, "- Total discovered src methods: ", length(targets))
        println(io, "- Audited methods: ", scounts[:auditable])
        println(io, "- Blocked methods: ", scounts[:blocked])
        println(io, "- Unsupported methods: ", scounts[:unsupported])
        println(
            io,
            "- Findings (high / medium / low): ",
            sev[:high],
            " / ",
            sev[:medium],
            " / ",
            sev[:low],
        )
        println(io)

        println(io, "## Findings By Module")
        println(io)
        if isempty(modcounts)
            println(io, "No actionable findings in audited methods.")
        else
            for (m, c) in modcounts
                println(io, "- `", m, "`: ", c)
            end
        end
        println(io)

        println(io, "## Top Hotspots")
        println(io)
        if isempty(findings)
            println(io, "No actionable hotspots detected.")
        else
            topn = min(length(findings), 12)
            for i in 1:topn
                f = findings[i]
                println(io, "### ", i, ". `", f.module_name, ".", f.function_name, "`")
                println(io)
                println(io, "- Severity: `", f.severity, "`")
                println(io, "- Signature: `", f.signature, "`")
                println(io, "- Source: `", f.source_file, ":", f.source_line, "`")
                println(io, "- Return type: `", f.return_type, "`")
                println(io, "- Body summary: `", f.body_summary, "`")
                println(
                    io, "- Covered by existing `@inferred` tests: ", f.covered_by_inferred_tests
                )
                println(io, "- Suggested fix: ", f.recommendation)
                println(io)
                println(io, "```text")
                println(io, f.evidence)
                println(io, "```")
                println(io)
            end
        end

        println(io, "## Full Findings Table")
        println(io)
        println(
            io,
            "| Severity | Module | Function | Source | Return Type | Body Summary | Inferred Tests |",
        )
        println(io, "|---|---|---|---|---|---|---|")
        for f in findings
            src = string("`", f.source_file, ":", f.source_line, "`")
            println(
                io,
                "| `",
                f.severity,
                "` | `",
                f.module_name,
                "` | `",
                f.function_name,
                "` | ",
                src,
                " | `",
                replace(f.return_type, "|" => "\\|"),
                "` | `",
                f.body_summary,
                "` | ",
                f.covered_by_inferred_tests,
                " |",
            )
        end
        println(io)

        println(io, "## Coverage And Exclusions")
        println(io)
        println(io, "### Blocked / Unsupported Reason Counts")
        println(io)
        for (reason, count) in reasons
            println(io, "- ", reason, ": ", count)
        end
        println(io)

        println(io, "### Notes")
        println(io)
        println(
            io, "- Actionable criterion: broad escaping return types or `Any` in hot paths."
        )
        println(
            io,
            "- Suppressions: finite concrete unions and intentional dynamic `getproperty(::Symbol)` cases.",
        )
        println(
            io,
            "- Existing inference suites cross-referenced from `test/solve_test.jl` and `test/type_stability_test.jl`.",
        )

        write(path, String(take!(io)))
        return path
    end

    function run_audit_once(; subset = nothing)
        targets = collect_method_targets()
        ctx = build_context()
        specs = default_fixture_specs()
        if subset !== nothing
            wanted = Set(subset)
            specs = filter(
                s ->
                (s.module_name in wanted) ||
                    (s.function_name in wanted) ||
                    any(w -> occursin(w, s.id), wanted),
                specs,
            )
        end
        records = audit_specs(ctx, specs)
        status_targets = _target_status_map(targets, records)
        inferred = _extract_inferred_targets()
        findings = build_findings(records, inferred)
        return status_targets, records, findings
    end

    function _parse_args(args)
        report_path = joinpath(_repo_root(), "test", "reports", "warntype_audit_report.md")
        verify_repeat = 2
        subset = nothing
        for a in args
            if startswith(a, "--report=")
                report_path = split(a, "="; limit = 2)[2]
                if !isabspath(report_path)
                    report_path = normpath(joinpath(_repo_root(), report_path))
                end
            elseif startswith(a, "--verify-repeat=")
                verify_repeat = parse(Int, split(a, "="; limit = 2)[2])
            elseif startswith(a, "--subset=")
                subset = split(split(a, "="; limit = 2)[2], ",")
            end
        end
        return report_path, verify_repeat, subset
    end

    function main(args = ARGS)
        report_path, verify_repeat, subset = _parse_args(args)

        snapshots = Tuple{Vector{MethodTarget}, Vector{AuditRecord}, Vector{AuditFinding}}[]
        digests = String[]

        for _ in 1:verify_repeat
            targets, records, findings = run_audit_once(; subset)
            push!(snapshots, (targets, records, findings))
            push!(digests, _audit_digest(targets, findings))
        end

        if length(unique(digests)) != 1
            error("Audit is not deterministic across repeated runs.")
        end

        targets, records, findings = snapshots[end]
        out = write_report(report_path, targets, records, findings)

        println("Warntype audit complete")
        println("report: ", out)
        println("discovered_methods: ", length(targets))
        println("audited_methods: ", count(t -> t.status == :auditable, targets))
        println("findings: ", length(findings))
        return nothing
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    WarntypeAudit.main()
end
