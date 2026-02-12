"""
    PolyhedralStartSolutionsIterator

An iterator providing start solutions for the polyhedral homotopy.
"""
struct PolyhedralStartSolutionsIterator{Iter}
    support::Vector{Matrix{Int32}}
    start_coefficients::Vector{Vector{ComplexF64}}
    lifting::Vector{Vector{Int32}}
    mixed_cells::Iter
    BSS::BinomialSystemSolver
    show_progress::Bool
end

function PolyhedralStartSolutionsIterator(
        support::AbstractVector{<:AbstractMatrix{<:Integer}},
        coeffs::AbstractVector{<:AbstractVector{<:Number}},
        lifting = map(c -> zeros(Int32, length(c)), coeffs),
        mixed_cells = MixedCell[],
        ;
        show_progress::Bool = true,
    )
    BSS = BinomialSystemSolver(length(support))
    return PolyhedralStartSolutionsIterator(
        convert(Vector{Matrix{Int32}}, support),
        coeffs,
        lifting,
        mixed_cells,
        BSS,
        show_progress,
    )
end

Base.show(io::IO, C::PolyhedralStartSolutionsIterator) = print(
    io, "PolyhedralStartSolutionsIterator()"
)
Base.show(
    io::IO, ::MIME"application/prs.juno.inline", x::PolyhedralStartSolutionsIterator
) = x
Base.IteratorSize(::Type{<:PolyhedralStartSolutionsIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:PolyhedralStartSolutionsIterator}) = Base.HasEltype()
Base.eltype(iter::PolyhedralStartSolutionsIterator) = Tuple{MixedCell, Vector{ComplexF64}}

function compute_mixed_cells!(iter::PolyhedralStartSolutionsIterator)
    first_cell = iterate(iter.mixed_cells)
    if isnothing(first_cell)
        res = MixedSubdivisions.fine_mixed_cells(
            iter.support; show_progress = iter.show_progress
        )
        if isnothing(res) || isempty(res[1])
            throw(OverflowError("Cannot compute a start system."))
        end
        mixed_cells, lifting = res
        empty!(iter.mixed_cells)
        append!(iter.mixed_cells, mixed_cells)

        for (i, w) in enumerate(lifting)
            empty!(iter.lifting[i])
            append!(iter.lifting[i], w)
        end
    end
    return iter
end


function Base.iterate(iter::PolyhedralStartSolutionsIterator)
    compute_mixed_cells!(iter)

    first_cell = iterate(iter.mixed_cells)
    isnothing(first_cell) && return nothing
    cell, inner_state = first_cell

    solve(iter.BSS, iter.support, iter.start_coefficients, cell)
    x = [iter.BSS.X[i, 1] for i in 1:length(iter.support)]

    # return the value _and_ the combined state:
    # (cell, inner_state, column_index)
    return (cell, x), (cell, inner_state, 1)
end

function Base.iterate(iter::PolyhedralStartSolutionsIterator, state::Tuple{Any, Any, Int})
    cell, inner_state, j = state
    ncol = size(iter.BSS.X, 2)

    if j < ncol
        # we still have more columns to emit for the current cell
        new_j = j + 1
        x = [iter.BSS.X[i, new_j] for i in 1:length(iter.support)]
        return (cell, x), (cell, inner_state, new_j)
    else
        # we finished the last column for `cell`, advance to the next cell
        next_cell = iterate(iter.mixed_cells, inner_state)
        next_cell === nothing && return nothing
        cell2, new_inner_state = next_cell

        solve(iter.BSS, iter.support, iter.start_coefficients, cell2)
        x = [iter.BSS.X[i, 1] for i in 1:length(iter.support)]

        return (cell2, x), (cell2, new_inner_state, 1)
    end
end

# Tracker
function polyhedral_system(support)
    m = 0
    p = Variable[]
    coeffs = map(support) do A
        c = variables(:c, (m + 1):(m + size(A, 2)))
        m += size(A, 2)
        append!(p, c)
        c
    end

    return System(support, coeffs; variables = variables(:x, 1:length(support)), parameters = p)
end

"""
    PolyhedralTracker <: AbstractPathTracker

This tracker realises the two step approach of the polyhedral homotopy.
See also [`polyhedral`].
"""
struct PolyhedralTracker{H1 <: ToricHomotopy, H2 <: AbstractHomotopy, M} <: AbstractPathTracker
    toric_tracker::Tracker{H1, M}
    generic_tracker::EndgameTracker{H2, M}
    support::Vector{Matrix{Int32}}
    lifting::Vector{Vector{Int32}}
end

"""
Internal polyhedral start-system kernel used by solve preparation.
Returns `(tracker, starts)` for a fixed `SystemProblem + PolyhedralAlgorithm` setup.
"""
_polyhedral_rng(alg::PolyhedralAlgorithm, rng::Union{Nothing, Random.AbstractRNG}) = isnothing(rng) ? Random.Xoshiro(
    alg.seed
) : rng

function _polyhedral_kernel(
        F::AbstractSystem,
        alg::PolyhedralAlgorithm;
        show_progress::Bool = true,
        rng::Union{Nothing, Random.AbstractRNG} = nothing,
    )
    _, n = size(F)
    @var x[1:n]
    return _polyhedral_kernel(System(F(x), x), alg; show_progress, rng)
end

function _polyhedral_kernel(
        f::System,
        alg::PolyhedralAlgorithm;
        show_progress::Bool = true,
        rng::Union{Nothing, Random.AbstractRNG} = nothing,
    )
    _validate_affine_square_system(f; check_square = false)
    _validate_affine_square_system(f)
    support, target_coeffs = support_coefficients(f)
    return _polyhedral_kernel(support, target_coeffs, alg; show_progress, rng)
end

function has_zero_col(A)
    for j in 1:size(A, 2)
        has_non_zero = false
        for i in 1:size(A, 1)
            if !iszero(A[i, j])
                has_non_zero = true
                break
            end
        end
        has_non_zero || return true
    end
    return false
end

_paths_to_track_polyhedral(f::AbstractSystem, alg::PolyhedralAlgorithm) = _paths_to_track_polyhedral(
    System(f), alg
)
function _paths_to_track_polyhedral(f::System, alg::PolyhedralAlgorithm)
    only_non_zero = alg.only_non_zero
    _validate_affine_square_system(f)
    supp, _ = support_coefficients(f)

    return if only_non_zero
        MixedSubdivisions.mixed_volume(supp)
    else
        supp′ = map(A -> has_zero_col(A) ? A : [A zeros(Int32, size(A, 1), 1)], supp)
        MixedSubdivisions.mixed_volume(supp′)
    end
end
MixedSubdivisions.mixed_volume(f::Union{System, AbstractSystem}) = _paths_to_track_polyhedral(
    f, PolyhedralAlgorithm(only_torus = true)
)

function _polyhedral_kernel(
        support::AbstractVector{<:AbstractMatrix},
        target_coeffs::AbstractVector,
        alg::PolyhedralAlgorithm;
        show_progress::Bool = true,
        rng::Union{Nothing, Random.AbstractRNG} = nothing,
    )
    effective_rng = _polyhedral_rng(alg, rng)
    start_coeffs = map(
        c -> rand_approx_unit(length(c); rng = effective_rng) .* LA.norm(c, Inf), target_coeffs
    )
    return _polyhedral_kernel(
        support, start_coeffs, target_coeffs, alg; show_progress, rng = effective_rng
    )
end

"""
    rand_approx_unit(n::Integer)::Vector{ComplexF64}

This samples uniformly from the rectangle ``[0.9,1.1] × [0,2π]`` and transforms the sampled
values with the map ``(r, φ) ↦ r * e^{i φ}``.
"""
rand_approx_unit(n; rng::Random.AbstractRNG = Random.default_rng()) = map(
    _ -> (0.9 + 0.2 * rand(rng)) * cis2pi(rand(rng)), 1:n
)

cis2pi(x) = complex(cospi(2x), sinpi(2x))

function _polyhedral_kernel(
        support::AbstractVector{<:AbstractMatrix},
        start_coeffs::AbstractVector,
        target_coeffs::AbstractVector,
        alg::PolyhedralAlgorithm;
        show_progress::Bool = true,
        rng::Union{Nothing, Random.AbstractRNG} = nothing,
    )
    rng = _polyhedral_rng(alg, rng)
    only_non_zero = alg.only_non_zero
    size.(support, 2) == length.(start_coeffs) == length.(target_coeffs) ||
        throw(ArgumentError("Number of terms do not coincide."))
    if only_non_zero
        min_vecs = minimum.(support; dims = 2)
        support = map((A, v) -> A .- v, support, min_vecs)
    elseif !all(has_zero_col, support)
        # Add 0 to each support
        starts = Vector{ComplexF64}[]
        targets = Vector{ComplexF64}[]
        supp = eltype(support)[]
        for (i, A) in enumerate(support)
            if has_zero_col(A)
                push!(starts, start_coeffs[i])
                push!(targets, target_coeffs[i])
                push!(supp, A)
            else
                push!(starts, vcat(start_coeffs[i], randn(rng, ComplexF64)))
                push!(targets, vcat(target_coeffs[i], 0.0))
                push!(supp, [A zeros(Int32, size(A, 1), 1)])
            end
        end
        support = supp
        start_coeffs = starts
        target_coeffs = targets
    end

    F = fixed(polyhedral_system(support); compile_mode = alg.compile_mode)
    H₁ = ToricHomotopy(F, start_coeffs)
    toric_tracker = Tracker(H₁; options = alg.tracker_options)

    H₂ = begin
        p = reduce(append!, start_coeffs; init = ComplexF64[])
        q = reduce(append!, target_coeffs; init = ComplexF64[])
        CoefficientHomotopy(F, p, q)
    end
    generic_tracker = EndgameTracker(
        Tracker(H₂; options = alg.tracker_options); options = alg.endgame_options
    )

    S = PolyhedralStartSolutionsIterator(support, start_coeffs; show_progress)
    tracker = PolyhedralTracker(toric_tracker, generic_tracker, S.support, S.lifting)

    return tracker, S
end

function track(
        PT::PolyhedralTracker,
        start_solution::Tuple{MixedCell, <:AbstractVector{ComplexF64}};
        path_number::Union{Nothing, Int} = nothing,
        debug::Bool = false,
    )
    cell, x∞ = start_solution
    H = PT.toric_tracker.homotopy
    # The tracker works in two stages
    # 1) Revert the totric degeneration by tracking 0 to 1
    #    Here we use the following strategy
    #    a) Reparameterize the path such that lowest power of t occuring is 1
    #    b)
    #      i) If maximal power is less than 10 track from 0 to 1
    #      ii) If maximal power is larger than 10 we can get some problems for values close to 1
    #          (since t^k is << 1 for t < 1 and k large)
    #          Therefore split in two stages:
    #           1) track from 0 to t₀ = clamp(0.1^(10 / max_weight), 0.9, 0.985)
    #           2) Second, reperamerize path s.t. max power of t is 10. This shifts the
    #              problem towards 0 again (where we have more accuracy available.)
    #              And then track to 1
    # 2) Perform coefficient homotopy with possible endgame
    min_weight, max_weight = update_weights!(
        H, PT.support, PT.lifting, cell; min_weight = 1.0
    )

    if debug
        println("Min-Max weight: ", min_weight, ", ", max_weight)
    end
    if max_weight < 10
        retcode = track!(
            PT.toric_tracker,
            x∞,
            0.0,
            1.0;
            ω = 20.0,
            μ = 1.0e-12,
            max_initial_step_size = 0.2,
            debug = debug,
        )
        @unpack μ, ω = PT.toric_tracker.state
    else
        t₀ = clamp(0.1^(10 / max_weight), 0.9, 1 - 1.0e-6)
        retcode = track!(
            PT.toric_tracker,
            x∞,
            0.0,
            t₀;
            ω = 20.0,
            μ = 1.0e-12,
            max_initial_step_size = 0.2,
            debug = debug,
        )

        @unpack μ, ω = PT.toric_tracker.state
        if is_success(retcode)
            min_weight, max_weight = update_weights!(
                H, PT.support, PT.lifting, cell; max_weight = 10.0
            )
            t_restart = t₀^(1 / min_weight)
            # set min_step_size to 0.0 to not accidentally loose a solution
            min_step_size = PT.toric_tracker.options.min_step_size
            PT.toric_tracker.options.min_step_size = 0.0

            retcode = track!(
                PT.toric_tracker,
                PT.toric_tracker.state.x,
                t_restart,
                1.0;
                ω = ω,
                μ = μ,
                τ = 0.1 * t_restart,
                keep_steps = true,
                debug = debug,
            )
            @unpack μ, ω = PT.toric_tracker.state
            PT.toric_tracker.options.min_step_size = min_step_size
        end
    end
    if !is_success(retcode)
        state = PT.toric_tracker.state
        return PathResult(;
            return_code = PathResultCode.polyhedral_failed,
            solution = copy(state.x),
            start_solution = start_solution,
            t = real(state.t),
            accuracy = state.accuracy,
            singular = false,
            condition_jacobian = NaN,
            residual = NaN,
            winding_number = nothing,
            last_path_point = (copy(state.x), real(state.t)),
            valuation = nothing,
            ω = state.ω,
            μ = state.μ,
            path_number = path_number,
            extended_precision = state.extended_prec,
            accepted_steps = state.accepted_steps,
            rejected_steps = state.rejected_steps,
            extended_precision_used = state.used_extended_prec,
        )
    end

    r = track(
        PT.generic_tracker,
        PT.toric_tracker.state.x;
        start_solution = start_solution,
        # Don't provide ω since this can be misleading a lead to a too large initial step
        # ω = ω,
        μ = μ,
        path_number = path_number,
        debug = debug,
    )
    # report accurate steps
    r.accepted_steps += PT.toric_tracker.state.accepted_steps
    r.rejected_steps += PT.toric_tracker.state.rejected_steps

    return r
end
