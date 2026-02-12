"""
Internal total-degree start-system kernel used by solve preparation.
Returns `(tracker, starts)` for a fixed `SystemProblem + TotalDegreeAlgorithm` setup.
"""

function _total_degree_impl(F::Union{System, AbstractSystem}, alg::TotalDegreeAlgorithm)
    _, n = size(F)

    @unique_var x[1:n] s[1:n]
    F₀ = System(F(x), x)
    _validate_affine_square_system(F; homogeneous_system = F₀)
    support, coeffs = support_coefficients(F₀)

    D = zeros(Int, length(support))
    for (k, A) in enumerate(support)
        d = 0
        for j in 1:size(A, 2)
            dⱼ = 0
            for i in 1:size(A, 1)
                dⱼ += A[i, j]
            end
            d = max(d, dⱼ)
        end
        D[k] = d
    end
    scaling = maximum.(abs ∘ float, coeffs)

    if F isa System
        F = fixed(F; compile_mode = alg.compile_mode)
    end
    G = fixed(System(s .* (x .^ D .- 1), x, s); compile_mode = alg.compile_mode)
    G = fix_parameters(G, scaling)

    H = StraightLineHomotopy(G, F; γ = alg.gamma)
    T = EndgameTracker(H; tracker_options = alg.tracker_options, options = alg.endgame_options)
    starts = total_degree_start_solutions(D)

    return T, starts
end

function _total_degree_kernel(F::Union{System, AbstractSystem}, alg::TotalDegreeAlgorithm)
    tracker, starts = _total_degree_impl(F, alg)
    # Keep historical behavior from solve preparation: zero starts are rejected.
    if iterate(starts) === nothing
        throw(
            ArgumentError("The number of start solutions is zero (total degree is zero).")
        )
    end
    return tracker, starts
end


struct TotalDegreeStartSolutionsIterator{Iter}
    degrees::Vector{Int}
    iterator::Iter
end
function TotalDegreeStartSolutionsIterator(degrees)
    iterator = Iterators.product(map(d -> 0:(d - 1), degrees)...)
    return TotalDegreeStartSolutionsIterator(degrees, iterator)
end
function Base.show(io::IO, iter::TotalDegreeStartSolutionsIterator)
    return print(
        io, "$(length(iter)) total degree start solutions for degrees $(iter.degrees)"
    )
end

function Base.iterate(iter::TotalDegreeStartSolutionsIterator)
    it = iterate(iter.iterator)
    it === nothing && return nothing
    return _value(iter, first(it)), last(it)
end
function Base.iterate(iter::TotalDegreeStartSolutionsIterator, state)
    it = iterate(iter.iterator, state)
    it === nothing && return nothing
    return _value(iter, first(it)), last(it)
end

function _value(iter::TotalDegreeStartSolutionsIterator, indices)
    return map((k, d) -> cis(2π * k / d), indices, iter.degrees)
end


Base.length(iter::TotalDegreeStartSolutionsIterator) = length(iter.iterator)
Base.eltype(iter::Type{<:TotalDegreeStartSolutionsIterator}) = Vector{Complex{Float64}}

"""
    total_degree_start_solutions(degrees)

Returns an iterator of the start solutions of the system
```math
\\begin{array}{rl}
    z_1^{d_1} - 1 &= 0 \\\\
    z_2^{d_2} - 1 &= 0 \\\\
    &\\vdots \\\\
    z_n^{d_n} - 1 &= 0 \\\\
\\end{array}
```
where ``d_i`` is `degrees[i]`.

## Example
```julia-repl
julia> iter = total_degree_start_solutions([2, 2])
4 total degree start solutions for degrees [2, 2]

julia> collect(iter)
4-element Array{Array{Complex{Float64},1},1}:
 [1.0 + 0.0im, 1.0 + 0.0im]
 [-1.0 + 1.2246467991473532e-16im, 1.0 + 0.0im]
 [1.0 + 0.0im, -1.0 + 1.2246467991473532e-16im]
 [-1.0 + 1.2246467991473532e-16im, -1.0 + 1.2246467991473532e-16im]
```
"""
total_degree_start_solutions(degrees::AbstractVector{<:Integer}) = TotalDegreeStartSolutionsIterator(
    degrees
)


function _paths_to_track_total_degree(
        f::Union{System, AbstractSystem}, ::TotalDegreeAlgorithm
    )::Int
    fixed_f = fix_parameters(f, zeros(ComplexF64, nparameters(f)))
    _, starts = _total_degree_kernel(fixed_f, TotalDegreeAlgorithm())
    return length(starts)
end
