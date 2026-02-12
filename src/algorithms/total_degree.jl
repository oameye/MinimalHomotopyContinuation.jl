"""
    total_degree(
        F::System;
        parameters = nothing,
        gamma = cis(2π * rand()),
        tracker_options = TrackerOptions(),
        endgame_options = EndgameOptions(),
    )

Solve the system `F` using a total degree homotopy.
This returns a path tracker ([`EndgameTracker`](@ref)) and an iterator to compute the start solutions.
"""
function total_degree(
        F::Union{System, AbstractSystem};
        compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE,
        target_parameters = nothing,
        gamma = cis(2π * rand()),
        tracker_options = TrackerOptions(),
        endgame_options = EndgameOptions(),
    )
    return total_degree_variables(
        F;
        target_parameters,
        gamma,
        tracker_options,
        endgame_options,
        compile_mode,
    )
end


function total_degree_variables(
        F::Union{System, AbstractSystem};
        target_parameters = nothing,
        gamma = cis(2π * rand()),
        tracker_options = TrackerOptions(),
        endgame_options = EndgameOptions(),
        compile_mode::AbstractCompileMode = DEFAULT_COMPILE_MODE,
    )
    _, n = size(F)

    @unique_var x[1:n] s[1:n]
    F₀ = System(F(x, target_parameters), x)
    _validate_affine_square_system(
        F;
        homogeneous_system = F₀,
    )
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
        F = fixed(F; compile_mode = compile_mode)
    end
    if target_parameters !== nothing
        F = fix_parameters(F, target_parameters; compile_mode = compile_mode)
    end
    G = fixed(System(s .* (x .^ D .- 1), x, s); compile_mode = compile_mode)
    G = fix_parameters(G, scaling)

    H = StraightLineHomotopy(G, F; γ = gamma)
    T = EndgameTracker(H, tracker_options = tracker_options, options = endgame_options)
    starts = total_degree_start_solutions(D)

    return T, starts
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
    return print(io, "$(length(iter)) total degree start solutions for degrees $(iter.degrees)")
end

function Base.iterate(iter::TotalDegreeStartSolutionsIterator)
    indices, state = iterate(iter.iterator)
    return _value(iter, indices), state
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
total_degree_start_solutions(
    degrees::AbstractVector{<:Integer},
) = TotalDegreeStartSolutionsIterator(degrees)


function paths_to_track(f::Union{System, AbstractSystem}, ::TotalDegreeAlgorithm)::Int
    target_parameters = zeros(ComplexF64, nparameters(f))
    _, starts = total_degree(f; target_parameters = target_parameters)
    return length(starts)
end
