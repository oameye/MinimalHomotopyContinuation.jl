function Base.show(io::IO, x::Result)
    s = statistics(x)
    total = s.nonsingular + s.singular
    header = "Result with $total $(plural("solution", total))"
    println(io, header)
    println(io, "="^(length(header)))
    println(io, "• $(ntracked(x)) $(plural("path", ntracked(x))) tracked")
    println(
        io,
        "• $(s.nonsingular) non-singular $(plural("solution", s.nonsingular)) ",
        "($(s.real_nonsingular) real)",
    )
    if s.singular > 0
        println(
            io,
            "• $(s.singular) singular $(plural("solution", s.singular)) ",
            "($(s.real_singular) real)",
        )
    end
    s.excess_solution > 0 && println(
        io,
        "• $(s.excess_solution) excess $(plural("solution", s.excess_solution))",
    )
    println(io, "• random_seed: ", sprint(show, seed(x)))
    if !(algorithm(x) isa UnknownAlgorithm)
        println(io, "• algorithm: ", nameof(typeof(algorithm(x))))
    end
    return nothing
end
