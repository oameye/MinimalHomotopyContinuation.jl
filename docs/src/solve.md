# Solving Systems

The canonical API is:

```julia
solve(problem, algorithm; kwargs...)
```

For polynomial systems with finitely many solutions, use `SystemProblem` together with
`PolyhedralAlgorithm` or `TotalDegreeAlgorithm`.

For parameter homotopies and sweeps, use `PathTrackingAlgorithm`.

```@docs
solve
init
Solver
PolyhedralAlgorithm
TotalDegreeAlgorithm
PathTrackingAlgorithm
CompileAll
CompileMixed
CompileNone
paths_to_track
IdentityReducer
MapReducer
FlatMapReducer
```
