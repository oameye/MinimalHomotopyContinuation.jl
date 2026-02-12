# Solving Systems

The canonical API is:

```julia
solve(problem, algorithm; kwargs...)
```

For polynomial systems with finitely many solutions, use `SystemProblem` together with
`PolyhedralAlgorithm` or `TotalDegreeAlgorithm`.

For parameter homotopies and sweeps, use `PathTrackingAlgorithm`.

## Supported API

The stable solver surface consists of:

* problem types: `SystemProblem`, `ParameterHomotopyProblem`, `ParameterSweepProblem`, `HomotopyProblem`
* algorithm types: `PolyhedralAlgorithm`, `TotalDegreeAlgorithm`, `PathTrackingAlgorithm`
* compile modes: `CompileAll()`, `CompileMixed()`, `CompileNone()`
* entry points: `init`, `solve`, `solve!`
* typed path counting: `paths_to_track(problem, algorithm)`
* result surface: `Result`, `ResultIterator`, and the result filtering/counting helpers

Low-level tracker/newton/start-system constructors are intentionally internal and unstable.

```@docs
solve
init
PolyhedralAlgorithm
TotalDegreeAlgorithm
PathTrackingAlgorithm
CompileAll
CompileMixed
CompileNone
SystemProblem
ParameterHomotopyProblem
ParameterSweepProblem
HomotopyProblem
paths_to_track
IdentityReducer
MapReducer
FlatMapReducer
```
