```@meta
CollapsedDocStrings = true
```

# API

This page documents the stable public API surface.

For current exports, see `src/public_api.jl`.
Low-level tracker/newton/start-system constructors are intentionally internal
and unstable.

## Contents

```@contents
Pages = ["API.md"]
Depth = 2
```

## Entry points

```@docs
init
solve
solve!
paths_to_track
```

## Problems

```@docs
AbstractHCProblem
SystemProblem
ParameterHomotopyProblem
ParameterSweepProblem
HomotopyProblem
```

## Algorithms

```@docs
AbstractHCAlgorithm
PolyhedralAlgorithm
TotalDegreeAlgorithm
PathTrackingAlgorithm
AbstractSweepReducer
IdentityReducer
MapReducer
FlatMapReducer
```

## Compile modes

```@docs
AbstractCompileMode
CompileAll
CompileMixed
CompileNone
DEFAULT_COMPILE_MODE
```

## Result types

```@docs
Result
ResultStatistics
statistics
algorithm
seed
path_results
```

## Result filtering

```@docs
results
solutions
real_solutions
nonsingular
singular
at_infinity
failed
```

## Result counting

```@docs
nresults
nfinite
nsolutions
nsingular
nat_infinity
nexcess_solutions
nfailed
nnonsingular
nreal
ntracked
```

## Result iterators

```@docs
ResultIterator
bitmask
bitmask_filter
solver
start_solutions
```

## Path results

```@docs
PathResult
PathResultCode
solution
accuracy
residual
steps
accepted_steps
rejected_steps
winding_number
path_number
start_solution
multiplicity
last_path_point
valuation
is_success
is_at_infinity
is_excess_solution
is_failed
is_finite
is_singular
is_nonsingular
is_real
```