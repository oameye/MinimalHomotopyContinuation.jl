# Exhaustive `@code_warntype` Audit Report

Deterministic entry command:

```bash
julia --project=test test/warntype_audit.jl
```

## Executive Summary

- Total discovered src methods: 1251
- Audited methods: 99
- Blocked methods: 952
- Unsupported methods: 200
- Findings (high / medium / low): 16 / 0 / 0

## Findings By Module

- `MinimalHomotopyContinuation`: 14
- `MinimalHomotopyContinuation.ModelKit`: 2

## Top Hotspots

### 1. `MinimalHomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(MinimalHomotopyContinuation.fixed), MinimalHomotopyContinuation.ModelKit.Homotopy, MinimalHomotopyContinuation.CompileAll}`
- Source: `src/homotopies/homotopies.jl:11`
- Return type: `MinimalHomotopyContinuation.ModelKit.CompiledHomotopy`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MODELKIT.COMPILEDHOMOTOPY
```

### 2. `MinimalHomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(MinimalHomotopyContinuation.fixed), MinimalHomotopyContinuation.ModelKit.Homotopy, MinimalHomotopyContinuation.CompileMixed}`
- Source: `src/homotopies/homotopies.jl:13`
- Return type: `MinimalHomotopyContinuation.MixedHomotopy`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MIXEDHOMOTOPY
```

### 3. `MinimalHomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(MinimalHomotopyContinuation.fixed), MinimalHomotopyContinuation.ModelKit.Homotopy}`
- Source: `src/homotopies/homotopies.jl:7`
- Return type: `MinimalHomotopyContinuation.MixedHomotopy`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MIXEDHOMOTOPY
```

### 4. `MinimalHomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(MinimalHomotopyContinuation.fixed), MinimalHomotopyContinuation.ModelKit.System, MinimalHomotopyContinuation.CompileAll}`
- Source: `src/systems/systems.jl:13`
- Return type: `MinimalHomotopyContinuation.ModelKit.CompiledSystem`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MODELKIT.COMPILEDSYSTEM
```

### 5. `MinimalHomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(MinimalHomotopyContinuation.fixed), MinimalHomotopyContinuation.ModelKit.System, MinimalHomotopyContinuation.CompileMixed}`
- Source: `src/systems/systems.jl:15`
- Return type: `MinimalHomotopyContinuation.MixedSystem`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MIXEDSYSTEM
```

### 6. `MinimalHomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(MinimalHomotopyContinuation.fixed), MinimalHomotopyContinuation.ModelKit.System}`
- Source: `src/systems/systems.jl:9`
- Return type: `MinimalHomotopyContinuation.MixedSystem`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MIXEDSYSTEM
```

### 7. `MinimalHomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), MinimalHomotopyContinuation.ParameterHomotopyProblem, MinimalHomotopyContinuation.PathTrackingAlgorithm, Type{MinimalHomotopyContinuation.ResultIterator}}`
- Source: `src/solve/api.jl:207`
- Return type: `MinimalHomotopyContinuation.PathIteratorSolveCache{SolverT, Vector{Vector{ComplexF64}}, Nothing} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.PathTrackingAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(MinimalHomotopyContinuation.ParameterHomotopy{T} where T<:MinimalHomotopyContinuation.MixedSystem)))`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.PATHITERATORSOLVECACHE{SOLVERT, VECTOR{VECTOR{COMPLEXF64}}, NOTHING} WHERE SOLVERT<:(HOMOTOPYCONTINUATION.SOLVER{T, HOMOTOPYCONTINUATION.PATHTRACKINGALGORITHM{HOMOTOPYCONTINUATION.COMPILEMIXED}} WHERE T<:(HOMOTOPYCONTINUATION.ENDGAMETRACKER{H, MATRIX{COMPLEXF64}} WHERE H<:(HOMOTOPYCONTINUATION.PARAMETERHOMOTOPY{T} WHERE T<:HOMOTOPYCONTINUATION.MIXEDSYSTEM)))
```

### 8. `MinimalHomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), MinimalHomotopyContinuation.ParameterHomotopyProblem, MinimalHomotopyContinuation.PathTrackingAlgorithm}`
- Source: `src/solve/api.jl:111`
- Return type: `MinimalHomotopyContinuation.PathSolveCache{SolverT, Vector{Vector{ComplexF64}}, typeof(MinimalHomotopyContinuation.always_false)} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.PathTrackingAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(MinimalHomotopyContinuation.ParameterHomotopy{T} where T<:MinimalHomotopyContinuation.MixedSystem)))`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.PATHSOLVECACHE{SOLVERT, VECTOR{VECTOR{COMPLEXF64}}, TYPEOF(HOMOTOPYCONTINUATION.ALWAYS_FALSE)} WHERE SOLVERT<:(HOMOTOPYCONTINUATION.SOLVER{T, HOMOTOPYCONTINUATION.PATHTRACKINGALGORITHM{HOMOTOPYCONTINUATION.COMPILEMIXED}} WHERE T<:(HOMOTOPYCONTINUATION.ENDGAMETRACKER{H, MATRIX{COMPLEXF64}} WHERE H<:(HOMOTOPYCONTINUATION.PARAMETERHOMOTOPY{T} WHERE T<:HOMOTOPYCONTINUATION.MIXEDSYSTEM)))
```

### 9. `MinimalHomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), MinimalHomotopyContinuation.ParameterSweepProblem, MinimalHomotopyContinuation.PathTrackingAlgorithm, Type{MinimalHomotopyContinuation.ResultIterator}}`
- Source: `src/solve/api.jl:240`
- Return type: `MinimalHomotopyContinuation.SweepIteratorSolveCache{SolversT, Vector{Vector{Float64}}, Nothing} where SolversT<:(Vector)`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.SWEEPITERATORSOLVECACHE{SOLVERST, VECTOR{VECTOR{FLOAT64}}, NOTHING} WHERE SOLVERST<:(VECTOR)
```

### 10. `MinimalHomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), MinimalHomotopyContinuation.ParameterSweepProblem, MinimalHomotopyContinuation.PathTrackingAlgorithm}`
- Source: `src/solve/api.jl:153`
- Return type: `MinimalHomotopyContinuation.SweepSolveCache{SolverT, Vector{Vector{Float64}}, Vector{Vector{Float64}}, MinimalHomotopyContinuation.MapReducer{Main.WarntypeAudit.var"#build_context##0#build_context##1"}} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.PathTrackingAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(MinimalHomotopyContinuation.ParameterHomotopy{T} where T<:MinimalHomotopyContinuation.MixedSystem)))`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.SWEEPSOLVECACHE{SOLVERT, VECTOR{VECTOR{FLOAT64}}, VECTOR{VECTOR{FLOAT64}}, HOMOTOPYCONTINUATION.MAPREDUCER{MAIN.WARNTYPEAUDIT.VAR"#BUILD_CONTEXT##0#BUILD_CONTEXT##1"}} WHERE SOLVERT<:(HOMOTOPYCONTINUATION.SOLVER{T, HOMOTOPYCONTINUATION.PATHTRACKINGALGORITHM{HOMOTOPYCONTINUATION.COMPILEMIXED}} WHERE T<:(HOMOTOPYCONTINUATION.ENDGAMETRACKER{H, MATRIX{COMPLEXF64}} WHERE H<:(HOMOTOPYCONTINUATION.PARAMETERHOMOTOPY{T} WHERE T<:HOMOTOPYCONTINUATION.MIXEDSYSTEM)))
```

### 11. `MinimalHomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), MinimalHomotopyContinuation.SystemProblem, MinimalHomotopyContinuation.PolyhedralAlgorithm, Type{MinimalHomotopyContinuation.ResultIterator}}`
- Source: `src/solve/api.jl:218`
- Return type: `MinimalHomotopyContinuation.PathIteratorSolveCache{SolverT, MinimalHomotopyContinuation.PolyhedralStartSolutionsIterator{Vector{MixedSubdivisions.MixedCell}}, Nothing} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.PolyhedralAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.PolyhedralTracker{H, H1, Matrix{ComplexF64}} where {H<:(MinimalHomotopyContinuation.ToricHomotopy{S} where S<:MinimalHomotopyContinuation.MixedSystem), H1<:(MinimalHomotopyContinuation.CoefficientHomotopy{T} where T<:MinimalHomotopyContinuation.MixedSystem)}))`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.PATHITERATORSOLVECACHE{SOLVERT, HOMOTOPYCONTINUATION.POLYHEDRALSTARTSOLUTIONSITERATOR{VECTOR{MIXEDSUBDIVISIONS.MIXEDCELL}}, NOTHING} WHERE SOLVERT<:(HOMOTOPYCONTINUATION.SOLVER{T, HOMOTOPYCONTINUATION.POLYHEDRALALGORITHM{HOMOTOPYCONTINUATION.COMPILEMIXED}} WHERE T<:(HOMOTOPYCONTINUATION.POLYHEDRALTRACKER{H, H1, MATRIX{COMPLEXF64}} WHERE {H<:(HOMOTOPYCONTINUATION.TORICHOMOTOPY{S} WHERE S<:HOMOTOPYCONTINUATION.MIXEDSYSTEM), H1<:(HOMOTOPYCONTINUATION.COEFFICIENTHOMOTOPY{T} WHERE T<:HOMOTOPYCONTINUATION.MIXEDSYSTEM)}))
```

### 12. `MinimalHomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), MinimalHomotopyContinuation.SystemProblem, MinimalHomotopyContinuation.PolyhedralAlgorithm}`
- Source: `src/solve/api.jl:125`
- Return type: `MinimalHomotopyContinuation.PathSolveCache{SolverT, MinimalHomotopyContinuation.PolyhedralStartSolutionsIterator{Vector{MixedSubdivisions.MixedCell}}, typeof(MinimalHomotopyContinuation.always_false)} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.PolyhedralAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.PolyhedralTracker{H, H1, Matrix{ComplexF64}} where {H<:(MinimalHomotopyContinuation.ToricHomotopy{S} where S<:MinimalHomotopyContinuation.MixedSystem), H1<:(MinimalHomotopyContinuation.CoefficientHomotopy{T} where T<:MinimalHomotopyContinuation.MixedSystem)}))`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.PATHSOLVECACHE{SOLVERT, HOMOTOPYCONTINUATION.POLYHEDRALSTARTSOLUTIONSITERATOR{VECTOR{MIXEDSUBDIVISIONS.MIXEDCELL}}, TYPEOF(HOMOTOPYCONTINUATION.ALWAYS_FALSE)} WHERE SOLVERT<:(HOMOTOPYCONTINUATION.SOLVER{T, HOMOTOPYCONTINUATION.POLYHEDRALALGORITHM{HOMOTOPYCONTINUATION.COMPILEMIXED}} WHERE T<:(HOMOTOPYCONTINUATION.POLYHEDRALTRACKER{H, H1, MATRIX{COMPLEXF64}} WHERE {H<:(HOMOTOPYCONTINUATION.TORICHOMOTOPY{S} WHERE S<:HOMOTOPYCONTINUATION.MIXEDSYSTEM), H1<:(HOMOTOPYCONTINUATION.COEFFICIENTHOMOTOPY{T} WHERE T<:HOMOTOPYCONTINUATION.MIXEDSYSTEM)}))
```

## Full Findings Table

| Severity | Module | Function | Source | Return Type | Body Summary | Inferred Tests |
|---|---|---|---|---|---|---|
| `high` | `MinimalHomotopyContinuation` | `fixed` | `src/homotopies/homotopies.jl:11` | `MinimalHomotopyContinuation.ModelKit.CompiledHomotopy` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `fixed` | `src/homotopies/homotopies.jl:13` | `MinimalHomotopyContinuation.MixedHomotopy` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `fixed` | `src/homotopies/homotopies.jl:7` | `MinimalHomotopyContinuation.MixedHomotopy` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `fixed` | `src/systems/systems.jl:13` | `MinimalHomotopyContinuation.ModelKit.CompiledSystem` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `fixed` | `src/systems/systems.jl:15` | `MinimalHomotopyContinuation.MixedSystem` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `fixed` | `src/systems/systems.jl:9` | `MinimalHomotopyContinuation.MixedSystem` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `init` | `src/solve/api.jl:207` | `MinimalHomotopyContinuation.PathIteratorSolveCache{SolverT, Vector{Vector{ComplexF64}}, Nothing} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.PathTrackingAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(MinimalHomotopyContinuation.ParameterHomotopy{T} where T<:MinimalHomotopyContinuation.MixedSystem)))` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `init` | `src/solve/api.jl:111` | `MinimalHomotopyContinuation.PathSolveCache{SolverT, Vector{Vector{ComplexF64}}, typeof(MinimalHomotopyContinuation.always_false)} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.PathTrackingAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(MinimalHomotopyContinuation.ParameterHomotopy{T} where T<:MinimalHomotopyContinuation.MixedSystem)))` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `init` | `src/solve/api.jl:240` | `MinimalHomotopyContinuation.SweepIteratorSolveCache{SolversT, Vector{Vector{Float64}}, Nothing} where SolversT<:(Vector)` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `init` | `src/solve/api.jl:153` | `MinimalHomotopyContinuation.SweepSolveCache{SolverT, Vector{Vector{Float64}}, Vector{Vector{Float64}}, MinimalHomotopyContinuation.MapReducer{Main.WarntypeAudit.var"#build_context##0#build_context##1"}} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.PathTrackingAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(MinimalHomotopyContinuation.ParameterHomotopy{T} where T<:MinimalHomotopyContinuation.MixedSystem)))` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `init` | `src/solve/api.jl:218` | `MinimalHomotopyContinuation.PathIteratorSolveCache{SolverT, MinimalHomotopyContinuation.PolyhedralStartSolutionsIterator{Vector{MixedSubdivisions.MixedCell}}, Nothing} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.PolyhedralAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.PolyhedralTracker{H, H1, Matrix{ComplexF64}} where {H<:(MinimalHomotopyContinuation.ToricHomotopy{S} where S<:MinimalHomotopyContinuation.MixedSystem), H1<:(MinimalHomotopyContinuation.CoefficientHomotopy{T} where T<:MinimalHomotopyContinuation.MixedSystem)}))` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `init` | `src/solve/api.jl:125` | `MinimalHomotopyContinuation.PathSolveCache{SolverT, MinimalHomotopyContinuation.PolyhedralStartSolutionsIterator{Vector{MixedSubdivisions.MixedCell}}, typeof(MinimalHomotopyContinuation.always_false)} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.PolyhedralAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.PolyhedralTracker{H, H1, Matrix{ComplexF64}} where {H<:(MinimalHomotopyContinuation.ToricHomotopy{S} where S<:MinimalHomotopyContinuation.MixedSystem), H1<:(MinimalHomotopyContinuation.CoefficientHomotopy{T} where T<:MinimalHomotopyContinuation.MixedSystem)}))` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `init` | `src/solve/api.jl:196` | `MinimalHomotopyContinuation.PathIteratorSolveCache{SolverT, MinimalHomotopyContinuation.TotalDegreeStartSolutionsIterator, Nothing} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.TotalDegreeAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(MinimalHomotopyContinuation.StraightLineHomotopy{S, T} where {S<:(MinimalHomotopyContinuation.FixedParameterSystem{S, Float64} where S<:MinimalHomotopyContinuation.MixedSystem), T<:MinimalHomotopyContinuation.MixedSystem})))` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation` | `init` | `src/solve/api.jl:97` | `MinimalHomotopyContinuation.PathSolveCache{SolverT, MinimalHomotopyContinuation.TotalDegreeStartSolutionsIterator, typeof(MinimalHomotopyContinuation.always_false)} where SolverT<:(MinimalHomotopyContinuation.Solver{T, MinimalHomotopyContinuation.TotalDegreeAlgorithm{MinimalHomotopyContinuation.CompileMixed}} where T<:(MinimalHomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(MinimalHomotopyContinuation.StraightLineHomotopy{S, T} where {S<:(MinimalHomotopyContinuation.FixedParameterSystem{S, Float64} where S<:MinimalHomotopyContinuation.MixedSystem), T<:MinimalHomotopyContinuation.MixedSystem})))` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation.ModelKit` | `evaluate` | `src/model_kit/symbolic.jl:1514` | `Any` | `body concrete` | false |
| `high` | `MinimalHomotopyContinuation.ModelKit` | `evaluate` | `src/model_kit/symbolic.jl:1193` | `Any` | `body concrete` | false |

## Coverage And Exclusions

### Blocked / Unsupported Reason Counts

- no fixture registered: 952
- generated/internal method name: 200

### Notes

- Actionable criterion: broad escaping return types or `Any` in hot paths.
- Suppressions: finite concrete unions and intentional dynamic `getproperty(::Symbol)` cases.
- Existing inference suites cross-referenced from `test/solve_test.jl` and `test/type_stability_test.jl`.
