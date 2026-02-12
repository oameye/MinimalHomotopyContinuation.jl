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

- `HomotopyContinuation`: 14
- `HomotopyContinuation.ModelKit`: 2

## Top Hotspots

### 1. `HomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(HomotopyContinuation.fixed), HomotopyContinuation.ModelKit.Homotopy, HomotopyContinuation.CompileAll}`
- Source: `src/homotopies/homotopies.jl:11`
- Return type: `HomotopyContinuation.ModelKit.CompiledHomotopy`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MODELKIT.COMPILEDHOMOTOPY
```

### 2. `HomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(HomotopyContinuation.fixed), HomotopyContinuation.ModelKit.Homotopy, HomotopyContinuation.CompileMixed}`
- Source: `src/homotopies/homotopies.jl:13`
- Return type: `HomotopyContinuation.MixedHomotopy`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MIXEDHOMOTOPY
```

### 3. `HomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(HomotopyContinuation.fixed), HomotopyContinuation.ModelKit.Homotopy}`
- Source: `src/homotopies/homotopies.jl:7`
- Return type: `HomotopyContinuation.MixedHomotopy`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MIXEDHOMOTOPY
```

### 4. `HomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(HomotopyContinuation.fixed), HomotopyContinuation.ModelKit.System, HomotopyContinuation.CompileAll}`
- Source: `src/systems/systems.jl:13`
- Return type: `HomotopyContinuation.ModelKit.CompiledSystem`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MODELKIT.COMPILEDSYSTEM
```

### 5. `HomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(HomotopyContinuation.fixed), HomotopyContinuation.ModelKit.System, HomotopyContinuation.CompileMixed}`
- Source: `src/systems/systems.jl:15`
- Return type: `HomotopyContinuation.MixedSystem`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MIXEDSYSTEM
```

### 6. `HomotopyContinuation.fixed`

- Severity: `high`
- Signature: `Tuple{typeof(HomotopyContinuation.fixed), HomotopyContinuation.ModelKit.System}`
- Source: `src/systems/systems.jl:9`
- Return type: `HomotopyContinuation.MixedSystem`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `fixed` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.MIXEDSYSTEM
```

### 7. `HomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), HomotopyContinuation.ParameterHomotopyProblem, HomotopyContinuation.PathTrackingAlgorithm, Type{HomotopyContinuation.ResultIterator}}`
- Source: `src/solve/api.jl:207`
- Return type: `HomotopyContinuation.PathIteratorSolveCache{SolverT, Vector{Vector{ComplexF64}}, Nothing} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.PathTrackingAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(HomotopyContinuation.ParameterHomotopy{T} where T<:HomotopyContinuation.MixedSystem)))`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.PATHITERATORSOLVECACHE{SOLVERT, VECTOR{VECTOR{COMPLEXF64}}, NOTHING} WHERE SOLVERT<:(HOMOTOPYCONTINUATION.SOLVER{T, HOMOTOPYCONTINUATION.PATHTRACKINGALGORITHM{HOMOTOPYCONTINUATION.COMPILEMIXED}} WHERE T<:(HOMOTOPYCONTINUATION.ENDGAMETRACKER{H, MATRIX{COMPLEXF64}} WHERE H<:(HOMOTOPYCONTINUATION.PARAMETERHOMOTOPY{T} WHERE T<:HOMOTOPYCONTINUATION.MIXEDSYSTEM)))
```

### 8. `HomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), HomotopyContinuation.ParameterHomotopyProblem, HomotopyContinuation.PathTrackingAlgorithm}`
- Source: `src/solve/api.jl:111`
- Return type: `HomotopyContinuation.PathSolveCache{SolverT, Vector{Vector{ComplexF64}}, typeof(HomotopyContinuation.always_false)} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.PathTrackingAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(HomotopyContinuation.ParameterHomotopy{T} where T<:HomotopyContinuation.MixedSystem)))`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.PATHSOLVECACHE{SOLVERT, VECTOR{VECTOR{COMPLEXF64}}, TYPEOF(HOMOTOPYCONTINUATION.ALWAYS_FALSE)} WHERE SOLVERT<:(HOMOTOPYCONTINUATION.SOLVER{T, HOMOTOPYCONTINUATION.PATHTRACKINGALGORITHM{HOMOTOPYCONTINUATION.COMPILEMIXED}} WHERE T<:(HOMOTOPYCONTINUATION.ENDGAMETRACKER{H, MATRIX{COMPLEXF64}} WHERE H<:(HOMOTOPYCONTINUATION.PARAMETERHOMOTOPY{T} WHERE T<:HOMOTOPYCONTINUATION.MIXEDSYSTEM)))
```

### 9. `HomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), HomotopyContinuation.ParameterSweepProblem, HomotopyContinuation.PathTrackingAlgorithm, Type{HomotopyContinuation.ResultIterator}}`
- Source: `src/solve/api.jl:240`
- Return type: `HomotopyContinuation.SweepIteratorSolveCache{SolversT, Vector{Vector{Float64}}, Nothing} where SolversT<:(Vector)`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.SWEEPITERATORSOLVECACHE{SOLVERST, VECTOR{VECTOR{FLOAT64}}, NOTHING} WHERE SOLVERST<:(VECTOR)
```

### 10. `HomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), HomotopyContinuation.ParameterSweepProblem, HomotopyContinuation.PathTrackingAlgorithm}`
- Source: `src/solve/api.jl:153`
- Return type: `HomotopyContinuation.SweepSolveCache{SolverT, Vector{Vector{Float64}}, Vector{Vector{Float64}}, HomotopyContinuation.MapReducer{Main.WarntypeAudit.var"#build_context##0#build_context##1"}} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.PathTrackingAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(HomotopyContinuation.ParameterHomotopy{T} where T<:HomotopyContinuation.MixedSystem)))`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.SWEEPSOLVECACHE{SOLVERT, VECTOR{VECTOR{FLOAT64}}, VECTOR{VECTOR{FLOAT64}}, HOMOTOPYCONTINUATION.MAPREDUCER{MAIN.WARNTYPEAUDIT.VAR"#BUILD_CONTEXT##0#BUILD_CONTEXT##1"}} WHERE SOLVERT<:(HOMOTOPYCONTINUATION.SOLVER{T, HOMOTOPYCONTINUATION.PATHTRACKINGALGORITHM{HOMOTOPYCONTINUATION.COMPILEMIXED}} WHERE T<:(HOMOTOPYCONTINUATION.ENDGAMETRACKER{H, MATRIX{COMPLEXF64}} WHERE H<:(HOMOTOPYCONTINUATION.PARAMETERHOMOTOPY{T} WHERE T<:HOMOTOPYCONTINUATION.MIXEDSYSTEM)))
```

### 11. `HomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), HomotopyContinuation.SystemProblem, HomotopyContinuation.PolyhedralAlgorithm, Type{HomotopyContinuation.ResultIterator}}`
- Source: `src/solve/api.jl:218`
- Return type: `HomotopyContinuation.PathIteratorSolveCache{SolverT, HomotopyContinuation.PolyhedralStartSolutionsIterator{Vector{MixedSubdivisions.MixedCell}}, Nothing} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.PolyhedralAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.PolyhedralTracker{H, H1, Matrix{ComplexF64}} where {H<:(HomotopyContinuation.ToricHomotopy{S} where S<:HomotopyContinuation.MixedSystem), H1<:(HomotopyContinuation.CoefficientHomotopy{T} where T<:HomotopyContinuation.MixedSystem)}))`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.PATHITERATORSOLVECACHE{SOLVERT, HOMOTOPYCONTINUATION.POLYHEDRALSTARTSOLUTIONSITERATOR{VECTOR{MIXEDSUBDIVISIONS.MIXEDCELL}}, NOTHING} WHERE SOLVERT<:(HOMOTOPYCONTINUATION.SOLVER{T, HOMOTOPYCONTINUATION.POLYHEDRALALGORITHM{HOMOTOPYCONTINUATION.COMPILEMIXED}} WHERE T<:(HOMOTOPYCONTINUATION.POLYHEDRALTRACKER{H, H1, MATRIX{COMPLEXF64}} WHERE {H<:(HOMOTOPYCONTINUATION.TORICHOMOTOPY{S} WHERE S<:HOMOTOPYCONTINUATION.MIXEDSYSTEM), H1<:(HOMOTOPYCONTINUATION.COEFFICIENTHOMOTOPY{T} WHERE T<:HOMOTOPYCONTINUATION.MIXEDSYSTEM)}))
```

### 12. `HomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), HomotopyContinuation.SystemProblem, HomotopyContinuation.PolyhedralAlgorithm}`
- Source: `src/solve/api.jl:125`
- Return type: `HomotopyContinuation.PathSolveCache{SolverT, HomotopyContinuation.PolyhedralStartSolutionsIterator{Vector{MixedSubdivisions.MixedCell}}, typeof(HomotopyContinuation.always_false)} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.PolyhedralAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.PolyhedralTracker{H, H1, Matrix{ComplexF64}} where {H<:(HomotopyContinuation.ToricHomotopy{S} where S<:HomotopyContinuation.MixedSystem), H1<:(HomotopyContinuation.CoefficientHomotopy{T} where T<:HomotopyContinuation.MixedSystem)}))`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.PATHSOLVECACHE{SOLVERT, HOMOTOPYCONTINUATION.POLYHEDRALSTARTSOLUTIONSITERATOR{VECTOR{MIXEDSUBDIVISIONS.MIXEDCELL}}, TYPEOF(HOMOTOPYCONTINUATION.ALWAYS_FALSE)} WHERE SOLVERT<:(HOMOTOPYCONTINUATION.SOLVER{T, HOMOTOPYCONTINUATION.POLYHEDRALALGORITHM{HOMOTOPYCONTINUATION.COMPILEMIXED}} WHERE T<:(HOMOTOPYCONTINUATION.POLYHEDRALTRACKER{H, H1, MATRIX{COMPLEXF64}} WHERE {H<:(HOMOTOPYCONTINUATION.TORICHOMOTOPY{S} WHERE S<:HOMOTOPYCONTINUATION.MIXEDSYSTEM), H1<:(HOMOTOPYCONTINUATION.COEFFICIENTHOMOTOPY{T} WHERE T<:HOMOTOPYCONTINUATION.MIXEDSYSTEM)}))
```

## Full Findings Table

| Severity | Module | Function | Source | Return Type | Body Summary | Inferred Tests |
|---|---|---|---|---|---|---|
| `high` | `HomotopyContinuation` | `fixed` | `src/homotopies/homotopies.jl:11` | `HomotopyContinuation.ModelKit.CompiledHomotopy` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `fixed` | `src/homotopies/homotopies.jl:13` | `HomotopyContinuation.MixedHomotopy` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `fixed` | `src/homotopies/homotopies.jl:7` | `HomotopyContinuation.MixedHomotopy` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `fixed` | `src/systems/systems.jl:13` | `HomotopyContinuation.ModelKit.CompiledSystem` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `fixed` | `src/systems/systems.jl:15` | `HomotopyContinuation.MixedSystem` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `fixed` | `src/systems/systems.jl:9` | `HomotopyContinuation.MixedSystem` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `init` | `src/solve/api.jl:207` | `HomotopyContinuation.PathIteratorSolveCache{SolverT, Vector{Vector{ComplexF64}}, Nothing} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.PathTrackingAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(HomotopyContinuation.ParameterHomotopy{T} where T<:HomotopyContinuation.MixedSystem)))` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `init` | `src/solve/api.jl:111` | `HomotopyContinuation.PathSolveCache{SolverT, Vector{Vector{ComplexF64}}, typeof(HomotopyContinuation.always_false)} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.PathTrackingAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(HomotopyContinuation.ParameterHomotopy{T} where T<:HomotopyContinuation.MixedSystem)))` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `init` | `src/solve/api.jl:240` | `HomotopyContinuation.SweepIteratorSolveCache{SolversT, Vector{Vector{Float64}}, Nothing} where SolversT<:(Vector)` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `init` | `src/solve/api.jl:153` | `HomotopyContinuation.SweepSolveCache{SolverT, Vector{Vector{Float64}}, Vector{Vector{Float64}}, HomotopyContinuation.MapReducer{Main.WarntypeAudit.var"#build_context##0#build_context##1"}} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.PathTrackingAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(HomotopyContinuation.ParameterHomotopy{T} where T<:HomotopyContinuation.MixedSystem)))` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `init` | `src/solve/api.jl:218` | `HomotopyContinuation.PathIteratorSolveCache{SolverT, HomotopyContinuation.PolyhedralStartSolutionsIterator{Vector{MixedSubdivisions.MixedCell}}, Nothing} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.PolyhedralAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.PolyhedralTracker{H, H1, Matrix{ComplexF64}} where {H<:(HomotopyContinuation.ToricHomotopy{S} where S<:HomotopyContinuation.MixedSystem), H1<:(HomotopyContinuation.CoefficientHomotopy{T} where T<:HomotopyContinuation.MixedSystem)}))` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `init` | `src/solve/api.jl:125` | `HomotopyContinuation.PathSolveCache{SolverT, HomotopyContinuation.PolyhedralStartSolutionsIterator{Vector{MixedSubdivisions.MixedCell}}, typeof(HomotopyContinuation.always_false)} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.PolyhedralAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.PolyhedralTracker{H, H1, Matrix{ComplexF64}} where {H<:(HomotopyContinuation.ToricHomotopy{S} where S<:HomotopyContinuation.MixedSystem), H1<:(HomotopyContinuation.CoefficientHomotopy{T} where T<:HomotopyContinuation.MixedSystem)}))` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `init` | `src/solve/api.jl:196` | `HomotopyContinuation.PathIteratorSolveCache{SolverT, HomotopyContinuation.TotalDegreeStartSolutionsIterator, Nothing} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.TotalDegreeAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(HomotopyContinuation.StraightLineHomotopy{S, T} where {S<:(HomotopyContinuation.FixedParameterSystem{S, Float64} where S<:HomotopyContinuation.MixedSystem), T<:HomotopyContinuation.MixedSystem})))` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `init` | `src/solve/api.jl:97` | `HomotopyContinuation.PathSolveCache{SolverT, HomotopyContinuation.TotalDegreeStartSolutionsIterator, typeof(HomotopyContinuation.always_false)} where SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.TotalDegreeAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(HomotopyContinuation.StraightLineHomotopy{S, T} where {S<:(HomotopyContinuation.FixedParameterSystem{S, Float64} where S<:HomotopyContinuation.MixedSystem), T<:HomotopyContinuation.MixedSystem})))` | `body concrete` | false |
| `high` | `HomotopyContinuation.ModelKit` | `evaluate` | `src/model_kit/symbolic.jl:1514` | `Any` | `body concrete` | false |
| `high` | `HomotopyContinuation.ModelKit` | `evaluate` | `src/model_kit/symbolic.jl:1193` | `Any` | `body concrete` | false |

## Coverage And Exclusions

### Blocked / Unsupported Reason Counts

- no fixture registered: 952
- generated/internal method name: 200

### Notes

- Actionable criterion: broad escaping return types or `Any` in hot paths.
- Suppressions: finite concrete unions and intentional dynamic `getproperty(::Symbol)` cases.
- Existing inference suites cross-referenced from `test/solve_test.jl` and `test/type_stability_test.jl`.
