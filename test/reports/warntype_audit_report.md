# Exhaustive `@code_warntype` Audit Report

Deterministic entry command:

```bash
julia --project=test test/warntype_audit.jl
```

## Executive Summary

- Total discovered src methods: 1231
- Audited methods: 51
- Blocked methods: 979
- Unsupported methods: 201
- Findings (high / medium / low): 2 / 0 / 0

## Findings By Module

- `HomotopyContinuation`: 2

## Top Hotspots

### 1. `HomotopyContinuation.fixed`

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

### 2. `HomotopyContinuation.init`

- Severity: `high`
- Signature: `Tuple{typeof(CommonSolve.init), HomotopyContinuation.SystemProblem, HomotopyContinuation.TotalDegreeAlgorithm}`
- Source: `src/solve/api.jl:100`
- Return type: `HomotopyContinuation.PathSolveCache{SolverT, StartsT, typeof(HomotopyContinuation.always_false)} where {SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.TotalDegreeAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(HomotopyContinuation.StraightLineHomotopy{S, T} where {S<:(HomotopyContinuation.FixedParameterSystem{S} where S<:HomotopyContinuation.MixedSystem), T<:HomotopyContinuation.MixedSystem}))), StartsT<:(HomotopyContinuation.TotalDegreeStartSolutionsIterator{Iter} where Iter<:Base.Iterators.ProductIterator)}`
- Body summary: `body concrete`
- Covered by existing `@inferred` tests: false
- Suggested fix: Stabilize return type for `init` by separating mixed-type branches or adding concretely-typed wrappers.

```text
Body::HOMOTOPYCONTINUATION.PATHSOLVECACHE{SOLVERT, STARTST, TYPEOF(HOMOTOPYCONTINUATION.ALWAYS_FALSE)} WHERE {SOLVERT<:(HOMOTOPYCONTINUATION.SOLVER{T, HOMOTOPYCONTINUATION.TOTALDEGREEALGORITHM{HOMOTOPYCONTINUATION.COMPILEMIXED}} WHERE T<:(HOMOTOPYCONTINUATION.ENDGAMETRACKER{H, MATRIX{COMPLEXF64}} WHERE H<:(HOMOTOPYCONTINUATION.STRAIGHTLINEHOMOTOPY{S, T} WHERE {S<:(HOMOTOPYCONTINUATION.FIXEDPARAMETERSYSTEM{S} WHERE S<:HOMOTOPYCONTINUATION.MIXEDSYSTEM), T<:HOMOTOPYCONTINUATION.MIXEDSYSTEM}))), STARTST<:(HOMOTOPYCONTINUATION.TOTALDEGREESTARTSOLUTIONSITERATOR{ITER} WHERE ITER<:BASE.ITERATORS.PRODUCTITERATOR)}
```

## Full Findings Table

| Severity | Module | Function | Source | Return Type | Body Summary | Inferred Tests |
|---|---|---|---|---|---|---|
| `high` | `HomotopyContinuation` | `fixed` | `src/systems/systems.jl:9` | `HomotopyContinuation.MixedSystem` | `body concrete` | false |
| `high` | `HomotopyContinuation` | `init` | `src/solve/api.jl:100` | `HomotopyContinuation.PathSolveCache{SolverT, StartsT, typeof(HomotopyContinuation.always_false)} where {SolverT<:(HomotopyContinuation.Solver{T, HomotopyContinuation.TotalDegreeAlgorithm{HomotopyContinuation.CompileMixed}} where T<:(HomotopyContinuation.EndgameTracker{H, Matrix{ComplexF64}} where H<:(HomotopyContinuation.StraightLineHomotopy{S, T} where {S<:(HomotopyContinuation.FixedParameterSystem{S} where S<:HomotopyContinuation.MixedSystem), T<:HomotopyContinuation.MixedSystem}))), StartsT<:(HomotopyContinuation.TotalDegreeStartSolutionsIterator{Iter} where Iter<:Base.Iterators.ProductIterator)}` | `body concrete` | false |

## Coverage And Exclusions

### Blocked / Unsupported Reason Counts

- no fixture registered: 979
- generated/internal method name: 201

### Notes

- Actionable criterion: broad escaping return types or `Any` in hot paths.
- Suppressions: finite concrete unions and intentional dynamic `getproperty(::Symbol)` cases.
- Existing inference suites cross-referenced from `test/solve_test.jl` and `test/type_stability_test.jl`.
