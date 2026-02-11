# Solve Examples


## A First Example (no start system)

We can solve the system ``F(x,y) = (x^2+y^2+1, 2x+3y-1)`` in the following way

```@examples $(simple_solve)
using HomotopyContinuation #hide
@var x y
F = System([x^2+y^2+1, 2x+3y-1], variables = [x, y])
prob = SystemProblem(F)
alg = PolyhedralAlgorithm()
solve(prob, alg; show_progress = false) #hide
solve(prob, alg)
```

Here, the call
```julia
F = System([x^2+y^2+1, 2x+3y-1], variables = [x, y])
```
also determines the ordering of the variables in the solution vectors.
By default, variables are ordered *lexciographically*. If this is okay, you can also
call `solve` without first constructing a system, i.e.,
```julia
solve(SystemProblem([x^2+y^2+1, 2x+3y-1]), PolyhedralAlgorithm())
```

## Parameter Homotopy

Using the problem and algorithm API
```julia
solve(problem, PathTrackingAlgorithm(); kwargs...)
```
We can track the given start solutions alogn the parameter homotopy
```math
H(x, t) = F(x, tp₁+(1-t)p₀),
```
where ``p₁`` (=`start_parameters`) and ``p₀`` (=`target_parameters`) are vectors of
parameter values for ``F`` where ``F`` is a [`System`](@ref) depending on parameters.

Assume we want to perform a parameter homotopy ``H(x,t) := F(x; t[1, 0]+(1-t)[2, 4])`` where
```math
F(x; a) := (x₁^2-a₁, x₁x₂-a₁+a₂)
```
and let's say we are only interested in tracking the solution ``[1,1]``.
This can be accomplished as follows
```@example
using HomotopyContinuation #hide
@var x[1:2] a[1:2]
F = System([x[1]^2-a[1], x[1]*x[2]-a[1]+a[2]], parameters = a)
start_solutions = [[1, 1]]
p₁ = [1, 0]
p₀ = [2, 4]
prob = ParameterHomotopyProblem(
    F,
    start_solutions;
    start_parameters = p₁,
    target_parameters = p₀,
)
solve(prob, PathTrackingAlgorithm(); show_progress=false) #hide
solve(prob, PathTrackingAlgorithm())
```


## Start Target Homotopy

Create a straight-line homotopy explicitly and solve it with a path-tracking algorithm.

This constructs the homotopy ``H(x,t) = tG(x)+(1-t)F(x)`` to compute solutions of the
system `F`.
`start_solutions` is a list of solutions of `G` which are tracked to solutions of `F`.
```@example
using HomotopyContinuation #hide
@var x y
G = System([x^2+1,y+1])
F = System([x^2+y^2+1, 2x+3y-1])
H = StraightLineHomotopy(G, F)
prob = HomotopyProblem(H, [[im, -1], [-im, -1]])
solve(prob, PathTrackingAlgorithm(); show_progress = false) #hide
solve(prob, PathTrackingAlgorithm())
```
