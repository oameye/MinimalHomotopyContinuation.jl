# Benchmarks

This file contains different scripts of polynomial systems and problems we encountered
during the development of this package.
These problems take too much time to include in our normal test suite, but nonetheless
these should be checked manually after major changes in the package.

For the typed solver refactor baseline, run:

```julia
julia --project -e 'include("benchmarks/refactor_baseline.jl")'
```
