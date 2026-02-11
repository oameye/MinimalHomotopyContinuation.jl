using PolynomialTestSystems

# see https://www.juliahomotopycontinuation.org/examples/tritangents/
f = equations(tritangents())
prob = SystemProblem(f)

res = solve(prob, TotalDegreeAlgorithm())
println("# nonsingular: ", nnonsingular(res)) # should be 720

res = solve(prob, PolyhedralAlgorithm())
println("# nonsingular: ", nnonsingular(res)) # should be 720
