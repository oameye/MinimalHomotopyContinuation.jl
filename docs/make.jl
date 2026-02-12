using Documenter, MinimalHomotopyContinuation

makedocs(;
    sitename = "MinimalHomotopyContinuation.jl",
    pages = [
        "Introduction" => "index.md",
        "API" => "API.md",
        "Example" => "solve_example.md",
    ],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    strict = false,
)

deploydocs(;
    repo = "github.com/JuliaHomotopyContinuation/MinimalHomotopyContinuation.jl.git",
    push_preview = false,
)
