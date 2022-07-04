using Documenter
using RealPolyhedralHomotopy

push!(LOAD_PATH,"../src/")

makedocs(
         sitename="RealPolyhedralHomotopy.jl Documentation",
         pages = [
            "Index" => "index.md",
         ],
         format = Documenter.HTML(prettyurls = false)
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/klee669/RealPolyhedralHomotopy.jl.git",
)