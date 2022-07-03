using Documenter
using RealPolyhedralHomotopy

makedocs(
    sitename = "RealPolyhedralHomotopy",
    format = Documenter.HTML(),
    modules = [RealPolyhedralHomotopy]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/klee669/RealPolyhedralHomotopy.jl.git",
    devbranch = "main"
)
