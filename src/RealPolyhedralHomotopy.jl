module RealPolyhedralHomotopy

  # internal
  using HomotopyContinuation
  using AbstractAlgebra
  using PolynomialRoots

  # exported
  using Reexport
  @reexport using HomotopyContinuation, AbstractAlgebra, PolynomialRoots

  include("certify_patchwork.jl")
  include("generate_binomials.jl")
  include("rph_track.jl")

end # module


