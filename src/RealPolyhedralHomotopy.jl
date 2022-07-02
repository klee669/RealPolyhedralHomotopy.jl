module RealPolyhedralHomotopy

  # internal
  using MixedSubdivisions
  using HomotopyContinuation
  using LinearAlgebra
  using AbstractAlgebra
  using PolynomialRoots

  # exported
  using Reexport
  @reexport using MixedSubdivisions, HomotopyContinuation, LinearAlgebra, AbstractAlgebra, PolynomialRoots

  include("certify_patchwork.jl")
  include("generate_binomials.jl")
  include("rph_track.jl")

end # module


