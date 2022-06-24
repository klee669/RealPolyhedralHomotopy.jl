module RealPolyhedralHomotopy

  # internal
  using MixedSubdivisions
  using HomotopyContinuation
  using DynamicPolynomials
  using LinearAlgebra
  using AbstractAlgebra
  using PolynomialRoots

  # exported
  using Reexport
  @reexport using MixedSubdivisions, HomotopyContinuation, DynamicPolynomials, LinearAlgebra, AbstractAlgebra, PolynomialRoots

  include("certify_patchwork.jl")
  include("generate_binomials.jl")
  include("rph_track.jl")

end # module


