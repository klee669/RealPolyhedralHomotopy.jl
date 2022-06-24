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
  @reexport using HomotopyContinuation

  include("certify_patchwork.jl")
  include("generate_binomials.jl")
  include("rph_track.jl")

end # module

