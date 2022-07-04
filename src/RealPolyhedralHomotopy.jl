"""
    RealPolyhedralHomotopy

A package for finding real roots of systems of polynomial equations using polyhedral homotopy.
"""
module RealPolyhedralHomotopy

  # internal
  using HomotopyContinuation
  using AbstractAlgebra
  using PolynomialRoots
  using MixedSubdivisions
  using LinearAlgebra

  # exported
  using Reexport
  @reexport using MixedSubdivisions, HomotopyContinuation, AbstractAlgebra, PolynomialRoots, LinearAlgebra

  include("certify_patchwork.jl")
  include("generate_binomials.jl")
  include("rph_track.jl")

end # module


