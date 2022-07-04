# RealPolyhedralHomotopy.jl

A package for finding real roots of systems of polynomial equations using polyhedral homotopy.
[RealPolyhedralHomotopy.jl](https://github.com/klee669/RealPolyhedralHomotopy.jl)
is a package for finding real roots of systems of polynomial equations using polyhedral homotopy.
The package implements the algorithm for the real polyhedral homotopy establised in [A Polyhedral Homotopy Algorithm For Real Zeros](https://arxiv.org/abs/1910.01957). The idea for the real polyhedral homotopy motivated from the celebrated article [A Polyhedral Method for Solving Sparse Polynomial Systems](https://www.jstor.org/stable/2153370?seq=1) to apply the polyhedral homotopy method for real roots finding.

## Installation

The package can be installed via the Julia package manager
```julia
pkg> add RealPolyhedralHomotopy
```

## Introduction

We support system input through the [HomotopyContinuation](https://www.juliahomotopycontinuation.org) package.
```julia
@var x y;
F = System([-1 - 24000*y + x^3, -9 + 50*x*y - 1*y^2]);
```

For finding real roots, the list of binomial systems corresponding to `F` is required as a start system.
```julia
B = generate_binomials(F);

```
```
2-element Vector{Any}:
 Expression[-24000*y + x^3, 50*x*y - y^2]
 Expression[-24000*y + x^3, -9 + 50*x*y]
```
Using the function `rph_track`, real roots are found by tracking the real polyhedral homotopy.
```julia
realSols = rph_track(B,F)
```
```
4-element Vector{Vector{Float64}}:
 [-1095.4451129504978, -54772.25548320812]
 [1095.4451137838312, 54772.255524874796]
 [8.111114476617955, 0.02219298606763958]
 [-8.103507635567631, -0.02221382112196499]
```

## Functions for the real polyhedral homotopy

```@docs
certify_patchwork
generate_binomials
rph_track
```