# RealPolyhedralHomotopy.jl

## Documents
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://klee669.github.io/RealPolyhedralHomotopy.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://klee669.github.io/RealPolyhedralHomotopy.jl/dev)

`RealPolyhedralHomotopy.jl` is a package for finding real roots of systems of polynomial equations using polyhedral homotopy.
The package implements the algorithm for the real polyhedral homotopy establised in [A Polyhedral Homotopy Algorithm For Real Zeros](https://arxiv.org/abs/1910.01957). The idea for the real polyhedral homotopy motivated from the celebrated article [A Polyhedral Method for Solving Sparse Polynomial Systems](https://www.jstor.org/stable/2153370?seq=1) to apply the polyhedral homotopy method for real roots finding.

## Installation

The package can be installed via the Julia package manager
```julia
pkg> add RealPolyhedralHomotopy
```

## Build status
[![Build status (Github Actions)](https://github.com/sylvaticus/MyAwesomePackage.jl/workflows/CI/badge.svg)](https://github.com/sylvaticus/MyAwesomePackage.jl/actions)
[![codecov.io](http://codecov.io/github/sylvaticus/MyAwesomePackage.jl/coverage.svg?branch=main)](http://codecov.io/github/sylvaticus/MyAwesomePackage.jl?branch=main)
