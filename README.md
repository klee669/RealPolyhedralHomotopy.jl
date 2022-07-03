# RealPolyhedralHomotopy.jl
A package for finding real roots of systems of polynomial equations using polyhedral homotopy.

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://klee669.github.io/RealPolyhedralHomotopy.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://klee669.github.io/RealPolyhedralHomotopy.jl/dev)

## Usage
To run RealPolyhedralHomotopy.jl, it is suggested to use Julia 1.6.1 or higher.
- Clone this repository.
    - Create a new julia notebook inside the `RealPolyhedralHomotopy.jl` folder and put the following inside the first cell.
    ```
    using Pkg
    Pkg.activate(".")
    using RealPolyhedralHomotopy
    ```
    
    - or run Julia in the `RealPolyhedralHomotopy.jl` folder and run the following commands:
    ```
    using Pkg
    Pkg.activate(".")
    using RealPolyhedralHomotopy
    ```

Polynomials can be constructed by declaring variables with `@var` and including them into a Julia expression.
The Package exports the functions `certify_patchwork`, `generate_binomials` and `rph_track`.

## Example
### Running the package.
Run the following commands in the folder `RealPolyhedralHomotopy.jl`:
```julia
using Pkg
Pkg.activate(".")
using RealPolyhedralHomotopy
```

### Constructing a system.
Input:
```julia
@var x y;
F = System([-1 - 24000*y + x^3, -9 + 50*x*y - 1*y^2]);
```
### Certifying patchworkedness.
Input:
```julia
result = certify_patchwork(F)
```

Output:
```
1
```


### Generating starting binomial systems.
Input:
```julia
B = generate_binomials(F)
```
Output:
```
2-element Vector{Any}:
 Expression[-24000*y + x^3, 50*x*y - y^2]
 Expression[-24000*y + x^3, -9 + 50*x*y]
```


### Tracking the the real solutions.
Input:
```julia
realSols = rph_track(B,F)
```
Output:
```
4-element Vector{Vector{Float64}}:
 [-1095.4451129504978, -54772.25548320812]
 [1095.4451137838312, 54772.255524874796]
 [8.111114476617955, 0.02219298606763958]
 [-8.103507635567631, -0.02221382112196499]
```




