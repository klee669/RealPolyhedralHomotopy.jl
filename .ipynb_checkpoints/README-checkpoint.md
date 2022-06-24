# RealPolyhedralHomotopy.jl
## A package for finding real roots of systems of polynomial equations using polyhedral homotopy.

## Usage
To run a notebook with RealPolyhedralHomotopy:
- Clone this repository.
    - Create a new julia notebook inside the `RealPolyhedralHomotopy` folder and put
    ```
    using Pkg
    Pkg.activate(".")
    using RealPolyhedralHomotopy
    ```
    inside the first cell.
    - or run Julia in the `RealPolyhedralHomotopy` folder and run the commands
    ```
    using Pkg
    Pkg.activate(".")
    using RealPolyhedralHomotopy
    ```

Polynomials can be constructed by declaring variables with `@polyvar` and including them into a Julia expression.
The Package exports the functions `generate_binomials`, `certify_patchwork` and `rph_track`.

### Example
Input:
```julia
using Pkg
Pkg.activate(".")
using RealPolyhedralHomotopy

@var x y;
F = System([-1 - 24000*y + x^3, -9 + 50*x*y - 1*y^2]);
#
B = generate_binomials(F)

# Tracks the the real solutions.
realSols = rph_track(B,F)

# Returns 1 if certified, 0 otherwise.
result = certify_patchwork(F)

```




