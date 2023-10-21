export generate_binomials, Binomial_system_data

"""
    generate_binomials(F::System)

Return a wrapper object Binomial_system_data from an input polynomial system. Binomial systems obtained from the mixed cells induced by the ``Log|C|``-lifting.
The object Binomial_system_data contains the binomial system, normal vectors, lifting vectors, and cells from the mixed subdivision computed.
The object Binomial_system_data is used as an input for the rph_track function.
# Arguments
* `F` : The target system for the real polyhedral homotopy. 
```julia
B = generate_binomials(F)
```
```
Binomial_system_data
```
```julia
B.binomial_system
```
```
2-element Vector{Any}:
 Expression[-24000*y + x^3, 50*x*y - y^2]
 Expression[-24000*y + x^3, -9 + 50*x*y]
```
"""

struct Binomial_system_data
  binomial_system::Vector{Any}
  normal_vectors::Vector{Any}
  lifts::Vector{Vector{Int64}}
  cells::Vector{MixedCell}
end

function Base.show(io::IO,x::Binomial_system_data)
  print(io,"Binomial_system_data")
end


function generate_binomials(F::System)

  x = variables(F);
  neqs = length(F);

  # Define matrices that are the monomial support and coefficients of F
  A = support_coefficients(F)[1];
  B = support_coefficients(F)[2];

  # Use Log(|C|) to define lift
  w1 = round.(-1*(10^6)*log.(abs.(support_coefficients(F)[2][1])));
  w1 = convert.(Int,w1);
  lifts = [w1];

  for i in 2:neqs
    w = round.(-1*(10^6)*log.(abs.(support_coefficients(F)[2][i])))
    w = convert.(Int,w)
    append!(lifts, [w])
  end

  # Compute mixed cells
  cells = mixed_cells(A, lifts);
  ncells = length(cells);


  # Define binomial systems from mixed cells
  binomial_systems = [];
  normal_vectors = [];
  for i in 1:ncells
    system = [];
    mons = indices(cells[i])
    for j in 1:neqs
      mons2 = mons[j]
      bi1 = transpose(A[j])[mons2[1]:mons2[1],1:end]
      bi2 = transpose(A[j])[mons2[2]:mons2[2],1:end]
      term1 = B[j][mons2[1]]*prod(x.^(transpose(bi1)))
      term2 = B[j][mons2[2]]*prod(x.^(transpose(bi2)))
      p = term1 + term2;
      append!(system,p)
    end
    append!(binomial_systems, [system])
    append!(normal_vectors, [normal(cells[i])])
  end

  # Convert binomial systems from type Any to type Expression
  for i in 1:length(binomial_systems)
    binomial_systems[i] = convert(Array{Expression,1}, binomial_systems[i])
  end

  return Binomial_system_data(binomial_systems, normal_vectors, lifts, cells)

end

