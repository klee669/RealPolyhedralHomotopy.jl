export generate_binomials

"""
    generate_binomials(poly_system::System)

Return a collection of binomial systems from an input polynomial system. Binomial systems obtained from the mixed cells induced by the ``\Log|C|``-lifting.
"""
function generate_binomials(poly_system)

  F = poly_system;
  neqs = length(F);
  varsF = variables(F);

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
  for i in 1:ncells
    system = [];
    mons = indices(cells[i])
    for j in 1:neqs
      mons2 = mons[j]
      bi1 = transpose(A[j])[mons2[1]:mons2[1],1:end]
      bi2 = transpose(A[j])[mons2[2]:mons2[2],1:end]
      term1 = B[j][mons2[1]]*prod(varsF.^(transpose(bi1)))
      term2 = B[j][mons2[2]]*prod(varsF.^(transpose(bi2)))
      p = term1 + term2;
      append!(system,p)
    end
    append!(binomial_systems, [system])
  end

  # Convert binomial systems from type Any to type Expression
  for i in 1:length(binomial_systems)
    binomial_systems[i] = convert(Array{Expression,1}, binomial_systems[i])
  end

  return binomial_systems

end
