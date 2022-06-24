using MixedSubdivisions, HomotopyContinuation, DynamicPolynomials, LinearAlgebra, AbstractAlgebra, PolynomialRoots

function rph_track(binomialSystem, targetSystem;Certification::Bool = false)
  F = targetSystem;
  neqs = length(F);
  n = neqs;
  varsF = variables(F);
  F_eqs = expressions(F);

  binomial_systems = binomialSystem;

  ncells = length(binomial_systems);

  # Find real solutions to binomial systems
  r_binomial_sols = [];


  for i in 1:ncells
    T = System(binomial_systems[i])
    D = zero_matrix(ZZ, n, n)
    B = zeros(n)
    T_mons = support_coefficients(T)[1]
    T_coeffs = support_coefficients(T)[2]
      for j in 1:n
        v = T_mons[j][1:end,1:1] - T_mons[j][1:end,2:2]
          for k in 1:length(v)
            D[k,j] = v[k]
          end
        B[j] = -1*T_coeffs[j][2]/T_coeffs[j][1]
      end
    D = transpose(D)
    H,U = hnf_with_transform(D)
    B_new = zeros(n);
    for i in 1:n
      v = Array(U[i:i, 1:end])
      v1 = transpose(B).^v
      B_new[i] = prod(v1)
    end

    #Create dictionary to store real solutions
    sols = Dict()

    #Initialize first entry of sols
    poly_root = zeros(H[n,n]+1)
    poly_root[end] = 1;
    poly_root[1] = -1*B_new[end]
    R = PolynomialRoots.roots(poly_root)

    real_R = findall(x->norm(x)<10^(-10), imag.(R));
    test = [];
    for i in 1:length(real_R)
      append!(test, [[real(R[real_R[i]])]])
    end
    sols[n] = test

    for k in 1:n-1
      poly_root = im*zeros(H[n-k,n-k]+1)
      test =[];

      for j in 1:length(sols[n-k+1])
        poly_root[end] = prod(sols[n-k+1][j][end-l]^H[n-k,n-l] for l in 0:k-1)
        poly_root[1] = -1*B_new[n-k]
        R = PolynomialRoots.roots(poly_root)
        real_R = findall(x->norm(x)<10^(-10), imag.(R));
        for l in 1:length(real_R)
          append!(test, [[real.(R[real_R[l]]),sols[n-k+1][j]]])
        end
      end

      for l in 1:length(test)
        test[l] = collect(Iterators.flatten(test[l]))
      end

      sols[n-k] = test
    end
    append!(r_binomial_sols, [sols[1]])
  end



  # G is an array of the polynomial systems that are the set of monomials not in each binomial system
  G = [];
  for i in 1:ncells
    append!(G, [F_eqs - binomial_systems[i]])
  end

  # Define our homotopy variable
  @var t

  # Create systems consisting of monomials not in binomial system multiplied by (1-t)^a for a positive integer a
  homotopies = [];
  for k in 1:ncells
  J = System(G[k])
  G_sup1 = support_coefficients(J)[1];
  G_sup2 = support_coefficients(J)[2];
    h1 = [];
    for i in 1:neqs
      t1, t2 = size(transpose(G_sup1[i]))
      monomial_list = [];
        for j in 1:t1
          a = abs.(rand(Int8))%10 + 1
          mon_vec = transpose(G_sup1[i])[j:j,1:end]
          mon1 = G_sup2[i][j]*prod(varsF.^(transpose(mon_vec)))
          append!(monomial_list, mon1*((1-t)^a))
        end
      append!(h1, [sum(monomial_list)])
      end
      append!(homotopies, [h1])
  end

  # Change homotopies from type Any to type Expression
  for i in 1:length(homotopies)
    homotopies[i] = convert(Array{Expression,1}, homotopies[i])
  end

  # Create our homotopy systems by adding each binomial system to other systems that have monomials with (1-t)^a attached
  homotopy_systems = [];
  for i in 1:ncells
    combined_sys = homotopies[i] + binomial_systems[i]
    append!(homotopy_systems,  [Homotopy(combined_sys, varsF, t)])
  end



  #This gives the real solutions after only tracking the real solutions from the binomial systems
  real_sols = [];

  # Track real solutions from binomial systems to target system
  for i in 1:length(r_binomial_sols)
    H1 = InterpretedHomotopy(homotopy_systems[i])
    R2 = HomotopyContinuation.solve(H1, r_binomial_sols[i],show_progress=false)
    append!(real_sols, real_solutions(R2))
  end

if Certification == true         
        for i in 1:ncells
        cr = certify(binomial_systems[i],convert(Array{Vector{Float64}},r_binomial_sols[i]));
            if nreal_certified(cr) < length(r_binomial_sols[i])
                return 0
            end
        end
        for i in 1:length(real_sols)
        cr = certify(F,real_sols[i]);
            if nreal_certified(cr) < 1
                return 0
            end
        end
        return (convert(Array{Vector{Float64}},real_sols),1)
else
    return convert(Array{Vector{Float64}},real_sols)
end
end
