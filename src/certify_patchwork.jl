export certify_patchwork

function certify_patchwork(polySystem;Number_Real_Solutions::Bool = false)
  F = polySystem
  neqs = length(F)
  n =neqs
  varsF = variables(F);

    
  # Define matrices that are the monomial support and coefficients of F
  A = support_coefficients(F)[1];
  B = support_coefficients(F)[2];




  vB = reduce(vcat, B);


  # Use Log(|C|) to define lift
  l1 = round.(-1*(10^6)*log.(abs.(support_coefficients(F)[2][1])));
  l1 = convert.(Int,l1);
  lifts = [l1];
  for i in 2:neqs
    l = round.(-1*(10^6)*log.(abs.(support_coefficients(F)[2][i])))
    l = convert.(Int,l)
    append!(lifts, [l])
  end

  # Compute mixed cells
  cells = mixed_cells(A, lifts);
  ncells = length(cells);



  #Construct the Cayley matrix
  mats = [];
  for i in 1:length(A)
    sz = size(A[i])[2]
    m1 = A[i];
    m2 = zeros(i-1, sz);
    m3 = ones(1, sz);
    m4 = zeros(neqs - i, sz)
    M = [m1 ; m2 ; m3 ; m4]
    append!(mats, [M])
  end

  M = reduce(hcat, mats);

  ## Make inequality for each patchworked system
  failed_cells = [];
  success_cells = [];
  scales = [];
  for i in 1:ncells
    in_cols = [];
    mixedCells = indices(cells[i])
    for j in 1:n
      offset = sum(size(A[k])[2] for k in 1:j) - size(A[j])[2]
      col1 = mixedCells[j][1] + offset
      col2 = mixedCells[j][2] + offset
      append!(in_cols, col1)
      append!(in_cols, col2)
    end

    out_cols = [];
    for j in 1:size(M)[2]
      t = findall(x->x==j, in_cols)
      if length(t) == 0
        append!(out_cols, j)
      end
    end

    fails = [];
    for j in 1:length(out_cols)
      cols = vcat(in_cols, out_cols[j]);
      sort!(cols);
      M_cells = M[1:end, cols];
      null = nullspace(M_cells);
      vBmod = vB[cols];
      lhs = abs.(dot(null, log.(abs.(vBmod))))
      rhs = log(size(M)[2])*norm(null,1)
      if lhs < rhs
        append!(fails, 1)
        append!(scales, rhs/lhs)
      end
    end

    if length(fails) == 0
      append!(success_cells, i)
    else
      append!(failed_cells, i)
      #println(failed)
    end
  end


        
##### Constructing binomial systems begins #####    

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

##### Constructing binomial systems ends #####    

    
##### Finding solutions for binomial systems begins #####    
  # Find real solutions to binomial systems
  r_binomial_sols = [];


  for i in 1:ncells
    T = System(binomial_systems[i])
    D = zero_matrix(ZZ, n, n)
    Bz = zeros(n)
    T_mons = support_coefficients(T)[1]
    T_coeffs = support_coefficients(T)[2]
      for j in 1:n
        v = T_mons[j][1:end,1:1] - T_mons[j][1:end,2:2]
          for k in 1:length(v)
            D[k,j] = v[k]
          end
        Bz[j] = -1*T_coeffs[j][2]/T_coeffs[j][1]
      end
    D = transpose(D)
    H,U = hnf_with_transform(D)
    Bz_new = zeros(n);
    for i in 1:n
      v = Array(U[i:i, 1:end])
      v1 = transpose(Bz).^v
      Bz_new[i] = prod(v1)
    end

    #Create dictionary to store real solutions
    sols = Dict()

    #Initialize first entry of sols
    poly_root = zeros(H[n,n]+1)
    poly_root[end] = 1;
    poly_root[1] = -1*Bz_new[end]
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
        poly_root[1] = -1*Bz_new[n-k]
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
##### Finding solutions for binomial systems ends #####    
        
        
if Number_Real_Solutions == true         
  if length(failed_cells)>0
    return 0
  else
    return (1, length(collect(Iterators.flatten(r_binomial_sols))));
  end
else
  if length(failed_cells)>0
    return 0
  else
    return 1
  end
end
            
end


