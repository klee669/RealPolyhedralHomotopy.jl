export rph_track


"""
    rph_track(B::Binomial_system_data, F::System; Certification = false)

Return the output of tracking the real solutions of a given list of binomial systems to the target system.
# Arguments
* `B` : The object `Binomial_system_data` obtained from `generate_binomials(F)`.
* `F` : The target system for the real polyhedral homotopy.
```julia
@var x y;
F = System([-1 - 24000*y + x^3, -9 + 50*x*y - 1*y^2]);
B = generate_binomials(F);
realSols = rph_track(B,F)
```
```
4-element Vector{Vector{Float64}}:
 [-1095.4451129504978, -54772.25548320812]
 [1095.4451137838312, 54772.255524874796]
 [8.111114476617955, 0.02219298606763958]
 [-8.103507635567631, -0.02221382112196499]
```
* The optional argument `Certification` certifies that all real solutions to a patchedworked system are found.
This is done by an a posteriori certification for numerical approximations obtained by the real polyhedral homotopy.
When the real polyhedral homotopy root-finding is certified, it returns a list of solutions to the target and `1`; otherwise, it returns `0`. The default value for the option is false.
```julia
realSols = rph_track(B,F; Certification = true)
```
```
([[-1095.4451129504978, -54772.25548320812], [1095.4451137838312, 54772.255524874796], [8.111114476617955, 0.02219298606763958], [-8.103507635567631, -0.022213821121964985]], 1)
```
"""
function rph_track(BData::Binomial_system_data,F::System;Certification::Bool = false)

  neqs = length(F);
  n = neqs;
  varsF = variables(F);

  binomial_systems = BData.binomial_system;
  normal_vectors = BData.normal_vectors;
  l = BData.lifts;

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
        v = T_mons[j][:,1] - T_mons[j][:,2]
          for k in 1:length(v)
            D[k,j] = v[k]
          end
        B[j] = -1*T_coeffs[j][2]/T_coeffs[j][1]
      end
    D = transpose(D)
    H,U = hnf_with_transform(D)
    B_new = zeros(n);

    for l in 1:n
      v = Array(U[l, :])
      v1 = transpose(B).^v
      B_new[l] = prod(v1)
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

    for l in 1:length(real_R)
      append!(test, [[real(R[real_R[l]])]])
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



### Define homotopies

  # Define our homotopy variable
  @var t


  # Define matrices that are the monomial support and coefficients of F
  A = support_coefficients(F)[1];
  B = support_coefficients(F)[2];






  #Create homotopy for each mixed cell

  homotopies = [];
  c = BData.cells;
  for k in 1:length(c)
    rn = normal_vectors[k]

    #Find exponents for homotopy variable t
    tvecs = [];
    A = support_coefficients(F)[1];
    B = support_coefficients(F)[2];
    Anew = A;
    for i in 1:n
      tlist = transpose(A[i])*rn + l[i]
      m = minimum(tlist)
      tlist = tlist - m*ones(length(tlist))
      tlist = round.(tlist);
      tlist = convert.(Int, tlist)
      append!(tvecs, [tlist])
      Anew[i] = [A[i] ; tlist']
    end

    varsFt = collect(Iterators.flatten([varsF, t]))
    eqs = [];
    for i in 1:n
      varmat = varsFt.^Anew[i]
      eq = sum(B[i][j]*prod(varmat[:, j]) for j in 1:size(varmat)[2])
      append!(eqs, [eq])
    end

    #send t->1-t since Homotopy goes from 1 to 0
    for i in 1:length(eqs)
      eqs[i] = subs(eqs[i], t=>(1-t))
    end


    append!(homotopies, [eqs])

  end


  # Change homotopies from type Any to type Expression
  for i in 1:length(homotopies)
    homotopies[i] = convert(Array{Expression,1}, homotopies[i])
  end

  # Create our homotopy systems

  homotopy_systems = [];
  for i in 1:ncells
    append!(homotopy_systems,  [Homotopy(homotopies[i], varsF, t)])
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
