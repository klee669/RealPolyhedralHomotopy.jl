using Test, RealPolyhedralHomotopy

@var x y;
F = System([-1 - 24000*y + x^3, -9 + 50*x*y - 1*y^2]);
result = certify_patchwork(F);

@test result == 1

result = certify_patchwork(F; Number_Real_Solutions = true);

@test result[2] == 4

B = generate_binomials(F);
@test length(B.binomial_system) == 2

realSols = rph_track(B,F)
@test length(realSols) == 4

println("tests completed successfully")
