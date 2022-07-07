var documenterSearchIndex = {"docs":
[{"location":"#RealPolyhedralHomotopy.jl","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy.jl","text":"","category":"section"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"RealPolyhedralHomotopy.jl is a package for finding real roots of systems of polynomial equations using polyhedral homotopy. The package implements the algorithm for the real polyhedral homotopy establised in A Polyhedral Homotopy Algorithm For Real Zeros. The idea for the real polyhedral homotopy motivated from the celebrated article A Polyhedral Method for Solving Sparse Polynomial Systems to apply the polyhedral homotopy method for real roots finding.","category":"page"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"The authors of this package are","category":"page"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"Kisun Lee\nJulia Lindberg\nJose Israel Rodriguez","category":"page"},{"location":"#Installation","page":"RealPolyhedralHomotopy","title":"Installation","text":"","category":"section"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"The package can be installed via the Julia package manager","category":"page"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"pkg> add RealPolyhedralHomotopy","category":"page"},{"location":"#Introduction","page":"RealPolyhedralHomotopy","title":"Introduction","text":"","category":"section"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"We support system input through the HomotopyContinuation package.","category":"page"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"using RealPolyhedralHomotopy\n\n@var x y;\nF = System([-1 - 24000*y + x^3, -9 + 50*x*y - 1*y^2]);","category":"page"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"For finding real roots, the list of binomial systems corresponding to F is required as a start system.","category":"page"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"B = generate_binomials(F);\n","category":"page"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"2-element Vector{Any}:\n Expression[-24000*y + x^3, 50*x*y - y^2]\n Expression[-24000*y + x^3, -9 + 50*x*y]","category":"page"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"Using the function rph_track, real roots are found by tracking the real polyhedral homotopy.","category":"page"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"realSols = rph_track(B,F)","category":"page"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"4-element Vector{Vector{Float64}}:\n [-1095.4451129504978, -54772.25548320812]\n [1095.4451137838312, 54772.255524874796]\n [8.111114476617955, 0.02219298606763958]\n [-8.103507635567631, -0.02221382112196499]","category":"page"},{"location":"#Functions-for-the-real-polyhedral-homotopy","page":"RealPolyhedralHomotopy","title":"Functions for the real polyhedral homotopy","text":"","category":"section"},{"location":"","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy","text":"certify_patchwork\ngenerate_binomials\nrph_track","category":"page"},{"location":"#RealPolyhedralHomotopy.certify_patchwork","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy.certify_patchwork","text":"certify_patchwork(F::System; Number_Real_Solutions = false)\n\nCertify if a given system is patchworked that all real solutions can be found using the real polyhedral homotopy. It returns the value 1 if the system F is certified to be patchworked according to the certification inequality.  Otherwise, 0 is returned. \n\nArguments\n\nF : The target system for the real polyhedral homotopy. \n\n@var x y;\nF = System([-1 - 24000*y + x^3, -9 + 50*x*y - 1*y^2]);\nresult = certify_patchwork(F)\n\n1\n\nThere is an option Number_Real_Solutions returning (1,k) where k is number of real solutions to the target system when the target system is patchedworked. The default value for the option is false.\n\nresult = certify_patchwork(F; Number_Real_Solutions = true)\n\n(1,4)\n\n\n\n\n\n","category":"function"},{"location":"#RealPolyhedralHomotopy.generate_binomials","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy.generate_binomials","text":"generate_binomials(F::System)\n\nReturn a collection of binomial systems from an input polynomial system. Binomial systems obtained from the mixed cells induced by the LogC-lifting.\n\nArguments\n\nF : The target system for the real polyhedral homotopy. \n\nB = generate_binomials(F)\n\n2-element Vector{Any}:\n Expression[-24000*y + x^3, 50*x*y - y^2]\n Expression[-24000*y + x^3, -9 + 50*x*y]\n\n\n\n\n\n","category":"function"},{"location":"#RealPolyhedralHomotopy.rph_track","page":"RealPolyhedralHomotopy","title":"RealPolyhedralHomotopy.rph_track","text":"rph_track(B::Vector{Any}, F::System; Certification = false)\n\nReturn the output of tracking the real solutions of a given list of binomial systems to the target system.\n\nArguments\n\nB : The list of binomial systems obtained from generate_binomials(F).\nF : The target system for the real polyhedral homotopy. \n\nrealSols = rph_track(B,F)\n\n4-element Vector{Vector{Float64}}:\n [-1095.4451129504978, -54772.25548320812]\n [1095.4451137838312, 54772.255524874796]\n [8.111114476617955, 0.02219298606763958]\n [-8.103507635567631, -0.02221382112196499]\n\nThe optional argument Certification certifies that all real solutions to a patchedworked system are found.  \n\nThis is done by an a posteriori certification for numerical approximations obtained by the real polyhedral homotopy.  When the real polyhedral homotopy root-finding is certified, it returns a list of solutions to the target and 1; otherwise, it returns 0. The default value for the option is false.\n\nrealSols = rph_track(B,F; Certification = true)\n\n([[-1095.4451129504978, -54772.25548320812], [1095.4451137838312, 54772.255524874796], [8.111114476617955, 0.02219298606763958], [-8.103507635567631, -0.022213821121964985]], 1)\n\n\n\n\n\n","category":"function"}]
}
