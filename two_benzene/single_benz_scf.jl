using QCBase
using ClusterMeanField
using ActiveSpaceSolvers
using RDM
using InCoreIntegrals
using PyCall
using JLD2
using NPZ

molecule = "
  H      1.2194     -0.1652      2.1600
  C      0.6825     -0.0924      1.2087
  C     -0.7075     -0.0352      1.1973
  H     -1.2644     -0.0630      2.1393
  C     -1.3898      0.0572     -0.0114
  H     -2.4836      0.1021     -0.0204
  C     -0.6824      0.0925     -1.2088
  H     -1.2194      0.1652     -2.1599
  C      0.7075      0.0352     -1.1973
  H      1.2641      0.0628     -2.1395
  C      1.3899     -0.0572      0.0114
  H      2.4836     -0.1022      0.0205
"
atoms = []
for (li,line) in enumerate(split(rstrip(lstrip(molecule)), "\n"))
    l = split(line)
    push!(atoms, Atom(li, l[1], parse.(Float64,l[2:4])))
end

basis = "sto-3g"
# Create FermiCG.Molecule type
mol     = Molecule(0, 1, atoms,basis);

pyscf = pyimport("pyscf")
tools = pyimport("pyscf.tools")
fcidump = pyimport("pyscf.tools.fcidump");

pymol = pyscf.gto.Mole(atom=molecule,
                            symmetry = false, spin =0,charge=0,
                            basis = basis)
pymol.build()
mf = pyscf.scf.RHF(pymol).run()
S = mf.get_ovlp(pymol)
C = mf.mo_coeff
#npzwrite("C_scf.npy", C)
#npzwrite("S_scf.npy", S)
println(mf.mo_energy)

println(size(C))

ndocc = pymol.nelec[2]
nsing = pymol.nelec[1] - ndocc
nvirt = size(C,1) - ndocc - nsing
println(ndocc)
println(nsing)
println(nvirt)



tools.molden.from_mo(pymol, "/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci/two_benzene/mo_coeffs_single.molden", C)
error("a")
# Just use alpha orbitals
Cdocc = mf.mo_coeff[:,collect(1:ndocc)]
Csing = mf.mo_coeff[:,collect(ndocc:ndocc+nsing)]
val = ndocc+nsing
val2 = ndocc+nsing+nvirt
Cvirt = mf.mo_coeff[:,collect(val:val2)]
println(size(Cdocc))
println(size(Csing))
println(size(Cvirt))

nbas = size(Cdocc,1)

@save "/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci/two_benzene/scf.jld2" C S pymol mf Cdocc Csing Cvirt 

