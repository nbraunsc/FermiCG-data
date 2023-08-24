using QCBase
using ClusterMeanField
using ActiveSpaceSolvers
using RDM
using InCoreIntegrals
using PyCall
using JLD2
using NPZ

molecule = "
C       -6.8604981940     -0.4376683883      0.0000000000
C       -7.0417938838      0.9473143671      0.0000000000
C       -5.5765818320     -0.9841959654      0.0000000000
C       -5.9160648559      1.7742824642      0.0000000000
C       -4.4255614400     -0.1700917199      0.0000000000
C       -4.6336867757      1.2242295397      0.0000000000
C       -3.0514687259     -0.7633813258      0.0000000000
C       -2.8520613255     -2.1570405955      0.0000000000
C       -1.5640941394     -2.6887327507      0.0000000000
C       -0.4451033014     -1.8584342989      0.0000000000
C       -0.5927492072     -0.4581981628      0.0000000000
C       -1.9041039939      0.0505338987      0.0000000000
H       -8.0463585128      1.3756739224      0.0000000000
H       -5.4821758879     -2.0696091478      0.0000000000
H       -6.0325833442      2.8604965477      0.0000000000
H       -3.7849879120      1.9076345469      0.0000000000
H       -3.6983932318     -2.8425398416      0.0000000000
H       -1.4298106521     -3.7726953568      0.0000000000
H        0.5426199813     -2.3168156796      0.0000000000
H       -2.0368540411      1.1280906363      0.0000000000
H       -7.7260669052     -1.1042415323      0.0000000000
H        0.22700           0.16873           0.00000
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
npzwrite("C_scf.npy", C)
npzwrite("S_scf.npy", S)


println(size(C))

ndocc = pymol.nelec[2]
nsing = pymol.nelec[1] - ndocc
nvirt = size(C,1) - ndocc - nsing
println(ndocc)
println(nsing)
println(nvirt)

println(mf.mo_energy)

tools.molden.from_mo(pymol, "/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci/two_benzene/mo_coeffs.molden", C)
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

