using QCBase
using ClusterMeanField
using ActiveSpaceSolvers
using RDM
using InCoreIntegrals
using PyCall
using JLD2
using NPZ
using LinearAlgebra
#using FermiCG

orb_part = pyimport("orbitalpartitioning")

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

basis = "cc-pvdz"
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

ndocc = pymol.nelec[2]
nsing = pymol.nelec[1] - ndocc
nvirt = size(C,2) - ndocc - nsing
println(ndocc)
println(nsing)
println(nvirt)
# Just use alpha orbitals
Cdocc = mf.mo_coeff[:,collect(1:ndocc)]
Csing = mf.mo_coeff[:,collect(ndocc:ndocc+nsing)]
Cvirt = mf.mo_coeff[:,collect(ndocc+nsing:ndocc+nsing+nvirt)]

# Find AO's corresponding to atoms
full = []
frag1 = []
frag2 = []
tmp = mf.mol.ao_labels(fmt=false)
tmp_ao = []
tmp2 = []
for i in tmp
    tmp2 = [convert(Int, i[1]), i[2], i[3], i[4]]
    push!(tmp_ao, tmp2)
end
println(tmp_ao)

for (ao_idx,ao) in enumerate(tmp_ao)
    if ao[1] in (0,1,2,3,4,5)
        if ao[3] == "2p"
            if ao[4] == "z"
                print(ao)
                push!(frag1, ao_idx)
                push!(full, ao_idx)
            end
        end
    end

    if ao[1] in (6,7,8,9,10,11)
        if ao[3] == "2p"
            if ao[4] == "z"
                push!(frag2, ao_idx)
                push!(full, ao_idx)
            end
        end
    end
end

frags = [frag1, frag2]
display(frags)

# Define projectors
nbas = size(Cdocc,1)
X = sqrt(S)
X = Matrix(1.0I, nbas, nbas)
Pfull = X[:,full]  # non-orthogonal
Pf = []

for f in frags
    push!(Pf, X[:,f])
end

# Project MOs onto all Fragments
(Oact, Sact, Vact), (Cenv, Cerr, _) = orb_part.svd_subspace_partitioning((Cdocc, Csing, Cvirt), Pfull, S)
Cact = hcat(Oact, Sact, Vact)
display(Cact)
pyscf.tools.molden.from_mo(mf.mol, "Cact.molden", Cact)
print(" Should be 1: ", det(Cact'*S*Cact))
error("a")

# Project active orbitals onto fragments
init_fspace = []
clusters = []
Cfrags = []
orb_index = 1

for (fi,f) in enumerate(frags)
    println(" \nFragment: ", f)
    (Of, Sf, Vf), (_, _, _) = orb_part.svd_subspace_partitioning((Oact, Sact, Vact), Pf[fi], S)

    push!(Cfrags, hcat(Of, Sf, Vf))
    ndocc_f = size(Of,2)
    push!(init_space, (ndocc_f+size(Sf,2), ndocc_f))
    nmof = size(Of,2) + size(Sf,2) + size(Vf,2)
    push!(clusters, [collect(orb_index:orb_index+nmof)])
    orb_index += nmof
end

# Orthogonalize Fragment orbitals
Cfrags = orb_part.sym_ortho(Cfrags, S)

# Pseudo canonicalize fragments
Cfrags = orb_part.canonicalize(Cfrags, F)

#np.save("Cfrags", Cfrags)
Cact = hcat(Cfrags)

npzwrite("Cact.npy", Cact)

print("nick: ", svd(Cact'*S*Cact)[2])
# Write Molden files for visualization
pyscf.tools.molden.from_mo(mf.mol, "Pfull.molden", Pfull)
pyscf.tools.molden.from_mo(mf.mol, "Cact.molden", Cact)
pyscf.tools.molden.from_mo(mf.mol, "Cenv.molden", Cenv)
for i in 1:length(frags)
    pyscf.tools.molden.from_mo(mf.mol, "Cfrag_sym_pc%i.molden"%i, Cfrags[i])
    #pyscf.tools.molden.from_mo(mf.mol, "Cfrag%i.molden"%i, Cfrags[i])
end
println(" init_fspace = ", init_fspace)
println(" clusters    = ", clusters)

@save "/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci/two_benzene/C.jld2" Cact Cfrags

# Make Integrals
d1_embed = 2 * Cenv*Cenv'

h0 = pyscf.gto.mole.energy_nuc(mf.mol)
h  = pyscf.scf.hf.get_hcore(mf.mol)
j, k = pyscf.scf.hf.get_jk(mf.mol, d1_embed, hermi=1)

h0 += tr(d1_embed*( h + .5*j - .25*k))

h = Cact'*h*Cact;
j = Cact'*j*Cact;
k = Cact'*k*Cact;
nact = size(h,1)

h2 = pyscf.ao2mo.kernel(pymol, Cact, aosym="s4", compact=False)
size(h2) = (nact, nact, nact, nact)
# The use of d1_embed only really makes sense if it has zero electrons in the
# active space. Let's warn the user if that's not true

S = pymol.intor("int1e_ovlp_sph")
n_act = tr(S*d1_embed*S*Cact*Cact')
if abs(n_act) > 1e-8 == False
    print(n_act)
    error(" I found embedded electrons in the active space?!")
end

h1 = h + j - .5*k;

npzwrite("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci/two_benzene/ints_h0", h0)
npzwrite("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci/two_benzene/ints_h1", h1)
npzwrite("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci/two_benzene/ints_h2", h2)
npzwrite("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci/two_benzene/mo_coeffs", Cact)
npzwrite("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci/two_benzene/overlap_mat", S)
ints = InCoreInts(h0, h1, h2)
@save "/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci/two_benzene/ints.jld2" ints


Pa = mf.make_rdm1()[1]
Pb = mf.make_rdm1()[2]
#np.save("Pa", Cact.T @ S @ Pa @ S @ Cact)
#np.save("Pb", Cact.T @ S @ Pb @ S @ Cact)

#tools.fcidump.from_mo(pymol, "fcidump.he07_oct", lo_ao)
