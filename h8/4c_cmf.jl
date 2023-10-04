using QCBase
using ClusterMeanField
using ActiveSpaceSolvers
using RDM
using InCoreIntegrals
using PyCall
using JLD2
using NPZ

molecule = "
H       -3.7181016563     -0.4088529731      0.0000000000
H       -4.7240997577      0.4284082328      0.0000000000
H       -3.7712553004      1.4028748967      0.0000000000
H       -3.9915160957     -1.7092506766     -0.0000000000
H       -1.9286222252      0.2081363722      0.0000000000
H       -1.4735643684      1.4641152794      0.0000000000
H       -1.0663746577     -0.5427171043      0.0000000000
H       -1.7599280378     -1.6484585439      0.0000000000
"
atoms = []
for (li,line) in enumerate(split(rstrip(lstrip(molecule)), "\n"))
    l = split(line)
    push!(atoms, Atom(li, l[1], parse.(Float64,l[2:4])))
end

basis = "sto-3g"
# Create FermiCG.Molecule type
pyscf = pyimport("pyscf")
tools = pyimport("pyscf.tools")
fcidump = pyimport("pyscf.tools.fcidump");

pymol = pyscf.gto.Mole(atom=molecule,
                       symmetry = false, spin =0,charge=0,
                       basis = basis)
pymol.build()

h0 = npzread("/Users/nicole/My Drive/code/FermiCG-data/h8/ints_h0.npy")
h1 = npzread("/Users/nicole/My Drive/code/FermiCG-data/h8/ints_h1.npy")
h2 = npzread("/Users/nicole/My Drive/code/FermiCG-data/h8/ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

Cact = npzread("/Users/nicole/My Drive/code/FermiCG-data/h8/mo_coeffs.npy")

n_clusters = 4
cluster_list = [[1,2], [3,4], [5,6], [7,8]]
clusters = [MOCluster(i,collect(cluster_list[i])) for i = 1:length(cluster_list)]
init_fspace = [ (2,2),(0,0), (2,2), (0,0)]
display(clusters)
ansatze = [FCIAnsatz(2, 2, 2), FCIAnsatz(2,0,0), FCIAnsatz(2,2,2), FCIAnsatz(2,0,0)]
display(ansatze)

rdm1 = zeros(size(ints.h1))
#run cmf_oo
e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, ansatze, RDM1(rdm1, rdm1), maxiter_oo = 100, tol_oo=1e-8, tol_ci=1e-10, diis_start=1);
ints_cmf = orbital_rotation(ints,U_cmf)
C_cmf = Cact*U_cmf
@save "/Users/nicole/My Drive/code/FermiCG-data/h8/cmf_fci_4c.jld2" ints_cmf U_cmf d1 e_cmf C_cmf ansatze init_fspace clusters
pyscf.tools.molden.from_mo(pymol, "/Users/nicole/My Drive/code/FermiCG-data/h8/C_cmf_fci_4c.molden", C_cmf)





