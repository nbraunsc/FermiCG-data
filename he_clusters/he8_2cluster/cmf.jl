using QCBase
using ClusterMeanField
using ActiveSpaceSolvers
using RDM
using InCoreIntegrals
using PyCall
using JLD2
using NPZ

#include("../rasci_natural_orbs.jl")
#@load "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/cmf_ras_cas.jld2"
#@load "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/cmf_ras_diis.jld2"
#@load "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/cmf_ras_newton.jld2"

molecule = "
He       0.0000000000      0.0000000000      0.0000000000
He       0.0000000000      2.0000000000      0.0000000000
He       2.0000000000      0.0000000000      0.0000000000
He       2.0000000000      2.0000000000      0.0000000000
He       0.0000000000      0.0000000000      5.0000000000
He       0.0000000000      2.0000000000      5.0000000000
He       2.0000000000      0.0000000000      5.0000000000
He       2.0000000000      2.0000000000      5.0000000000
"
atoms = []
for (li,line) in enumerate(split(rstrip(lstrip(molecule)), "\n"))
    l = split(line)
    push!(atoms, Atom(li, l[1], parse.(Float64,l[2:4])))
end

basis = "6-31g*"
# Create FermiCG.Molecule type
pyscf = pyimport("pyscf")
tools = pyimport("pyscf.tools")
fcidump = pyimport("pyscf.tools.fcidump");

pymol = pyscf.gto.Mole(atom=molecule,
                       symmetry = false, spin =0,charge=0,
                       basis = basis)
pymol.build()

h0 = npzread("/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/ints_h0.npy")
h1 = npzread("/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/ints_h1.npy")
h2 = npzread("/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/ints_h2.npy")
ints = InCoreInts(h0, h1, h2)
#
Cact = npzread("/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/mo_coeffs.npy")

clusters_in    = [(1:8),(9:16)]
n_clusters = 2
cluster_list = [collect(1:8), collect(9:16)]
clusters = [MOCluster(i,collect(cluster_list[i])) for i = 1:length(cluster_list)]
init_fspace = [ (4,4) for i in 1:n_clusters]
display(clusters)
#ansatze = [RASCIAnsatz(8, 4, 4, (3,2,3), max_h=0, max_p=0), RASCIAnsatz(8,4,4,(3,2,3), max_h=0, max_p=0)] #CASCI
#ansatze = [RASCIAnsatz(8, 4, 4, (3,2,3), max_h=1, max_p=1), RASCIAnsatz(8,4,4,(3,2,3), max_h=1, max_p=1)] #Singles
ansatze = [RASCIAnsatz(8, 4, 4, (3,2,3), max_h=2, max_p=2), RASCIAnsatz(8,4,4,(3,2,3), max_h=2, max_p=2)] #Doubles
#ansatze = [RASCIAnsatz(8, 4, 4, (3,2,3), max_h=6, max_p=6), RASCIAnsatz(8,4,4,(3,2,3), max_h=6, max_p=6)] #FCI
display(ansatze)


#Cfrags = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci/two_benzene/Cfrags.npy")
#rdms = ras_cluster_rdms(ints, clusters, cluster_list, init_cluster_ansatz)

#for i in 1:n_clusters
#    Cfrags[i,:,:] .= Cfrags[i,:,:]*rdms[i]
#end

#CMF has default of delta_elec=0 for all clusters
#delta_elec = [1,1,1,1]
#ansatze = ActiveSpaceSolvers.generate_cluster_fock_ansatze(init_fspace, clusters, init_cluster_ansatz, delta_elec)

#run cmf_oo

#DIIS
rdm1 = zeros(size(ints.h1))
e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, ansatze, RDM1(rdm1, rdm1), maxiter_oo = 100, tol_oo=1e-8, tol_ci=1e-10, verbose=0, diis_start=1);

#Newton Step from with no guess
#e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_gd(ints, clusters, init_fspace, ansatze, RDM1(rdm1, rdm1), maxiter_oo = 100, tol_oo=1e-8, tol_ci=1e-10, verbose=0, zero_intra_rots = false);

#Newton Step
#e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_gd(ints_cmf, clusters, init_fspace, ansatze, RDM1(d1.a, d1.b), maxiter_oo = 600, tol_oo=1e-8, tol_ci=1e-10, verbose=0, zero_intra_rots_g = false);

#Gradient descent
#e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_gd(ints_cmf, clusters, init_fspace, ansatze, RDM1(d1.a, d1.b), maxiter_oo = 500, tol_oo=1e-8, tol_ci=1e-10, verbose=0, zero_intra_rots_g = true, orb_hess=false);
#
#ints_cmf = orbital_rotation(ints_cmf,U_cmf)
ints_cmf = orbital_rotation(ints,U_cmf)
C_cmf = Cact*U_cmf
#@save "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/cmf_ras_newton.jld2" ints_cmf U_cmf d1 e_cmf C_cmf ansatze init_fspace clusters
@save "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/cmf_ras_doubles_diis.jld2" ints_cmf U_cmf d1 e_cmf C_cmf ansatze init_fspace clusters

#@save "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/cmf_rascas_diis.jld2" ints_cmf U_cmf d1 e_cmf C_cmf ansatze init_fspace clusters
#pyscf.tools.molden.from_mo(pymol, "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/C_cmf_diis.molden", C_cmf)






