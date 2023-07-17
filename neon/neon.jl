using QCBase
using ClusterMeanField
using ActiveSpaceSolvers
using RDM
using InCoreIntegrals
using PyCall
using JLD2
using NPZ
    
## Integrals
pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("/Users/nicole/My Drive/code/FermiCG-data/neon/fcidump_nat.neon");
#ctx = fcidump.read("/Users/nicole/My Drive/code/FermiCG-data/neon/fcidump.neon");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, size(h,2))
ints = InCoreInts(ecore,h,g)
#@save "/Users/nicole/My Drive/code/FermiCG-data/neon/ints_nat.jld2" ints
#error("a")
    
clusters_in    = [(1:9),(10:18)]
n_clusters = 2
cluster_list = [collect(1:9), collect(10:18)]
clusters = [MOCluster(i,collect(cluster_list[i])) for i = 1:length(cluster_list)]
init_fspace = [ (5,5) for i in 1:n_clusters]
display(clusters)
init_cluster_ansatz = [RASCIAnsatz(9, 5, 5, (4,2,3), max_h=1, max_p=1), RASCIAnsatz(9,5,5,(4,2,3), max_h=1, max_p=1)]
display(init_cluster_ansatz)
    
ansatze = ActiveSpaceSolvers.generate_cluster_fock_ansatze(init_fspace, clusters, init_cluster_ansatz)
        
rdm1 = zeros(size(ints.h1))

#include("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci_natural_orbs.jl")
#rdms = ras_cluster_rdms(ints, clusters, RDM1(rdm1, rdm1), init_cluster_ansatz)
#@save "/Users/nicole/My Drive/code/FermiCG-data/neon/rdms.jld2" rdms
#error("hi")
        
#e_cmf, U_cmf, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, ansatze, RDM1(rdm1, rdm1), max_iter_oo=200, verbose=0, gconv=1e-6, method="cg");
e_cmf, U_cmf, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, ansatze, RDM1(rdm1, rdm1), max_iter_oo=200, verbose=3, gconv=1e-6, method="bfgs");
#e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, ansatze, RDM1(rdm1, rdm1), maxiter_oo = 200, verbose=3, zero_intra_rots = false, diis_start=1);






    
