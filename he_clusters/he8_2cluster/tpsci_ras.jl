using FermiCG, NPZ, JLD2
using LinearAlgebra
using Printf
using QCBase
using RDM
using ActiveSpaceSolvers

#@load "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/cmf_ras_doubles_diis.jld2"
#@load "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/cmf_rascas_diis.jld2"
#@load "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/cmf_fci_diis.jld2"
#@load "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/cmf_ras_cas_3A.jld2"
@load "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/cmf_fci_100A.jld2"

ints = deepcopy(ints_cmf)
C = deepcopy(C_cmf);
ecore = ints.h0

M = 118
#M = 16

n_clusters = 2
cluster_list = [collect(1:8), collect(9:16)]
clusters = [MOCluster(i,collect(cluster_list[i])) for i = 1:length(cluster_list)]
init_fspace = [ (4,4) for i in 1:n_clusters]

ref_fock = FockConfig(init_fspace)
#ansatze = [RASCIAnsatz(8, 4, 4, (3,2,3), max_h=1, max_p=0), RASCIAnsatz(8,4,4,(3,2,3), max_h=1, max_p=0)] #doubles
ansatze = [RASCIAnsatz(8, 4, 4, (3,2,3), max_h=1, max_p=1), RASCIAnsatz(8,4,4,(3,2,3), max_h=1, max_p=1)] #Singles
#
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [0,0], ref_fock, ansatze, max_roots=M, verbose=1);
#

nroots = 256
nroots = 3
#nroots = 472
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1])] = zeros(Float64,nroots)
FermiCG.expand_each_fock_space!(ci_vector, cluster_bases)
FermiCG.eye!(ci_vector)

# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

#
# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

#
# Add cmf hamiltonians for doing MP-style PT2 
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b, verbose=0);

e0b, v0b = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham, solver="arpack", max_iter=200, conv_thresh=1e-8);

#e0b, v0b = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham,
#                            incremental  = true,
#                            thresh_cipsi = 0.0001,
#                            thresh_foi   = 1e-5,
#                            thresh_asci  = -1,
#                            ci_max_iter = 200,
#                            max_mem_ci = 50.0);
#
#@time e2 = FermiCG.compute_pt2_energy(v0b, cluster_ops, clustered_ham, thresh_foi=1e-8, verbose=0)

#println()
#println("   *======TPSCI results======*")
#@printf("TCI Thresh: %8.6f  Dim:%8d\n",0.0001,size(v0b)[1])
#println()
#@printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)")
#for r in 1:nroots
#    @printf("TCI %5s %12.8f %12.8f\n",r, e0b[r] + ecore, e0b[r] + ecore)
#end
#
clustered_S2 = FermiCG.extract_S2(ci_vector.clusters)
#
#println()
#println("   *======TPSCI S2 results======*")
#@printf(" %-50s", "Compute FINAL S2 expectation values: ")
@time s2 = FermiCG.compute_expectation_value_parallel(v0b, cluster_ops, clustered_S2)
#
#@printf(" %5s %12s %12s\n", "Root", "Energy", "S2")
#for r in 1:nroots
#    @printf(" %5s %12.8f %12.8f\n",r, e0b[r]+ecore, abs(s2[r]))
#end
#@save "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/ras_tpsci_3A_1root.jld2" cluster_bases e0b v0b s2 ecore
#@save "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/ras_tpsci_100A_1root.jld2" e0b v0b s2
@save "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/ras_tpsci_100A_3roots.jld2" e0b v0b s2

