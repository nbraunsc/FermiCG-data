using FermiCG, NPZ, JLD2
using LinearAlgebra
using Printf
using QCBase
using RDM
using ActiveSpaceSolvers
#using ClusterMeanField

@load "/Users/nicole/My Drive/code/FermiCG-data/h8/cmf_fci_4c.jld2" 

ints = deepcopy(ints_cmf)
ecore = ints.h0

M = 500

n_clusters = 4
cluster_list = [[1,2], [3,4], [5,6], [7,8]]
clusters = [MOCluster(i,collect(cluster_list[i])) for i = 1:length(cluster_list)]
init_fspace = [ (2,2),(0,0), (2,2), (0,0)]
display(clusters)
ansatze = [FCIAnsatz(2, 2, 2), FCIAnsatz(2,0,0), FCIAnsatz(2,2,2), FCIAnsatz(2,0,0)]
display(ansatze)


ref_fock = FockConfig(init_fspace)
#
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [4,4,4,4], ref_fock, max_roots=M, verbose=1)
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [4,4,4,4], ref_fock, ansatze, max_roots=M, verbose=1)
#cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, d1, [4,4,4,4], ref_fock, ansatze, max_roots=M, verbose=1)
#
#cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, init_fspace=init_fspace, delta_elec=8, rdm1a=d1.a, rdm1b=d1.b, max_roots=M, verbose=1);

nroots = 10
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
#FermiCG.add_fock_configs_for_rasci(ref_fock, ci_vector, ex_level="hh", n_clusters=6)
#FermiCG.expand_each_fock_space!(ci_vector, cluster_bases)
FermiCG.expand_to_full_space!(ci_vector, cluster_bases, 4,4)
FermiCG.eye!(ci_vector)

# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
#@save "/Users/nicole/My Drive/code/FermiCG-data/h8/tps_fci_4c_ham_woutremoval.jld2" cluster_bases clustered_ham clusters

#
# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

#
# Add cmf hamiltonians for doing MP-style PT2 
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b, verbose=0);

e0b, v0b = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham, max_iter=200, conv_thresh=1e-8);


println()
println("   *======TPSCI results======*")
@printf("TCI Thresh: %8.6f  Dim:%8d\n",0.0001,size(v0b)[1])
println()
@printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)")
for r in 1:nroots
    @printf("TCI %5s %12.8f %12.8f\n",r, e0b[r] + ecore, e0b[r] + ecore)
    #@printf("TCI %5s %12.8f %12.8f\n",r, e0b[r] + ecore, e0b[r] + e2[r] + ecore)
end

clustered_S2 = FermiCG.extract_S2(ci_vector.clusters)

println()
println("   *======TPSCI S2 results======*")
@printf(" %-50s", "Compute FINAL S2 expectation values: ")
@time s2 = FermiCG.compute_expectation_value_parallel(v0b, cluster_ops, clustered_S2)

@printf(" %5s %12s %12s\n", "Root", "Energy", "S2")
for r in 1:nroots
    @printf(" %5s %12.8f %12.8f\n",r, e0b[r]+ecore, abs(s2[r]))
end
@save "/Users/nicole/My Drive/code/FermiCG-data/h8/tps_fci_4c_new2.jld2" cluster_bases e0b v0b s2 ecore
#@save "/Users/nicole/My Drive/code/FermiCG-data/h8/tps_fci_4c_no4body.jld2" cluster_bases e0b v0b s2 ecore clustered_ham
