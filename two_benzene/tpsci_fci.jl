using FermiCG, NPZ, JLD2
using LinearAlgebra
using Printf
using QCBase
using RDM
using ActiveSpaceSolvers
#using ClusterMeanField

#@load "/Users/nicole/My Drive/code/FermiCG-data/two_benzene/cmf_fci.jld2"
@load "/Users/nicole/My Drive/code/FermiCG-data/two_benzene/cmf_fci_100A.jld2" 

ints = deepcopy(ints_cmf)
#C = deepcopy(C_cmf);
ecore = ints.h0

M = 400

init_fspace = [ (1,1), (2,2), (0,0), (1,1), (2,2), (0,0) ]
ref_fock = FockConfig(init_fspace)
ansatze = [FCIAnsatz(1,1,1), FCIAnsatz(4,2,2), FCIAnsatz(1,0,0), FCIAnsatz(1,1,1), FCIAnsatz(4,2,2), FCIAnsatz(1,0,0)]
#ansatze = [FCIAnsatz(6, 3, 3), FCIAnsatz(6,3,3)]
#
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [6,6,6,6,6,6], ref_fock, ansatze, max_roots=M, verbose=0);

nroots = 6
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
FermiCG.add_fock_configs_for_rasci(ref_fock, ci_vector, ex_level=1)
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

#ci_vector = FermiCG.add_spin_focksectors(ci_vector)

# Add the lowest energy single exciton to basis
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([2,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,2])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([3,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,3])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([4,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,4])] = zeros(Float64,nroots)

e0b, v0b = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham, max_iter=200, conv_thresh=1e-8);

#@time e2 = FermiCG.compute_pt2_energy(v0b, cluster_ops, clustered_ham, thresh_foi=1e-8, verbose=0)

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
@save "/Users/nicole/My Drive/code/FermiCG-data/two_benzene/fci_tpsci_100A.jld2" cluster_bases e0b v0b s2 ecore
