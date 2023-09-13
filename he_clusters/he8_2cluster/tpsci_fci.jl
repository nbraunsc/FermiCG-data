using FermiCG, NPZ, JLD2
using LinearAlgebra
using Printf
using QCBase
using RDM
using ActiveSpaceSolvers

include("/Users/nicole/My Drive/code/FermiCG/src/fci_ras_fockspaces.jl")

@load "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/cmf_fci_3A.jld2"

ints = deepcopy(ints_cmf)
C = deepcopy(C_cmf);
ecore = ints.h0

M = 118

ref_fock = FockConfig(init_fspace)
ansatze = [FCIAnsatz(3,3,3), FCIAnsatz(2,1,1), FCIAnsatz(3,0,0), FCIAnsatz(3,3,3), FCIAnsatz(2,1,1), FCIAnsatz(3,0,0)]
#
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [8,8,8,8,8,8], ref_fock, ansatze, max_roots=M, verbose=0);

nroots = 19
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);

# Add the lowest energy single exciton to basis
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,1,1])] = zeros(Float64,nroots)

add_fock_configs_for_rasci(ref_fock, ci_vector)
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

e0b, v0b = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham,
                            incremental  = true,
                            thresh_cipsi = 0.0001,
                            thresh_foi   = 1e-6,
                            thresh_asci  = -1,
                            #threaded = false,
                            ci_max_iter = 200,
                            max_mem_ci = 30.0);

#@time e2 = FermiCG.compute_pt2_energy(v0b, cluster_ops, clustered_ham, thresh_foi=1e-8)
#@time e2 = FermiCG.compute_pt2_energy(v0b, cluster_ops, clustered_ham, thresh_foi=1e-8, verbose=0)
e2 = 0

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
#@save "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/fci_tpsci_notspincompleted_de4.jld2" cluster_bases e0b v0b s2 ecore
@save "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/fci_tpsci_19roots_3A.jld2" cluster_bases e0b v0b s2 e2 ecore
#@save "/Users/nicole/My Drive/code/FermiCG-data/he_clusters/he8_2cluster/fci_tpsci.jld2" cluster_bases e0b v0b s2 ecore
