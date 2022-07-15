using FermiCG
using JLD2

@load "/home/nbraunsc/FermiCG-data/he_clusters/06_oct/after_cmf.jld2"

#Build Clustered Operator
cluster_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

#Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, Da, Db);

#Need to find reference state 
ref_fock = FermiCG.FockConfig(init_fspace)
nroots = 3
#ci_vector = FermiCG.TPSCIstate(clusters, ref_fock, R=nroots)
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
#ci_vector = FermiCG.ClusteredState(clusters, ref_fock, R=nroots);
#Need to find the automated way to define these other excited configs away from ref state, example is to large
#to do by hand
#probably something to do with building p spaces and q spaces
ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,1])] = [1,0,0]
ci_vector[ref_fock][ClusterConfig([1,2,1,1,1,1])] = [0,1,0]
ci_vector[ref_fock][ClusterConfig([1,1,3,1,1,1])] = [0,0,1]
#ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,1,1])] = [1,0,0,0,0,0,0,0]
#ci_vector[ref_fock][ClusterConfig([2,1,1,1,1,1,1])] = [0,1,0,0,0,0,0,0]
#ci_vector[ref_fock][ClusterConfig([1,2,1,1,1,1,1])] = [0,0,1,0,0,0,0,0]
#ci_vector[ref_fock][ClusterConfig([1,1,2,1,1,1,1])] = [0,0,0,1,0,0,0,0]
#ci_vector[ref_fock][ClusterConfig([1,1,1,2,1,1,1])] = [0,0,0,0,1,0,0,0]
#ci_vector[ref_fock][ClusterConfig([1,1,1,1,2,1,1])] = [0,0,0,0,0,1,0,0]
#ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,2,1])] = [0,0,0,0,0,0,1,0]
#ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,1,2])] = [0,0,0,0,0,0,0,1]

e0, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, cluster_ham,
                          thresh_cipsi=0.0001, # Threshold for adding to P-space
                          thresh_foi=1e-2,    # Threshold for keeping terms when defining FOIS
                          thresh_asci=0,     # Threshold of P-space configs to search from
                          max_iter=10);

@time e2 = FermiCG.compute_pt2_energy(v0, cluster_ops, cluster_ham, thresh_foi=1e-8)

e0_001, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, cluster_ham,
                          thresh_cipsi=0.001, # Threshold for adding to P-space
                          thresh_foi=1e-2,    # Threshold for keeping terms when defining FOIS
                          thresh_asci=0,     # Threshold of P-space configs to search from
                          max_iter=10);

@time e2_001 = FermiCG.compute_pt2_energy(v0, cluster_ops, cluster_ham, thresh_foi=1e-8)


e0_01, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, cluster_ham,
                          thresh_cipsi=0.01, # Threshold for adding to P-space
                          thresh_foi=1e-2,    # Threshold for keeping terms when defining FOIS
                          thresh_asci=0,     # Threshold of P-space configs to search from
                          max_iter=10);

@time e2_01 = FermiCG.compute_pt2_energy(v0, cluster_ops, cluster_ham, thresh_foi=1e-8)

println("0.0001")
println(e0)
println(e2)

println("0.001")
println(e0_001)
println(e2_001)

println("0.01")
println(e0_01)
println(e2_01)

@save "tpsci_scan.jld2" e0 e2 e0_001 e2_001 e0_01 e2_01
