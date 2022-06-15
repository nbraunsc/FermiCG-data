using FermiCG
using JLD2

@load "/home/nbraunsc/FermiCG-data/he_clusters/07_oct/after_cmf.jld2"

#Build Clustered Operator
cluster_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

#Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, Da, Db);

#Need to find reference state 
ref_fock = FermiCG.FockConfig(init_fspace)
nroots = 2
#ci_vector = FermiCG.TPSCIstate(clusters, ref_fock, R=nroots)
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
#ci_vector = FermiCG.ClusteredState(clusters, ref_fock, R=nroots);
#Need to find the automated way to define these other excited configs away from ref state, example is to large
#to do by hand
#probably something to do with building p spaces and q spaces
ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,1,1])] = [1,0,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([2,1,1,1,1,1,1])] = [0,1,0,0,0,0,0,0]
#ci_vector[ref_fock][ClusterConfig([1,2,1,1,1,1,1])] = [0,0,1,0,0,0,0,0]
#ci_vector[ref_fock][ClusterConfig([1,1,2,1,1,1,1])] = [0,0,0,1,0,0,0,0]
#ci_vector[ref_fock][ClusterConfig([1,1,1,2,1,1,1])] = [0,0,0,0,1,0,0,0]
#ci_vector[ref_fock][ClusterConfig([1,1,1,1,2,1,1])] = [0,0,0,0,0,1,0,0]
#ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,2,1])] = [0,0,0,0,0,0,1,0]
#ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,1,2])] = [0,0,0,0,0,0,0,1]

display(ci_vector)

e0, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, cluster_ham,
                          thresh_cipsi=1e-3, # Threshold for adding to P-space
                          thresh_foi=1e-7,    # Threshold for keeping terms when defining FOIS
                          thresh_asci=0.001,     # Threshold of P-space configs to search from
                          max_iter=10,
                          matvec=3);
println(e0)

