using FermiCG
using JLD2

@load "/home/nbraunsc/FermiCG-data/he_clusters/08_01/after_cmf.jld2"

#Build Clustered Operator
cluster_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

#Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, Da, Db);

#Need to find reference state 
ref_fock = FermiCG.FockConfig(init_fspace)
nroots = 9
#ci_vector = FermiCG.TPSCIstate(clusters, ref_fock, R=nroots)
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
#ci_vector = FermiCG.ClusteredState(clusters, ref_fock, R=nroots);
#Need to find the automated way to define these other excited configs away from ref state, example is to large
#to do by hand
#probably something to do with building p spaces and q spaces
ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,1,1,1])] = [1,0,0,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([4,1,1,1,1,1,1,1])] = [0,1,0,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,4,1,1,1,1,1,1])] = [0,0,1,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,1,4,1,1,1,1,1])] = [0,0,0,1,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,1,1,4,1,1,1,1])] = [0,0,0,0,1,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,1,1,1,4,1,1,1])] = [0,0,0,0,0,1,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,4,1,1])] = [0,0,0,0,0,0,1,0,0]
ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,1,4,1])] = [0,0,0,0,0,0,0,1,0]
ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,1,1,4])] = [0,0,0,0,0,0,0,0,1]

display(ci_vector)

e0, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, cluster_ham,
                          thresh_cipsi=1e-3, # Threshold for adding to P-space
                          thresh_foi=1e-7,    # Threshold for keeping terms when defining FOIS
                          thresh_asci=0.001,     # Threshold of P-space configs to search from
                          max_iter=10,
                          matvec=3);
println(e0)

@time e2 = FermiCG.compute_pt2_energy(v0, cluster_ops, cluster_ham, thresh_foi=1e-8)
@save "eq_tpsci_results.jld2" e0 e2 v0
    
println()
println("	*======TPSCI results======*")
@printf("TCI Thresh: %8.6f  Dim:%8d\n",1e-3,size(v0)[1])
println()
@printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
for r in 1:nroots
    @printf("TCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2a[r] + ecore)
end

#global ci_vector = v0

#ci_vector1 = FermiCG.ClusteredState(clusters, ref_fock, R=1)

#vtmp = FermiCG.get_vector(clusters, root=4)

#print(v0)

for r in 1:nroots
    @printf("TPSCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2a[r] + ecore)
    display(v0,thresh=1e-4,root=r)
end


