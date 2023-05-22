using FermiCG
using JLD2
using PyCall
#using Plots
using LinearAlgebra
using Printf
using QCBase
using RDM
using ClusterMeanField


pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("/home/nbraunsc/FermiCG-data/excited_paper/p5/benz-6-fcidump");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, size(h,2))
ints = InCoreInts(ecore,h,g);

rdm1 = zeros(size(ints.h1))


na = 18
nb = 18

#clusters_in    = [(1:6),(7:12),(13:18),(19:24),(25:30),(31:36)]
#init_fspace = [(3,3),(3,3),(3,3),(3,3),(3,3),(3,3)]
n_clusters = 6

# define clusters
#cluster_list = [collect(1:6), collect(7:12), collect(13:18), collect(19:24), collect(25:30), collect(31:36)]
#clusters = [MOCluster(i,collect(cluster_list[i])) for i = 1:length(cluster_list)]
#init_fspace = [ (3,3) for i in 1:n_clusters]
#display(clusters)

##run cmf_oo
#e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, RDM1(rdm1, rdm1), verbose=0, diis_start=3);
#
#ints = orbital_rotation(ints,U_cmf)
#@save  "/home/nbraunsc/FermiCG-data/benz6/cmf_diis.jld2" ints d1 clusters init_fspace
@load  "/home/nbraunsc/FermiCG-data/excited_paper/p5/cmf_diis.jld2"

M = 100

ref_fock = FockConfig(init_fspace)

# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [4,4,4,4,4,4], ref_fock, max_roots=M, verbose=1);
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3,3,3], ref_fock, max_roots=M, verbose=1);
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [6,6,6,6,6,6], ref_fock, max_roots=M, verbose=1);
#
# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
#
# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
#
# Add cmf hamiltonians for doing MP-style PT2 
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b, verbose=0);

nroots = 25
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);

ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,1,1])] = zeros(Float64,nroots)

#Triplets
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([2,1,1,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,2,1,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,2,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,2,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,2,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,1,2])] = zeros(Float64,nroots)

#Singlets
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([3,1,1,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,3,1,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,3,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,3,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,3,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,1,3])] = zeros(Float64,nroots)                       

#2nd Triplets
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([4,1,1,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,4,1,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,4,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,4,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,4,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,1,4])] = zeros(Float64,nroots)

#2nd Singlets
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([5,1,1,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,5,1,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,5,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,5,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,5,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,1,5])] = zeros(Float64,nroots)

#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([6,1,1,1,1,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,6,1,1,1,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,6,1,1,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,6,1,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,6,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,1,6])] = zeros(Float64,nroots)
#
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([7,1,1,1,1,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,7,1,1,1,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,7,1,1,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,7,1,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,7,1])] = zeros(Float64,nroots)
#ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,1,7])] = zeros(Float64,nroots)

#display(ci_vector)
FermiCG.eye!(ci_vector)


#thresh_list = [0.005,0.002,0.001,0.0007,0.0005,0.0003,0.0001]
#thresh = [0.01, 0.005, 0.001, 0.0007, 0.0005, 0.0001]
#thresh = [0.001, 0.0007, 0.0005, 0.0002, 0.0001]
thresh = [0.0005]

for i in thresh
    e0, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham,
                              thresh_asci =-1,     # Threshold of P-space configs to search from
                              thresh_foi  =1e-5,    # Threshold for keeping terms when defining FOIS
                              thresh_cipsi=i, # Threshold for adding to P-space
                              max_iter=10, 
                              max_mem_ci=100.0);

    @time e2 = FermiCG.compute_pt2_energy(v0, cluster_ops, clustered_ham, thresh_foi=1e-8);
    name = "25roots_thresh_"*string(i)*".jld2"
    
    println()
    println("	*======TPSCI results======*")
    @printf("TCI Thresh: %8.6f  Dim:%8d\n",i,size(v0)[1])
    println()
    @printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
    for r in 1:nroots
        @printf("TCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2[r] + ecore)
    end
    
    clustered_S2 = FermiCG.extract_S2(ci_vector.clusters)
    
    println()
    println("	*======TPSCI S2 results======*")
    @printf(" %-50s", "Compute FINAL S2 expectation values: ")
    @time s2 = FermiCG.compute_expectation_value_parallel(v0, cluster_ops, clustered_S2)
    
    @printf(" %5s %12s %12s\n", "Root", "Energy", "S2") 
    for r in 1:nroots
        @printf(" %5s %12.8f %12.8f\n",r, e0[r]+ecore, abs(s2[r]))
    end
    @save "/home/nbraunsc/FermiCG-data/excited_paper/p5/"*string(name) clusters d1 ints cluster_bases ci_vector e0 v0 e2 ecore s2
end
