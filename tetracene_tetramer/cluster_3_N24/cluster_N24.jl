using Pkg; Pkg.activate("../../../FermiCG/")
using FermiCG
using JLD2
using NPZ

#load integrals from disk
ints = InCoreInts(
    npzread("integrals_h0.npy"), 
    npzread("integrals_h1.npy"), 
    npzread("integrals_h2.npy") 
);
C = npzread("mo_coeffs_act.npy")
S = npzread("../overlap_mat.npy");
D = npzread("../density_mat.npy");

@load "cmf_data.jld2"


max_roots = 100

#
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, verbose=0, max_roots=max_roots,
        init_fspace=init_fspace, rdm1a=Da, rdm1b=Db, delta_elec=4);
#
# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

#
# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

#
# Add cmf hamiltonians for doing MP-style PT2 
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, Da, Db, verbose=0);

nroots = 31
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
# Add the lowest energy single exciton to basis
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([2,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,2,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,2,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,2])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([3,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,3,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,3,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,3])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([4,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,4,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,4,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,4])] = zeros(Float64,nroots)

# TT states ms=0
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([2,2,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([2,1,2,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([2,1,1,2])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,2,2,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,2,1,2])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,2,2])] = zeros(Float64,nroots)

# Spin-flip states
fspace_0 = FermiCG.FockConfig(init_fspace)

## ba
tmp_fspace = FermiCG.replace(fspace_0, (1,2), ([2,4],[4,2]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
tmp_fspace = FermiCG.replace(fspace_0, (1,3), ([2,4],[4,2]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
tmp_fspace = FermiCG.replace(fspace_0, (1,4), ([2,4],[4,2]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
tmp_fspace = FermiCG.replace(fspace_0, (2,3), ([2,4],[4,2]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
tmp_fspace = FermiCG.replace(fspace_0, (2,4), ([2,4],[4,2]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
tmp_fspace = FermiCG.replace(fspace_0, (3,4), ([2,4],[4,2]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)

## ab
tmp_fspace = FermiCG.replace(fspace_0, (1,2), ([4,2],[2,4]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
tmp_fspace = FermiCG.replace(fspace_0, (1,3), ([4,2],[2,4]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
tmp_fspace = FermiCG.replace(fspace_0, (1,4), ([4,2],[2,4]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
tmp_fspace = FermiCG.replace(fspace_0, (2,3), ([4,2],[2,4]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
tmp_fspace = FermiCG.replace(fspace_0, (2,4), ([4,2],[2,4]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
tmp_fspace = FermiCG.replace(fspace_0, (3,4), ([4,2],[2,4]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)

FermiCG.eye!(ci_vector)

e0, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham,
                            thresh_asci =1e-2,     # Threshold of P-space configs to search from
                            thresh_foi  =1e-5,    # Threshold for keeping terms when defining FOIS
                            thresh_cipsi=1e-4, # Threshold for adding to P-space
                            max_iter=10);
@time e2 = FermiCG.compute_pt2_energy(v0, cluster_ops, clustered_ham, thresh_foi=1e-8);

@save "M100.jld2" clusters Da Db ints C cluster_bases ci_vector e0 v0 e2
