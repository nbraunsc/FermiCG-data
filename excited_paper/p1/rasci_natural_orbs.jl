using QCBase
using ClusterMeanField
using ActiveSpaceSolvers
using RDM
using InCoreIntegrals
using PyCall
using JLD2
using NPZ
using Printf
using LinearAlgebra


function ras_cluster_rdms(ints::InCoreInts{T}, clusters::Vector{MOCluster}, rdm1::RDM1{T}, ansatze; verbose = 0) where T
    
    #cluster_list = [collect(1:6), collect(7:12), collect(13:18), collect(19:24)]
    #cluster_list = [collect(1:9), collect(10:18)]
    no_transform = []
    for i in 1:length(clusters)
        ci = clusters[i]
        verbose == 0 || display(ci)
        #ints_i = subset(ints, ci, rdm1) 
        ints_i = subset(ints, cluster_list[i]) 
        display(ints_i.h1)

        # Do sparse build 
        sol = solve(ints_i, ansatze[i], SolverSettings(nroots=1))
        rdm = ActiveSpaceSolvers.compute_1rdm(sol)
        rdm = rdm[1]+rdm[2]
        e, v = eigen(rdm)
        println(e)
        push!(no_transform, v)
    end
    return no_transform
end





        
