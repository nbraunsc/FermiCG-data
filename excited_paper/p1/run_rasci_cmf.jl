using QCBase
using ClusterMeanField
using ActiveSpaceSolvers
using RDM
using InCoreIntegrals
using PyCall
using JLD2
using NPZ

function run_cmf(d_guess=nothing)
   
    ## Integrals with sym_orth and pseudo canoncalization
    #h0 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h0.npy")
    #h1 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h1.npy")
    #h2 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h2.npy")
    
    ## Integrals with sym_orth but no pseudo canoncalization
    #h0 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h0_sym_nopc.npy")
    #h1 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h1_sym_nopc.npy")
    #h2 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h2_sym_nopc.npy")
    
    ## Integrals without sym_orth and without pseudo canoncalization
    #h0 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h0_no_pc.npy")
    #h1 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h1_no_pc.npy")
    #h2 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h2_no_pc.npy")
    
    ## Integrals without sym_orth and but do have pseudo canoncalization
    #h0 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h0_no_ortho.npy")
    #h1 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h1_no_ortho.npy")
    #h2 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h2_no_ortho.npy")
    
    ## Integrals from vibin, fermicluster 2pz and Boys localization (used these in the paper)
    #pyscf = pyimport("pyscf");
    #fcidump = pyimport("pyscf.tools.fcidump");
    #ctx = fcidump.read("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/fcidump_4mer");
    #h = ctx["H1"];
    #g = ctx["H2"];
    #ecore = ctx["ECORE"];
    #g = pyscf.ao2mo.restore("1", g, size(h,2))
    #ints = InCoreInts(ecore,h,g)
    
    ## RHF Natural Orbital integrals
    #h0 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h0_no.npy")
    #h1 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h1_no.npy")
    #h2 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h2_no.npy")
   
    # 2pz natural orbs
    #h0 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h0_ras_no_old.npy")
    #h1 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h1_ras_no_old.npy")
    #h2 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h2_ras_no_old.npy")
    
    #h0 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h0_ras_no_sympc.npy")
    #h1 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h1_ras_no_sympc.npy")
    #h2 = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_h2_ras_no_sympc.npy")
    #ints = InCoreInts(h0,h1,h2)

    clusters_in    = [(1:6),(7:12),(13:18),(19:24)]
    n_clusters = 4
    cluster_list = [collect(1:6), collect(7:12), collect(13:18), collect(19:24)]
    clusters = [MOCluster(i,collect(cluster_list[i])) for i = 1:length(cluster_list)]
    init_fspace = [ (3,3) for i in 1:n_clusters]
    display(clusters)
    #init_cluster_ansatz = [FCIAnsatz(6, 3, 3), FCIAnsatz(6, 3, 3), FCIAnsatz(6, 3, 3), FCIAnsatz(6, 3, 3)]
    #init_cluster_ansatz = [RASCIAnsatz(6, 3, 3, (2,2,2), max_h=2, max_p=2), RASCIAnsatz(6,3,3,(2,2,2), max_h=2, max_p=2), RASCIAnsatz(6,3,3,(2,2,2), max_h=2, max_p=2), RASCIAnsatz(6,3,3,(2,2,2), max_h=2, max_p=2)]
    init_cluster_ansatz = [RASCIAnsatz(6, 3, 3, (2,2,2), max_h=1, max_p=1), RASCIAnsatz(6,3,3,(2,2,2), max_h=1, max_p=1), RASCIAnsatz(6,3,3,(2,2,2), max_h=1, max_p=1), RASCIAnsatz(6,3,3,(2,2,2), max_h=1, max_p=1)]
    #init_cluster_ansatz = [RASCIAnsatz(6, 3, 3, (2,2,2), max_h=4, max_p=4), RASCIAnsatz(6,3,3,(2,2,2), max_h=4, max_p=4), RASCIAnsatz(6,3,3,(2,2,2), max_h=4, max_p=4), RASCIAnsatz(6,3,3,(2,2,2), max_h=4, max_p=4)]
    display(init_cluster_ansatz)
    
    ### FCI-CMF converged integrals
    @load "/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_fci.jld2" 
    ints = ints_fci
    #Cfrags = npzread("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/Cfrags.npy")
    #rdm1 = zeros(size(ints.h1))
    #
    #include("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/rasci_natural_orbs.jl")
    #print(typeof(init_cluster_ansatz))

    #rdms = ras_cluster_rdms(ints, clusters, RDM1(rdm1, rdm1), init_cluster_ansatz)

    #for i in 1:n_clusters
    #    Cfrags[i,:,:] .= Cfrags[i,:,:]*rdms[i]
    #end

    #npzwrite("/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/C_nat_orbs_sym_pc.npy", Cfrags)
    #error("hi")

    #CMF has default of delta_elec=0 for all clusters
    #delta_elec = [1,1,1,1]
    #ansatze = ActiveSpaceSolvers.generate_cluster_fock_ansatze(init_fspace, clusters, init_cluster_ansatz, delta_elec)
    ansatze = ActiveSpaceSolvers.generate_cluster_fock_ansatze(init_fspace, clusters, init_cluster_ansatz)
    display(ansatze)
    
    if d_guess == nothing
        rdm1 = zeros(size(ints.h1))
        #run cmf_oo
        #e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, ansatze, RDM1(rdm1, rdm1), maxiter_oo = 300, verbose=0, zero_intra_rots = true, diis_start=1);
        #e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, ansatze, RDM1(rdm1, rdm1), maxiter_oo = 1000, verbose=0, zero_intra_rots = false, diis_start=1);
        e_cmf, U_cmf, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, ansatze, RDM1(rdm1, rdm1), max_iter_oo=200, verbose=0, gconv=1e-6, method="cg");
        ints_cg = orbital_rotation(ints,U_cmf)
        @save "/Users/nicole/My Drive/code/FermiCG-data/excited_paper/p1/ints_cg.jld2" ints_bfgs U_cmf
        return e_cmf, U_cmf, d1
    else
        println("here")
        e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, ansatze, RDM1(d_guess.a, d_guess.b), maxiter_oo = 300, verbose=0, diis_start=3);
        return e_cmf, U_cmf, d1
    end

    #e_cmf, U_cmf, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, ansatze, RDM1(rdm1, rdm1), max_iter_oo=200, verbose=0, gconv=1e-6, method="bfgs");
    #cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3], init_fspace, max_roots=20, verbose=1);
    #@save "testing_ansatz.jld2" ints clusters d1 init_fspace rdm1 
    #return d1
end

