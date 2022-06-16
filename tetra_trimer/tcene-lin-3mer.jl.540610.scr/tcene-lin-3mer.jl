using FermiCG
using PyCall
using Plots
using LinearAlgebra
using Printf


pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("/home/nbraunsc/FermiCG-data/tetra_trimer/fcidump_tcene_lin_3mer");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, size(h,2))
ints = InCoreInts(ecore,h,g)

rdm1 = zeros(size(ints.h1))


na = 9
nb = 9

clusters_in    = [(1:6),(7:12),(13:18)]
init_fspace = [(3,3),(3,3),(3,3)]

init_fspace_tts1 = [(4,2),(2,4),(3,3)]
init_fspace_tts2 = [(2,4),(4,2),(3,3)]
init_fspace_tst1 = [(4,2),(3,3),(2,4)]
init_fspace_tst2 = [(2,4),(3,3),(4,2)]
init_fspace_stt1 = [(3,3),(4,2),(2,4)]
init_fspace_stt2 = [(3,3),(2,4),(4,2)]



# define clusters
clusters = [Cluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]
display(clusters)

#FermiCG.cmf_ci(ints, clusters, init_fspace, rdm1, rdm1, verbose=2,sequential=false)

e_cmf, U, Da, Db  = FermiCG.cmf_oo(ints, clusters, init_fspace, rdm1,rdm1,
                                        max_iter_oo=100, verbose=0, gconv=1e-7, method="gd",sequential=false);
ints = FermiCG.orbital_rotation(ints,U)

max_roots = 100
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, verbose=0, max_roots=max_roots,
        init_fspace=init_fspace, rdm1a=Da, rdm1b=Db);

#cluster_bases = FermiCG.compute_cluster_est_basis(ints, clusters, Da, Db, thresh_schmidt=5e-5, init_fspace=init_fspace)

#
# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

# Add cmf hamiltonians for doing MP-style PT2 
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, Da, Db, verbose=0);



nroots = 19

ref_fock = FermiCG.FockConfig(init_fspace)
TTS1_fock = FermiCG.FockConfig(init_fspace_tts1)
TTS2_fock = FermiCG.FockConfig(init_fspace_tts2)
TST1_fock = FermiCG.FockConfig(init_fspace_tst1)
TST2_fock = FermiCG.FockConfig(init_fspace_tst2)
STT1_fock = FermiCG.FockConfig(init_fspace_stt1)
STT2_fock = FermiCG.FockConfig(init_fspace_stt2)

ci_vector = FermiCG.TPSCIstate(clusters, ref_fock, R=nroots)
ci_vector[ref_fock][ClusterConfig([2,1,1])] = [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,2,1])] = [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,1,2])] = [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

ci_vector[ref_fock][ClusterConfig([2,2,1])] = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([2,1,2])] = [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,2,2])] = [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]

ci_vector[ref_fock][ClusterConfig([3,1,1])] = [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,3,1])] = [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,1,3])] = [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]

ci_vector[ref_fock][ClusterConfig([4,1,1])] = [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,4,1])] = [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]
ci_vector[ref_fock][ClusterConfig([1,1,4])] = [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0]

FermiCG.add_fockconfig!(ci_vector,TTS1_fock)
FermiCG.add_fockconfig!(ci_vector,TTS2_fock)
FermiCG.add_fockconfig!(ci_vector,TST1_fock)
FermiCG.add_fockconfig!(ci_vector,TST2_fock)
FermiCG.add_fockconfig!(ci_vector,STT1_fock)
FermiCG.add_fockconfig!(ci_vector,STT2_fock)

ci_vector[TTS1_fock][ClusterConfig([1,1,1])] = [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0]
ci_vector[TTS2_fock][ClusterConfig([1,1,1])] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0]
ci_vector[TST1_fock][ClusterConfig([1,1,1])] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0]
ci_vector[TST2_fock][ClusterConfig([1,1,1])] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0]
ci_vector[STT1_fock][ClusterConfig([1,1,1])] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0]
ci_vector[STT2_fock][ClusterConfig([1,1,1])] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]


e0, vec_out = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham)
#e0, vec_out = FermiCG.tps_ci_davidson(ci_vector, cluster_ops, clustered_ham)

display(e0)
display(ci_vector)
ci_vector = vec_out
display(ci_vector)



display(ci_vector,thresh=0.1)

thresh_list = [0.005,0.002,0.001,0.0007,0.0005,0.0003,0.0001,0.00005]
thresh_list = [0.0005,0.0003,0.0001,0.00005]


for thresh_cipsi in thresh_list
    e0, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham,
    			    thresh_cipsi=thresh_cipsi, # Threshold for adding to P-space
    			    thresh_foi=1e-6,    # Threshold for keeping terms when defining FOIS    
            		    ci_max_ss_vecs = 2,
    			    thresh_asci=0.001,     # Threshold of P-space configs to search from
    			    max_iter=10,
    			    matvec=3);

    e2a = FermiCG.compute_pt2_energy(v0, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=0)

    println()
    println("	*======TPSCI results======*")
    @printf("TCI Thresh: %8.6f  Dim:%8d\n",thresh_cipsi,size(v0)[1])
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
end

