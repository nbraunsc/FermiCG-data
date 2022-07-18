using FermiCG
using PyCall
using Plots
using LinearAlgebra
using Printf
using JLD2

molecule = "
He         0.00000        0.00000        0.00000
He         4.24264        0.00000        0.00000
He         0.00000        4.24264        0.00000
He         4.24264        4.24264        0.00000
He         2.12132        2.12132        3.00000
He         2.12132        2.12132       -3.00000
"

atoms = []
for (li,line) in enumerate(split(rstrip(lstrip(molecule)), "\n"))
    l = split(line)
    push!(atoms, Atom(li, l[1], parse.(Float64,l[2:4])))
end

basis = "cc-pvdz"

# Create FermiCG.Molecule type
mol     = Molecule(0, 1, atoms,basis);

pyscf = pyimport("pyscf")

n_steps = 40
step_size = .05

pymol_init = pyscf.gto.Mole(atom=molecule,
                            symmetry = false, spin =0,charge=0,
                            basis = basis)
pymol_init.build()

io = open("traj.xyz", "w");
energies_ground = []
energies_t1 = []
energies_t2 = []
energies_t3 = []
energies_t4 = []
energies_t5 = []
energies_t6 = []

energies_s1 = []
energies_s2 = []
energies_s3 = []
energies_s4 = []
energies_s5 = []
energies_s6 = []

energies_t7 = []
energies_t8 = []
energies_t9 = []
energies_t10 = []
energies_t11 = []
energies_t12 = []

lo = pyimport("pyscf.lo.orth")
tools = pyimport("pyscf.tools")
fcidump = pyimport("pyscf.tools.fcidump");

for R in 1:n_steps

    pymol = deepcopy(pymol_init)
    scale = 1+R*step_size

    xyz = @sprintf("%5i\n\n", length(mol.atoms))
    tmp = []
    for a in mol.atoms
        push!(tmp, ["He", (a.xyz[1]*scale, a.xyz[2]*scale, a.xyz[3]*scale)])
        xyz = xyz * @sprintf("%6s %24.16f %24.16f %24.16f \n", a.symbol, a.xyz[1]*scale, a.xyz[2]*scale, a.xyz[3]*scale)
    end
    pymol.atom = tmp
    pymol.build()


    println(xyz)
    write(io, xyz);


    #println(pymol.format_atom(1))

    #     mol_R = Molecule(0, 1, [a[0]pymol.atom, pymol.basis)


    mf = pyscf.scf.RHF(pymol).run()
    
    s = mf.get_ovlp(pymol)

    lo_ao = lo.lowdin(s)
    println("size of Lowdin ortho AO's:", size(lo_ao))

    #write fci dump file from the modified mo coefficients
    tools.fcidump.from_mo(pymol, "fcidump.he06_oct", lo_ao)

    #Can just read in pyscf dump file for integrals (once you have already run an scf calculation)
    ctx = fcidump.read("fcidump.he06_oct");
    h = ctx["H1"];
    g = ctx["H2"];
    ecore = ctx["ECORE"];
    g = pyscf.ao2mo.restore("1", g, size(h,2))

    #This one below was not working. Error: setfield! immutable struct of type InCoreInts cannot be changed
    ints = InCoreInts(ecore,h,g);

    #Run cmf
    #Define clusters and intial Fock space for inital CMF calc for 14 orbs each He
    #clusters_in = [(1:14),(15:28), (29:42), (43:56), (57:70), (71:84), (85:98)]

    #Define clusters and intial Fock space for inital CMF calc for 9 orbs each He
    #clusters_in = [(1:9),(10:18), (19:27), (28:36), (37:45), (46:54), (55:63)]

    #Define clusters and intial Fock space for inital CMF calc for 5 orbs each He
    clusters_in = [(1:5),(6:10), (11:15), (16:20), (21:25), (26:30)]
    init_fspace = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
    rdm1 = zeros(size(ints.h1))
    na=6
    nb=6

    #Define clusters now using FermiCG code
    clusters = [Cluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]
    display(clusters)

    #@save "before_cmf.jld2" ints clusters init_fspace

    print(size(ints.h1))
    rdm1 = zeros(size(ints.h1))

    #do a CMF calculation to optimize cluster orbitals

    e_cmf, U, Da, Db = FermiCG.cmf_oo(ints, clusters, init_fspace, rdm1, rdm1, max_iter_oo=200, verbose=0, gconv=1e-6, method="bfgs");


    #rotate the integrals by the cmf calculation
    ints = FermiCG.orbital_rotation(ints, U);
    max_roots = 100

    #Build Cluster Basis (delta n is here)
    cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, verbose=1, max_roots=max_roots, init_fspace=init_fspace, rdm1a=Da, rdm1b=Db);

    #@save "scan_after_cmf.jld2" ints Da Db e_cmf cluster_bases clusters init_fspace

    #Build Clustered Operator
    cluster_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
    
    #Build Cluster Operators
    cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
    FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, Da, Db);

    #Need to find reference state 
    ref_fock = FermiCG.FockConfig(init_fspace)
    nroots = 19
    #ci_vector = FermiCG.TPSCIstate(clusters, ref_fock, R=nroots)
    ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
    #ci_vector = FermiCG.ClusteredState(clusters, ref_fock, R=nroots);
    #Need to find the automated way to define these other excited configs away from ref state, example is to large
    #to do by hand
    #probably something to do with building p spaces and q spaces
    ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,1])] = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([2,1,1,1,1,1])] = [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,2,1,1,1,1])] = [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,1,2,1,1,1])] = [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,1,1,2,1,1])] = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,1,1,1,2,1])] = [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,2])] = [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]
    
    ci_vector[ref_fock][ClusterConfig([3,1,1,1,1,1])] = [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,3,1,1,1,1])] = [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,1,3,1,1,1])] = [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,1,1,3,1,1])] = [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,1,1,1,3,1])] = [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,3])] = [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0]
    
    ci_vector[ref_fock][ClusterConfig([4,1,1,1,1,1])] = [0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,4,1,1,1,1])] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,1,4,1,1,1])] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0]
    ci_vector[ref_fock][ClusterConfig([1,1,1,4,1,1])] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0]
    ci_vector[ref_fock][ClusterConfig([1,1,1,1,4,1])] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0]
    ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,4])] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    display(ci_vector)


    #thresh_list = [0.01, 0.001, 0.003, 0.005, 0.0001]

    #for thresh_cipsi in thresh_list
    e0, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, cluster_ham,
                              thresh_cipsi=1e-3, # Threshold for adding to P-space
                              #thresh_cipsi=thresh_cipsi, # Threshold for adding to P-space
                              thresh_foi=1e-3,    # Threshold for keeping terms when defining FOIS
                              thresh_asci=0,     # Threshold of P-space configs to search from
                              max_iter=10);

    @time e2 = FermiCG.compute_pt2_energy(v0, cluster_ops, cluster_ham, thresh_foi=1e-8)
    #name = "eq_tpsci_results"*string(0.001)*".jld2"
    #@save name e0 e2 v0 ecore
    
    println()
    println("	*======TPSCI results======*")
    @printf("TCI Thresh: %8.6f  Dim:%8d\n",1e-2,size(v0)[1])
    println()
    @printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
    for r in 1:nroots
        @printf("TCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2[r] + ecore)
    end
    #end

    push!(energies_ground, e0[1]+e2[1]+ecore)
    push!(energies_t1, e0[2]+e2[2]+ecore)
    push!(energies_t2, e0[3]+e2[3]+ecore)
    push!(energies_t3, e0[4]+e2[4]+ecore)
    push!(energies_t4, e0[5]+e2[5]+ecore)
    push!(energies_t5, e0[6]+e2[6]+ecore)
    push!(energies_t6, e0[7]+e2[7]+ecore)
    
    push!(energies_s1, e0[8]+e2[8]+ecore)
    push!(energies_s2, e0[9]+e2[9]+ecore)
    push!(energies_s3, e0[10]+e2[10]+ecore)
    push!(energies_s4, e0[11]+e2[11]+ecore)
    push!(energies_s5, e0[12]+e2[12]+ecore)
    push!(energies_s6, e0[13]+e2[13]+ecore)    
    
    push!(energies_t7, e0[14]+e2[14]+ecore)
    push!(energies_t8, e0[15]+e2[15]+ecore)
    push!(energies_t9, e0[16]+e2[16]+ecore)
    push!(energies_t10, e0[17]+e2[17]+ecore)
    push!(energies_t11, e0[18]+e2[18]+ecore)
    push!(energies_t12, e0[19]+e2[19]+ecore)
end

@save "scan_energies_triplet.jld2" energies_ground energies_t1 energies_t2 energies_t3 energies_t4 energies_t5 energies_t6
@save "scan_energies_singlet.jld2" energies_s1 energies_s2 energies_s3 energies_s4 energies_s5 energies_s6
@save "scan_energies_second_triple.jld2" energies_t7 energies_t8 energies_t9 energies_t10 energies_t11 energies_t12
