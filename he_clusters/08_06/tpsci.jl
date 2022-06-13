using FermiCG
using PyCall
using Plots
using LinearAlgebra
using Printf
using JLD2

molecule = "
He         0.00000        0.00000       -5.93406
He         0.00000        0.00000       -4.23728
He         0.00000        0.00000       -2.54161
He         0.00000        0.00000       -0.84720
He         0.00000        0.00000        0.84720
He         0.00000        0.00000        2.54161
He         0.00000        0.00000        4.23728
He         0.00000        0.00000        5.93406
"
atoms = []
for (li,line) in enumerate(split(rstrip(lstrip(molecule)), "\n"))
    l = split(line)
    push!(atoms, Atom(li, l[1], parse.(Float64,l[2:4])))
end

#basis = "aug-cc-pvdz" #9 orbs on each He
basis = "cc-pvdz" #5 orbs on each He
#basis = "cc-pvtz" # 6 orbs on each He
#basis = "sto-3g"

# Create FermiCG.Molecule type
mol = Molecule(0,1,atoms,basis)

pyscf = pyimport("pyscf")
pymol = pyscf.gto.Mole(atom=molecule, spin=0, charge=0, basis=basis)
pymol.build()
mf = pyscf.scf.RHF(pymol).run()
s = mf.get_ovlp(pymol)

lo = pyimport("pyscf.lo.orth")
lo_ao = lo.lowdin(s)
println("size of Lowdin ortho AO's:", size(lo_ao))

FermiCG.pyscf_write_molden(mol, lo_ao, filename="lowdin_ao_ccpvdz.molden")

#write fci dump file from the modified mo coefficients
tools = pyimport("pyscf.tools")
tools.fcidump.from_mo(pymol, "fcidump.he08_06", lo_ao)

#Can just read in pyscf dump file for integrals (once you have already run an scf calculation)
pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("fcidump.he08_06");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, size(h,2))
#This one below was not working. Error: setfield! immutable struct of type InCoreInts cannot be changed
ints = InCoreInts(ecore,h,g);


#Run cmf
#Define clusters and intial Fock space for inital CMF calc for 5 orbs each He
clusters_in = [(1:5),(6:10), (11:15), (16:20), (21:25), (26:30), (31:35), (36:40)]
init_fspace = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
rdm1 = zeros(size(ints.h1))
na=8
nb=8

#Define clusters now using FermiCG code
clusters = [Cluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]
display(clusters)

@save "before_cmf.jld2" ints clusters init_fspace

print(size(ints.h1))
rdm1 = zeros(size(ints.h1))

#do a CMF calculation to optimize cluster orbitals

e_cmf, U, Da, Db = FermiCG.cmf_oo(ints, clusters, init_fspace, rdm1, rdm1, max_iter_oo=100, verbose=0, gconv=1e-6, method="bfgs");

#rotate the integrals by the cmf calculation
ints = FermiCG.orbital_rotation(ints, U);
max_roots = 100

#Build Cluster Basis (delta n is here)
cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, verbose=1, max_roots=max_roots, init_fspace=init_fspace, rdm1a=Da, rdm1b=Db);

@save "after_cmf.jld2" ints Da Db e_cmf cluster_bases clusters init_fspace

FermiCG.pyscf_write_molden(mol,lo_ao*U, filename="cmf_he08_06.molden");

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
ci_vector[ref_fock][ClusterConfig([1,1,1,1,1,1,1,1])] = [1,0,0]
ci_vector[ref_fock][ClusterConfig([2,1,1,1,1,1,1,1])] = [0,1,0]
ci_vector[ref_fock][ClusterConfig([1,2,1,1,1,1,1,1])] = [0,0,1]

display(ci_vector)

e0, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, cluster_ham,
                          thresh_cipsi=1e-3, # Threshold for adding to P-space
                          thresh_foi=1e-7,    # Threshold for keeping terms when defining FOIS
                          thresh_asci=0.001,     # Threshold of P-space configs to search from
                          max_iter=10,
                          matvec=3);

@time e2 = FermiCG.compute_pt2_energy(v0, cluster_ops, cluster_ham, thresh_foi=1e-8)
