using Pkg; Pkg.activate("../../../FermiCG/")
using FermiCG
using PyCall
using Plots
using LinearAlgebra
using Printf
using JLD2

molecule = "
He         0.83846       -0.83846        0.83846
He         0.83846       -0.83846       -0.83846
He         0.83846        0.83846       -0.83846
He         0.83846        0.83846        0.83846
He        -0.83846        0.83846        0.83846
He        -0.83846        0.83846       -0.83846
He        -0.83846       -0.83846       -0.83846
He        -0.83846       -0.83846        0.83846
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

#FermiCG.pyscf_write_molden(mol, lo_ao, filename="he08_02/lowdin_ao_augccpvdz.molden")

#write fci dump file from the modified mo coefficients
tools = pyimport("pyscf.tools")
tools.fcidump.from_mo(pymol, "fcidump.he08_02", lo_ao)

#Can just read in pyscf dump file for integrals (once you have already run an scf calculation)
pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("fcidump.he08_02");
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



