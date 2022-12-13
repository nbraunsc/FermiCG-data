using FermiCG
using JLD2
using PyCall
using LinearAlgebra
using Printf
using QCBase
using RDM
using ClusterMeanField

molecule = "
C       -1.9966725193      0.4342546555      0.0000000000
C        3.0351564726      0.9167598283      0.0000000000
C       -0.3000353437      5.8124300241      1.5796949175
C       -1.3746008454      3.9325988906     -2.4741342248
"
atoms = []
for (li,line) in enumerate(split(rstrip(lstrip(molecule)), "\n"))
    l = split(line)
    push!(atoms, Atom(li, l[1], parse.(Float64,l[2:4])))
end

basis = "cc-pvdz"

# Create FermiCG.Molecule type
mol = Molecule(0,1,atoms,basis)

pyscf = pyimport("pyscf")
pymol = pyscf.gto.Mole(atom=molecule, spin=0, charge=0, basis=basis)
pymol.build()
println(pymol.nelec)
mf = pyscf.scf.RHF(pymol).run()

s = mf.get_ovlp(pymol)

lo = pyimport("pyscf.lo.orth")
lo_ao = lo.lowdin(s)
println("size of Lowdin ortho AO's:", size(lo_ao))

FermiCG.pyscf_write_molden(mol, lo_ao, filename="four_lo.molden")

#write fci dump file from the modified mo coefficients
tools = pyimport("pyscf.tools")
tools.fcidump.from_mo(pymol, "fcidump.four", lo_ao)

#Can just read in pyscf dump file for integrals (once you have already run an scf calculation)
pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("fcidump.four");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, size(h,2))
#This one below was not working. Error: setfield! immutable struct of type InCoreInts cannot be changed
ints = InCoreInts(ecore,h,g);

rdm1 = zeros(size(ints.h1))

na = 12
nb = 12

#clusters_in    = [(1:6),(7:12),(13:18),(19:24)]
n_clusters = 4

# define clusters
cluster_list = [collect(1:14), collect(15:28), collect(29:42), collect(43:56)]
clusters = [MOCluster(i,collect(cluster_list[i])) for i = 1:length(cluster_list)]
init_fspace = [ (3,3) for i in 1:n_clusters]
display(clusters)

#run cmf_oo
e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, RDM1(rdm1, rdm1), verbose=0, diis_start=3);

ints = orbital_rotation(ints,U_cmf)
ecore = ints.h0

M = 150
ref_fock = FockConfig(init_fspace)

# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [1,1,1,1], ref_fock, max_roots=M, verbose=1);
#
# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
display(clustered_ham)

@save "four_ham.jld2" clusters clustered_ham
