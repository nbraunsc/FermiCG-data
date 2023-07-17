using FermiCG
using PyCall
using LinearAlgebra
using Printf
using JLD2

molecule = "
Ne 0.00000000 0.00000000 0.00000000
Ne 3.41421356 0.00000000 0.0000000
"
atoms = []
for (li,line) in enumerate(split(rstrip(lstrip(molecule)), "\n"))
    l = split(line)
    push!(atoms, Atom(li, l[1], parse.(Float64,l[2:4])))
end

#basis = "aug-cc-pvdz" #9 orbs on each He

#basis = "cc-pvdz" #5 orbs on each He
#basis = "cc-pvtz" # 6 orbs on each He
basis = "6-31g"
#basis = "sto-3g"

# Create FermiCG.Molecule type
mol = Molecule(0,1,atoms,basis)

pyscf = pyimport("pyscf")
pymol = pyscf.gto.Mole(atom=molecule, spin=0, charge=0, basis=basis)
pymol.build()
mf = pyscf.scf.RHF(pymol).run()
s = mf.get_ovlp(pymol)

println("nelec:", pymol.nelec)

lo = pyimport("pyscf.lo.orth")
lo_ao = lo.lowdin(s)
println("size of Lowdin ortho AO's:", size(lo_ao))

#write fci dump file from the modified mo coefficients
tools = pyimport("pyscf.tools")
tools.molden.from_mo(pymol, "lo_ao.molden", lo_ao)
tools.fcidump.from_mo(pymol, "fcidump.neon", lo_ao)


@load "/Users/nicole/My Drive/code/FermiCG-data/neon/rdms.jld2"

frag1 = lo_ao[:,1:9]
frag2 = lo_ao[:, 10:18]

frag1_no = frag1*rdms[1]
frag2_no = frag2*rdms[2]

C_no = zeros(size(lo_ao))
C_no[:,1:9] .=frag1_no
C_no[:,10:18] .=frag2_no
println(size(C_no))

tools.molden.from_mo(pymol, "/Users/nicole/My Drive/code/FermiCG-data/neon/C_no_ras.molden", C_no)
tools.fcidump.from_mo(pymol, "/Users/nicole/My Drive/code/FermiCG-data/neon/fcidump_nat.neon", C_no)




