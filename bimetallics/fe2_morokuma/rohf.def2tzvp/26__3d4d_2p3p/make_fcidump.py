import pyscf
import pyscf.tools
from pyscf.tools import fcidump

import numpy as np

molecule = """
 Fe 1.67785607 0.00052233 0.06475932
 O 0.00000000 0.00000000 -0.47099074
 Fe -1.67785607 -0.00052233 0.06475932
 Cl 1.87002704 -1.09796437 1.99091682
 Cl 2.93244917 -0.98210488 -1.47467288
 Cl 2.37160936 2.07954091 -0.50446591
 Cl -1.87002704 1.09796437 1.99091682
 Cl -2.93244917 0.98210488 -1.47467288
 Cl -2.37160936 -2.07954091 -0.50446591
 """

basis = "def2-tzvp"
pymol = pyscf.gto.Mole(
        atom    =   molecule,
        symmetry=   True,
<<<<<<< HEAD
        spin    =   10, # number of unpaired electrons
=======
        spin    =   0, # number of unpaired electrons
>>>>>>> ef26398330eb8def75e3d4261ae2679e9ddd8a06
        charge  =   0,
        basis   =   basis)


pymol.build()
print(pymol.nelectron)

nelec = 16
nmo = 26
ms = 0

h0 = np.load("ints_h0.npy")
h1 = np.load("ints_h1.npy")
h2 = np.load("ints_h2.npy")
C = np.load("mo_coeffs.npy")

fcidump.from_integrals('FCIDUMP', h1, h2, nmo, nelec, nuc=h0, ms=0)

ctx = fcidump.read("FCIDUMP");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, h.shape[1])

print("Ecore:")
print(ecore, " ", h0)
 
print(g.shape)

