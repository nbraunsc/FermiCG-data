import sys, os
import numpy as np
import scipy
import itertools
import time
from math import factorial
import copy as cp
import sys

from fermicluster import *
from pyscf_helper import *
import pyscf
from pyscf import gto, scf, ao2mo, molden, lo, mo_mapping, mcscf
import tools 

pyscf.lib.num_threads(1) #with degenerate states and multiple processors there can be issues
np.set_printoptions(suppress=True, precision=3, linewidth=1500)


molecule = '''
C     4.93354    -2.84833     0.00000
C     4.93321    -1.45861     0.00000
C     3.73075    -0.73168     0.00000
C     2.49032    -1.43783     0.00000
C     1.23513    -0.71311     0.00000
C     1.23513     0.71311     0.00000
C     2.49032     1.43783     0.00000
C     3.73075     0.73168     0.00000
C     4.93321     1.45861     0.00000
C     4.93354     2.84833     0.00000
C     3.72991     3.54302     0.00000
C     2.49908     2.86514     0.00000
C     1.23174     3.59686    -0.00000
C     1.20340     5.00170    -0.00000
C     0.00000     5.69682    -0.00000
C    -1.20340     5.00170    -0.00000
C    -1.23174     3.59686    -0.00000
C     0.00000     2.87564    -0.00000
C     0.00000     1.42626    -0.00000
C    -1.23513     0.71311    -0.00000
C    -2.49032     1.43783     0.00000
C    -2.49908     2.86514    -0.00000
C    -3.72991     3.54302    -0.00000
C    -4.93354     2.84833     0.00000
C    -4.93321     1.45861     0.00000
C    -3.73075     0.73168     0.00000
C    -3.73075    -0.73168     0.00000
C    -4.93321    -1.45861     0.00000
C    -4.93354    -2.84833     0.00000
C    -3.72991    -3.54302     0.00000
C    -2.49908    -2.86514    -0.00000
C    -2.49032    -1.43783    -0.00000
C    -1.23513    -0.71311    -0.00000
C    -0.00000    -1.42626    -0.00000
C    -0.00000    -2.87564    -0.00000
C    -1.23174    -3.59686    -0.00000
C    -1.20340    -5.00170    -0.00000
C     0.00000    -5.69682    -0.00000
C     1.20340    -5.00170    -0.00000
C     1.23174    -3.59686    -0.00000
C     2.49908    -2.86514    -0.00000
C     3.72991    -3.54302     0.00000
H     5.87949    -3.39439     0.00000
H     5.89025    -0.93984     0.00000
H     5.89025     0.93984     0.00000
H     5.87949     3.39439     0.00000
H     3.75925     4.63123     0.00000
H     2.13117     5.57119    -0.00000
H     0.00000     6.78906    -0.00000
H    -2.13117     5.57119    -0.00000
H    -3.75925     4.63123    -0.00000
H    -5.87949     3.39439     0.00000
H    -5.89025     0.93984     0.00000
H    -5.89025    -0.93984     0.00000
H    -5.87949    -3.39439     0.00000
H    -3.75925    -4.63123     0.00000
H    -2.13117    -5.57119    -0.00000
H     0.00000    -6.78906    -0.00000
H     2.13117    -5.57119    -0.00000
H     3.75925    -4.63123    -0.00000
'''

cas_nel = 42
cas_norb = 42



#PYSCF inputs
mol = gto.Mole(atom=molecule,
    symmetry = True,basis = 'ccpvdz' )
mol.build()
print("symmertry: ",mol.topgroup)

#SCF
mf = scf.RHF(mol)
mf.verbose = 4
mf.conv_tol = 1e-12
mf.conv_tol_grad = 1e-9
mf.run(max_cycle=200)

## Active space selection


h,ecore,g,C = get_pi_space(mol,mf,cas_norb,cas_nel,local=True)
print("core %16.8f"%ecore)



mc = mulliken_ordering(mol,h.shape[0],C)
idx = np.where(mc>.9)[1]  #gives index map from atom to local orbital corresponding to that orbital

# Reorder
h,g = reorder_integrals(idx,h,g)
print(h)
C = C[:,idx] # make sure u reorder this too
molden.from_mo(mol, 'cas.molden', C)


from pyscf import tools
tools.fcidump.from_integrals('hexabenzo-fcidump', h, g, cas_norb,
                                 cas_nel, nuc=ecore, ms=0)


