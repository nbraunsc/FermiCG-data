import sys, os
import numpy as np
import scipy
import itertools
import time
from math import factorial
import copy as cp
import sys
import tools 

from fermicluster import *
from pyscf_helper import *
import pyscf
from pyscf import gto, scf, ao2mo, molden, lo, mo_mapping, mcscf

sys.path.append('/Users/nicole/My Drive/code/OrbitalPartitioning/orbitalpartitioning')
import orbitalpartitioning
from orbitalpartitioning import *

pyscf.lib.num_threads(1) #with degenerate states and multiple processors there can be issues
np.set_printoptions(suppress=True, precision=3, linewidth=1500)


molecule = '''
C       4.6336867757    -1.2242295397     0.0000000000
C       4.4255614400     0.1700917199     0.0000000000
C       5.9160648559    -1.7742824642     0.0000000000
C       5.5765818320     0.9841959654     0.0000000000
C       7.0417938838    -0.9473143671     0.0000000000
C       6.8604981940     0.4376683883     0.0000000000
C       0.5927492072     0.4581981628     0.0000000000
C       0.4451033014     1.8584342989     0.0000000000
C       1.9041039939    -0.0505338987     0.0000000000
C       1.5640941394     2.6887327507     0.0000000000
C       3.0514687259     0.7633813258     0.0000000000
C       2.8520613255     2.1570405955     0.0000000000
C      -6.8604981940    -0.4376683883     0.0000000000
C      -7.0417938838     0.9473143671     0.0000000000
C      -5.5765818320    -0.9841959654     0.0000000000
C      -5.9160648559     1.7742824642     0.0000000000
C      -4.4255614400    -0.1700917199     0.0000000000
C      -4.6336867757     1.2242295397     0.0000000000
C      -3.0514687259    -0.7633813258     0.0000000000
C      -2.8520613255    -2.1570405955     0.0000000000
C      -1.5640941394    -2.6887327507     0.0000000000
C      -0.4451033014    -1.8584342989     0.0000000000
C      -0.5927492072    -0.4581981628     0.0000000000
C      -1.9041039939     0.0505338987     0.0000000000
H       3.7849879120    -1.9076345469     0.0000000000
H       6.0325833442    -2.8604965477     0.0000000000
H       5.4821758879     2.0696091478     0.0000000000
H       8.0463585128    -1.3756739224     0.0000000000
H      -0.5426199813     2.3168156796     0.0000000000
H       2.0368540411    -1.1280906363     0.0000000000
H       1.4298106521     3.7726953568     0.0000000000
H       3.6983932318     2.8425398416     0.0000000000
H      -8.0463585128     1.3756739224     0.0000000000
H      -5.4821758879    -2.0696091478     0.0000000000
H      -6.0325833442     2.8604965477     0.0000000000
H      -3.7849879120     1.9076345469     0.0000000000
H      -3.6983932318    -2.8425398416     0.0000000000
H      -1.4298106521    -3.7726953568     0.0000000000
H       0.5426199813    -2.3168156796     0.0000000000
H      -2.0368540411     1.1280906363     0.0000000000
H      -7.7260669052    -1.1042415323     0.0000000000
H       7.7260669052     1.1042415323     0.0000000000
'''
cas_nel = 24
cas_norb = 24


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

F = mf.get_fock()

## Active space selection

h,ecore,g,C = get_pi_space(mol,mf,cas_norb,cas_nel,local=True)

# Run a CAS-CI calculation for comparison
#from pyscf import fci
#cisolver = fci.direct_spin1.FCI()
#ecas, vcas = cisolver.kernel(h, g, cas_norb, nelec=cas_nel, ecore=ecore,nroots =1,verbose=100)
#print("CAS-CI:%10.8f"%(ecas))
#print(" CASCI           %12.8f      Dim:%6d" % (ecas,vcas.shape[0]*vcas.shape[1]))
print("core %16.8f"%ecore)


mc = mulliken_ordering(mol,h.shape[0],C)
idx = np.where(mc>.9)[1]  #gives index map from atom to local orbital corresponding to that orbital

# Reorder
h,g = reorder_integrals(idx,h,g)
print(h)
C = C[:,idx] # make sure u reorder this too
molden.from_mo(mol, 'cas.molden', C)

print(C.shape)
Cfrags = []
Cfrags.append(C[:,0:5])
Cfrags.append(C[:,6:11])
Cfrags.append(C[:,12:17])
Cfrags.append(C[:,18:23])

Cf = canonicalize(Cfrags, F)

#make integrals with pseudo-canonicalized coeffs
mycas = mcscf.CASCI(mf, cas_norb, cas_nel)

Cfinal = np.hstack(Cf)
print(Cfinal.shape)
# Get the active space integrals and the frozen core energy
h, ecore = mycas.get_h1eff(Cfinal)
g = ao2mo.kernel(mol,Cfinal, aosym = 's4', compact = False).reshape(4*((cas_norb),))
#g = ao2mo.kernel(mol,Cfinal[:,focc:focc+cas_norb], aosym = 's4', compact = False).reshape(4*((cas_norb),))
print(h.shape)
print(g.shape)

from pyscf import tools
tools.fcidump.from_integrals('FCIDUMP_canonalized', h, g, cas_norb,
                                 cas_nel, nuc=ecore, ms=0)

