import json
import sys
from pyscf import gto, scf, ao2mo, cc, fci, mcscf,lib 
from functools import reduce
import numpy

name=sys.argv[1]
with open(name) as json_file:
    pyscf_ins=json.load(json_file)

if(pyscf_ins['units']=='atom'):
        mol = gto.M(
        atom = pyscf_ins['atoms'],
        basis = pyscf_ins['bs'],
        verbose = 4,#pyscf_ins['verbosity'],
        spin=pyscf_ins['spin'],
        charge=pyscf_ins['charge'],
        # symmetry_subgroup = pyscf_ins['symmetry_subgroup'], #0 is code for A1 point group
        )
        myhf=scf.RHF(mol)
        myhf.kernel()
        eri = mol.intor('int2e')
        
        c = myhf.mo_coeff
        # Get 1-electron integrals and convert to MO basis
        h1e = reduce(numpy.dot, (c.T, myhf.get_hcore(), c))
        # Get 2-electron integrals and transform them
        eri = ao2mo.kernel(mol, c)
        # Ignore all permutational symmetry, and write as four-index tensor, in chemical notation
        eri_full = ao2mo.restore(1, eri, c.shape[1])
        # Scalar nuclear repulsion energy
        Hnuc = myhf.energy_nuc()    

       
else:
    mol = gto.M(
    unit = pyscf_ins['units'],
    atom = pyscf_ins['atoms'],
    basis = pyscf_ins['bs'],
    verbose = 4,#pyscf_ins['verbosity'],
    symmetry = pyscf_ins['symmetry'],
    spin=pyscf_ins['spin'],
    charge=pyscf_ins['charge'],
    # symmetry_subgroup = pyscf_ins['symmetry_subgroup'], #0 is code for A1 point group
    )
    myhf=scf.RHF(mol)
    myhf.kernel()
    """Obtaining one and two electron integrals from pyscf calculation
    Code adapted from George Booth"""
    # Extract AO->MO transformation matrix
    c = myhf.mo_coeff
    # Get 1-electron integrals and convert to MO basis
    h1e = reduce(numpy.dot, (c.T, myhf.get_hcore(), c))
    # Get 2-electron integrals and transform them
    eri = ao2mo.kernel(mol, c)
    # Ignore all permutational symmetry, and write as four-index tensor, in chemical notation
    eri_full = ao2mo.restore(1, eri, c.shape[1])
    # Scalar nuclear repulsion energy
    Hnuc = myhf.energy_nuc()

mycc = cc.CCSD(myhf).run()
print('CCSD total energy', mycc.e_tot)
et = mycc.ccsd_t()
print('CCSD(T) total energy', mycc.e_tot + et)

fci_solver = fci.FCI(mol, myhf.mo_coeff)
e_fci, myci = fci_solver.kernel(h1e=h1e, eri=eri)
print(e_fci)
print((myci.shape))

num_determinants = len(fci_solver.ci)
print("Number of Slater determinants:", num_determinants)