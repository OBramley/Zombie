import inputs
from pyscf import gto, scf, ao2mo, fci, mcscf,lib 
from functools import reduce
import numpy
from pyscf.fci import cistring,fci_slow
from pyscf.fci import fci_slow


if(inputs.pyscf['units']=='atom'):
        mol = gto.M(
        atom = inputs.pyscf['atoms'],
        basis = inputs.pyscf['bs'],
        verbose = inputs.pyscf['verbosity'],
        spin=inputs.pyscf['spin'],
        charge=inputs.pyscf['charge'],
        # symmetry_subgroup = inputs.pyscf['symmetry_subgroup'], #0 is code for A1 point group
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
    unit = inputs.pyscf['units'],
    atom = inputs.pyscf['atoms'],
    basis = inputs.pyscf['bs'],
    verbose = inputs.pyscf['verbosity'],
    symmetry = inputs.pyscf['symmetry'],
    spin=inputs.pyscf['spin'],
    charge=inputs.pyscf['charge'],
    # symmetry_subgroup = inputs.pyscf['symmetry_subgroup'], #0 is code for A1 point group
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


fci_solver = fci.FCI(mol, myhf.mo_coeff)
e_fci, myci = fci_solver.kernel(h1e=h1e, eri=eri)
print(e_fci)
print((myci.shape))
count=0
for i in range(len(myci)):
    for j in range(len(myci[0])):
        if(myci[i][j]!=0):
            count+=1
            
print(count)
