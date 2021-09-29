from pyscf import gto, scf, ao2mo
import math, numpy
from functools import reduce
import hell

hell.hi()

def spatospin1(H1ea,norb):
    """Converting H1ea from spatial to spin"""
    if norb%2 != 0:
        raise ValueError('norb must be even')
    nspao = int(norb/2)
    H1ei = numpy.zeros((norb,norb))
    for i in range(nspao):
        for j in range(nspao):
            ii = i*2
            jj = j*2
            H1ei[ii,jj] = H1ea[i,j] # alpha spin
            H1ei[ii+1,jj+1] = H1ea[i,j] # beta spin
    return H1ei


# Program to generate one and two electron integrals from PyScf
# Some of this code has been adapted from George Booth

def pyscf_gen(norb):
    mol = gto.M(
    unit = 'Bohr',
    atom = 'Li 0 0 0; Li 0 0 6', #[['H', 0, 0, i] for i in range(6)],
    basis = '6-31g**',
    verbose = 1,
    symmetry = True,
    # spin=1,
    # charge=-1
    symmetry_subgroup = 0, #0 is code for A1 point group
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
    # Now convert from sppatial to spin orbitals
    print(numpy.shape(h1e))
    h1e=numpy.asarray(h1e)
    H1ei = spatospin1(h1e,norb)
    test = numpy.zeros((norb,norb))
    test2 = elec.spatospin(h1e,norb,test)
    # print(H1ei)
    # H2ei = spatospin2(eri_full,norb)
    # H1ei, H2ei = spatospin(h1e,eri_full,norb)
    return mol, myhf 
    # return Hnuc, H1ei, H2ei


pyscf_gen(10)










def spatospin2(H2ea,norb):
    """Converting H2ea from spatial orbitals to spin orbitals
    where the spatial orbitals are in the chemist's notation"""
    if norb%2 != 0:
        raise ValueError('norb must be even')
    nspao = int(norb/2)
    H2ei =  numpy.zeros((norb,norb,norb,norb))
    for i in range(nspao):
        for j in range(nspao):
            ii = i*2
            jj = j*2
            for k in range(nspao):
                for l in range(nspao):
                    kk = k*2
                    ll = l*2
                    Ht = H2ea[i,j,k,l]
                    H2ei[ii,kk,jj,ll] = Ht
                    H2ei[ii+1,kk,jj+1,ll] = Ht
                    H2ei[ii,kk+1,jj,ll+1] = Ht
                    H2ei[ii+1,kk+1,jj+1,ll+1] = Ht
    return H2ei






   
    



def molpro_read(filename,norb):
    file = open(filename,'r')
    if norb%2 != 0:
        raise ValueError('norb must be even')
    nspao = int(norb/2)
    H1e = numpy.zeros((nspao,nspao))
    H2e = numpy.zeros((nspao,nspao,nspao,nspao))
    Hnr = 0.0
    for line in file:
        bits = line.split()
        #print(len(bits))
        en = float(bits[0])
        i = int(bits[1]) - 1
        j = int(bits[2]) - 1
        k = int(bits[3]) - 1
        l = int(bits[4]) - 1
        #print(i,j,k,l,en)
        if k != -1 and l != -1:
            H2e[i,j,k,l] = en
            H2e[j,i,k,l] = en
            H2e[i,j,l,k] = en
            H2e[j,i,l,k] = en
            H2e[k,l,i,j] = en
            H2e[l,k,i,j] = en
            H2e[k,l,j,i] = en
            H2e[l,k,j,i] = en
            # print('2e',i,j,k,l,en)
        elif k == -1 and l == -1 and i != -1 and j != -1:
            H1e[i,j] = en
            if i != j:
                H1e[j,i] = en
        elif i == -1 and j == -1 and k == -1 and l == -1:
            Hnr = en
            # print('nucr',en)
        else:
            print('error',bits)
    H1ei = numpy.zeros((norb,norb))
    H2ei = numpy.zeros((norb,norb,norb,norb))
    for i in range(nspao):
        for j in range(nspao):
            ii = i*2
            jj = j*2
            H1ei[ii,jj] = H1e[i,j] # alpha spin
            H1ei[ii+1,jj+1] = H1e[i,j] # beta spin
            for k in range(nspao):
                for l in range(nspao):
                    kk = k*2
                    ll = l*2
                    Ht = H2e[i,j,k,l]
                    H2ei[ii,kk,jj,ll] = Ht
                    H2ei[ii+1,kk,jj+1,ll] = Ht
                    H2ei[ii,kk+1,jj,ll+1] = Ht
                    H2ei[ii+1,kk+1,jj+1,ll+1] = Ht
    return Hnr, H1ei, H2ei