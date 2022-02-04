from pyscf import gto, scf, ao2mo
import math, numpy
from functools import reduce
import inputs
import csv

def write(Hnuc, H1ei, H2ei,norb):
    with open('electron_integrals_1.csv','w',newline='')as csvfile:
        spamwriter=csv.writer(csvfile, delimiter=',')
        spamwriter.writerow([Hnuc,0,0,0,0,0,0,0,0,0])
        spamwriter.writerows(H1ei)
        for i in range(norb):
            for j in range(norb):
                spamwriter.writerows(H2ei[i,j,:,:])


    return

def read_in(norb, filename):
    H2ei =  numpy.zeros((norb,norb,norb,norb))
    file = open(filename)
    numpy_array=numpy.loadtxt(file,delimiter=',')
    Hnuc=numpy_array[0,0]
    H1ei=numpy_array[1:norb+1,:]
    counter=norb+1

    for i in range(norb):
        for j in range(norb):
            H2ei[i,j,0:norb,0:norb]=numpy_array[counter:(counter+norb),:]
            counter=+norb

    return Hnuc,H1ei,H2ei
        

                



# Program to generate one and two electron integrals from PyScf
# Some of this code has been adapted from George Booth
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

def pyscf_gen(norb):
    mol = gto.M(
    unit = inputs.pyscf['units'],
    atom = inputs.pyscf['atoms'],
    basis = inputs.pyscf['bs'],
    verbose = inputs.pyscf['verbosity'],
    symmetry = inputs.pyscf['symmetry'],
    spin=inputs.pyscf['spin'],
    charge=inputs.pyscf['charge'],
    symmetry_subgroup = inputs.pyscf['symmetry_subgroup'], #0 is code for A1 point group
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
    h1e=numpy.asarray(h1e)
    H1ei = spatospin1(h1e,norb)
    H2ei = spatospin2(eri_full,norb)
    write(Hnuc, H1ei, H2ei,norb)
    return Hnuc, H1ei, H2ei

# Routine to read in MOLPRO data
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
    
    write(Hnr, H1ei, H2ei,norb)
    return Hnr, H1ei, H2ei