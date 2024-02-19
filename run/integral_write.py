import numpy 
import csv 
from pyscf import gto, scf, ao2mo
from functools import reduce

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


def one_elec_setup(norb,H1ei,EXDIR1):

    h1count=0
    g=numpy.zeros(((norb*norb),3))

    for i in range(norb): # i is annihilation 
        for j in range(norb): #j is creation
            if(H1ei[i,j]!=0.0):
                h1count=h1count+1
                g[h1count,0]=H1ei[i,j]
                g[h1count,1]=int(i+1)
                g[h1count,2]=int(j+1)
               

    g[0,0]=h1count
    with open(EXDIR1+"/integrals/h1e.csv",'w', newline='')as csvfile:
            spamwriter=csv.writer(csvfile, delimiter=',')
            spamwriter.writerows(g[0:h1count+1,:])
           
        
    
    return


def two_elec_setup(norb,H2ei,EXDIR1):

    
    h2count=0
    g=numpy.zeros(((norb*norb*norb*norb),5))
    for i in range(norb):
        for j in range(norb):
            if(i!=j):
                for k in range(norb):
                    for l in range(norb):
                        if(k!=l):
                            if(H2ei[i,j,k,l]!=0.0):
                                h2count=h2count+1  
                                g[h2count,0]=H2ei[i,j,k,l]
                                g[h2count,1]=i+1
                                g[h2count,2]=j+1
                                g[h2count,3]=k+1
                                g[h2count,4]=l+1
    g[0,0]=h2count
   
    with open(EXDIR1+"/integrals/h2e.csv",'w', newline='')as csvfile:
        spamwriter=csv.writer(csvfile, delimiter=',')
        spamwriter.writerows(g[0:h2count+1,:])


    return 



def elec_writer(h1e,eri_full,norb,EXDIR1):

    H1ei = spatospin1(h1e,norb)
    H2ei = spatospin2(eri_full,norb)

    one_elec_setup(norb,H1ei,EXDIR1)
    two_elec_setup(norb,H2ei,EXDIR1)

def pyscf_do(pyscf_ins,norbs,EXDIR1):
    if(pyscf_ins['units']=='atom'):
        mol = gto.M(
        atom = pyscf_ins['atoms'],
        basis = pyscf_ins['bs'],
        verbose = pyscf_ins['verbosity'],
        spin=pyscf_ins['spin'],
        charge=pyscf_ins['charge'],
        # symmetry_subgroup = pyscf_ins['symmetry_subgroup'], #0 is code for A1 point group
        )
        myhf=scf.RHF(mol)
        myhf.kernel()
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
        verbose = pyscf_ins['verbosity'],
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

    elec_writer(h1e,eri_full,norbs,EXDIR1)

    with open(EXDIR1+"/integrals/hnuc.csv",'w', newline='')as csvfile:
        spamwriter=csv.writer(csvfile)
        spamwriter.writerow([Hnuc,0])


def molpro_read(filename,norbs,EXDIR1):
    file = open(filename,'r')
    if norbs%2 != 0:
        raise ValueError('norb must be even')
    nspao = int(norbs/2)
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
            Hnuc = en
            # print('nucr',en)
        else:
            print('error',bits)
    file.close()

    elec_writer(H1e,H2e,norbs,EXDIR1)

    with open(EXDIR1+"/integrals/hnuc.csv",'w', newline='')as csvfile:
        spamwriter=csv.writer(csvfile)
        spamwriter.writerow([Hnuc,0])
   
    
    
  
    
   



  