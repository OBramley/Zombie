from pyscf import gto, scf, ao2mo
import numpy
from functools import reduce
import op
import in_outputs

        # Program to generate one and two electron integrals from PyScf
# Some of this code has been adapted from George Booth
def spatospin1(H1ea,norb):
    print(norb)
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
    # in_outputs.write_ham(H1ei,"H1ei.csv")
    
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

    # for i in range(norb):
    #     for j in range(norb):
    #         obj=H2ei[i,j,:,:]
    #         in_outputs.write_ham(obj,"H2ei_"+str(i)+"_"+str(j)+".csv")
    return H2ei

class system:
    """Class holding the system parameters for zombie state calculation"""
    def __init__(self, norb, Hnr, H1ei, H2ei):
        self.Hnr = Hnr
        self.H1ei = H1ei
        self.H2ei = H2ei
        self.norb = norb
        if norb%2 != 0:
            raise ValueError('norb must be even')
        self.nspao = int(norb/2)
    def Ham1z(self,zom1,zom2):
        Ht1 = 0.0
        for ii in range(self.norb):
            for jj in range(self.norb):
                zomt = numpy.copy(zom2)
                zomt = op.an(zomt,ii)
                zomt = op.cr(zomt,jj)
                ov = op.overlap(zom1,zomt)
                # print(ii,jj,ov)
                Ht1 += ov*self.H1ei[ii,jj]
        return Ht1
    def Ham2z_v5(self,zom1,zom2):
        Ht2 = 0.0
        if type(zom1[0,0].item()) is complex:
            Z1ij = numpy.zeros((self.norb,self.norb,self.norb,2),dtype=complex)
        else:
            Z1ij = numpy.zeros((self.norb,self.norb,self.norb,2),dtype=float)
        if type(zom2[0,0].item()) is complex:
            Z2k = numpy.zeros((self.norb,self.norb,2),dtype=complex)
        else:
            Z2k = numpy.zeros((self.norb,self.norb,2),dtype=float)
        for ii in range(self.norb):
            for jj in range(self.norb):
                zomt = numpy.copy(zom1)
                zomt = op.an(zomt,ii)
                zomt = op.an(zomt,jj)
                Z1ij[ii,jj,:,:] = zomt[:,:]
        for kk in range(self.norb):
            zomt = numpy.copy(zom2)
            zomt = op.an(zomt,kk)
            Z2k[kk,:,:] = zomt[:,:]
        for ii in range(self.norb):
            if zom1[ii,1] == 0.0:
                continue
            ispin = ii%2
            for jj in range(self.norb):
                if op.iszero(Z1ij[ii,jj,:,:]):
                    continue
                jspin = jj%2
                for kk in range(ispin,self.norb,2):
                    if zom2[kk,1] == 0.0:
                        continue
                    Ht2 += op.z_an_z3(Z1ij[ii,jj,:,:],Z2k[kk,:,:], \
                                        self.norb,self.H2ei[ii,jj,kk,:])
        return 0.5*Ht2
    def HTot(self,zom1,zom2):
        H1et = self.Ham1z(zom1,zom2)
        H2et = self.Ham2z_v5(zom1,zom2)
        HH = H1et + H2et + self.Hnr*op.overlap_f(zom1,zom2)
        return HH


def pyscf_gen(norb,filenamer,pyscf):
    mol = gto.M(
    unit = pyscf['units'],
    atom = pyscf['atoms'],
    basis = pyscf['bs'],
    verbose = pyscf['verbosity'],
    symmetry = pyscf['symmetry'],
    spin=pyscf['spin'],
    charge=pyscf['charge'],
    symmetry_subgroup = pyscf['symmetry_subgroup'], #0 is code for A1 point group
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
    print(Hnuc)
    # in_outputs.write_integrals(Hnuc, H1ei, H2ei,norb)
    Ham=system(norb, Hnuc, H1ei, H2ei)
    filename=filenamer+'_inegrals.pkl'
    in_outputs.save_object(Ham,filename)
    return Ham

# Routine to read in MOLPRO data
def molpro_read(filename,norb,filenamer):
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
    
    # in_outputs.write_integrals(Hnr, H1ei, H2ei,norb)
    Ham=system(norb, Hnr, H1ei, H2ei)
    filename1=filenamer+'_inegrals.pkl'
    in_outputs.save_object(Ham,filename1)
    return Ham

def hamiltonian(ndet,Ham,zstore,filenamer):
    Kover = numpy.zeros((ndet,ndet))
    Bigham = numpy.zeros((ndet,ndet))
    for idet in range(ndet):
        for jdet in range(idet,ndet):
            Kover[idet,jdet] = op.overlap_f(zstore[idet].zs, zstore[jdet].zs)
            Kover[jdet,idet] = Kover[idet,jdet]
            Bigham[idet,jdet] = Ham.HTot(zstore[idet].zs, zstore[jdet].zs)
            Bigham[jdet,idet] = Bigham[idet,jdet]
        print('Hamiltonian row '+ str(idet)+' completed')
    hamname=filenamer+'_hamiltonian.csv'
    ovrlname=filenamer+'_overlap.csv'
    in_outputs.write_ham(Kover,ovrlname)
    in_outputs.write_ham(Bigham,hamname)
    return Bigham, Kover