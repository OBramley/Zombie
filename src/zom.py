import math, numpy, scipy
from functools import reduce
from pyscf import gto, scf, ao2mo
from pyscf import tools
from pyscf import symm

def make(norb):
    # Returns a completely empty zombie state (all zeros)
    zz = numpy.zeros((norb,2))
    return zz

def populate(zz,occ):
    zz[:,:] = 0.
    for iorb in range(len(zz[:,0])):
        zz[iorb,occ[iorb]] = 1.
    return zz

def new(norb,occ):
    zz = make(norb)
    zz = populate(zz,occ)
    return zz

def new_ran(norb):
    rantemp = 2.0*math.pi*numpy.random.random((norb))
    zz = make(norb)
    zz[:,0] = numpy.cos(rantemp)
    zz[:,1] = numpy.sin(rantemp)
    return zz

# Trying a zombie class
class zom:
    """Class of a single zombie state"""
    def __init__(self,norb,typ='empty',occ=None,coefs=None, \
                 ib=None,nel=None,thetas=None):
        self.norb = norb
        # Possible types if initialisation
        if typ == 'empty':
            self.zs = make(self.norb)
            self.zs[:,0] = 0.0
        elif typ == 'occ':
            self.zs = new(self.norb,occ)
        elif typ == 'ran':
            self.zs = new_ran(self.norb)
        elif typ == 'coef':
            self.zs = coefs
        elif typ == 'binary':
            bini = numtodet(ib,self.norb)
            self.zs = new(self.norb,bini)
        elif typ == 'aufbau':
            occ = numpy.zeros((self.norb),dtype=int)
            occ[:nel] = 1
            self.zs = new(self.norb,occ)
        elif typ == 'theta':
            self.zs = make(self.norb)
            self.zs[:,0] = numpy.cos(2.0*math.pi*thetas)
            self.zs[:,1] = numpy.sin(2.0*math.pi*thetas)
        else:
            raise ValueError('Invalid type')
    def ov(self):
        return overlap_f(self.zs,self.zs)
    def cr(self,iorb):
        self.zs[iorb,1] = self.zs[iorb,0]
        self.zs[iorb,0] = 0.0
        self.zs[:iorb,1] *= -1.0
    def an(self,iorb):
        self.zs[iorb,0] = self.zs[iorb,1]
        self.zs[iorb,1] = 0.0
        self.zs[:iorb,1] *= -1.0
    def num(self):
        "Number of electrons in the zombie state"
        return numf(self.zs,self.zs)
    def isdet(self):
        return isdet(self.zs,self.norb)
    def sz(self):
        return szf(self.zs,self.zs,self.norb)
# Make a class for the Hamiltonian parameters and functions?
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
                zomt = an(zomt,ii)
                zomt = cr(zomt,jj)
                ov = overlap(zom1,zomt)
                # print(ii,jj,ov)
                Ht1 += ov*self.H1ei[ii,jj]
        return Ht1
    def Ham2z_v3(self,zom1,zom2):
        """Algorithm to calculate the two-electron Hamiltonian bra-ket
        Between two zombie states
        Breaking up annihilation and creation calculation
        Scales as O(M^5) but with a low prefactor
        Needs some attention to deal with complex numbers"""
        Ht2 = 0.0
        Z1ij = numpy.zeros((self.norb,self.norb,self.norb,2),dtype=complex)
        Z2lk = numpy.zeros((self.norb,self.norb,self.norb,2),dtype=complex)
        for ii in range(self.norb):
            for jj in range(self.norb):
                zomt = numpy.copy(zom1)
                zomt = an(zomt,ii)
                zomt = an(zomt,jj)
                Z1ij[ii,jj,:,:] = zomt[:,:]
        for kk in range(self.norb):
            for ll in range(self.norb):
                zomt = numpy.copy(zom2)
                zomt = an(zomt,kk)
                zomt = an(zomt,ll)
                Z2lk[ll,kk,:,:] = zomt[:,:]
        for ii in range(self.norb):
            if zom1[ii,1] == 0.0:
                continue
            for jj in range(self.norb):
                if zom1[jj,1] == 0.0:
                    continue
                for kk in range(self.norb):
                    if zom2[kk,1] == 0.0:
                        continue 
                    for ll in range(self.norb):
                        if self.H2ei[ii,jj,kk,ll] == 0.0 or zom2[ll,1] == 0.0:
                            continue
                        ov = overlap_f(Z1ij[ii,jj,:,:],Z2lk[ll,kk,:,:])
                        if ov == 0.0:
                            continue
                        Ht2 += 0.5*ov*self.H2ei[ii,jj,kk,ll]
        return Ht2
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
                zomt = an(zomt,ii)
                zomt = an(zomt,jj)
                Z1ij[ii,jj,:,:] = zomt[:,:]
        for kk in range(self.norb):
            zomt = numpy.copy(zom2)
            zomt = an(zomt,kk)
            Z2k[kk,:,:] = zomt[:,:]
        for ii in range(self.norb):
            if zom1[ii,1] == 0.0:
                continue
            ispin = ii%2
            for jj in range(self.norb):
                if iszero(Z1ij[ii,jj,:,:]):
                    continue
                jspin = jj%2
                for kk in range(ispin,self.norb,2):
                    if zom2[kk,1] == 0.0:
                        continue
                    Ht2 += z_an_z3(Z1ij[ii,jj,:,:],Z2k[kk,:,:], \
                                        self.norb,self.H2ei[ii,jj,kk,:])
        return 0.5*Ht2
    def HTot(self,zom1,zom2):
        H1et = self.Ham1z(zom1,zom2)
        H2et = self.Ham2z_v5(zom1,zom2)
        HH = H1et + H2et + self.Hnr*overlap_f(zom1,zom2)
        return HH