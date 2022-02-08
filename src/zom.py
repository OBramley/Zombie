import op
import math, numpy
import in_outputs


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


def zom_gen(norb,ndet,type,filename):
    zstore=[]
    if(type == 'ran'):
        for i in range(ndet):
            zstore.append(zom(norb,typ='ran'))
    elif(type=='HF'):
        for i in range(ndet):
            zstore.append(zom(norb, typ='binary', ib=i))
    elif(type=='bb'):
        print('need to write biasing module')
    in_outputs.save_object(zstore,filename)

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
            bini = op.numtodet(ib,self.norb)
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
        return op.overlap_f(self.zs,self.zs)
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
        return op.numf(self.zs,self.zs)
    def isdet(self):
        return op.isdet(self.zs,self.norb)
    def sz(self):
        return op.szf(self.zs,self.zs,self.norb)
# Make a class for the Hamiltonian parameters and functions?

