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

def spf(norb,thetas):
    zz=make(norb)
    for iorb in range(norb):
        if iorb in thetas:
            zz[iorb,0]=1
        else:
            zz[iorb,1]=1
    return zz


def biased_basis(norb,zom_bias,ndet,zstore):
    orbitals=int(norb/2)
    alive=zom_bias['alive']
    astart=zom_bias['alive_start']
    aend=zom_bias['alive_end']
    mid=zom_bias['mid']
    dstart=zom_bias['dead_start']
    dend=zom_bias['dead_end']

    def bform(x):
        return 1/2*(1-numpy.tanh((x)))

    def bforminv(y):
        return numpy.arctanh(-1*(2*y-1))

    values=numpy.zeros(norb)
    astartx=bforminv(astart)
    aendx=bforminv(aend)
    dead=orbitals-alive-mid
    dstartx=bforminv(dstart)
    dendx=bforminv(dend)
    astep=abs((astartx-aendx)/(alive-1))
    dstep=abs((dstartx-dendx)/(dead-1))

    for i in range(alive):
        val=bform((astartx+(i*astep)))
        j=2*i
        values[j]=val
        values[j+1]=val
    
    values[(2*alive):2*(alive+mid)]=0.5

    for i in range(dead):
        val=bform(dstartx+(i*dstep))
        j= 2*(i+alive+mid)
        values[j]=val
        values[j+1]=val
    
    randoms=(numpy.random.uniform(0.0,0.25,size=(ndet,norb)))
    randoms=numpy.exp(-randoms)
    for i in range(ndet):
        theta=numpy.multiply(values,randoms[i,:])
        zstore.append(zom(norb,typ='theta',thetas=theta))
   
    return zstore
    

def zom_gen(norb,ndet,type,filenamer,seed,zom_bias):
    zstore=[]
    numpy.random.seed(seed)
    if(type == 'ran'):
        for i in range(ndet):
            zstore.append(zom(norb,typ='ran'))
    elif(type=='HF'):
        for i in range(ndet):
            zstore.append(zom(norb, typ='binary', ib=i))
    elif(type=='bb'):
        zstore=biased_basis(norb,zom_bias,ndet,zstore)
    filename=filenamer+'_zombie_states.pkl'
    in_outputs.save_object(zstore,filename)
    return zstore

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
        elif typ == 'clean':
            self.zs = spf(norb,thetas)
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

