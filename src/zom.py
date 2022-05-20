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
            zz[iorb,1]=1
        else:
            zz[iorb,0]=1
    return zz




def biased_basis(norb,zom_bias,ndet,zstore):
    orbitals=int(norb/2)
    alive=zom_bias['alive']
    astart=zom_bias['alive_start']
    aend=zom_bias['alive_end']
    mid=zom_bias['mid']
    dstart=zom_bias['dead_start']
    dend=zom_bias['dead_end']

    def bform2(x):
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
    # print(astartx)
    # print(aendx)
    # print(dead)
    # print(dstartx)
    # print(dendx)
    # print(astep)
    # print(dstep)

    values[0]=astartx
    values[1]=astartx
    for i in range(1,alive):
        val=(astartx+(i*astep))
        j=2*i
        values[j]=val
        values[j+1]=val
        # print(values)
    
    values[2*(alive+mid)]=dstartx
    values[2*(alive+mid)+1]=dstartx

    for i in range(1,dead):
        val=dstartx+(i*dstep)
        j= 2*(i+alive+mid)
        values[j]=val
        values[j+1]=val
    
    # values=bform2(values)

    # print(values)
    # print('complete')
    # randoms=(numpy.random.uniform(0.0,0.25,size=(ndet,norb)))
    # randoms=numpy.exp(-randoms)
    randoms=numpy.zeros((ndet,norb))
    for i in range(norb):
        randoms[:,i]=numpy.random.normal(loc=values[i],scale=0.01,size=ndet)

    randoms=bform2(randoms)
    
    for i in range(ndet):
        theta=randoms[i,:]
        zstore.append(zom(norb,typ='biased',thetas=theta))
    
    for i in range(ndet):
        in_outputs.write_ham(zstore[i].zs,'zombie_'+str(i)+'.csv')

    return zstore

def biased_basis2(norb,ndet, zstore):
    for idet in range(ndet):
        num = numpy.zeros(norb)
        for i in range(norb):
            if(i<4):
                num[i]=0.25
            elif((i>3)and(i<6)):
                mu=0.25
                sig=0.175
                num[i]=numpy.random.normal(mu,sig)
            elif((i>5)and(i<7)):
                mu=0
                sig=0.351
                num[i]=numpy.random.normal(mu,sig)
            elif((i>6)and(i<8)):
                mu=0
                sig=0.351
                num[i]=numpy.random.normal(mu,sig)
            else:
                mu=0
                sig=0.120
                num[i]=numpy.random.normal(mu,sig)
        iteration=zom(norb,typ='theta',thetas=num)
        zstore.append(iteration)
    for i in range(ndet):
        in_outputs.write_ham(zstore[i].zs,'zombie_'+str(i)+'.csv')

    return zstore

def biased_basis3(norb, ndet, zstore):
    test=numpy.array((0.25,0.213632469,0.193380738,0.001262455,0.000505343,0.00062495,0.000530594,9.57371E-06,0.000169358,3.27753E-05,0.004644281,0.000396432,0.000387224,5.15685E-05,0.004644276,0.000396434,0.000387213,5.16551E-05,9.58165E-06))
    random=numpy.exp(-numpy.random.uniform(0.0,0.02021,size=(ndet,norb)))
    for idet in range(ndet):
        randompick=random[idet,:]
        num=numpy.zeros(norb)
        for j in range(0,norb,2):
            centre=test[int(j/2)]
            num[j]=centre*randompick[j]
            num[j+1]=centre*randompick[j+1]
        zstore.append(zom(norb,typ='theta',thetas=num))

    for i in range(ndet):
        in_outputs.write_ham(zstore[i].zs,'zombie_'+str(i)+'.csv')
    
    return zstore

def zom_gen(norb,ndet,zombs,filenamer,seed,zom_bias):
    type=zombs['zomtyp']
    rhf=zombs['rhf_1']
    nelc=zombs['nel']
    zstore=[]
    numpy.random.seed(seed)
    if(type == 'ran'):
        # Random basis
        for i in range(ndet):
            zstore.append(zom(norb,typ='ran'))
    # HArtree Fock determinant basis
    elif(type=='HF'):
        for i in range(ndet):
            zstore.append(zom(norb, typ='binary', ib=i))
    # Biased Basis
    elif(type=='bb'):
        # zstore=biased_basis(norb,zom_bias,ndet,zstore)
        # zstore=biased_basis2(norb,ndet,zstore)
        zstore=biased_basis3(norb,ndet,zstore)
    if((rhf=='y')):
        zstore[0]=zom(norb,typ='aufbau',nel=nelc)
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
        elif typ=='biased':
            self.zs=make(self.norb)
            self.zs[:,0]=numpy.sqrt((1-numpy.square(thetas)))
            self.zs[:,1]=thetas
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



def biased_check(zstore,norb, ndet, filenamer):
    dvec=in_outputs.read_dvec(filenamer+"_dvect.csv" )
    checker=numpy.zeros((norb,2))
    # print(ndet)
    for i in range(ndet):
        in_outputs.write_ham(zstore[i].zs,'zombie_'+str(i)+'.csv')
        checker=checker+zstore[i].zs
    
    # print(checker)
    # print(checker/ndet)
    # print(dvec)
    return
