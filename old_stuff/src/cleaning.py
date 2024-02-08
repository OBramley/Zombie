from xml.dom import NoDataAllowedErr
import numpy
import in_outputs
import zom
from itertools import combinations
import ham
import op

def clean_setup(norb, nel,ndet, spin, Ham, zstore, filenamer):
    filenamer=filenamer+"_clean"
    cstore=[]
    combs=numpy.asarray(list(combinations(range(norb),nel)))
    rawdet=combs.shape[0]
    ndew=0
    spinchecker=combs%2
    for i in range(rawdet):
        if(numpy.sum((spinchecker[i]))==(((nel/2)-spin))):
            cstore.append(zom.zom(norb,typ='clean',thetas=combs[i]))
            ndew=ndew+1
    
    magnitude=numpy.zeros((ndew))
    cstore_f=[]
    ndew2=0
  
    for i in range(ndew):
        for j in range(ndet):
            magnitude[i]=magnitude[i]+op.overlap_f(cstore[i].zs,zstore[j].zs)
        if(magnitude[i]>0.00001):
            cstore_f.append(cstore[i])
            ndew2=ndew2+1
    cstore_f=cstore
    del cstore

    # in_outputs.save_object(cstore,filenamer+'_zombie_states.pkl')
    # cleanham, cleanovrl = ham.hamiltonian(ndew,Ham,cstore,filenamer)

    # output zombie states
    in_outputs.save_object(cstore_f,filenamer+'_zombie_states.pkl')

    # Generate the clean Hamiltonian
    cleanham, cleanovrl = ham.hamiltonian(ndew2,Ham,cstore_f,filenamer)

    return

def cleaner(ndet, zstore, ham_name, zom_name,filenamer):
    clean_ham = in_outputs.read_ham(ham_name)
    sdstore = in_outputs.read_object(zom_name)
    dvec = in_outputs.read_dvec(filenamer+'_dvect.csv')
    ndet_clean=len(sdstore)
    dvec_new=numpy.zeros((ndet_clean))
    norm=0

    for kdet in range(ndet_clean):
        hfz=sdstore[kdet].zs
        dvec_c=numpy.conjugate(dvec)
        for idet in range(ndet):
            izom=zstore[idet].zs
            ioverlap=op.overlap_f(hfz,izom)
            dvec_new[kdet]=dvec_new[kdet]+(dvec[idet]*ioverlap)
            for jdet in range(ndet):
                jzom=zstore[jdet].zs
                val=dvec_c[jdet]*dvec[idet]*op.overlap_f(jzom,hfz)*ioverlap
                norm=norm+val
    
    in_outputs.write_dvec(dvec_new,filenamer+'_clean_dvect.csv')
    energy=numpy.einsum('i,ij,j',dvec_new,clean_ham,dvec_new)

    return energy, norm