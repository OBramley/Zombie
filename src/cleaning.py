import numpy
import in_outputs
import zom
from itertools import combinations
import ham
import op

def clean_setup(norb, nel, spin, Ham, filenamer):
    filenamer=filenamer+"_clean"
    zstore=[]
    combs=list(combinations(range(norb),nel))
    ndet=len(combs)
    for i in range(ndet):
        zstore.append(zom.zom(norb,typ='clean',thetas=combs[i]))

    # output zombie states
    in_outputs.save_object(zstore,filenamer+'_zombie_states.pkl')

    # Generate the clean Hamiltonian
    cleanham, cleanovrl = ham.hamiltonian(ndet,Ham,zstore,filenamer)

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