import numpy
import zom
import op
import math
import in_outputs

# Routine to set he starting vector to a restricted Hartree Fock energy of a deterimant of 
# a specified number
def hf_state(norb,hfnum,zstore,ndet,Ki):
    zrhf=zom.zom(norb,typ='aufbau',nel=hfnum)
    dtcr=numpy.zeros((ndet))
    for i in range(ndet):
        dtcr[i]=op.overlap_f(zrhf.zs,zstore[i].zs)
    dvec=numpy.dot(Ki,dtcr)
    return dvec

# Imaginar time propagation routine
def itime_prop(Bigham, Kover, beta, steps, norb, hfflg, hfnum, zstore, ndet,filenamer):
    Ki = numpy.linalg.inv(Kover)
    in_outputs.write_ham(Ki,filenamer+'_inverse.csv')
    if(hfflg=='y'):
        dvec=hf_state(norb, hfnum, zstore, ndet, Ki)
    elif(hfflg=='n'):
        dvec=numpy.zeros((ndet))
        dvec[0]=1
    KinvH=numpy.matmul(Ki,Bigham)
    in_outputs.write_ham(Ki,filenamer+'_kinvh.csv')
    db=beta/steps
    eb=numpy.zeros((steps+1,2))
    for i in range(steps+1):
        den=numpy.einsum('i,ij,j',dvec,Bigham,dvec)
        eb[i,0]=i*db
        eb[i,1]=den
        ddot = -numpy.dot(KinvH,dvec)
        dvec = dvec + db*ddot
        norm = abs(numpy.einsum('i,ij,j',dvec,Kover,dvec))
        dvec /= math.sqrt(norm)
    filename=filenamer+'_energy.csv'
    filename1=filenamer+"_dvect.csv" 
    in_outputs.write_ham(eb,filename)
    in_outputs.write_dvec(dvec,filename1)
    return eb


def gs(dvecs,Kr,nst):
    # Gram-Schmidt orthogonalization
    uecs = numpy.copy(dvecs)
    for ist in range(1,nst):
        for jst in range(0,ist):
            numer = numpy.dot(dvecs[:,ist],numpy.dot(Kr,uecs[:,jst]))
            den = numpy.dot(uecs[:,jst],numpy.dot(Kr,uecs[:,jst]))
            uecs[:,ist] -= uecs[:,jst]*numer/den
    for ist in range(nst):
        fac = numpy.dot(uecs[:,ist],numpy.dot(Kr,uecs[:,ist]))
        uecs[:,ist] /= math.sqrt(fac)
    return uecs

# Imaginar time propagation routine
def itime_prop_gs(Bigham, Kover, beta, steps, norb, hfflg, hfnum, zstore, ndet, gramnum, filenamer):
    Ki = numpy.linalg.inv(Kover)
    if(hfflg=='y'):
        dvec=hf_state(norb, hfnum, zstore, ndet, Ki)
    elif(hfflg=='n'):
        dvec=numpy.zeros((ndet))
        dvec[0]=1

    dvecs = numpy.zeros((ndet,gramnum))
    dvecs[:,0]=dvec
    for i in range(1,gramnum):
        dvecs[i:,i]=1.0

    dvecs=gs(dvecs,Kover,gramnum)
    KinvH=numpy.matmul(Ki,Bigham)
    db=beta/steps
    eb=numpy.zeros((steps+1,gramnum+1))
    for i in range(steps+1):
        den=numpy.einsum('ik,ij,jk->k',dvecs,Bigham,dvecs)
        eb[i,0]=i*db
        eb[i,1:]=den
        ddot = -numpy.dot(KinvH,dvecs)
        dvecs = dvecs + db*ddot
        dvecs = gs(dvecs,Kover,gramnum)

    filename=filenamer+'_energy.csv'
    filename1=filenamer+"_dvect.csv"   
    in_outputs.write_ham(eb,filename)
    in_outputs.write_ham(dvecs,filename1)
    return eb

