import numpy

def overlap(z1,z2):
    # coefficients can be real or complex
    temp = 1.0
    for iorb in range(len(z1[:,0])):
        tt = z1[iorb,0].conjugate()*z2[iorb,0] \
             + z1[iorb,1].conjugate()*z2[iorb,1]
        temp = temp*tt
    return temp

def overlap_f(z1,z2):
    # coefficients can be real or complex
    temp = 1.0
    for iorb in range(len(z1[:,0])):
        tt = z1[iorb,0].conjugate()*z2[iorb,0] \
             + z1[iorb,1].conjugate()*z2[iorb,1]
        if tt == 0.0:
            return 0.0
        temp = temp*tt
    return temp

def cr(z1,iorb):
    # Creation operator on orbital iorb
    z1[iorb,1] = z1[iorb,0]
    z1[iorb,0] = 0.0
    z1[:iorb,1] *= -1.0
    return z1

def an(z1,iorb):
    # Annihilation operator on orbital iorb
    z1[iorb,0] = z1[iorb,1]
    z1[iorb,1] = 0.0
    z1[:iorb,1] *= -1.0
    return z1

def numf(z1,z2):
    # faster algorithm for application of number operator
    # needs adaptation for complex calculation
    norb = len(z1[:,0])
    cc = numpy.zeros((norb))
    dd = numpy.zeros((norb))
    mult = numpy.zeros((norb))
    multb = numpy.zeros((norb))
    for iorb in range(norb):
        mult[iorb] = z1[iorb,1].conjugate()*z2[iorb,1]
        multb[iorb] = mult[iorb] + z1[iorb,0].conjugate()*z2[iorb,0]
    cc[0] = multb[0]
    dd[-1] = multb[-1]
    for iorb in range(1,norb):
        cc[iorb] = cc[iorb-1]*multb[iorb]
    for iorb in range(norb-2,-1,-1):# I think indices are correct
        dd[iorb] = dd[iorb+1]*multb[iorb]
    temp = mult[0]*dd[1]
    #print('mult',mult)
    #print('cc',cc)
    #print('dd',dd)
    for iorb in range(1,norb-1):
        temp += cc[iorb-1]*mult[iorb]*dd[iorb+1]
        #print(iorb,temp,mult[iorb])
    temp += cc[norb-2]*mult[-1]
    return temp

def num(z1,iorb):
    # Number operator on a specific orbital iorb
    z1[iorb,0] = 0.0
    return z1

def nsq(z1,z2,norb):
    "N squared, used numo to cope with complex"
    temp = 0.0
    for iorb in range(norb):
        zt = numpy.copy(z2)
        zt = num(zt,iorb)
        temp += numf(z1,zt)
    return temp

def isdet(zom,norb):
    """Checks if a given zombie state
    is equal to a single determinant
    and the entries for each spinorbital are either 0 and 1 or 
    0 and -1"""
    for iorb in range(norb):
        if zom[iorb,0] == 0.0:
            if zom[iorb,1] in [1.0,-1.0]:
                pass
            else:
                return False
        elif zom[iorb,0] in [1.0,-1.0]:
            if zom[iorb,1] == 0.0:
                pass
            else:
                return False
        else:
            return False
        #if (zom[iorb,0] == 0.0 and zom[iorb,1] == 1.0) or \
        #   (zom[iorb,0] == 1.0 and zom[iorb,1] == 0.0):
        #    pass
        #else:
        #    return False
    return True

def iszero(z1):
    """Determines if a given zombie state vanishes
    Returns True if state vanishes
    False otherwise"""
    # coefficients can be real or complex
    for iorb in range(len(z1[:,0])):
        tt = z1[iorb,0].conjugate()*z1[iorb,0] \
             + z1[iorb,1].conjugate()*z1[iorb,1]
        if tt == 0.0:
            return True
    return False

def numtodet(i,norb):
    """Turns an integer number 0 <= j <= 2**norb-1
    Into an arrange length norb with orbital occupancy
    by converting the integer into binary"""
    if i >= 2**norb:
        raise ValueError('i too big')
    bini = numpy.zeros((norb),dtype = int)
    it = i
    for j in range(norb):
        bini[j] = it%2
        it -= bini[j]
        it /= 2
    return bini

def dettonum(bini,norb):
    """Turns an integer array of 0s and 1s length norb
    into its corresponding binary number
    The reverse of numtodet"""
    twop = numpy.zeros((norb),dtype = int)
    for i in range(norb):
        twop[i] = 2**i
    return numpy.dot(bini,twop)

def szf(z1,z2,norb):
    """Fast algorithm for application of sz operator
    O(N) steps"""
    if norb%2 != 0:
        raise ValueError('Even norb required')
    cc = numpy.zeros((norb))
    dd = numpy.zeros((norb))
    mult = numpy.zeros((norb))
    multb = numpy.zeros((norb))
    for iorb in range(norb):
        mult[iorb] = z1[iorb,1].conjugate()*z2[iorb,1]
        multb[iorb] = mult[iorb] + z1[iorb,0].conjugate()*z2[iorb,0]
    cc[0] = multb[0]
    dd[-1] = multb[-1]
    for iorb in range(1,norb):
        cc[iorb] = cc[iorb-1]*multb[iorb]
    for iorb in range(norb-2,-1,-1):
        dd[iorb] = dd[iorb+1]*multb[iorb]
    temp = mult[0]*dd[1]
    for iorb in range(1,norb-1):
        temp += cc[iorb-1]*mult[iorb]*dd[iorb+1]*(-1)**iorb
    temp -= cc[norb-2]*mult[-1]
    return 0.5*temp

def sz2f(zs1,zs2,norb):
    """Computing <zs1 | S_z^2 | zs2 >
    O(M^2) not O(M^3)"""
    temp = 0.0
    for iorb in range(norb):
        zs2t = numpy.copy(zs2)
        zs2t = num(zs2t,iorb)
        temp += szf(zs1,zs2t,norb)*(-1)**iorb
    return temp*0.5

def spsmfast(zs1,zs2,norb):
    """Fastest calcualtion of <zs1 |S_+S_- |zs2>"""
    kmax=int(norb/2)
    cc = numpy.zeros((kmax))
    dd = numpy.zeros((kmax))
    ss = numpy.zeros((kmax))
    tt= numpy.zeros((kmax))
    for i in range(kmax):
        a=2*i
        b=(2*i)+1
        cc[i]=(zs1[a,1].conjugate()*zs2[a,1]+zs1[a,0].conjugate()*zs2[a,0])*(zs1[b,1].conjugate()*zs2[b,1]+zs1[b,0].conjugate()*zs2[b,0])
        dd[i]=(zs1[a,0].conjugate()*zs2[a,1])*(zs1[b,1].conjugate()*zs2[b,0])
        ss[i]=(zs1[a,1].conjugate()*zs2[a,0])*(zs1[b,0].conjugate()*zs2[b,1])
        tt[i]=(zs1[a,1].conjugate()*zs2[a,1])*(zs1[b,0].conjugate()*zs2[b,0])
    tot=0.0 
    for i in range(kmax):
        for j in range(i, kmax):
            p1=0
            p2=0
            if(i==0):
                    p1=1
            elif(i==1):
                p1=cc[0]
            else:
                p1=numpy.prod(cc[:i])
            if(j==kmax-1):
                p2=1
            elif(j==kmax-2):
                p2=cc[kmax-1]
            else:
                p2=numpy.prod(cc[j+1:])
            if(j==i):
                tot=tot+(p1*p2*tt[i])
            elif(j==i+1):
                tot=tot+(p1*p2*ss[i]*dd[j])+(p1*p2*ss[j]*dd[i])
            elif(j==i+2):
                tot=tot+(p1*p2*ss[i]*dd[j]*cc[i+1])+(p1*p2*ss[j]*dd[i]*cc[i+1])
            else:
                p3=numpy.prod(cc[i+1:j])
                tot=tot+(p1*p2*ss[i]*dd[j]*p3)+(p1*p2*ss[j]*dd[i]*p3)

    return tot

def Stotfast(zs1,zs2,norb):
    return spsmfast(zs1,zs2,norb) -szf(zs1,zs2,norb)+sz2f(zs1,zs2,norb)

def z_an_z3(zom1,zom2,norb,vec):
    """Finding sum_k <zom1 | b_k | zom2> vec_k for all k in norb
    Dynamically handles real vs complex input data
    Faster than the first two versions"""
    vmult = numpy.multiply(zom1.conjugate(),zom2)
    if type(vmult[0,0].item()) is complex:
        gg = numpy.zeros((norb),dtype=complex)
        hh = numpy.zeros((norb),dtype=complex)
    else:
        gg = numpy.zeros((norb),dtype=float)
        hh = numpy.zeros((norb),dtype=float)        
    gg[0] = vmult[0,0] - vmult[0,1]
    gmax = norb
    for ii in range(1,norb):
        gg[ii] = gg[ii-1]*(vmult[ii,0]-vmult[ii,1])
        if gg[ii] == 0.0:
            gmax = ii
            break
    hh[-1] = vmult[-1,0] + vmult[-1,1]   
    hmin = 0
    for ii in range(norb-2,-1,-1):
        hh[ii] = hh[ii+1]*(vmult[ii,0]+vmult[ii,1])
        if hh[ii] == 0.0:
            hmin = ii
            break
    #an_out = numpy.zeros(norb)
    an = 0.0
    if gmax < hmin:
        return 0.0
    if vec[0] != 0:
        an += zom1[0,0].conjugate()*zom2[0,1]*hh[1]*vec[0]
    for ii in range(1,norb-1,1):
        if vec[ii]!=0.0:
            an += gg[ii-1]*zom1[ii,0].conjugate()*zom2[ii,1]*hh[ii+1]*vec[ii]
    if vec[-1] !=0: 
        an += gg[-2]*zom1[-1,0].conjugate()*zom2[-1,1]*vec[-1]
    return an