import inputs
import numpy 
import csv 



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



def elec_writer(h1e,eri_full,EXDIR1):

    norb=int(inputs.zombs['norb']*2)
    H1ei = spatospin1(h1e,norb)
    H2ei = spatospin2(eri_full,norb)

    one_elec_setup(norb,H1ei,EXDIR1)
    two_elec_setup(norb,H2ei,EXDIR1)
   
    

   



  