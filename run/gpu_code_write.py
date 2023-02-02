import inputs
import numpy 
import csv 

# choice=("j","k")
# alive="alive{0}({1})"
# dead="dead{0}({1})"
# aa=alive.format("1","{0}")+"*"+alive.format("2","{0}")
# dd=dead.format("1","{0}")+"*"+dead.format("2","{0}")
# ad=alive.format("1","{0}")+"*"+dead.format("2","{0}")
# da=dead.format("1","{0}")+"*"+alive.format("2","{0}")
# standard="("+aa+"+"+dd+")"
# negative="("+dd+"-"+aa+")"
# body=(" ",standard,negative)
# cr="("+ad+")"
# cr_neg="(-"+ad+")"
# cr_choice=(" ",cr,cr_neg)
# an="("+da+")"
# an_neg="(-"+da+")"
# an_choice=(" ",an,an_neg)
# cran="("+aa+")"
# cran_neg="(-"+aa+")"
# cran_choice=(" ",cran,cran_neg)
# omp="!$omp task"
# ompend="!$omp end task"
# end="return"
# el%h1ei("+str(i+1)+"


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


def one_elec_setup(norb,H1ei,occupancy_an_cr,occupancy_an_cr_v,EXDIR1):

    h1=""
    h1count=0
    g=numpy.zeros(((norb*norb),3))

    for i in range(norb): # i is annihilation 
        for j in range(norb): #j is creation
            if(H1ei[i,j]!=0.0):
                h1count=h1count+1
                g[h1count,0]=H1ei[i,j]
                g[h1count,1]=int(i+1)
                g[h1count,2]=int(j+1)
                # a_string=""
                # d_string=""
                # negatives=("","","-")
                # choice=("","a2(","d2(")
                # ends=("","{}","{}")               
                # values=("0.0","","")
                # for p in range(norb):
                #     oan=int(occupancy_an_cr[i,j,0,p])
                #     oanv=int(occupancy_an_cr_v[i,j,0,p])
                #     a_string=a_string+"{}".format(negatives[oan])+values[oanv]+"{}".format(choice[oanv])+"{}".format(ends[oanv]).format(str(p+1)+")")
                #     oan=int(occupancy_an_cr[i,j,1,p])
                #     oanv=int(occupancy_an_cr_v[i,j,1,p])
                #     d_string=d_string+"{}".format(negatives[oan])+values[oanv]+"{}".format(choice[oanv])+"{}".format(ends[oanv]).format(str(p+1)+")")
        
                #     if(p!=norb-1):
                #         a_string=a_string+","
                #         d_string=d_string+","
                #     else:
                #         a_string=a_string+"]"
                #         d_string=d_string+"]"
                # h1=h1+"        a2_op(:,"+str(h1count)+") = [real(kind=8) ::"+a_string+"\n        d2_op(:,"+str(h1count)+") = [real(kind=8) ::"+d_string+"\n"

    g[0,0]=h1count
    with open(EXDIR1+"/integrals/h1e.csv",'w', newline='')as csvfile:
            spamwriter=csv.writer(csvfile, delimiter=',')
            spamwriter.writerows(g[0:h1count+1,:])
           
        
    
    return h1


def two_elec_setup(norb,H2ei,occupancy_2cr_2an,occupancy_2cr_2an_v,EXDIR1):

    h2=""
    h2count=0
    g=numpy.zeros(((norb*norb*norb*norb),5))

    for i in range(norb):
        for j in range(norb):
            if(i!=j):
                if(i%2==0):
                    jspin=0
                else: 
                    jspin=1
                for k in range(jspin,norb,2):
                    for l in range(norb):
                        if(k!=l):
                            if(H2ei[i,j,k,l]!=0.0):
                                h2count=h2count+1  
                                g[h2count,0]=H2ei[i,j,k,l]
                                g[h2count,1]=i+1
                                g[h2count,2]=j+1
                                g[h2count,4]=k+1
                                g[h2count,3]=l+1
                                # a_string=""
                                # d_string=""
                                # negatives=("","","-")
                                # choice=("","a2(","d2(")
                                # ends=("","{}","{}")               
                                # values=("0.0","","")
                                # for p in range(norb):
                                #     oan=int(occupancy_2cr_2an[i,j,k,l,0,p])
                                #     oanv=int(occupancy_2cr_2an_v[i,j,k,l,0,p])
                                #     a_string=a_string+"{}".format(negatives[oan])+values[oanv]+"{}".format(choice[oanv])+"{}".format(ends[oanv]).format(str(p+1)+")")
                                #     oan=int(occupancy_2cr_2an[i,j,k,l,1,p])
                                #     oanv=int(occupancy_2cr_2an_v[i,j,k,l,1,p])
                                #     d_string=d_string+"{}".format(negatives[oan])+values[oanv]+"{}".format(choice[oanv])+"{}".format(ends[oanv]).format(str(p+1)+")")
                        
                                #     if(p!=norb-1):
                                #         a_string=a_string+","
                                #         d_string=d_string+","
                                #     else:
                                #         a_string=a_string+"]"
                                #         d_string=d_string+"]"
                                # h2=h2+"        a2_op(:,"+str(h2count)+") = [real(kind=8) ::"+a_string+"\n        d2_op(:,"+str(h2count)+") = [real(kind=8) ::"+d_string+"\n"


    g[0,0]=h2count
    # with open(EXDIR1+"h2e.csv",'w', newline='')as csvfile:
    with open(EXDIR1+"/integrals/h2e.csv",'w', newline='')as csvfile:
        spamwriter=csv.writer(csvfile, delimiter=',')
        spamwriter.writerows(g[0:h2count+1,:])
        # spamwriter.writerow([h2count,' '])
        # for j in range(h2count):
        #     spamwriter.writerow([g[j],' '])


    return h2





def ham_writer(h1e,eri_full,EXDIR1):

    norb=inputs.zombs['norb']*2
    
    occupancy_an_cr=numpy.ones((norb,norb,2,norb))
    occupancy_an=numpy.ones((norb,2,norb))
    occupancy_2an=numpy.ones((norb,norb,2,norb))
    occupancy_2cr_2an=numpy.ones((norb,norb,norb,norb,2,norb))
    occupancy_cr_2an=numpy.ones((norb,norb,norb,2,norb))
    occupancy_an_cr_v=numpy.ones((norb,norb,2,norb))
    occupancy_an_v=numpy.ones((norb,2,norb))
    occupancy_2an_v=numpy.ones((norb,norb,2,norb))
    occupancy_2cr_2an_v=numpy.ones((norb,norb,norb,norb,2,norb))
    occupancy_cr_2an_v=numpy.ones((norb,norb,norb,2,norb))

    occupancy_an_cr_v[:,:,1,:]=2
    occupancy_an_v[:,1,:]=2
    occupancy_2an_v[:,:,1,:]=2
    occupancy_2cr_2an_v[:,:,:,:,1,:]=2
    occupancy_cr_2an_v[:,:,:,1,:]=2

    # for i in range(norb):
    #     occupancy_an[i,0,i]=0
    #     occupancy_an_v[i,0,i]=0
    #     for p in range(i-1,-1,-1):
    #         occupancy_an[i,0,p]=-1
    #     for j in range(norb):
    #         occupancy_2an[i,j,:,:]=occupancy_an[i,:,:]
    #         occupancy_2an_v[i,j,:,:]=occupancy_an_v[i,:,:]
    #         occupancy_2an[i,j,1,j]=occupancy_2an[i,j,0,j]
    #         occupancy_2an_v[i,j,1,j]=occupancy_2an_v[i,j,0,j]
    #         occupancy_2an[i,j,0,j]=0
    #         occupancy_2an_v[i,j,0,j]=0
    #         occupancy_an_cr[i,j,:,:]=occupancy_an[i,:,:]
    #         occupancy_an_cr[i,j,0,j]=occupancy_an_cr[i,j,1,j]
    #         occupancy_an_cr[i,j,1,j]=0
    #         occupancy_an_cr_v[i,j,:,:]=occupancy_an_v[i,:,:]
    #         occupancy_an_cr_v[i,j,0,j]=occupancy_an_cr_v[i,j,1,j]
    #         occupancy_an_cr_v[i,j,1,j]=0
    #         for p in range(j-1,-1,-1):
    #             occupancy_2an[i,j,0,p]=occupancy_2an[i,j,0,p]*(-1)
    #             occupancy_an_cr[i,j,0,p]=occupancy_an_cr[i,j,0,p]*(-1)
    #         for k in range(norb):
    #             occupancy_cr_2an[i,j,k,:,:]=occupancy_2an[i,j,:,:]
    #             occupancy_cr_2an[i,j,k,0,k]=occupancy_cr_2an[i,j,k,1,k]
    #             occupancy_cr_2an[i,j,k,1,k]=0
    #             occupancy_cr_2an_v[i,j,k,:,:]=occupancy_2an_v[i,j,:,:]
    #             occupancy_cr_2an_v[i,j,k,0,k]=occupancy_cr_2an_v[i,j,k,1,k]
    #             occupancy_cr_2an_v[i,j,k,1,k]=0
    #             for p in range(k-1,-1,-1):
    #                 occupancy_cr_2an[i,j,k,0,p]=occupancy_cr_2an[i,j,k,0,p]*(-1)
    #             for l in range(norb):
    #                 occupancy_2cr_2an[i,j,k,l,:,:]=occupancy_cr_2an[i,j,k,:,:]
    #                 occupancy_2cr_2an[i,j,k,l,0,l]=occupancy_2cr_2an[i,j,k,l,1,l]
    #                 occupancy_2cr_2an[i,j,k,l,1,l]=0
    #                 occupancy_2cr_2an_v[i,j,k,l,:,:]=occupancy_cr_2an_v[i,j,k,:,:]
    #                 occupancy_2cr_2an_v[i,j,k,l,0,l]=occupancy_2cr_2an_v[i,j,k,l,1,l]
    #                 occupancy_2cr_2an_v[i,j,k,l,1,l]=0
    #                 for l in range(j-1,-1,-1):
    #                     occupancy_2cr_2an[i,j,k,l,0,p]=occupancy_2cr_2an[i,j,k,l,0,p]*(-1)

  
    # occupancy_an_cr[occupancy_an_cr<0]=2
    # occupancy_2cr_2an[occupancy_2cr_2an<0]=2
   

    H1ei = spatospin1(h1e,norb)
    H2ei = spatospin2(eri_full,norb)

    h1=one_elec_setup(norb,H1ei,occupancy_an_cr,occupancy_an_cr_v,EXDIR1)
    h2=two_elec_setup(norb,H2ei,occupancy_2cr_2an,occupancy_2cr_2an_v,EXDIR1)
        
    

    with open('ham_2.f90','w',encoding="utf-8") as f:
        f.write("MODULE ham_2 \n")
        f.write(" \n")
        f.write("   use globvars \n")
        f.write(" \n")
        f.write("   contains \n")
        f.write("\n")
        f.write("   subroutine cr_an(a2,d2,a2_op,d2_op)\n\n")
        f.write("       implicit none\n\n")
        f.write("       real(kind=8), dimension(:),intent(in)::a2,d2\n")
        f.write("       real(kind=8),dimension(:,:),intent(inout)::a2_op,d2_op\n")
        f.write("       integer::j\n\n")
        f.write(" "+h1+"\n\n")
        f.write("       return\n")
        f.write("   end subroutine cr_an\n")
        f.write("\n\n")
        f.write("   subroutine cr2_an2(a2,d2,a2_op,d2_op)\n\n")
        f.write("       implicit none\n\n")
        f.write("       real(kind=8), dimension(:),intent(in)::a2,d2\n")
        f.write("       real(kind=8),dimension(:,:),intent(inout)::a2_op,d2_op\n")
        f.write("       integer::j\n\n")
        f.write(" "+h2+"\n\n")
        f.write("       return\n")
        f.write("   end subroutine cr2_an2\n")



        f.write(" \n")
        f.write("END MODULE ham_2 \n")




# def one_elec_setup(norb,H1ei,occupancy_an_cr,EXDIR1):

#     h1=""
#     h1count=0
#     g=numpy.zeros(((norb*norb)))

#     for i in range(norb): # i is annihilation 
#         for j in range(norb): #j is creation
#             if(H1ei[i,j]!=0.0):
#                 g[h1count]=H1ei[i,j]
#                 h1count=h1count+1
#                 newline=0
#                 if(h1count==1):
#                     srt="if(j.eq.1)then"
#                 else:
#                     srt="else if(j.eq."+str(h1count)+") then"
#                 if(i!=j):
#                     h1=h1+""+srt+" \n               "+"h1e=el*"+an.format(str(i+1),choice[0],choice[1])+"*{0}".format(cr_choice[int(occupancy_an_cr[i,j,0,j])]).format(str(j+1),choice[0],choice[1])+"&\n               " 
#                     orb=2
#                     for p in range(norb):
                        
#                         if((p==i)or(p==j)):
#                             continue 
#                         else:
#                             orb=orb+1
#                             h1=h1+ "*{0}".format(body[int(occupancy_an_cr[i,j,0,p])]).format(str(p+1),choice[0],choice[1])
#                             newline=newline+1
#                             if(orb!=norb):
#                                 if(newline==2):
#                                     h1=h1+"&\n               "
#                                     newline=0
#                             else:
#                                 h1=h1+"\n               "+end+"\n       "
#                 else:
#                     h1=h1+""+srt+" \n               "+"h1e=el*"+cran.format(str(j+1),choice[0],choice[1])
#                     newline=1
#                     orb=1
#                     for p in range(norb):
#                         if((p==i)):
#                             continue 
#                         else:
#                             orb=orb+1
#                             h1=h1+"*"+standard.format(str(p+1),choice[0],choice[1]) 
#                             newline=newline+1
#                             if(orb!=norb):
#                                 if(newline==2):
#                                     h1=h1+"&\n               "
#                                     newline=0
#                             else:
#                                 h1=h1+"\n               "+end+"\n       "
#             else:
#                 continue
    
#     with open(EXDIR1+"/integrals/h1e.csv",'w', newline='')as csvfile:
#             spamwriter=csv.writer(csvfile, delimiter=',')
#             spamwriter.writerow([h1count,' '])
#             for j in range(h1count):
#                 spamwriter.writerow([g[j],' '])
        
    
#     return h1



# def two_elec_setup(norb,H2ei,occupancy_2cr_2an,EXDIR1):
   
#     h2count=0
#     h2=" "
#     g=numpy.zeros((norb*norb*norb*norb))
#     for i in range(norb):
#         for j in range(norb):
#             if(i!=j):
#                 for k in range(norb):
#                     for l in range(norb):
#                         if(k!=l):
#                             if(H2ei[i,j,k,l]!=0.0):
#                                 g[h2count]=H2ei[i,j,k,l]
#                                 h2count=h2count+1
#                                 newline=0
#                                 orb=0
#                                 if(h2count==1):
#                                     srt="if(j.eq.1)then"
#                                 else:
#                                     srt="else if(j.eq."+str(h2count)+") then"
#                                 if(i!=k): #i≠k
#                                     if(i!=l): #i≠k and i≠l
#                                         if(j!=k): #i≠k and i≠l and j≠k
#                                             if(j!=l): #i≠k and i≠l and j≠k j≠l i,j,k,l
#                                                 h2=h2+""+srt+" \n               "+"h2e=el*"+"{0}".format(an_choice[int(occupancy_2cr_2an[i,j,k,l,1,i])]).format(str(i+1),choice[0],choice[1])+\
#                                                     "*{0}".format(an_choice[int(occupancy_2cr_2an[i,j,k,l,1,j])]).format(str(j+1),choice[0],choice[1])+"&\n               "
#                                                 h2=h2+"*{0}".format(cr_choice[int(occupancy_2cr_2an[i,j,k,l,0,k])]).format(str(k+1),choice[0],choice[1])+"*{0}".format(cr_choice[int(occupancy_2cr_2an[i,j,k,l,0,l])]).format(str(l+1),choice[0],choice[1])+"&\n               "
#                                                 orb=4
#                                                 newline=0
#                                                 for p in range(norb):
#                                                     if((p==i)or(p==j)or(p==k)or(p==l)):
#                                                         continue 
#                                                     else:
#                                                         orb=orb+1
#                                                         h2=h2+ "*{0}".format(body[int(occupancy_2cr_2an[i,j,k,l,0,p])]).format(str(p+1),choice[0],choice[1])
#                                                         newline=newline+1
#                                                         if(orb!=norb):
#                                                             if(newline==2):
#                                                                 h2=h2+"&\n               "
#                                                                 newline=0
#                                                         else:
#                                                             h2=h2+"\n               "+end+"\n       "

#                                             else:#i≠k and i≠l and j≠k j=l i,jcran,k
#                                                 h2=h2+""+srt+" \n               "+"h2e=el*"+"{0}".format(an_choice[int(occupancy_2cr_2an[i,j,k,l,1,i])]).format(str(i+1),choice[0],choice[1])+\
#                                                     "*{0}".format(cran_choice[int(occupancy_2cr_2an[i,j,k,l,0,j])]).format(str(j+1),choice[0],choice[1])+"&\n               "
#                                                 h2=h2+"*{0}".format(cr_choice[int(occupancy_2cr_2an[i,j,k,l,0,k])]).format(str(k+1),choice[0],choice[1])
#                                                 orb=3
#                                                 newline=1
#                                                 for p in range(norb):
                                                    
#                                                     if((p==i)or(p==j)or(p==k)):
#                                                         continue 
#                                                     else:
#                                                         orb=orb+1
#                                                         h2=h2+ "*{0}".format(body[int(occupancy_2cr_2an[i,j,k,l,0,p])]).format(str(p+1),choice[0],choice[1])
#                                                         newline=newline+1
#                                                         if(orb!=norb):
#                                                             if(newline==2):
#                                                                 h2=h2+"&\n               "
#                                                                 newline=0
#                                                         else:
#                                                             h2=h2+"\n               "+end+"\n       "
#                                         else: #i≠l and i≠k and j=k j≠l i, jcran, l
#                                             h2=h2+""+srt+" \n               "+"h2e=el*"+"{0}".format(an_choice[int(occupancy_2cr_2an[i,j,k,l,1,i])]).format(str(i+1),choice[0],choice[1])+\
#                                                 "*{0}".format(cran_choice[int(occupancy_2cr_2an[i,j,k,l,0,j])]).format(str(j+1),choice[0],choice[1])+"&\n               "
#                                             h2=h2+"*{0}".format(cr_choice[int(occupancy_2cr_2an[i,j,k,l,0,l])]).format(str(l+1),choice[0],choice[1])
#                                             newline=1
#                                             orb=3
#                                             for p in range(norb):
                                                
#                                                 if((p==i)or(p==j)or(p==l)):
#                                                     continue 
#                                                 else:
#                                                     orb=orb+1
#                                                     h2=h2+ "*{0}".format(body[int(occupancy_2cr_2an[i,j,k,l,0,p])]).format(str(p+1),choice[0],choice[1])
#                                                     newline=newline+1
#                                                     if(orb!=norb):
#                                                         if(newline==2):
#                                                             h2=h2+"&\n               "
#                                                             newline=0
#                                                     else:
#                                                         h2=h2+"\n               "+end+"\n       "
#                                     else:  #i≠k i=l
#                                         if(j!=k): #i≠k and i=l and j≠k j≠l icran, j,k
#                                             h2=h2+""+srt+" \n               "+"h2e=el*"+"{0}".format(an_choice[int(occupancy_2cr_2an[i,j,k,l,1,j])]).format(str(j+1),choice[0],choice[1])+\
#                                                 "*{0}".format(cran_choice[int(occupancy_2cr_2an[i,j,k,l,0,i])]).format(str(i+1),choice[0],choice[1])+"&\n               "
#                                             h2=h2+"*{0}".format(cr_choice[int(occupancy_2cr_2an[i,j,k,l,0,k])]).format(str(k+1),choice[0],choice[1])
#                                             newline=1
#                                             orb=3
#                                             for p in range(norb):
#                                                 if((p==i)or(p==j)or(p==k)):
#                                                     continue 
#                                                 else:
#                                                     orb=orb+1
#                                                     h2=h2+ "*{0}".format(body[int(occupancy_2cr_2an[i,j,k,l,0,p])]).format(str(p+1),choice[0],choice[1])
#                                                     newline=newline+1
#                                                     if(orb!=norb):
#                                                         if(newline==2):
#                                                             h2=h2+"&\n               "
#                                                             newline=0
#                                                     else:
#                                                         h2=h2+"\n               "+end+"\n       "
#                                         else: #i≠k and i=l and j=k j≠l icran,jcran
#                                             h2=h2+""+srt+" \n               "+"h2e=el*"+\
#                                                 "*{0}".format(cran_choice[int(occupancy_2cr_2an[i,j,k,l,0,i])]).format(str(i+1),choice[0],choice[1])+"*{0}".format(cran_choice[int(occupancy_2cr_2an[i,j,k,l,0,j])]).format(str(j+1),choice[0],choice[1])+"&\n               "
#                                             newline=0
#                                             orb=2
#                                             for p in range(norb):
#                                                 if((p==i)or(p==j)):
#                                                     continue 
#                                                 else:
#                                                     orb=orb+1
#                                                     h2=h2+ "*{0}".format(body[int(occupancy_2cr_2an[i,j,k,l,0,p])]).format(str(p+1),choice[0],choice[1])
#                                                     newline=newline+1
#                                                     if(orb!=norb):
#                                                         if(newline==2):
#                                                             h2=h2+"&\n               "
#                                                             newline=0
#                                                     else:
#                                                         h2=h2+"\n               "+end+"\n       "
#                                 else: #i=k i≠l
#                                     if(j!=l): #i=k and i≠l j≠l and j≠k icran,j,l
#                                         h2=h2+""+srt+" \n               "+"h2e=el*"+"{0}".format(an_choice[int(occupancy_2cr_2an[i,j,k,l,1,j])]).format(str(j+1),choice[0],choice[1])+\
#                                             "*{0}".format(cran_choice[int(occupancy_2cr_2an[i,j,k,l,0,i])]).format(str(i+1),choice[0],choice[1])+"&\n               "+"*{0}".format(cr_choice[int(occupancy_2cr_2an[i,j,k,l,0,l])]).format(str(l+1),choice[0],choice[1])
#                                         newline=1
#                                         orb=3
#                                         for p in range(norb):
#                                             if((p==i)or(p==j)or(p==l)):
#                                                 continue 
#                                             else:
#                                                 orb=orb+1
#                                                 h2=h2+ "*{0}".format(body[int(occupancy_2cr_2an[i,j,k,l,0,p])]).format(str(p+1),choice[0],choice[1])
#                                                 newline=newline+1
#                                                 if(orb!=norb):
#                                                     if(newline==2):
#                                                         h2=h2+"&\n               "
#                                                         newline=0
#                                                 else:
#                                                     h2=h2+"\n               "+end+"\n       "
                                    

#                                     else: #i=k and i≠l j=l and j≠k icran,jcran
#                                         h2=h2+""+srt+" \n               "+"h2e=el*"+\
#                                             "*{0}".format(cran_choice[int(occupancy_2cr_2an[i,j,k,l,0,i])]).format(str(i+1),choice[0],choice[1])+"*{0}".format(cran_choice[int(occupancy_2cr_2an[i,j,k,l,0,j])]).format(str(j+1),choice[0],choice[1])+"&\n               "
#                                         orb=2
#                                         newline=0
#                                         for p in range(norb):
#                                             if((p==i)or(p==j)):
#                                                 continue 
#                                             else:
#                                                 orb=orb+1
#                                                 h2=h2+ "*{0}".format(body[int(occupancy_2cr_2an[i,j,k,l,0,p])]).format(str(p+1),choice[0],choice[1])
#                                                 newline=newline+1
#                                                 if(orb!=norb):
#                                                     if(newline==2):
#                                                         h2=h2+"&\n               "
#                                                         newline=0
#                                                 else:
#                                                     h2=h2+"\n               "+end+"\n       "
#                             else:
#                                 continue
#                         else:
#                             continue
#             else:
#                 continue

#     with open(EXDIR1+"/integrals/h2e.csv",'w', newline='')as csvfile:
#             spamwriter=csv.writer(csvfile, delimiter=',')
#             spamwriter.writerow([h2count,' '])
#             for j in range(h2count):
#                 spamwriter.writerow([g[j],' '])
    
    
#     return h2