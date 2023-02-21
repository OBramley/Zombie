MODULE gradient_descent

    use globvars
    use alarrays
    use ham
    use grad_d
    use imgtp
    use outputs
    use infnan_mod
  
    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators

    ! Subroutine that controls and calcualtes all of the hamiltonian variables 
    subroutine he_full_row(haml,zstore,elecs,size,an_cr,an2_cr2,diff_state)

        implicit none 

        type(hamiltonian), intent(inout)::haml
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::size,diff_state
        integer, allocatable,dimension(:)::IPIV1
        real(kind=8),allocatable,dimension(:)::WORK1
        
        integer::ierr


        if (errorflag .ne. 0) return
        ierr=0
        
        haml%ovrlp(:,diff_state)=ovrlp_column(zstore(diff_state)%val,zstore,diff_state)
        haml%ovrlp(diff_state,:)=haml%ovrlp(:,diff_state)
       
        haml%hjk(:,diff_state)=haml%ovrlp(:,diff_state)*elecs%hnuc 
        call haml_column(haml%hjk(:,diff_state),zstore(diff_state)%val,zstore,an_cr,an2_cr2,elecs,1)
        haml%hjk(diff_state,:)=haml%hjk(:,diff_state)
        haml%inv=haml%ovrlp
        allocate(WORK1(size),IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV or WORK1 vector allocation . ierr had value ", ierr
            errorflag=1
        end if 

        Call dgetrf(size, size, haml%inv, size, IPIV1, ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if
        if (ierr==0) call dgetri(size,haml%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if

        deallocate(WORK1,IPIV1)

        call DGEMM("N","N",size,size,size,1.d0,haml%inv,size,haml%hjk,size,0.d0,haml%kinvh,size)


    end subroutine he_full_row


    subroutine grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvec,grad_fin,en,orb)

        implicit none 

        type(zombiest),dimension(:),intent(in)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvec
        type(hamiltonian),intent(inout)::haml
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::pick,orb
        type(energy),intent(inout)::en

        if (errorflag .ne. 0) return
 
        call gradient_zs(haml,zstore,elect,an_cr,an2_cr2,pick,orb)
        
        call imgtime_prop(dvec,en,haml,pick,0)
       
        call final_grad(dvec(1),haml,grad_fin,pick,orb)
        
        en%erg=0
        en%t=0
      
        grad_fin%grad_avlb(pick)=1
    
        return

    end subroutine grad_calc

    subroutine orbital_gd(zstore,grad_fin,elect,dvecs,an_cr,an2_cr2,en,haml,temp_ham,temp_zom,&
        epoc_cnt,alphain,b,picker,maxloop,rjct_cnt,epoc_max) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore,temp_zom
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml,temp_ham
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(inout)::epoc_cnt,rjct_cnt
        integer,intent(in)::maxloop,epoc_max
        real(kind=8),intent(in)::b,alphain
        integer,dimension(ndet-1),intent(inout)::picker
        integer::next,acpt_cnt,pick,pickorb,acpt_cnt2,bloops,lralt,lralt_temp,loop_max
        real(kind=8)::t,fxtdk,l2_rglrstn,gradient_norm
        integer,dimension(norb)::pickerorb
        integer,dimension(norb)::chng_trk2
        integer::j,k,l,n
        
    
        lralt=0    ! power alpha is raised to 
        chng_trk2=0 !stores which orbitals in the ZS have changed 
        l2_rglrstn=0.1! !L2 regularisation lambda paramter
        acpt_cnt2=0
        loop_max=20
        bloops=0
       
        pickerorb=scramble_norb(norb)
        do while(rjct_cnt.lt.(ndet-1)*2+4)
            bloops=bloops+1
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent steps &
                        &   | Orbitals altered '
            do j=1,(ndet-1)
               
                pick=picker(j)
                chng_trk2=0
                acpt_cnt=0
                grad_fin%current_erg=grad_fin%prev_erg
                do n=1,norb
                    pickorb=pickerorb(n)
                    lralt_temp=lralt
                  
                    do while(lralt_temp.lt.loop_max)

                        t=b*(alphain**lralt_temp)
                    
                        ! Setup temporary zombie state
                        
                        temp_zom=zstore
                        temp_zom(pick)%phi(pickorb)=zstore(pick)%phi(pickorb)-(t*(grad_fin%vars(pick,pickorb)))
                        temp_zom(pick)%phi(pickorb)=temp_zom(pick)%phi(pickorb)+&
                                            l2_rglrstn*((grad_fin%vars(pick,pickorb))*grad_fin%vars(pick,pickorb))
                       
                        temp_zom(pick)%sin(pickorb)=sin(temp_zom(pick)%phi(pickorb))
                        temp_zom(pick)%cos(pickorb)=cos(temp_zom(pick)%phi(pickorb))
                        temp_zom(j)%val(1:norb)=temp_zom(j)%sin
                        temp_zom(j)%val(norb+1:)=temp_zom(j)%cos
                       
                        temp_ham=haml
                        call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
        
                        ! Imaginary time propagation for back tracing
                        en%erg=0
                        en%t=0
                        call imgtime_prop(dvecs,en,temp_ham,0,0)
                        fxtdk=en%erg(1,timesteps+1)
                        ! print*,fxtdk
                        ! Check if energy is lower and accept or reject
                        if((fxtdk.lt.grad_fin%prev_erg).and.&
                        (abs((grad_fin%prev_erg)-fxtdk).lt.0.05).and.&
                        (abs(grad_fin%prev_erg-fxtdk)/abs(grad_fin%prev_erg).lt.0.0001))then
                            acpt_cnt=acpt_cnt+1
                            acpt_cnt2=acpt_cnt2+1
                            zstore(pick)=temp_zom(pick)
                            haml=temp_ham
                            zstore(pick)=temp_zom(pick)
                            chng_trk2(acpt_cnt)=pickorb
                            haml=temp_ham
                            grad_fin%grad_avlb=0
                            grad_fin%vars=0.0
                            haml%diff_hjk=0
                            haml%diff_invh=0
                            haml%diff_ovrlp=0
                            haml%diff_ov_dov=0
                            haml%diff_in_dhjk=0
                            grad_fin%prev_mmntm(pick,pickorb)=zstore(pick)%phi(pickorb)
                            grad_fin%prev_erg=fxtdk
                            Exit      
                        else 
                            lralt_temp=lralt_temp+1
                        end if
                        
                    end do
                    
                    if((n.ne.norb))then
                        call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,pickerorb(n+1))
                    end if
                end do

                if(acpt_cnt.gt.0)then
                    call epoc_writer(grad_fin%prev_erg,epoc_cnt,pick,t,0)
                    call zombiewriter(zstore(pick),pick,0)
                    zstore(pick)%update_num=zstore(pick)%update_num+1
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,*(i0:','))") '       ', pick,'              ', &
                    grad_fin%current_erg,'               ',grad_fin%prev_erg,'            ',chng_trk2(1:acpt_cnt) 
                    
                else
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,i0)") '       ', pick,'              ', &
                    grad_fin%current_erg,'               ',grad_fin%prev_erg,'            ',0
                end if 
          
                
                if(j.eq.(ndet-1))then 
                    epoc_cnt=epoc_cnt+1
                    picker=scramble(ndet-1)
                    next=picker(1)
                    !Every 100 epoc brings phi values back within the normal 0-2pi range
                    if(modulo(epoc_cnt,100).eq.0)then 
                        do k=2,ndet 
                            do l=1,norb 
                                if((zstore(k)%phi(l).gt.2*pirl).or.(zstore(k)%phi(l).lt.0))then 
                                    zstore(k)%phi(l)=asin(real(zstore(k)%sin(l)))
                                    if((zstore(k)%phi(l).gt.2*pirl))then
                                        zstore(k)%phi(l)=zstore(k)%phi(l)-2*pirl
                                    else if((zstore(k)%phi(l).lt.0))then
                                        zstore(k)%phi(l)=zstore(k)%phi(l)+2*pirl
                                    end if
                                end if
                            end do
                        end do
                    end if
                else 
                    pickerorb=scramble_norb(norb)
                    next=picker(j+1)
                    if(grad_fin%grad_avlb(next).eq.0)then
                        call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,pickerorb(1))
                    end if 
                end if
                
            end do
    
            write(6,"(a,i0,a,f21.16)") "Energy after epoch no. ",epoc_cnt,": ",grad_fin%prev_erg
            if(acpt_cnt2.gt.0)then
                rjct_cnt=0
            else 
                rjct_cnt=rjct_cnt+1
            end if 

            if(bloops.ge.maxloop)then 
                call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,0)
                exit 
            else 
                pickerorb=scramble_norb(norb)
                call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,pickerorb(1))
            end if
            if(epoc_cnt.gt.epoc_max)then 
                exit 
            end if
        end do

        return 

    end subroutine orbital_gd


    subroutine full_zs_gd(zstore,grad_fin,elect,dvecs,an_cr,an2_cr2,en,&
        haml,temp_ham,temp_zom,epoc_cnt,alphain,b,picker,epoc_max) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore,temp_zom
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml,temp_ham
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(inout)::epoc_cnt,epoc_max
        real(kind=8),intent(in)::b,alphain
        integer,dimension(ndet-1),intent(inout)::picker
        integer::lralt,rjct_cnt,next,acpt_cnt,pick,lralt_temp,mmntmflg,loop_max,d_cnt
        real(kind=8)::newb,t,fxtdk,l2_rglrstn,alpha,mmnmtb,btl,d_avrg_chk,d_avrg,d_chck_num
        real(kind=8),dimension(norb)::gradient_norm
        real(kind=8),dimension(ndet)::mmntm,mmntma
        integer,dimension(ndet-1)::rsrtpass
        DOUBLE PRECISION, external::ZBQLU01
        integer::j,k,l

    
        if (errorflag .ne. 0) return

        if(rstrtflg.eq.'y')then
            rsrtpass=1
            rstrtflg='n'
        else
            rsrtpass=0
        end if 

        alpha=alphain  ! learning rate reduction
        lralt=0    ! power alpha is raised to  
        newb=b
        rjct_cnt=0 !tracks how many rejections 
        t=newb*(alpha**lralt) !learning rate
        acpt_cnt=0  !counts how many ZS have been changed
        l2_rglrstn=0.1! !L2 regularisation lambda paramter
        mmntm=0!0.4
        mmntma=1
        mmntmflg=0
        loop_max=30
        d_cnt=0 
        btl=0.000001
        d_avrg_chk=9.9d-8
        d_chck_num=5
        
        do while(rjct_cnt.lt.(ndet-1)*2)
           
            acpt_cnt=0
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent step &
            &   | Learning rate | Acceptance count | Rejection count   |     Difference'
            
            do j=1,(ndet-1)
                
                pick=picker(j)
                lralt_temp=lralt
                
                gradient_norm=((grad_fin%vars(pick,:))*grad_fin%vars(pick,:))
                do while(lralt_temp.lt.(loop_max))
                
                    t=newb*(alpha**lralt_temp)
                    temp_zom=zstore
                   
                    temp_zom(pick)%phi(:)=zstore(pick)%phi(:)-(t*(grad_fin%vars(pick,:)))
                    ! temp_zom(pick)%phi(:)=temp_zom(pick)%phi(:)+l2_rglrstn*gradient_norm
                    ! if((mmntmflg.ne.0).and.(mmntm(pick).ne.0))then
                    !     mmnmtb=(t*mmntm(pick))/mmntma(pick)
                    !     temp_zom(pick)%phi(:)=temp_zom(pick)%phi(:)+mmnmtb*(zstore(pick)%phi(:)-grad_fin%prev_mmntm(pick,:))
                    ! end if
                   
                   
                    temp_zom(pick)%sin=sin(temp_zom(pick)%phi)
                    temp_zom(pick)%cos=cos(temp_zom(pick)%phi)
                    temp_zom(pick)%val(1:)=temp_zom(pick)%sin
                    temp_zom(pick)%val(norb+1:)=temp_zom(pick)%cos

                    temp_ham=haml
                   
                    call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
        
                    ! Imaginary time propagation for back tracing
                    en%erg=0
                    en%t=0
                    call imgtime_prop(dvecs,en,temp_ham,0,0)
                
                    fxtdk=en%erg(1,timesteps+1)
                    ! print*,t*sum(gradient_norm)i
                    ! if(((fxtdk.lt.grad_fin%prev_erg)).and.&
                    if(((fxtdk.lt.grad_fin%prev_erg+((btl*t*sum(gradient_norm))))).and.&
                    (abs(grad_fin%prev_erg-fxtdk).lt.0.0025).and.&
                    (abs(grad_fin%prev_erg-fxtdk)/abs(grad_fin%prev_erg).lt.0.001))then!.and.&
                    !(grad_fin%prev_erg-fxtdk.lt.0.0002))then
                        d_avrg=d_avrg+(grad_fin%prev_erg-fxtdk)
                        acpt_cnt=acpt_cnt+1
                        zstore=temp_zom
                        haml=temp_ham
                        zstore(pick)=temp_zom(pick)
                        zstore(pick)%update_num=zstore(pick)%update_num+1
                        call zombiewriter(zstore(pick),pick,rsrtpass(pick))
                        rsrtpass(pick)=0
                        haml=temp_ham
                        rjct_cnt=0
                        grad_fin%grad_avlb=0
                        grad_fin%vars=0.0
                        haml%diff_hjk=0
                        haml%diff_invh=0
                        haml%diff_ovrlp=0
                        haml%diff_ov_dov=0
                        haml%diff_in_dhjk=0
                        grad_fin%prev_mmntm(pick,:)=zstore(pick)%phi(:)
                        mmntma(pick)=t
                        if((mmntm(pick).ne.0).and.(mmntmflg.ne.0))then
                            mmntm(pick)=mmnmtb
                        else 
                            mmntm(pick)=0.4
                        end if
                       
                        write(6,"(a,i0,a,f21.16,a,f21.16,a,f12.10,a,i0,a,i0,a,f21.16)") '       ', pick,'              ', &
                    grad_fin%prev_erg,'               ',fxtdk,'            ',t,'          ',acpt_cnt,'                 ',rjct_cnt, &
            '          ',(fxtdk-grad_fin%prev_erg)
                        grad_fin%prev_erg=fxtdk
                        call epoc_writer(grad_fin%prev_erg,epoc_cnt,pick,t,0)
                        Exit 
                    else 
                        lralt_temp=lralt_temp+1
                    end if
                    
        
                end do

                if(lralt_temp.ge.loop_max)then
                    rjct_cnt=rjct_cnt+1
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,f12.10,a,i0,a,i0)") '       ', pick,'              ', &
                grad_fin%prev_erg,'               ',fxtdk,'            ',0.0,'          ',acpt_cnt,'                 ',rjct_cnt
                end if
            
                
                if(j.eq.(ndet-1))then 
                    
                    if(mmntmflg.eq.0)then
                        mmntmflg=1
                    end if
                    if(acpt_cnt.gt.0)then 
                        picker=scramble(ndet-1)
                    else  
                        mmntmflg=0
                    end if
                    next=picker(1)
    
                    !Every 100 epoc brings phi values back within the normal 0-2pi range
                    if(modulo(epoc_cnt,100).eq.0)then 
                        do k=2,ndet 
                            do l=1,norb 
                                if((zstore(k)%phi(l).gt.2*pirl).or.(zstore(k)%phi(l).lt.0))then 
                                    zstore(k)%phi(l)=asin(real(zstore(k)%sin(l)))
                                    if((zstore(k)%phi(l).gt.2*pirl))then
                                        zstore(k)%phi(l)=zstore(k)%phi(l)-2*pirl
                                    else if((zstore(k)%phi(l).lt.0))then
                                        zstore(k)%phi(l)=zstore(k)%phi(l)+2*pirl
                                    end if
                                end if
                            end do
                        end do
                    end if
                else 
                    next=picker(j+1)
                end if
                
                !Set up gradients for next pass
                if(grad_fin%grad_avlb(next).eq.0)then
                    call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,0)
                end if

            end do

            ! call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,0.0d0,0)
            write(6,"(a,i0,a,f21.16,a,i0,a)") "Energy after epoc no. ",epoc_cnt,": ", &
                grad_fin%prev_erg, ".  ", acpt_cnt, " Zombie state(s) altered."
            epoc_cnt=epoc_cnt+1
            ! if(acpt_cnt.lt.1)then 
            !     call orbital_gd(zstore,grad_fin,elect,dvecs,an_cr,an2_cr2,en,haml,temp_ham,temp_zom,&
            !     epoc_cnt,alphain,b,picker,1,rjct_cnt,epoc_max) 
            ! end if 

            if(epoc_cnt.gt.epoc_max)then
                exit
            end if 

            ! if(acpt_cnt.gt.0)then
            !     d_avrg=d_avrg/acpt_cnt
            !     if(abs(d_avrg).lt.d_avrg_chk)then
            !         d_cnt=d_cnt+1
            !         print*,'Diiference low'
            !     else 
            !         d_cnt=0
            !     end if 
            !     if((d_cnt.ge.d_chck_num).and.(epoc_cnt.gt.1000))then
            !         ! if(btl.lt.9.9d-10)then
            !         !     d_chck_num=d_chck_num-1
            !         ! end if 
            !         newb=newb*alpha
            !         d_cnt=0
            !         print*,'REDUCING THE BTL'
            !         btl=btl*0.1
            !         d_avrg_chk=d_avrg_chk*0.1
            !     end if
            ! end if 

            if((epoc_cnt.gt.2000).and.(modulo(epoc_cnt,500)==0))then
                newb=newb*alpha
            end if 
            
        end do

        return 

    end subroutine full_zs_gd


    subroutine zombie_alter(zstore,grad_fin,haml,elect,en,dvecs,an_cr,an2_cr2)

        implicit none

        type(zombiest),dimension(:),intent(inout)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(oprts),intent(in)::an_cr,an2_cr2
        type(dvector),dimension(:),intent(inout)::dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml
        type(zombiest),dimension(:),allocatable::temp_zom
        type(hamiltonian)::temp_ham
        integer::epoc_cnt,epoc_max
        real(kind=8)::alpha,b
        integer,dimension(ndet-1)::picker

    
        integer::ierr,j,k,l
        
        if (errorflag .ne. 0) return

        epoc_cnt=0 !epoc counter
        if(rstrtflg.eq.'y')then 
            open(unit=450,file='epoc.csv',status="old",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening epoc file to read in. ierr had value ", ierr
                errorflag=1
                return
            end if
            do 
                read(450,*,iostat=ierr)
                if(ierr<0)then
                    close(450)
                    exit
                else if (ierr/=0) then
                    write(0,"(a,i0)") "Error in counting epocs. ierr had value ", ierr
                    errorflag=1
                    return
                end if
                epoc_cnt=epoc_cnt+1
            end do
            close(450) 

            ! call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,0.0d0,1)
            call epoc_writer(grad_fin%prev_erg,epoc_cnt,0,0.0d0,1)
            ierr=0
        else
            epoc_cnt=1
        end if

        alpha=0.5  ! learning rate reduction
        b=8.0D0 !starting learning rate
        
        epoc_max=10000

        do j=1, ndet-1
            picker(j)=j+1
        end do

        call alloczs(temp_zom,ndet)
        call allocham(temp_ham,ndet,norb)
       

    
        if(epoc_cnt.lt.epoc_max)then
            call full_zs_gd(zstore,grad_fin,elect,dvecs,an_cr,an2_cr2,en,haml,temp_ham,&
            temp_zom,epoc_cnt,alpha,b,picker,epoc_max) 
        
            ! call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,an_cr,an2_cr2,en,&
            ! haml,temp_ham,chng_trk,temp_zom,epoc_cnt,alpha,b,picker,(epoc_max-epoc_cnt))
        end if 

        !Brings phi values back within the normal 0-2pi range
        do k=2,ndet 
            do l=1,norb 
                if((zstore(k)%phi(l).gt.2*pirl).or.(zstore(k)%phi(l).lt.0))then 
                    zstore(k)%phi(l)=asin(real(zstore(k)%sin(l)))
                    if((zstore(k)%phi(l).gt.2*pirl))then
                        zstore(k)%phi(l)=zstore(k)%phi(l)-2*pirl
                    else if((zstore(k)%phi(l).lt.0))then
                        zstore(k)%phi(l)=zstore(k)%phi(l)+2*pirl
                    end if
                end if
            end do
        end do

        call deallocham(temp_ham) 
        call dealloczs(temp_zom)
        
        return

    end subroutine zombie_alter

    ! Produces a random order for the ZS to be posisbly changed
    function scramble( number_of_values ) result(out)
        
        implicit none

        integer,intent(in)::number_of_values
        integer::out(number_of_values),array(number_of_values+1)
        integer::n,m,k,j,l,jtemp
        !real::r
        DOUBLE PRECISION, external::ZBQLU01

        out=[(j,j=1,number_of_values)]
        array=[(j,j=1,number_of_values+1)]
        n=1; m=number_of_values
        do k=1,2
            do j=1,m+1
                l = n + FLOOR((m+1-n)*ZBQLU01(1))
                jtemp=array(l)
                array(l)=array(j)
                array(j)=jtemp
            end do
        end do
        n=1
        do j=1,m+1
            if(array(j).eq.1)then 
                cycle
            end if
            out(n)=array(j)
            n=n+1
        end do
        return
     end function scramble

    ! Produces a random order for the ZS to be posisbly changed
    function scramble_norb( number_of_values ) result(out)
    
        implicit none

        integer,intent(in)::number_of_values
        integer::out(number_of_values)
        integer::n,m,k,j,l,jtemp
        DOUBLE PRECISION, external::ZBQLU01
        !real::r
        out=[(j,j=1,number_of_values)]
        
        n=1; m=number_of_values
        do k=1,2
            do j=1,m
                l = n + FLOOR((m-n)*ZBQLU01(1))
                jtemp=out(l)
                out(l)=out(j)
                out(j)=jtemp

            end do
        end do

        return
    end function scramble_norb

END MODULE gradient_descent