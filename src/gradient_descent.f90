MODULE gradient_descent

    use globvars
    use alarrays
    use ham
    use grad_d
    use imgtp
    use outputs
    use infnan_mod
    ! use ieee_arithmetic
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
        call haml_column(haml%hjk(:,diff_state),zstore(diff_state)%val,zstore,an_cr%ham,an2_cr2%ham,elecs,1)
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

        return

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
        DOUBLE PRECISION, external::ZBQLU01
        type(energy),intent(inout)::en
        integer::typ

        if (errorflag .ne. 0) return
        if(grad_fin%grad_avlb(0,pick).eq.2) return

        typ=0
    
        if(grad_fin%grad_avlb(0,pick).eq.0)then
            dvec(1)%d_diff(:,pick,:)=0
            call gradient_zs(haml,zstore,elect,an_cr,an2_cr2,pick,orb,grad_fin%grad_avlb(1:ndet,pick))
        else if(grad_fin%grad_avlb(0,pick).eq.1)then
            dvec(1)%d_diff(:,pick,:)=0
            call sub_matrices(haml,pick)
            call imgtime_prop(dvec,en,haml,pick,0)
        else if(grad_fin%grad_avlb(0,pick).eq.3)then
            dvec(1)%d_diff(:,pick,:)=0
        else if(grad_fin%grad_avlb(0,pick).eq.4)then
            dvec(1)%d_diff(:,pick,:)=0
            call imgtime_prop(dvec,en,haml,pick,0)
        end if
      
        
        call final_grad(dvec(1),haml,grad_fin,pick,orb)
      
        en%erg=0
        en%t=0

    
        if( grad_fin%grad_avlb(0,pick).eq.1)then 
            grad_fin%grad_avlb(0:,pick)=2
        else if( grad_fin%grad_avlb(0,pick).eq.0)then 
            grad_fin%grad_avlb(0:,pick)=1
        else if( grad_fin%grad_avlb(0,pick).eq.3)then
            grad_fin%grad_avlb(0:,pick)=4
        else if( grad_fin%grad_avlb(0,pick).eq.4)then
            grad_fin%grad_avlb(0:,pick)=3
        end if 
       

        return

    end subroutine grad_calc


    subroutine orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,&
        epoc_cnt,alphain,b,picker,maxloop,an_cr,an2_cr2,rjct_cnt_in) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvecs,temp_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml,temp_ham
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(inout)::epoc_cnt,rjct_cnt_in
        integer,intent(in)::maxloop
        real(kind=8),intent(in)::b,alphain
        integer,dimension(:),intent(inout)::picker
        type(zombiest),dimension(:),allocatable::temp_zom
        integer::rjct_cnt,acpt_cnt,pick,pickorb,rjct_cnt2,loops,lralt_zs,acpt_cnt_2,ierr
        real(kind=8)::t,fxtdk
        integer::j,k,l,n,p
        integer,dimension(:),allocatable::chng_trk,fibs,chng_trk2,pickerorb
 
        if (errorflag .ne. 0) return
        ierr=0

        call alloczs(temp_zom,int(ndet,kind=16))
        allocate(pickerorb(norb),stat=ierr)
        if(ierr==0) allocate(chng_trk(ndet-1),stat=ierr)
        if(ierr==0) allocate(chng_trk2(norb),stat=ierr)
        if(ierr==0) allocate(fibs(12),stat=ierr)
        ! if(ierr==0) allocate(fibs(5),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent allocations . ierr had value ", ierr
            errorflag=1
            return
        end if 

        fibs=[0,1,2,3,5,8,13,21,34,55,89,144]
        ! fibs=[0,1,5,34,144]

        lralt_zs=0    ! power alpha is raised to 
        chng_trk2=0 !stores which orbitals in the ZS have changed 
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        acpt_cnt_2=0
        rjct_cnt2=0
        loops=0

        p=70-norb
       
        do while(rjct_cnt2.lt.(norb*100))
            loops=loops+1
            do p=1, ((norb-8)/2)
                write(6,'(1a)',advance='no') ' '
            end do 
            write(6,"(a)",advance='no') 'Progress'
            do p=1, ((norb-8)/2)
                write(6,'(1a)',advance='no') ' '
            end do 
            write(6,"(a)") ' | Zombie state | Previous Energy     | Energy after Gradient Descent steps   | Orbitals altered '

            chng_trk=0
            acpt_cnt_2=0  

            do j=1,ndet-1
               
                grad_fin%current_erg=grad_fin%prev_erg
                pick=picker(j)
                chng_trk2=0
                acpt_cnt=0
                pickerorb=scramble_norb(norb)

                do n=1,norb
                    rjct_cnt=0
                    pickorb=pickerorb(n)
                    grad_fin%grad_avlb=0
                    call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,pickorb)

                    do lralt_zs=1,12!5

                        t=b*(alphain**fibs(lralt_zs))
                        
                        ! Setup temporary zombie state
                        temp_zom=zstore
                        temp_zom(pick)%phi(pickorb)=zstore(pick)%phi(pickorb)-(t*grad_fin%vars(pick,pickorb))
                        temp_zom(pick)%sin(pickorb)=sin(temp_zom(pick)%phi(pickorb))
                        temp_zom(pick)%cos(pickorb)=cos(temp_zom(pick)%phi(pickorb))
                        temp_zom(pick)%val(pickorb)=temp_zom(pick)%sin(pickorb)
                        temp_zom(pick)%val(norb+pickorb)=temp_zom(pick)%cos(pickorb)
                        temp_ham%hjk=haml%hjk
                        temp_ham%ovrlp=haml%ovrlp
                        call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
                      
                        ! Imaginary time propagation for back tracing
                        temp_dvecs(1)%d=0
                        en%erg=0
                        en%t=0
                        call imgtime_prop(temp_dvecs,en,temp_ham,0,0)
                   
                        fxtdk=en%erg(1,timesteps+1)
                        if((is_nan(fxtdk).eqv..true.))then 
                            dvecs(1)%d_diff=0.0d0
                            grad_fin%grad_avlb(:,0)=0
                            grad_fin%vars=0
                            fxtdk=0.0
                            EXIT
                        end if 
                        
                        
                        if((fxtdk.lt.grad_fin%prev_erg))then
                            ! Check if energy is lower and accept or reject                          
                            acpt_cnt=acpt_cnt+1
                            zstore(pick)=temp_zom(pick)
                            dvecs(1)%d=temp_dvecs(1)%d
                            chng_trk2(acpt_cnt)=pickorb
                            rjct_cnt=0
                            rjct_cnt2=0
                            haml%ovrlp=temp_ham%ovrlp
                            haml%hjk=temp_ham%hjk
                            haml%kinvh=temp_ham%kinvh
                            dvecs(1)%d=temp_dvecs(1)%d
                            grad_fin%grad_avlb(:,0)=0
                            grad_fin%grad_avlb(pick,:)=0
                            grad_fin%grad_avlb(:,pick)=0
                            grad_fin%vars=0.0
                            haml%diff_hjk(pick,:,:)=0
                            haml%diff_hjk(:,:,pick)=0
                            haml%diff_ovrlp(pick,:,:)=0
                            haml%diff_ovrlp(:,:,pick)=0
                            dvecs(1)%d_diff=0
                            ! grad_fin%prev_mmntm(pick,pickerorb)=(t*grad_fin%vars(pick,pickorb))
                            grad_fin%prev_erg=fxtdk
                            rjct_cnt_in=0
                            EXIT 
                        end if
                    
                        rjct_cnt=rjct_cnt+1
                       
                    end do
                    write(6,'(1a)',advance='no') '|'
                
                end do
               
                if(acpt_cnt.gt.0)then
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,*(i0:','))")'  ', pick,'          ', &
                    grad_fin%current_erg,'             ',grad_fin%prev_erg,'          ',chng_trk2(1:acpt_cnt) 
                    zstore(pick)%update_num=zstore(pick)%update_num+1
                    call zombiewriter(zstore(pick),pick,0)
                    acpt_cnt_2=acpt_cnt_2+1
                    chng_trk(acpt_cnt_2)=pick
                else 
                    ! grad_fin%prev_mmntm=0
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,i0)")' ',pick,'          ', &
                    grad_fin%current_erg,'             ',grad_fin%prev_erg,'          ',0
                    rjct_cnt2=rjct_cnt2+1
                end if

            end do
            write(6,"(a,i0,a,f21.16)") "Energy after epoch no. ",epoc_cnt,": ",grad_fin%prev_erg
            call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,0)
            epoc_cnt=epoc_cnt+1

            picker=scramble(ndet-1)
          
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

            if(rjct_cnt_in.ne.0)then
                grad_fin%grad_avlb=0
                haml%diff_hjk=0
                haml%diff_ovrlp=0
                dvecs(1)%d_diff=0
                grad_fin%prev_mmntm=0
            end if

            if(loops.ge.maxloop)then
                grad_fin%grad_avlb=0
                exit
            end if
            acpt_cnt_2=0
        end do

        call dealloczs(temp_zom)
        deallocate(pickerorb,stat=ierr)
        if(ierr==0) deallocate(chng_trk,stat=ierr)
        if(ierr==0) deallocate(chng_trk2,stat=ierr)
        if(ierr==0) deallocate(fibs,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent deallocations . ierr had value ", ierr
            errorflag=1
            return
        end if 

        return

    end subroutine orbital_gd


    subroutine full_zs_gd(zstore,grad_fin,elect,dvecs,en,haml,temp_ham,&
      epoc_cnt,alphain,b,an_cr,an2_cr2,epoc_max) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvecs
        type(dvector), dimension(:), allocatable:: temp_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml,temp_ham
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(inout)::epoc_cnt
        integer,intent(in)::epoc_max
        real(kind=8),intent(in)::b,alphain
        type(zombiest),dimension(:),allocatable::temp_zom
        integer::lralt,rjct_cnt,rjct_cnt2,acpt_cnt,pick,lralt_temp,loop_max,orb_cnt,ierr
        real(kind=8)::newb,t,fxtdk,alpha,mmntma
        integer::j,k,l
        integer,dimension(:),allocatable::chng_trk,fibs,picker
        real(kind=8),dimension(:),allocatable::mmntm,gradient_norm,lr_chng_trk,erg_chng_trk
        DOUBLE PRECISION, external::ZBQLU01

        ! integer(kind=8)::beginning,rate,end
        if (errorflag .ne. 0) return
        
        ierr=0
    
        allocate(mmntm(norb),stat=ierr)
        if(ierr==0) allocate(gradient_norm(norb),stat=ierr)
        if(ierr==0) allocate(chng_trk(ndet-1),stat=ierr)
        if(ierr==0) allocate(lr_chng_trk(ndet-1),stat=ierr)
        if(ierr==0) allocate(erg_chng_trk(ndet-1),stat=ierr)
        allocate(picker(ndet-1),stat=ierr)
        if(ierr==0) allocate(fibs(12),stat=ierr)
        ! if(ierr==0) allocate(fibs(5),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent allocations . ierr had value ", ierr
            errorflag=1
            return
        end if 
       
        call alloczs(temp_zom,int(ndet,kind=16))

        fibs=[0,1,2,3,5,8,13,21,34,55,89,144]
        ! fibs=[0,1,5,34,144]
        


        alpha=alphain  ! learning rate reduction
        lralt=1    ! power alpha is raised to  
        newb=b
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        loop_max=12!5!
        mmntma=0.9
        orb_cnt=100  

        call allocdv(temp_dvecs,1,ndet,norb)
        picker=scramble(ndet-1)


        do while(rjct_cnt.lt.(ndet-1)*50)
            ! call system_clock(beginning, rate)
    
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent step &
            &   |       Learning rate      | Acceptance count | Rejection count'
            
            chng_trk=0
            lr_chng_trk=0
            erg_chng_trk=0
            rjct_cnt2=0
            do j=1,(ndet-1)
                
               
                pick=picker(j)
                call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,0)

                mmntm=((grad_fin%vars(pick,:)))!+(mmntma*grad_fin%prev_mmntm(pick,:)))
                do lralt_temp=1,loop_max

                    t=newb*(alpha**fibs(lralt_temp))
                   
                    ! Setup temporary zombie state
                    temp_zom=zstore
                    temp_zom(pick)%phi(:)=zstore(pick)%phi(:)-t*mmntm
                    temp_zom(pick)%sin(:)=sin(temp_zom(pick)%phi(:))
                    temp_zom(pick)%cos(:)=cos(temp_zom(pick)%phi(:))
                    temp_zom(pick)%val(1:)=temp_zom(pick)%sin
                    temp_zom(pick)%val(norb+1:)=temp_zom(pick)%cos

                    temp_ham%hjk=haml%hjk
                    temp_ham%ovrlp=haml%ovrlp
                    
                    call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
                    
                    ! Imaginary time propagation for back tracing
                    en%erg=0
                    en%t=0
                  
                    call imgtime_prop(temp_dvecs,en,temp_ham,0,0)
                
                    fxtdk=en%erg(1,timesteps+1)
                   
                    if((is_nan(fxtdk).eqv..true.))then 
                        dvecs(1)%d_diff=0.0d0
                        dvecs(1)%d=0.0d0
                        call imgtime_prop(dvecs,en,haml,0,0)
                        grad_fin%grad_avlb(:,0)=0
                        grad_fin%vars=0
                        haml%diff_hjk=0
                        haml%diff_ovrlp=0
                        fxtdk=0.0
                        rjct_cnt2=loop_max
                        EXIT
                        ! call emergency(haml,dvecs,temp_dvecs,en)
                        ! call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,0)
                        ! temp_zom(pick)%sin(:)=sin(temp_zom(pick)%phi(:))
                        ! temp_zom(pick)%cos(:)=cos(temp_zom(pick)%phi(:))
                        ! temp_zom(pick)%val(1:)=temp_zom(pick)%sin
                        ! temp_zom(pick)%val(norb+1:)=temp_zom(pick)%cos
                        ! temp_ham%hjk=haml%hjk
                        ! temp_ham%ovrlp=haml%ovrlp
                        ! call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
                        ! fxtdk=en%erg(1,timesteps+1)
                        ! if(is_nan(fxtdk).eqv..true.)then
                        !     rjct_cnt=rjct_cnt+1
                        !     EXIT
                        ! end if 
                    end if 

                    if((fxtdk.lt.grad_fin%prev_erg))then
                        acpt_cnt=acpt_cnt+1
                        chng_trk(acpt_cnt)=pick
                        lr_chng_trk(acpt_cnt)=t
                        erg_chng_trk(acpt_cnt)=fxtdk
                        zstore(pick)=temp_zom(pick)
                        zstore(pick)%update_num=zstore(pick)%update_num+1
                        call zombiewriter(zstore(pick),pick,0)
                        haml%ovrlp=temp_ham%ovrlp
                        haml%hjk=temp_ham%hjk
                        haml%kinvh=temp_ham%kinvh
                        dvecs(1)%d=temp_dvecs(1)%d
                        rjct_cnt=0
                        rjct_cnt2=0
                        grad_fin%grad_avlb(:,0)=0
                        grad_fin%grad_avlb(pick,:)=0
                        grad_fin%grad_avlb(:,pick)=0
                        grad_fin%vars=0.0
                        haml%diff_hjk(pick,:,:)=0
                        haml%diff_hjk(:,:,pick)=0
                        haml%diff_ovrlp(pick,:,:)=0
                        haml%diff_ovrlp(:,:,pick)=0
                        dvecs(1)%d_diff=0
                        ! grad_fin%prev_mmntm(pick,:)=t*mmntm
                        write(6,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
                        grad_fin%prev_erg,'               ',fxtdk,'             ',t,'        ',acpt_cnt,'          ',rjct_cnt
                        grad_fin%prev_erg=fxtdk
                        Exit

                    end if 
                    rjct_cnt2=rjct_cnt2+1
                end do

                if(rjct_cnt2.eq.loop_max)then
                    rjct_cnt=rjct_cnt+1
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
                grad_fin%prev_erg,'               ',fxtdk,'             ',0.0,'        ',acpt_cnt,'          ',rjct_cnt
                end if
            
              
            end do
          
            ! call system_clock(end)
            ! print *, "elapsed time: ", real(end - beginning) / real(rate)
            
            write(6,"(a,i0,a,f21.16,a,i0,a)") "Energy after epoc no. ",epoc_cnt,": ", &
            grad_fin%prev_erg, ". ", acpt_cnt, " Zombie state(s) altered."
            
            if(acpt_cnt.gt.0)then
                picker=scramble(ndet-1)
                call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,erg_chng_trk,lr_chng_trk,0) 
                epoc_cnt=epoc_cnt+1
            end if 
           
            orb_cnt=orb_cnt-1
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
           
            if((rjct_cnt.ge.((ndet-1)*2)))then
                call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,&
                epoc_cnt,alphain,b,picker,1,an_cr,an2_cr2,rjct_cnt)
                orb_cnt=150
            else if((orb_cnt.le.0))then !(epoc_cnt.eq.2).or.
                call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,&
                epoc_cnt,alphain,b,picker,1,an_cr,an2_cr2,rjct_cnt)
                orb_cnt=150
            end if
 
            acpt_cnt=0
            if(epoc_cnt.gt.epoc_max)then 
                exit 
            end if
        end do

        rjct_cnt=0
        call dealloczs(temp_zom)
        deallocate(mmntm,stat=ierr)
        if(ierr==0) deallocate(gradient_norm,stat=ierr)
        if(ierr==0) deallocate(chng_trk,stat=ierr)
        if(ierr==0) deallocate(lr_chng_trk,stat=ierr)
        if(ierr==0) deallocate(erg_chng_trk,stat=ierr)
        if(ierr==0) deallocate(fibs,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent deallocations . ierr had value ", ierr
            errorflag=1
            return
        end if 

        call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,&
        epoc_cnt,alphain,b,picker,epoc_max-epoc_cnt,an_cr,an2_cr2,rjct_cnt)

        deallocate(picker,stat=ierr)

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

        type(hamiltonian)::temp_ham
        integer::epoc_cnt,epoc_max,cnt
        real(kind=8)::alpha,b
        integer::ierr,j,k,l
        
        if (errorflag .ne. 0) return

       
        epoc_cnt=1 !epoc counter
        if(rstrtflg.eq.'y')then 
            ierr=0
            open(unit=450,file='epoc.csv',status="old",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening epoc file to read in. ierr had value ", ierr
                errorflag=1
                return
            end if
            cnt=0
            do 
                read(450,*,iostat=ierr)
                if(ierr<0)then
                    exit
                else if(ierr>0)then
                    write(0,"(a,i0)") "Error in counting epocs. ierr had value ", ierr
                    errorflag=1
                    return
                else 
                    cnt=cnt+1
                end if
                
            end do
            close(450) 
            open(unit=450,file='epoc.csv',status="old",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening epoc file to read in. ierr had value ", ierr
                errorflag=1
                return
            end if
            do j=1,cnt-1
                read(450,*,iostat=ierr)
            end do 
            read(450,*,iostat=ierr)epoc_cnt
    
            close(450) 

            write(0,"(a,i0)") "Epoc read in as ", epoc_cnt
        
            call epoc_writer(grad_fin%prev_erg,epoc_cnt,0,0.0d0,1)
            epoc_cnt=epoc_cnt+1
        end if

        alpha=0.8  ! learning rate reduction
        b=1.0D0 !starting learning rate
        
        epoc_max=10000
      
        call allocham(temp_ham,ndet,norb)
    
        if(epoc_cnt.lt.epoc_max)then
            call full_zs_gd(zstore,grad_fin,elect,dvecs,en,haml,temp_ham,&
            epoc_cnt,alpha,b,an_cr,an2_cr2,epoc_max)
            
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
       
       
        
        return

    end subroutine zombie_alter


    ! 

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

    subroutine emergency(haml,dvecs,dvecs_temp,en)

            implicit none
            type(dvector),dimension(:),intent(inout)::dvecs,dvecs_temp
            type(energy),intent(inout)::en
            type(hamiltonian),intent(inout)::haml
            integer::j,k,l
            logical::checker
    
            ! do j=2,ndet
            !     do k=1,norb 
            !         if(is_nan(zstore(j)%phi(k)).eqv..true.)then
            !             call random_number(r)
            !             zstore(j)%phi(k)=2*pirl*r 
            !             zstore(j)%sin=sin(zstore(j)%phi(k))
            !             zstore(j)%cos=cos(zstore(j)%phi(k))
            !             zstore(j)%val(k)=zstore(j)%sin(k)
            !             zstore(j)%val(norb+k)=zstore(j)%cos(k)
            !             write(0,"(a,i0,i0)") "Error in ZS value number and orbital ", j,k
            !         end if
            !     end do 
            ! end do
            
            checker=.true.
            do j=1,ndet 
                if(is_nan(dvecs(1)%d(j)).eqv..true.)then
                    dvecs(1)%d=0.0d0
                    dvecs(1)%d_diff=0.0d0
                    dvecs(1)%norm=0.0d0
                    dvecs_temp(1)%d=0.0d0
                    dvecs_temp(1)%d_diff=0.0d0
                    dvecs_temp(1)%norm=0.0d0
                    checker=.false.
                    exit
                end if 
                do k=1,norb
                    do l=1,ndet 
                        if(is_nan(dvecs(1)%d_diff(j,k,l)).eqv..true.)then
                            dvecs(1)%d=0.0d0
                            dvecs(1)%d_diff=0.0d0
                            dvecs(1)%norm=0.0d0
                            dvecs_temp(1)%d=0.0d0
                            dvecs_temp(1)%d_diff=0.0d0
                            dvecs_temp(1)%norm=0.0d0
                            checker=.false.
                            exit
                        end if
                    end do 
                    if(checker.eqv..false.)then
                        exit 
                    end if 
                end do 
                if(checker.eqv..false.)then
                    exit 
                end if 
            end do 
    
            call imgtime_prop(dvecs,en,haml,0,0)
    
            ! call  dealloc(temp_zom)
    
            ! call  deallocham(haml)
            ! call  deallocham(haml_temp)
            ! call  deallocgrad(grad_fin)
            ! call  deallocdv(dvecs)
            ! call  deallocdv(temp_dvecs)
            ! call  deallocerg(en)
    
            ! call allocham(haml,ndet,norb)
            ! call alloc(dvecs,1)
            ! call allocdv(dvecs,1,ndet,norb)
            ! call allocerg(en,1)
            ! call  allocgrad(grad_fin,ndet,norb)
    
            ! call hamgen(haml,zstore,elect,ndet,an_cr,an2_cr2,1)
            ! call imgtime_prop(dvecs,en,haml,diff_state,0)
    
    
            return 
    
        end subroutine emergency



END MODULE gradient_descent