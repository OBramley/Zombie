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
    subroutine he_full_row(haml,zstore,zs_diff,elecs,size,an_cr,an2_cr2,diff_state)

        implicit none 

        type(hamiltonian), intent(inout)::haml
        type(zombiest),dimension(:),intent(in)::zstore
        type(zombiest),intent(in)::zs_diff
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::size,diff_state
        integer, allocatable,dimension(:)::IPIV1
        real(kind=8),allocatable,dimension(:)::WORK1
        
        integer::ierr


        if (errorflag .ne. 0) return
        ierr=0
       
        call haml_ovrlp_column(haml,zs_diff%val,zstore,ndet,an_cr%ham,an2_cr2%ham,elecs,diff_state)
        
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


    subroutine grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvec,grad_fin,en,orb,epoc_cnt)

        implicit none 

        type(zombiest),dimension(:),intent(in)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvec
        type(hamiltonian),intent(inout)::haml
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::pick,orb,epoc_cnt
        DOUBLE PRECISION, external::ZBQLU01
        type(energy),intent(inout)::en
        integer::j,strt
       
        if (errorflag .ne. 0) return

        strt=0
        if(epoc_cnt.ne.0)then
            if(modulo(epoc_cnt,2).eq.0)then
                strt=1
            else 
                strt=2
            end if 
        end if

       
        grad_fin%vars(pick,:)=0

        if(grad_fin%grad_avlb(0,pick).eq.0)then
            dvec(1)%d_diff(:,pick,:)=0
            call gradient_zs(haml,zstore,elect,an_cr,an2_cr2,pick,orb,grad_fin%grad_avlb(1:ndet,pick),strt)
            if(ZBQLU01(1).lt.0.4.and.(orb.ne.0))then
                grad_fin%grad_avlb(0,pick)=1
            end if
        end if 
       
        if(grad_fin%grad_avlb(0,pick).eq.1)then
            dvec(1)%d_diff(:,pick,:)=0
            call sub_matrices(haml,pick,strt)
            call imaginary_time_prop2(dvec,en,haml,ndet,pick,0)
        else if(grad_fin%grad_avlb(0,pick).eq.2)then
            dvec(1)%d_diff(:,pick,:)=0
        end if
        
    
        call final_grad(dvec(1),haml,grad_fin,pick,orb,strt)

        en%erg=0
        en%t=0
        ! grad_fin%vars(pick,:)=grad_fin%vars(pick,:)+ (8.0 * pirl)
       
        ! grad_fin%vars(pick,:)=modulo(grad_fin%vars(pick,:), 2.0 * pirl)
        if( grad_fin%grad_avlb(0,pick).eq.0)then 
            grad_fin%grad_avlb(0:,pick)=1
        else if( grad_fin%grad_avlb(0,pick).eq.1)then 
            grad_fin%grad_avlb(0:,pick)=2
        else if( grad_fin%grad_avlb(0,pick).eq.2)then 
            grad_fin%grad_avlb(0:,pick)=1
        end if 

        if(orb.eq.0)then
            do j=1,norb
                if(is_nan(grad_fin%vars(pick,j)).eqv..true.)then
                    ! write(0,"(a,i0,a,i0)")"nan detected in zombie state ",pick, " orbtial ",j
                    grad_fin%grad_avlb(j,pick)=0
                    grad_fin%vars(pick,j)=0
                    grad_fin%grad_avlb(0,pick)=0
                end if
            end do
        else
            if(is_nan(grad_fin%vars(pick,orb)).eqv..true.)then
                ! write(0,"(a,i0,a,i0)")"nan detected in zombie state ",pick, " orbtial ",orb
                grad_fin%grad_avlb(orb,pick)=0
                grad_fin%vars(pick,orb)=0
                grad_fin%grad_avlb(0,pick)=0
            end if
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
        type(dvector),dimension(:),allocatable::thread_d,global_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml,temp_ham
        type(hamiltonian)::thread_ham,global_ham
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(inout)::epoc_cnt,rjct_cnt_in
        integer,intent(in)::maxloop
        real(kind=8),intent(in)::alphain,b
        integer,dimension(:),intent(inout)::picker
        type(zombiest)::temp_zom,thread_zom,global_zom
        integer::rjct_cnt,acpt_cnt,pick,pickorb,rjct_cnt2,loops,lralt_zs,acpt_cnt_2,ierr
        real(kind=8)::t,fxtdk,erg_str
        integer::j,k,l,n,p,loop_max
        integer,dimension(:),allocatable::chng_trk,fibs,chng_trk2,pickerorb
        real(kind=8)::global_min_fxtdk,min_fxtdk,g_grad
        integer::global_min_fxtdk_idx,min_fxtdk_idx

        if (errorflag .ne. 0) return
        ierr=0

        call alloczf(temp_zom)
        call alloczf(thread_zom)
        call alloczf(global_zom)
        call allocdv(thread_d,1,ndet,norb)
        call allocdv(global_dvecs,1,ndet,norb)
        call allocham(thread_ham,ndet,norb)
        call allocham(global_ham,ndet,norb)

        allocate(pickerorb(norb),stat=ierr)
        if(ierr==0) allocate(chng_trk(ndet-1),stat=ierr)
        if(ierr==0) allocate(chng_trk2(norb),stat=ierr)
        ! if(ierr==0) allocate(fibs(13),stat=ierr)
        ! if(ierr==0) allocate(fibs(5),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent allocations . ierr had value ", ierr
            errorflag=1
            return
        end if 

        ! fibs=[0,1,2,3,5,8,13,21,34,55,89,144,160]
        ! fibs=[0,1,5,34,144]

        lralt_zs=0    ! power alpha is raised to 
        chng_trk2=0 !stores which orbitals in the ZS have changed 
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        acpt_cnt_2=0
        rjct_cnt2=0
        loops=0
        loop_max=13!5
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
               
                erg_str=grad_fin%prev_erg
                pick=picker(j)
                chng_trk2=0
                acpt_cnt=0
                pickerorb=scramble_norb(norb)

                do n=1,norb
                    rjct_cnt=0
                    pickorb=pickerorb(n)
                    grad_fin%grad_avlb=0
                    haml%diff_hjk=0
                    haml%diff_ovrlp=0
                    call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,pickorb,0)
                    global_min_fxtdk=grad_fin%prev_erg
                    global_min_fxtdk_idx=-1
                    !$OMP PARALLEL DEFAULT(NONE) SHARED(loop_max, b, alphain, zstore, grad_fin, haml, elect, ndet, &
                    !$OMP & an_cr, an2_cr2, timesteps,global_min_fxtdk,global_min_fxtdk_idx,global_zom,global_ham,&
                    !$omp & global_dvecs,pick,norb,pickorb) &
                    !$OMP & PRIVATE(lralt_zs, temp_zom, temp_ham, en, fxtdk, min_fxtdk, min_fxtdk_idx, thread_zom, thread_ham, &
                    !$OMP & thread_d,temp_dvecs,t)
                    min_fxtdk = grad_fin%prev_erg !0
                    min_fxtdk_idx = -1
                    !$omp do 
                    do lralt_zs=1,45!24!loop_max !45
                        !$omp cancellation point do
                        ! t=b*(alphain**fibs(lralt_zs))
                        t=b*0.5**(lralt_zs-1)
                        temp_zom=zstore(pick)
                        ! Setup temporary zombie state
                        temp_zom%phi(pickorb)=zstore(pick)%phi(pickorb)-(t*grad_fin%vars(pick,pickorb))
                        temp_zom%sin(pickorb)=sin(temp_zom%phi(pickorb))
                        temp_zom%cos(pickorb)=cos(temp_zom%phi(pickorb))
                        temp_zom%val(pickorb)=temp_zom%sin(pickorb)
                        temp_zom%val(norb+pickorb)=temp_zom%cos(pickorb)

                       
                        temp_ham%hjk=haml%hjk
                        temp_ham%ovrlp=haml%ovrlp
                        call he_full_row(temp_ham,zstore,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
                       
                        ! Imaginary time propagation for back tracing
                        en%erg=0
                        en%t=0

                        call imaginary_time_prop2(temp_dvecs,en,temp_ham,ndet,0,0)
                        fxtdk=en%erg(1,timesteps+1)
                        ! print*,fxtdk,t
                        if((fxtdk .lt. min_fxtdk))then
                            min_fxtdk = fxtdk
                            min_fxtdk_idx = lralt_zs
                            thread_ham=temp_ham
                            ! thread_ham%hjk=temp_ham%hjk
                            ! thread_ham%ovrlp=temp_ham%ovrlp
                            thread_zom = temp_zom
                            thread_d(1)%d = temp_dvecs(1)%d
                             !$omp cancel do
                        end if
                    end do  
                    !$omp end do
                    if(min_fxtdk_idx.ne.-1)then
                        !$OMP CRITICAL
                        ! Check the minimum fxtdk value across all threads
                        if ((min_fxtdk .lt. global_min_fxtdk)) then
                            ! Update the global minimum fxtdk and corresponding temp_zom and temp_ham
                            global_min_fxtdk = min_fxtdk
                            global_min_fxtdk_idx = min_fxtdk_idx
                            global_zom = thread_zom
                            global_ham = thread_ham
                            global_dvecs(1)%d = thread_d(1)%d
                        end if
                        !$OMP END CRITICAL
                    end if
                    !$OMP END PARALLEL
                    ! print*,global_min_fxtdk,global_min_fxtdk_idx
                    
                    if(global_min_fxtdk_idx .ne. -1)then
                        ! t=b*(alphain**fibs(global_min_fxtdk_idx))
                        t=b*0.5**(global_min_fxtdk_idx-1)
                        acpt_cnt=acpt_cnt+1
                        zstore(pick)%sin(pickorb)=global_zom%sin(pickorb)
                        zstore(pick)%cos(pickorb)=global_zom%cos(pickorb)
                        zstore(pick)%val(pickorb)=zstore(pick)%sin(pickorb)
                        zstore(pick)%val(norb+pickorb)=zstore(pick)%cos(pickorb)
                        dvecs(1)%d=global_dvecs(1)%d
                        chng_trk2(acpt_cnt)=pickorb
                        rjct_cnt2=0
                        haml%ovrlp(:,pick)=global_ham%ovrlp(:,pick); haml%ovrlp(pick,:)=haml%ovrlp(:,pick)
                        haml%hjk(:,pick)=global_ham%hjk(:,pick); haml%hjk(pick,:)=haml%hjk(:,pick)
                        haml%kinvh=global_ham%kinvh
                        grad_fin%grad_avlb=0
                        haml%diff_hjk=0
                        haml%diff_ovrlp=0
                        ! grad_fin%grad_avlb(:,0)=0
                        ! grad_fin%grad_avlb(pick,:)=0
                        ! grad_fin%grad_avlb(:,pick)=0
                        grad_fin%vars=0.0
                        ! haml%diff_hjk(pick,:,:)=0
                        ! haml%diff_hjk(:,:,pick)=0
                        ! haml%diff_ovrlp(pick,:,:)=0
                        ! haml%diff_ovrlp(:,:,pick)=0
                        dvecs(1)%d_diff=0.0d0
                        grad_fin%prev_erg=global_min_fxtdk
                        rjct_cnt_in=0
                    end if 
                      
                        ! if((fxtdk.le.grad_fin%prev_erg))then
                        !     grad_fin%current_erg=grad_fin%prev_erg
                        !     ! Check if energy is lower and accept or reject                          
                        !     acpt_cnt=acpt_cnt+1
                        !     zstore(pick)%sin(pickorb)=temp_zom%sin(pickorb)
                        !     zstore(pick)%cos(pickorb)=temp_zom%cos(pickorb)
                        !     zstore(pick)%val(pickorb)=zstore(pick)%sin(pickorb)
                        !     zstore(pick)%val(norb+pickorb)=zstore(pick)%cos(pickorb)
                        !     dvecs(1)%d=temp_dvecs(1)%d
                        !     chng_trk2(acpt_cnt)=pickorb
                        !     rjct_cnt=0
                        !     rjct_cnt2=0
                        !     haml%ovrlp(:,pick)=temp_ham%ovrlp(:,pick); haml%ovrlp(pick,:)=haml%ovrlp(:,pick)
                        !     haml%hjk(:,pick)=temp_ham%hjk(:,pick); haml%hjk(pick,:)=haml%hjk(:,pick)
                        !     haml%kinvh=temp_ham%kinvh
                        !     dvecs(1)%d=temp_dvecs(1)%d
                        !     grad_fin%grad_avlb(:,0)=0
                        !     grad_fin%grad_avlb(pick,:)=0
                        !     grad_fin%grad_avlb(:,pick)=0
                        !     grad_fin%vars=0.0
                        !     haml%diff_hjk(pick,:,:)=0
                        !     haml%diff_hjk(:,:,pick)=0
                        !     haml%diff_ovrlp(pick,:,:)=0
                        !     haml%diff_ovrlp(:,:,pick)=0
                        !     dvecs(1)%d_diff=0.0d0
                        !     grad_fin%prev_erg=fxtdk
                        !     rjct_cnt_in=0
                           
                        !     EXIT 
                        ! end if
                   
                    write(6,'(1a)',advance='no') '|'
                
                end do
               
                if(acpt_cnt.gt.0)then
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,*(i0:','))")'  ', pick,'          ', &
                    erg_str,'             ',grad_fin%prev_erg,'          ',chng_trk2(1:acpt_cnt) 
                    zstore(pick)%update_num=zstore(pick)%update_num+1
                    call zombiewriter(zstore(pick),pick,0)
                    acpt_cnt_2=acpt_cnt_2+1
                    chng_trk(acpt_cnt_2)=pick
                else 
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,i0)")' ',pick,'          ', &
                    erg_str,'             ',grad_fin%prev_erg,'          ',0
                    rjct_cnt2=rjct_cnt2+1
                end if

            end do
            
            write(6,"(a,i0,a,f21.16)") "Energy after epoch no. ",epoc_cnt,": ",grad_fin%prev_erg
            if(acpt_cnt_2.gt.0)then
                call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,0)
                epoc_cnt=epoc_cnt+1
                ! picker=scramble(ndet-1)
            else 
                grad_fin%prev_erg=grad_fin%prev_erg+1.0d-8
            end if

            picker=scramble(ndet-1)
            
            !Every 100 epoc brings phi values back within the normal 0-2pi range
            ! if(modulo(epoc_cnt,100).eq.0)then 
            !     do k=2,ndet 
            !         do l=1,norb 
            !             if((zstore(k)%phi(l).gt.2*pirl).or.(zstore(k)%phi(l).lt.0))then 
            !                 zstore(k)%phi(l)=asin(real(zstore(k)%sin(l)))
            !                 if((zstore(k)%phi(l).gt.2*pirl))then
            !                     zstore(k)%phi(l)=zstore(k)%phi(l)-2*pirl
            !                 else if((zstore(k)%phi(l).lt.0))then
            !                     zstore(k)%phi(l)=zstore(k)%phi(l)+2*pirl
            !                 end if
            !             end if
            !         end do
            !     end do
            ! end if

            ! if(rjct_cnt_in.ne.0)then
            grad_fin%grad_avlb=0
            ! haml%diff_hjk=0
            ! haml%diff_ovrlp=0
            ! dvecs(1)%d_diff=0
            ! end if

            if(loops.ge.maxloop)then
                grad_fin%grad_avlb=0
                exit
            end if
            acpt_cnt_2=0
        end do

        grad_fin%grad_avlb=0
        call dealloczf(temp_zom)
        call dealloczf(thread_zom)
        call dealloczf(global_zom)
        call deallocham(global_ham)
        call deallocdv(thread_d) 

        deallocate(pickerorb,stat=ierr)
        if(ierr==0) deallocate(chng_trk,stat=ierr)
        if(ierr==0) deallocate(chng_trk2,stat=ierr)
        ! if(ierr==0) deallocate(fibs,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent deallocations . ierr had value ", ierr
            errorflag=1
            return
        end if 

        return

    end subroutine orbital_gd


    subroutine full_zs_gd(zstore,grad_fin,elect,dvecs,en,haml,&
      epoc_cnt,alphain,b,an_cr,an2_cr2,epoc_max) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvecs
        type(dvector), dimension(:), allocatable:: temp_dvecs,thread_d,global_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml
        type(hamiltonian)::temp_ham,thread_ham,global_ham
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(inout)::epoc_cnt
        integer,intent(in)::epoc_max
        real(kind=8),intent(in)::b,alphain
        type(zombiest)::temp_zom,thread_zom,global_zom
        integer::lralt,rjct_cnt,rjct_cnt2,acpt_cnt,pick,lralt_temp,loop_max,orb_cnt,ierr
        real(kind=8)::newb,t,fxtdk,alpha
        integer::j,k,l
        integer,dimension(:),allocatable::chng_trk,fibs,picker
        real(kind=8),dimension(:),allocatable::lr_chng_trk,erg_chng_trk
        DOUBLE PRECISION, external::ZBQLU01

        integer(kind=8)::beginning,rate,end

        real(kind=8)::global_min_fxtdk,min_fxtdk,g_grad
        integer::global_min_fxtdk_idx,min_fxtdk_idx

        if (errorflag .ne. 0) return
        
        ierr=0
    
       
        allocate(chng_trk(ndet-1),stat=ierr)
        if(ierr==0) allocate(lr_chng_trk(ndet-1),stat=ierr)
        if(ierr==0) allocate(erg_chng_trk(ndet-1),stat=ierr)
        allocate(picker(ndet-1),stat=ierr)
        ! if(ierr==0) allocate(fibs(13),stat=ierr)
        ! if(ierr==0) allocate(fibs(5),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent allocations . ierr had value ", ierr
            errorflag=1
            return
        end if 
       
        

        ! fibs=[0,1,2,3,5,8,13,21,34,55,89,144,150]
        ! fibs=[0,1,5,34,144]
        


        alpha=alphain  ! learning rate reduction
        lralt=1    ! power alpha is raised to  
        newb=b
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        loop_max=13!5
        if(epoc_cnt.eq.1)then
            orb_cnt=10
        else
            orb_cnt=0
        end if 
        
        call alloczf(temp_zom)
        call alloczf(thread_zom)
        call alloczf(global_zom)
        call allocdv(temp_dvecs,1,ndet,norb)
        call allocdv(thread_d,1,ndet,norb)
        call allocdv(global_dvecs,1,ndet,norb)
        call allocham(temp_ham,ndet,norb)
        call allocham(thread_ham,ndet,norb)
        call allocham(global_ham,ndet,norb)

        picker=scramble(ndet-1)

        temp_ham=haml
        ! temp_dvecs=dvecs
          
        ! call system_clock(beginning, rate)
      
        do while(rjct_cnt.lt.(ndet-1)*30)
            ! call system_clock(beginning, rate)
    
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent step &
            &   |       Learning rate      | Acceptance count | Rejection count'
            
            chng_trk=0
            lr_chng_trk=0
            erg_chng_trk=0
            rjct_cnt2=0
            ! call system_clock(beginning, rate)
            do j=1,(ndet-1)
                
               
                pick=picker(j)
                rjct_cnt2=0
                call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,0,epoc_cnt)
                global_min_fxtdk= grad_fin%prev_erg
                ! g_grad=dot_product(grad_fin%vars(pick,:),grad_fin%vars(pick,:))
                global_min_fxtdk_idx=-1
                !$OMP PARALLEL DEFAULT(NONE) SHARED(loop_max, newb, alpha, zstore, grad_fin, haml, elect, ndet, &
                !$OMP & an_cr, an2_cr2, timesteps,global_min_fxtdk,global_min_fxtdk_idx,global_zom,global_ham,&
                !$omp & global_dvecs,pick,norb) &
                !$OMP & PRIVATE(lralt_temp, temp_zom, temp_ham, en, fxtdk, min_fxtdk, min_fxtdk_idx, thread_zom, thread_ham, &
                !$OMP & thread_d,temp_dvecs,t,g_grad)
                min_fxtdk = grad_fin%prev_erg !0
                min_fxtdk_idx = -1
                !$omp do
                do lralt_temp=1,45!24!loop_max
                    !$omp cancellation point do
                    ! t=newb*(alpha**fibs(lralt_temp))
                    t=newb*(0.5**(lralt_temp-1))
                   
                    temp_zom%phi=zstore(pick)%phi(:)-(t*grad_fin%vars(pick,:))
                    temp_zom%sin=sin(temp_zom%phi)
                    temp_zom%cos=cos(temp_zom%phi)
                 
                    temp_zom%val(1:norb)=temp_zom%sin
                    temp_zom%val(norb+1:)=temp_zom%cos
                  
                    
                    temp_ham%hjk=haml%hjk
                    temp_ham%ovrlp=haml%ovrlp
              
                    call he_full_row(temp_ham,zstore,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
                    
                    ! Imaginary time propagation for back tracing
                    en%erg=0
                    en%t=0

                    call imaginary_time_prop2(temp_dvecs,en,temp_ham,ndet,0,0)

                    fxtdk=en%erg(1,timesteps+1)
            
                    ! print*,'erg',fxtdk,t, 'inequality', (grad_fin%prev_erg-(t*g_grad*1.0d-6))
                    ! if((fxtdk .lt. grad_fin%prev_erg-(t*g_grad*1.0d-10)).and.(fxtdk .lt. min_fxtdk))then
                    if((fxtdk .lt. min_fxtdk))then
                        min_fxtdk = fxtdk
                        min_fxtdk_idx = lralt_temp
                        thread_ham=temp_ham
                        ! thread_ham%hjk=temp_ham%hjk
                        ! thread_ham%ovrlp=temp_ham%ovrlp
                        thread_zom = temp_zom
                        thread_d(1)%d = temp_dvecs(1)%d 
                        !$omp cancel do
                    end if
                    
                end do 
                !$omp end do
                if(min_fxtdk_idx.ne.-1)then
                !$OMP CRITICAL
                ! Check the minimum fxtdk value across all threads
                if ((min_fxtdk .lt. global_min_fxtdk)) then
                    ! Update the global minimum fxtdk and corresponding temp_zom and temp_ham
                    global_min_fxtdk = min_fxtdk
                    global_min_fxtdk_idx = min_fxtdk_idx
                    global_zom = thread_zom
                    global_ham = thread_ham
                    global_dvecs(1)%d = thread_d(1)%d
                end if
                !$OMP END CRITICAL
                end if
                !$OMP END PARALLEL
               

                if(global_min_fxtdk_idx .eq. -1)then
                    rjct_cnt=rjct_cnt+1
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
                    grad_fin%prev_erg,'               ',0.0,'             ',0.0,'        ',acpt_cnt,'          ',rjct_cnt
                else 
                 
                    ! t=newb*(alpha**fibs(global_min_fxtdk_idx))
                    t=newb*(0.5**(global_min_fxtdk_idx-1))
                    acpt_cnt=acpt_cnt+1
                    chng_trk(acpt_cnt)=pick
                    lr_chng_trk(acpt_cnt)=t
                    erg_chng_trk(acpt_cnt)=global_min_fxtdk
                    zstore(pick)=global_zom
                    zstore(pick)%update_num=zstore(pick)%update_num+1
                    call zombiewriter(zstore(pick),pick,0)
                    dvecs(1)%d=global_dvecs(1)%d
                    haml%ovrlp(:,pick)=global_ham%ovrlp(:,pick); haml%ovrlp(pick,:)=haml%ovrlp(:,pick)
                    haml%hjk(:,pick)=global_ham%hjk(:,pick); haml%hjk(pick,:)=haml%hjk(:,pick)
                    haml%kinvh=global_ham%kinvh
                    rjct_cnt=0
                    grad_fin%grad_avlb(:,0)=0
                    grad_fin%grad_avlb(pick,:)=0
                    grad_fin%grad_avlb(:,pick)=0
                    grad_fin%vars=0.0
                    haml%diff_hjk(pick,:,:)=0
                    haml%diff_hjk(:,:,pick)=0
                    haml%diff_ovrlp(pick,:,:)=0
                    haml%diff_ovrlp(:,:,pick)=0
                    dvecs(1)%d_diff=0

                    write(6,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
                grad_fin%prev_erg,'               ',global_min_fxtdk,'             ',t,'        ',acpt_cnt,'          ',rjct_cnt
                    
                    grad_fin%prev_erg=global_min_fxtdk
                end if 
               
                   
                !    print*, 'fxtdk', fxtdk
                !     if((fxtdk.lt.grad_fin%prev_erg))then
                !         grad_fin%current_erg=grad_fin%prev_erg
                !         acpt_cnt=acpt_cnt+1
                !         chng_trk(acpt_cnt)=pick
                !         lr_chng_trk(acpt_cnt)=t
                !         erg_chng_trk(acpt_cnt)=fxtdk
                !         zstore(pick)=temp_zom
                !         zstore(pick)%update_num=zstore(pick)%update_num+1
                !         call zombiewriter(zstore(pick),pick,0)
                !         haml%ovrlp(:,pick)=temp_ham%ovrlp(:,pick); haml%ovrlp(pick,:)=haml%ovrlp(:,pick)
                !         haml%hjk(:,pick)=temp_ham%hjk(:,pick); haml%hjk(pick,:)=haml%hjk(:,pick)
                !         haml%kinvh=temp_ham%kinvh
                !         dvecs(1)%d=temp_dvecs(1)%d
                !         rjct_cnt=0
                !         rjct_cnt2=0
                !         grad_fin%grad_avlb(:,0)=0
                !         grad_fin%grad_avlb(pick,:)=0
                !         grad_fin%grad_avlb(:,pick)=0
                !         grad_fin%vars=0.0
                !         haml%diff_hjk(pick,:,:)=0
                !         haml%diff_hjk(:,:,pick)=0
                !         haml%diff_ovrlp(pick,:,:)=0
                !         haml%diff_ovrlp(:,:,pick)=0
                !         dvecs(1)%d_diff=0
                !         write(6,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
                !         grad_fin%prev_erg,'               ',fxtdk,'             ',t,'        ',acpt_cnt,'          ',rjct_cnt
                !         grad_fin%prev_erg=fxtdk
                     
                !         Exit

                !     end if 
                !     rjct_cnt2=rjct_cnt2+1
                ! end do
               
                ! if(rjct_cnt2.eq.loop_max)then
                !     rjct_cnt=rjct_cnt+1
                !     write(6,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
                ! grad_fin%prev_erg,'               ',fxtdk,'             ',0.0,'        ',acpt_cnt,'          ',rjct_cnt
                ! end if
            
              
            end do
          
        
            write(6,"(a,i0,a,f21.16,a,i0,a)") "Energy after epoc no. ",epoc_cnt,": ", &
            grad_fin%prev_erg, ". ", acpt_cnt, " Zombie state(s) altered."

           
            picker=scramble(ndet-1)
            if(acpt_cnt.gt.0)then
                ! picker=scramble(ndet-1)
                call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,erg_chng_trk,lr_chng_trk,0) 
                epoc_cnt=epoc_cnt+1
            else 
                grad_fin%prev_erg=grad_fin%prev_erg+1.0d-10
            end if
           
            orb_cnt=orb_cnt-1
            !Every 100 epoc brings phi values back within the normal 0-2pi range
            ! if(modulo(epoc_cnt,100).eq.0)then 
            !     do k=2,ndet 
            !         do l=1,norb 
            !             if((zstore(k)%phi(l).gt.2*pirl).or.(zstore(k)%phi(l).lt.0))then 
            !                 zstore(k)%phi(l)=asin(real(zstore(k)%sin(l)))
            !                 if((zstore(k)%phi(l).gt.2*pirl))then
            !                     zstore(k)%phi(l)=zstore(k)%phi(l)-2*pirl
            !                 else if((zstore(k)%phi(l).lt.0))then
            !                     zstore(k)%phi(l)=zstore(k)%phi(l)+2*pirl
            !                 end if
            !             end if
            !         end do
            !     end do
            ! end if
           
            ! if(acpt_cnt.eq.0)then
            if(rjct_cnt.ge.(ndet*4))then
                call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,&
                epoc_cnt,alphain,newb,picker,1,an_cr,an2_cr2,rjct_cnt)
                orb_cnt=orb_cnt+1
            else if((orb_cnt.le.0))then
                call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,&
                epoc_cnt,alphain,newb,picker,20,an_cr,an2_cr2,rjct_cnt)
                orb_cnt=20
            end if
 
            acpt_cnt=0
            if(epoc_cnt.gt.epoc_max)then 
                exit 
            end if
        end do

        rjct_cnt=0
        
    
      
        deallocate(chng_trk,stat=ierr)
        if(ierr==0) deallocate(lr_chng_trk,stat=ierr)
        if(ierr==0) deallocate(erg_chng_trk,stat=ierr)
        ! if(ierr==0) deallocate(fibs,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent deallocations . ierr had value ", ierr
            errorflag=1
            return
        end if 
        call dealloczf(temp_zom)
        call dealloczf(thread_zom)
        call dealloczf(global_zom)

        call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,&
        epoc_cnt,alphain,b,picker,epoc_max-epoc_cnt,an_cr,an2_cr2,rjct_cnt)
      
       
        call deallocham(temp_ham)
        call deallocham(global_ham)
        call deallocdv(temp_dvecs)
        call deallocdv(thread_d) 
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
        b=1.D1 !starting learning rate
        
        epoc_max=1000
      
    
        if(epoc_cnt.lt.epoc_max)then
            call full_zs_gd(zstore,grad_fin,elect,dvecs,en,haml,epoc_cnt,alpha,b,an_cr,an2_cr2,epoc_max)
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

    subroutine best_start_zom(num_tries,zstore,elect,an_cr,an2_cr2)
        
        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        real(kind=8),dimension(:,:),allocatable::test_randoms
        integer,intent(in) :: num_tries
        real(kind=8)::erg_min
        type(zombiest),dimension(:),allocatable::zstore_try
        type(hamiltonian)::haml
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),allocatable::dvecs
        type(oprts),intent(in)::an_cr,an2_cr2
        type(energy)::en
        DOUBLE PRECISION, external::ZBQLU01
        integer :: l, j, k,ierr,choice
        

        ierr=0
        allocate(test_randoms(norb, num_tries),stat=ierr)
        
        call alloczs(zstore_try,int(2,kind=16))

       
        
        call allocham(haml,2,norb)
        call allocdv(dvecs,1,2,norb)
        call allocerg(en,1)
        !$omp parallel private(j,k,l,choice,test_randoms,zstore_try,erg_min,haml,dvecs,en) &
        !$omp & shared(zstore,elect,an_cr,an2_cr2,num_tries)
        zstore_try(1)%phi(1:nel)=0.5*pirl
        zstore_try(1)%phi(nel+1:)=0
        zstore_try(1)%sin=0
        zstore_try(1)%cos=1
        zstore_try(1)%sin(1:nel)=1
        zstore_try(1)%cos(1:nel)=0
        zstore_try(1)%val(1:norb)=zstore_try(1)%sin
        zstore_try(1)%val(norb+1:)=zstore_try(1)%cos
        !$omp do
        do l=2,ndet 
            do j=1,num_tries 
                do k=1,nel
                    test_randoms(k,j)=0.25*pirl+(0.75*pirl*(ZBQLU01(1)))
                end do 
                do k=nel+1,norb
                    test_randoms(k,j)=(ZBQLU01(1)) 
                end do 
            end do
            choice=1
            do j=1,num_tries
        
                zstore_try(2)%phi(:)=test_randoms(:,j)
                zstore_try(2)%sin=sin(zstore_try(2)%phi)
                zstore_try(2)%cos=cos(zstore_try(2)%phi)
                zstore_try(2)%val(1:)=zstore_try(2)%sin
                zstore_try(2)%val(norb+1:)=zstore_try(2)%cos
            
                call hamgen(haml,zstore_try,elect,2,an_cr,an2_cr2,0)
                call imaginary_time_prop2(dvecs,en,haml,2,0,0)
                if(j.eq.1)then
                    erg_min=en%erg(1,timesteps+1)
                else 
                    if(en%erg(1,timesteps+1).lt.erg_min)then
                        erg_min=en%erg(1,timesteps+1)
                        choice=j
                    end if
                end if
            end do
           
            zstore(l)%phi(:)=test_randoms(:,choice)
            zstore(l)%sin=sin(zstore(l)%phi)
            zstore(l)%cos=cos(zstore(l)%phi)
            zstore(l)%val(1:)=zstore(l)%sin
            zstore(l)%val(norb+1:)=zstore(l)%cos
        end do 
        !$omp end do
        !$omp end parallel 
        
        zstore(1)=zstore_try(1)

        call dealloczs(zstore_try)
        call deallocham(haml)
        call deallocdv(dvecs)
        call deallocerg(en)
        deallocate(test_randoms)

    end subroutine best_start_zom


    ! subroutine partial_correlation(num_observations,elect,an_cr,an2_cr2)
    !         implicit none
            
    !         real(kind=8),dimension(:),allocatable::dependent_variable
    !         real(kind=8),dimension(:,:),allocatable::independent_variables
    !         integer,intent(in) :: num_observations
    !         real(kind=8):: partial_correlations(norb)
    !         type(zombiest),dimension(:),allocatable::zstore
    !         type(hamiltonian)::haml
    !         type(elecintrgl),intent(in)::elect
    !         type(dvector),dimension(:),allocatable::dvecs
    !         type(oprts),intent(in)::an_cr,an2_cr2
    !         type(energy)::en
    !         DOUBLE PRECISION, external::ZBQLU01
    !         real(kind=8):: r(norb+1,norb+1), r_inv(norb+1,norb+1)
    !         real(kind=8):: r_lj, r_lk, r_jk, r_pk, numerator, denominator
    !         integer :: l, j, k, p,ierr
    !         integer, allocatable,dimension(:)::IPIV1
    !         real(kind=8),allocatable,dimension(:)::WORK1

    !         ierr=0
    !         allocate(dependent_variable(num_observations),independent_variables(norb, num_observations))
    !         allocate(WORK1(norb+1),IPIV1(norb+1),stat=ierr)

    !         do j=1,num_observations 
    !             do k=1,norb
    !                 independent_variables(k,j)=2*pirl*(ZBQLU01(1)) 
    !             end do 
    !         end do
            
    !         call alloczs(zstore,int(2,kind=16))
    !         zstore(1)%phi(1:nel)=0.5*pirl
    !         zstore(1)%phi(nel+1:)=0
    !         zstore(1)%sin=0
    !         zstore(1)%cos=1
    !         zstore(1)%sin(1:nel)=1
    !         zstore(1)%cos(1:nel)=0
    !         zstore(1)%val(1:norb)=zstore(1)%sin
    !         zstore(1)%val(norb+1:)=zstore(1)%cos
           
    !         call allocham(haml,2,norb)
    !         call allocdv(dvecs,1,2,norb)
    !         call allocerg(en,1)
            
    !         do j=1,num_observations
               
    !             zstore(2)%phi(:)=independent_variables(:,j)
              
              
    !             zstore(2)%sin=sin(zstore(2)%phi)
    !             zstore(2)%cos=cos(zstore(2)%phi)
    !             zstore(2)%val(1:)=zstore(2)%sin
    !             zstore(2)%val(norb+1:)=zstore(2)%cos
              
    !             call hamgen(haml,zstore,elect,2,an_cr,an2_cr2,0)
    !             call imaginary_time_prop2(dvecs,en,haml,2,0,0)

    !             dependent_variable(j)=abs(en%erg(1,timesteps+1))
    !         end do
            
    !         ! print*,dependent_variable
    !         ! print*,independent_variables(:,1)
    !         ! print*,independent_variables(:,2)
    !         ! Compute correlation matrix
    !         do l = 1, norb
    !           do j = 1, norb
    !             r(l,j) = corr(independent_variables(l,:), independent_variables(j,:), num_observations)
    !             ! r(j,l) = r(l,j)
    !           end do
    !         end do

    !         do j = 1, norb
    !             r(norb+1,j) = corr(dependent_variable, independent_variables(j,:), num_observations)
    !             r(j,norb+1) = corr(independent_variables(j,:),dependent_variable,  num_observations)
    !         end do
           
    !         r_inv=r
    !         ! print*,r
            
    !         ! Compute inverse of correlation matrix
    !         Call dgetrf(norb+1, norb+1, r_inv, norb+1, IPIV1, ierr)
    !         if (ierr/=0) then
    !             write(0,"(a,i0)")"Error in DGETRF",ierr
    !         end if
    !         call dgetri(norb+1,r_inv,norb+1,IPIV1,WORK1,norb+1,ierr)
    !         if (ierr/=0) then
    !             write(0,"(a,i0)")"Error in DGETRF",ierr
    !         end if
    !         ! print*,'************'
    !         ! print*,r_inv
        
    !         do j=1,norb
    !             print*,r_inv(j,norb+1),r_inv(j,j),r_inv(norb+1,norb+1)
    !             ! partial_correlations(j)=(-1)*r_inv(j,norb+1)/sqrt(r_inv(j,j)*r_inv(norb+1,norb+1))
    !             partial_correlations(j)=(-1)*r_inv(j,norb+1)/((r_inv(j,j)*r_inv(norb+1,norb+1))-(r_inv(j,norb+1)**2))
    !         end do
    !         ! ! Compute partial correlations for each independent variable
    !         ! do l = 1, norb
    !         !   numerator = 0.0
    !         !   denominator = 1.0
    !         !   do j = 1, 10
    !         !     if (l /= j) then
    !         !       do k = 1, 10
    !         !         if (l /= k .and. j /= k) then
    !         !           numerator = numerator + r(l,j) * r(j,k) * r(k,l)
    !         !           denominator = denominator * sqrt(1.0 - r(k,j)**2)
    !         !         end if
    !         !       end do
    !         !     end if
    !         !   end do
    !         !   partial_correlations(l) = numerator / denominator
    !         ! end do

    !         print*, partial_correlations
    !         deallocate(WORK1,IPIV1)
    !         deallocate(dependent_variable,independent_variables)

    !         return 
            
    ! end subroutine partial_correlation
          
    ! function corr(x, y, n) result(r)
    !     implicit none
        
    !     real(kind=8), intent(in) :: x(n), y(n)
    !     integer, intent(in) :: n
    !     real(kind=8):: x_mean, y_mean
    !     real(kind=8) :: r, r_num, r_denom, x_sum, y_sum, x_sq_sum, y_sq_sum
    !     integer :: j
        
    !     ! Compute means if not provided
       
    !     x_sum = sum(x)
    !     x_mean = x_sum / n
    !     y_sum = sum(y)
    !     y_mean = y_sum / n
       
        
    !     ! Compute numerator and denominator of correlation coefficient
    !     r_num = 0.0
    !     x_sq_sum = 0.0
    !     y_sq_sum = 0.0
    !     do j = 1, n
    !         r_num = r_num + (x(j) - x_mean) * (y(j) - y_mean)
    !         x_sq_sum = x_sq_sum + (x(j) - x_mean)**2
    !         y_sq_sum = y_sq_sum + (y(j) - y_mean)**2
    !     end do
    !     r_denom = sqrt(x_sq_sum * y_sq_sum)
        
    !     ! Compute correlation coefficient
    !     r = r_num / r_denom
    !     ! print*,'final',x_sq_sum,y_sq_sum,r_denom,r_num,r
    !     return 
    
    ! end function corr

    ! subroutine predict_dependent_variable(x, coeffs, n, y_pred)
    !     implicit none
        
    !     real, intent(in) :: x(n,norb), coeffs(norb,norb), y_mean, y_sd
    !     integer, intent(in) :: n
    !     real, intent(out) :: y_pred(n)
    !     real :: x_mean(norb), x_sd(norb), x_partial(norb)
    !     integer :: l, j
        
    !     ! Compute means and standard deviations of independent variables
    !     do j = 1, 10
    !         x_mean(j) = mean(x(:,j))
    !         x_sd(j) = std_dev(x(:,j))
    !     end do
        
    !     ! Compute mean and standard deviation of dependent variable
    !     y_mean = mean(x(:,11))
    !     y_sd = std_dev(x(:,11))
        
    !     ! Compute predicted values of dependent variable
    !     do l = 1, n
    !         do j = 1, 10
    !         x_partial(j) = (x(l,j) - x_mean(j)) / x_sd(j)
    !         end do
    !         y_pred(l) = y_mean + y_sd * dot_product(x_partial, coeffs(:,11))
    !     end do
        
    

    ! end subroutine predict_dependent_variable
    
    ! ! Compute the mean of an array
    ! function mean(a) result(m)
    ! implicit none
    ! real(kind=8), intent(in) :: a(:)
    ! real(kind=8):: m
    ! integer :: n
    
    ! n = size(a)
    ! m = sum(a) / real(n)
    ! end function mean
          
    !       ! Compute the standard deviation of an array
    ! function std_dev(a) result(sd)
    ! implicit none
    ! real, intent(in) :: a(:)
    ! real :: sd, m
    ! integer :: n
    
    ! n = size(a)
    ! m = mean(a)
    ! sd = sqrt(sum((a - m)**2) / real(n))
    ! end function std_dev



END MODULE gradient_descent