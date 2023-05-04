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
       
        call haml_ovrlp_column(haml,zs_diff%val,zstore,an_cr%ham,an2_cr2%ham,elecs,diff_state)
        
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
        integer::j

       
        if (errorflag .ne. 0) return
       
        grad_fin%vars(pick,:)=0

        if(grad_fin%grad_avlb(0,pick).eq.0)then
            dvec(1)%d_diff(:,pick,:)=0
            call gradient_zs(haml,zstore,elect,an_cr,an2_cr2,pick,orb,grad_fin%grad_avlb(1:ndet,pick))
            if(ZBQLU01(1).lt.0.4.and.(orb.ne.0))then
                grad_fin%grad_avlb(0,pick)=1
            end if
        end if 
       
        if(grad_fin%grad_avlb(0,pick).eq.1)then
            dvec(1)%d_diff(:,pick,:)=0
            call sub_matrices(haml,pick)
            call imaginary_time_prop2(dvec,en,haml,pick,0)
        else if(grad_fin%grad_avlb(0,pick).eq.2)then
            dvec(1)%d_diff(:,pick,:)=0
        end if
        
    
        call final_grad(dvec(1),haml,grad_fin,pick,orb)

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
        if(ierr==0) allocate(fibs(13),stat=ierr)
        ! if(ierr==0) allocate(fibs(5),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent allocations . ierr had value ", ierr
            errorflag=1
            return
        end if 

        fibs=[0,1,2,3,5,8,13,21,34,55,89,144,160]
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
                    call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,pickorb)
                    global_min_fxtdk=grad_fin%prev_erg
                    global_min_fxtdk_idx=-1
                    !$OMP PARALLEL DEFAULT(NONE) SHARED(loop_max, b, alphain, fibs, zstore, grad_fin, haml, elect, ndet, &
                    !$OMP & an_cr, an2_cr2, timesteps,global_min_fxtdk,global_min_fxtdk_idx,global_zom,global_ham,&
                    !$omp & global_dvecs,pick,norb,pickorb) &
                    !$OMP & PRIVATE(lralt_zs, temp_zom, temp_ham, en, fxtdk, min_fxtdk, min_fxtdk_idx, thread_zom, thread_ham, &
                    !$OMP & thread_d,temp_dvecs,t)
                    min_fxtdk = grad_fin%prev_erg !0
                    min_fxtdk_idx = -1
                    !$omp do
                    do lralt_zs=1,45!24!loop_max !45

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

                        call imaginary_time_prop2(temp_dvecs,en,temp_ham,0,0)
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
                picker=scramble(ndet-1)
            end if
            
            
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
        if(ierr==0) deallocate(fibs,stat=ierr)
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
        if(ierr==0) allocate(fibs(13),stat=ierr)
        ! if(ierr==0) allocate(fibs(5),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent allocations . ierr had value ", ierr
            errorflag=1
            return
        end if 
       
        

        fibs=[0,1,2,3,5,8,13,21,34,55,89,144,150]
        ! fibs=[0,1,5,34,144]
        


        alpha=alphain  ! learning rate reduction
        lralt=1    ! power alpha is raised to  
        newb=b
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        loop_max=13!5
        ! if(epoc_cnt.eq.1)then
        !     orb_cnt=100
        ! else
            orb_cnt=0
        ! end if 
        
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
                call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,0)
                global_min_fxtdk= grad_fin%prev_erg
                ! g_grad=dot_product(grad_fin%vars(pick,:),grad_fin%vars(pick,:))
                global_min_fxtdk_idx=-1
                !$OMP PARALLEL DEFAULT(NONE) SHARED(loop_max, newb, alpha, fibs, zstore, grad_fin, haml, elect, ndet, &
                !$OMP & an_cr, an2_cr2, timesteps,global_min_fxtdk,global_min_fxtdk_idx,global_zom,global_ham,&
                !$omp & global_dvecs,pick,norb) &
                !$OMP & PRIVATE(lralt_temp, temp_zom, temp_ham, en, fxtdk, min_fxtdk, min_fxtdk_idx, thread_zom, thread_ham, &
                !$OMP & thread_d,temp_dvecs,t,g_grad)
                min_fxtdk = grad_fin%prev_erg !0
                min_fxtdk_idx = -1
                !$omp do
                do lralt_temp=1,45!24!loop_max

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

                    call imaginary_time_prop2(temp_dvecs,en,temp_ham,0,0)

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
                 
                    t=newb*(alpha**fibs(global_min_fxtdk_idx))
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

           
    
            if(acpt_cnt.gt.0)then
                picker=scramble(ndet-1)
                call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,erg_chng_trk,lr_chng_trk,0) 
                epoc_cnt=epoc_cnt+1
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
            if(rjct_cnt.ge.(ndet*2)+1)then
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
        if(ierr==0) deallocate(fibs,stat=ierr)
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



END MODULE gradient_descent