MODULE gradient_descent

    use mod_types
    use dnad
    use globvars
    use alarrays
    use ham
    use imgtp
    use outputs
    use infnan_mod
    use zom 

    implicit none 
    integer::d_grad_flg=1
    
    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators

    ! Subroutine that controls and calcualtes all of the hamiltonian variables 
    subroutine he_full_row(temp,zstore,elecs,size,diff_state)

        implicit none 

        type(grad_do),intent(inout)::temp
        ! type(hamiltonian), intent(inout)::haml
        type(zombiest),dimension(:),intent(in)::zstore
        ! type(dual2),dimension(0:),intent(in)::zs_diff
        ! type(zombiest),intent(in)::zs_diff
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::size,diff_state
        integer, allocatable,dimension(:)::IPIV1
        real(dp),allocatable,dimension(:)::WORK1
        
        integer::ierr


        if (errorflag .ne. 0) return
        ierr=0
       
        call haml_ovrlp_column(temp,zstore,ndet,elecs,diff_state)
        
        temp%inv=temp%ovrlp
       
        allocate(WORK1(size),IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV or WORK1 vector allocation . ierr had value ", ierr
            errorflag=1
        end if 

      
        Call dgetrf(size, size, temp%inv%x, size, IPIV1, ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if
        if (ierr==0) call dgetri(size,temp%inv%x,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if

        deallocate(WORK1,IPIV1)
        ! if(d_grad_flg==2)then 
        !     call kinvh_grad(temp%inv,temp%hjk,temp%ovrlp,temp%kinvh)
        ! else
        call DGEMM("N","N",size,size,size,1.d0,temp%inv%x,size,temp%hjk%x,size,0.d0,temp%kinvh%x,size)
        ! end if
        return

    end subroutine he_full_row

    
    subroutine grad_calculate(haml,pick,dvec,grad_fin,erg)

        implicit none 

        type(grad),intent(inout)::grad_fin
        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(inout)::haml
        integer,intent(in)::pick
        type(dual), dimension(:),intent(inout)::erg
        
        integer::j
       
        if (errorflag .ne. 0) return

        call dx_zero(haml%hjk)
        call dx_zero(haml%ovrlp)
        call dx_zero(haml%kinvh)
        call dx_zero(dvec%d)

        do j=1,ndet
            haml%hjk(pick,j)%dx=haml%diff_hjk(:,j,pick); haml%hjk(j,pick)%dx=haml%diff_hjk(:,j,pick)
            haml%ovrlp(pick,j)%dx=haml%diff_ovrlp(:,j,pick); haml%ovrlp(j,pick)%dx=haml%diff_ovrlp(:,j,pick)
        end do
       
        ! if(d_grad_flg==2)then 
        !     call kinvh_grad(haml%inv,haml%hjk,haml%ovrlp,haml%kinvh)
        ! end if
       
        if(grad_fin%grad_avlb(pick)==d_grad_flg)then 
            return
        else
            ! if(d_grad_flg==2)then
            !     call imaginary_time_prop2(dvec,erg,haml,ndet)
            !     grad_fin%vars(pick,:)=erg(1+timesteps)%dx
            !     grad_fin%grad_avlb(pick)=2
            ! else
                erg(timesteps+1)=ergcalc(haml%hjk,dvec%d)
                grad_fin%vars(pick,:)=erg(1+timesteps)%dx
                grad_fin%grad_avlb(pick)=1
            ! end if 

        end if

        ! print*,grad_fin%vars(pick,:)
        ! grad_fin%vars(pick,:)=modulo(grad_fin%vars(pick,:), 2.0 * pirl)
       
        call var_check(grad_fin%vars(pick,:))
      
        return

    end subroutine grad_calculate

    elemental subroutine var_check(var)
    
            implicit none 
    
            real(dp),intent(inout)::var
    
            if(is_nan(var).eqv..true.)then
               
                var=0
            end if
    
            return
    
    end subroutine var_check

    subroutine kinvh_grad(inv,ham,ovrlp,kinvh)

        implicit none 
        type(dual),dimension(:,:),intent(in)::inv,ham,ovrlp
        type(dual),dimension(:,:),intent(inout)::kinvh
        type(dual),dimension(:,:),allocatable::temp_inv_grad,temp_inv_ham
        type(dual)::temp_val1,temp_val2
        real(wp),dimension(:,:),allocatable::temp_kinvh
        integer::j,k,l
        integer::ierr=0

        if (errorflag .ne. 0) return
        
        allocate(temp_inv_grad(ndet,ndet),temp_inv_ham(ndet,ndet),stat=ierr)
        allocate(temp_kinvh(ndet,ndet),stat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in temp_inv_grad, temp_inv_ham or ovrlp_temp, allocation . ierr had value ", ierr
            errorflag=1
            return
        end if
  
        do j=1,ndet
            do k=1,ndet 
                temp_val1=0; temp_val2=0
                do l=1,ndet
                    temp_val1%dx=temp_val1+inv(j,l)%x*ovrlp(l,k)%dx
                    temp_val2=temp_val2+inv(j,l)%x*ham(l,k)
                end do
                temp_inv_grad(j,k)%dx=temp_val1
                temp_inv_ham(j,k)=temp_val2
            end do
        end do

        temp_kinvh=temp_inv_ham%x
        kinvh=matmul(temp_inv_grad,temp_inv_ham)
        kinvh%x=temp_kinvh

        deallocate(temp_inv_grad,temp_inv_ham,temp_kinvh,stat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in temp_inv_grad, temp_inv_ham or ovrlp_temp, deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if
            

    end subroutine kinvh_grad


    subroutine orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,epoc_cnt,alpha,b,picker,maxloop,&
        rjct_cnt_in,epoc_max) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(dual),dimension(:), intent(inout)::erg
        type(grad),intent(inout)::grad_fin
        real(dp),intent(in)::alpha,b
        integer,intent(inout)::epoc_cnt,rjct_cnt_in
        integer,intent(in)::maxloop,epoc_max
        integer,dimension(:),intent(inout)::picker

        type(grad_do)::temp,thread,global

        ! type(dvector):: temp_dvecs,thread_d,global_dvecs       
        ! type(hamiltonian)::temp_ham,thread_ham,global_ham
        ! type(zombiest)::temp_zom,thread_zom,global_zom

        ! type(dual2),dimension(norb)::temp_zom_phi,thread_zom_phi,global_zom_phi
        ! type(dual2),dimension(0:2*norb)::temp_zom_val,thread_zom_val,global_zom_val
        real(wp)::global_min_erg,min_erg
        integer::global_min_idx,min_idx


        integer::rjct_cnt,acpt_cnt,pick,pickorb,rjct_cnt2,loops,lralt_zs,acpt_cnt_2,ierr
        real(dp)::t,erg_str
        integer::j,n,p,loop_max
        integer,dimension(:),allocatable::chng_trk,chng_trk2,pickerorb
       
        
        DOUBLE PRECISION, external::ZBQLU01

        if (errorflag .ne. 0) return
       
        ierr=0

        call alloc_grad_do(temp,ndet,norb)
        call alloc_grad_do(thread,ndet,norb)
        call alloc_grad_do(global,ndet,norb)

        ! call allocdv(temp_dvecs,ndet)
        ! call allocdv(thread_d,ndet)
        ! call allocdv(global_dvecs,ndet)
        ! call allocham(thread_ham,ndet,norb)
        ! call allocham(temp_ham,ndet,norb)
        ! call allocham(global_ham,ndet,norb)

        ! call alloczf(global_zom)
        ! call alloczf(thread_zom)
        ! call alloczf(temp_zom)
       

        allocate(pickerorb(norb),stat=ierr)

        if(ierr==0) allocate(chng_trk(ndet-1),stat=ierr)
        if(ierr==0) allocate(chng_trk2(norb),stat=ierr)
       
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent allocations . ierr had value ", ierr
            errorflag=1
            return
        end if 



        lralt_zs=0    ! power alpha is raised to 
        chng_trk2=0 !stores which orbitals in the ZS have changed 
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        acpt_cnt_2=0
        rjct_cnt2=0
        loops=0
        loop_max=13!5
        p=70-norb
        grad_fin%grad_avlb=0

 
        temp%hjk=haml%hjk; thread%hjk=haml%hjk; global%hjk=haml%hjk
        temp%ovrlp=haml%ovrlp; thread%ovrlp=haml%ovrlp; global%ovrlp=haml%ovrlp   
       
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

                call grad_calculate(haml,pick,dvecs,grad_fin,erg)
                temp%zom=zstore(pick); thread%zom=zstore(pick); global%zom=zstore(pick)
                call dx_zero(temp%hjk);call dx_zero(temp%ovrlp);call dx_zero(temp%kinvh);call dx_zero(temp%dvec%d)
                call dx_zero(thread%hjk);call dx_zero(thread%ovrlp);call dx_zero(thread%kinvh);call dx_zero(thread%dvec%d)
                call dx_zero(global%hjk);call dx_zero(global%ovrlp);call dx_zero(global%kinvh);call dx_zero(global%dvec%d)
                ! global_zom=zstore(pick)
               
                do n=1,norb
                    rjct_cnt=0
                    pickorb=pickerorb(n)
                   
                    global_min_erg=grad_fin%prev_erg
                    global_min_idx=-1
                  
                    ! temp_zom=global_zom
                    !$OMP PARALLEL DEFAULT(NONE) SHARED(loop_max, b, alpha, zstore, grad_fin, haml, elect, ndet, &
                    !$OMP &  timesteps,global_min_erg,global_min_idx,global_zom_phi,global_zom_val,&
                    !$omp & global_ham,global_dvecs,pick,norb,pickorb) &
                    !$OMP & PRIVATE(lralt_zs, temp_zom_phi,temp_zom_val, temp_ham, erg, fxtdk, min_fxtdk, min_fxtdk_idx,& 
                    !$OMP & thread_zom_phi, thread_zom_val, thread_ham,thread_d,temp_dvecs,t)
                   
                    min_erg = grad_fin%prev_erg !0
                    min_idx = -1

                    !$omp do !ordered schedule(static,1)
                    do lralt_zs=1,50
                        t=b*alpha**(lralt_zs-1)
                       
                        temp%zom%phi(pickorb)=global%zom%phi(pickorb)
                        temp%zom%phi(pickorb)%x = temp%zom%phi(pickorb)%x-(t*grad_fin%vars(pick,pickorb))
                        call val_set(temp%zom,pickorb)

                       
                        ! temp_ham=haml
                        
                        call he_full_row(temp,zstore,elect,ndet,pick)
                     
                        ! Imaginary time propagation for back tracing

                        call imaginary_time_erg(temp,ndet)
                     
                        ! call imaginary_time_prop2(temp_dvecs,erg,temp_ham,ndet)
                        ! fxtdk=temp%erg(timesteps+1)
            
                       
                        if((temp%erg%x .lt. min_erg))then
                            min_erg = temp%erg%x
                            min_idx = lralt_zs
                            call grad_do_transfer(temp,thread,pick)
                            ! thread_ham=temp_ham
                            ! thread_zom = temp_zom
                            ! thread_d = temp_dvecs
                            exit
                            !$omp cancel do
                        end if
                       !$omp cancellation point do
                    end do
                    
                    !$omp end do
                    if(min_idx.ne.-1)then
                        !$OMP CRITICAL
                        ! Check the minimum fxtdk value across all threads
                        if ((min_erg .lt. global_min_erg)) then
                            ! Update the global minimum fxtdk and corresponding temp_zom and temp_ham
                            global_min_erg = min_erg
                            global_min_idx = min_idx
                            call grad_do_transfer(thread,global,pick)
                            
                            ! global_zom = thread_zom
                            ! global_ham = thread_ham
                            ! global_dvecs = thread_d
                        end if
                        !$OMP END CRITICAL
                    end if
                    !$OMP END PARALLEL
                    
                    if(global_min_idx .ne. -1)then
                        t=b*0.5**(global_min_idx-1)
                        acpt_cnt=acpt_cnt+1
                        ! dvecs=global_dvecs
                        chng_trk2(acpt_cnt)=pickorb
                        rjct_cnt2=0
                        ! haml=global_ham
                        grad_fin%grad_avlb=0
                        grad_fin%grad_avlb(pick)=d_grad_flg
                        grad_fin%prev_erg=global_min_erg
                        grad_fin%vars(pick,:)=global%erg%dx
                        rjct_cnt_in=0
                    end if

                    call grad_do_transfer(global,temp,pick)
                    
                    write(6,'(1a)',advance='no') '|'
                    flush(6)
                 
                end do
               
                if(acpt_cnt.gt.0)then
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,*(i0:','))")'  ', pick,'          ', &
                    erg_str,'             ',grad_fin%prev_erg,'          ',chng_trk2(1:acpt_cnt) 
                    
                    call grad_do_haml_transfer(global,haml,zstore(pick),dvecs,grad_fin%vars(pick,:),pick)
                    ! zstore(pick)=global_zom
                    call zombiewriter(zstore(pick),pick,0)
                    acpt_cnt_2=acpt_cnt_2+1
                    chng_trk(acpt_cnt_2)=pick
                else 
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,i0)")'  ',pick,'          ', &
                    erg_str,'             ',grad_fin%prev_erg,'          ',0
                   
                    rjct_cnt2=rjct_cnt2+1

                    temp%hjk(:,pick)=haml%hjk(:,pick); thread%hjk(:,pick)=haml%hjk(:,pick); global%hjk(:,pick)=haml%hjk(:,pick)
                    temp%hjk(pick,:)=haml%hjk(pick,:); thread%hjk(pick,:)=haml%hjk(pick,:); global%hjk(pick,:)=haml%hjk(pick,:)
                    temp%ovrlp(:,pick)=haml%ovrlp(:,pick); thread%ovrlp(:,pick)=haml%ovrlp(:,pick)
                    global%ovrlp(:,pick)=haml%ovrlp(:,pick)
                    temp%ovrlp(pick,:)=haml%ovrlp(pick,:); thread%ovrlp(pick,:)=haml%ovrlp(pick,:)
                    global%ovrlp(pick,:)=haml%ovrlp(pick,:)      
                end if
  
            end do
           
            write(6,"(a,i0,a,f21.16)") "Energy after epoch no. ",epoc_cnt,": ",grad_fin%prev_erg
          
            if(acpt_cnt_2.gt.0)then
                call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,0)
                epoc_cnt=epoc_cnt+1
            else
                loops=loops-1
            end if 

            picker=scramble(ndet-1)
           
            ! if(d_grad_flg==1)then
            !     d_grad_flg=2
            ! else if(d_grad_flg==2)then
            !     d_grad_flg=1
            ! end if   

            if(loops.ge.maxloop)then
                grad_fin%grad_avlb=0
                exit
            end if
            acpt_cnt_2=0
            if(epoc_cnt.gt.epoc_max)then 
                exit 
            end if
        end do
        
        call dealloc_grad_do(temp)
        call dealloc_grad_do(thread)
        call dealloc_grad_do(global)
        ! call deallocham(global_ham)
        ! call deallocdv(thread_d)
        ! call deallocham(temp_ham)
        ! call deallocdv(temp_dvecs)

        ! call dealloczf(global_zom)
        ! call dealloczf(thread_zom)
        ! call dealloczf(temp_zom)

        deallocate(pickerorb,stat=ierr)
        if(ierr==0) deallocate(chng_trk,stat=ierr)
        if(ierr==0) deallocate(chng_trk2,stat=ierr)
       
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent deallocations . ierr had value ", ierr
            errorflag=1
            return
        end if 
       
        
        return

    end subroutine orbital_gd


    subroutine full_zs_gd(zstore,elect,dvecs,haml,erg,epoc_cnt,alpha,b,grad_fin,epoc_max,picker) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(dual),dimension(:), intent(inout)::erg
        type(grad),intent(inout)::grad_fin
        real(wp),intent(in)::b,alpha
        integer,intent(inout)::epoc_cnt
        integer,intent(in)::epoc_max
        integer,dimension(:),intent(inout)::picker
        type(grad_do)::temp,thread,global

        ! type(dvector):: temp_dvecs,thread_d,global_dvecs       
        ! type(hamiltonian)::temp_ham,thread_ham,global_ham
        ! type(zombiest)::temp_zom,thread_zom,global_zom

        ! type(dual),dimension(norb)::temp_zom_phi,thread_zom_phi,global_zom_phi
        ! type(dual2),dimension(0:2*norb)::temp_zom_val,thread_zom_val,global_zom_val

        real(wp)::global_min_erg,min_erg
        integer::global_min_idx,min_idx

        integer::lralt,rjct_cnt,rjct_cnt2,acpt_cnt,pick,lralt_temp,loop_max,orb_cnt,ierr
        real(wp)::t
        integer::j
        integer,dimension(:),allocatable::chng_trk
        real(wp),dimension(:),allocatable::lr_chng_trk,erg_chng_trk

        DOUBLE PRECISION, external::ZBQLU01
       
        

        if (errorflag .ne. 0) return
        
        ierr=0
    
      
        allocate(chng_trk(ndet-1),stat=ierr)
        if(ierr==0) allocate(lr_chng_trk(ndet-1),stat=ierr)
        if(ierr==0) allocate(erg_chng_trk(ndet-1),stat=ierr)
       
       
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent allocations . ierr had value ", ierr
            errorflag=1
            return
        end if 
       
        lralt=1    ! power alpha is raised to  
       
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        loop_max=50

        if(epoc_cnt.eq.1)then
            orb_cnt=2
        else
            orb_cnt=2
        end if 

        call alloc_grad_do(temp,ndet,norb)
        call alloc_grad_do(thread,ndet,norb)
        call alloc_grad_do(global,ndet,norb)
        
        ! call allocdv(temp_dvecs,ndet)
        ! call allocdv(thread_d,ndet)
        ! call allocdv(global_dvecs,ndet)
        ! call allocham(temp_ham,ndet,norb)
        ! call allocham(thread_ham,ndet,norb)
        ! call allocham(global_ham,ndet,norb)

        ! call alloczf(global_zom)
        ! call alloczf(thread_zom)
        ! call alloczf(temp_zom)
       
       
        ! temp_ham=haml
        temp%hjk=haml%hjk; thread%hjk=haml%hjk; global%hjk=haml%hjk
        temp%ovrlp=haml%ovrlp; thread%ovrlp=haml%ovrlp; global%ovrlp=haml%ovrlp   
        
        do while(rjct_cnt.lt.(ndet-1)*30)

            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent step &
            &   |       Learning rate      | Acceptance count | Rejection count'
            chng_trk=0
            lr_chng_trk=0
            erg_chng_trk=0
            rjct_cnt2=0

            do j=1,(ndet-1)
                
                pick=picker(j)
                rjct_cnt2=0
                
                
   
                call grad_calculate(haml,pick,dvecs,grad_fin,erg)
                call dx_zero(temp%hjk);call dx_zero(temp%ovrlp);call dx_zero(temp%kinvh);call dx_zero(temp%dvec%d)
                global_min_erg=grad_fin%prev_erg
                global_min_idx=-1
                !$OMP PARALLEL DEFAULT(NONE) SHARED(loop_max, b, alpha, zstore, grad_fin, haml, elect, ndet, &
                !$OMP &  timesteps,global_min_erg,global_min_idx,global_zom_phi,global_zom_val,global_ham,&
                !$omp & global_dvecs,pick,norb) &
                !$OMP & PRIVATE(lralt_temp, temp_zom_phi,temp_zom_val, temp_ham, erg, fxtdk, min_erg, min_fxtdk_idx,&
                !$OMP & thread_zom_phi, thread_zom_val, thread_ham, thread_d,temp_dvecs,t)
                min_erg = grad_fin%prev_erg !0
                min_idx = -1
                !$omp do !ordered schedule(static,1)
                do lralt_temp=1,loop_max
                    !$omp cancellation point do
                
                    t=b*(alpha**(lralt_temp-1))
                    temp%zom=zstore(pick)
                    temp%zom%phi=zstore(pick)%phi
                    temp%zom%phi%x=zstore(pick)%phi%x-(t*grad_fin%vars(pick,:))
                    call val_set(temp%zom)
                 
                  
                    call he_full_row(temp,zstore,elect,ndet,pick)
                    call imaginary_time_erg(temp,ndet)
                
            
                   
                    if((temp%erg%x .lt. min_erg))then
                        min_erg = temp%erg%x
                        min_idx = lralt_temp
                        call grad_do_transfer(temp,thread,pick)
                        ! thread_ham=temp_ham
                        ! thread_zom=temp_zom
                        ! thread_d = temp_dvecs
                        exit
                       !$omp cancel do
                    end if
                    !$omp cancellation point do
                end do 
                !$omp end do
                if(min_idx.ne.-1)then
                !$OMP CRITICAL
                ! Check the minimum fxtdk value across all threads
                if ((min_erg .lt. global_min_erg)) then
                    ! Update the global minimum fxtdk and corresponding temp_zom and temp_ham
                    global_min_erg = min_erg
                    global_min_idx = min_idx
                    call grad_do_transfer(thread,global,pick)
                    ! global_zom=thread_zom
                    ! global_ham = thread_ham
                    ! global_dvecs = thread_d
                end if
                !$OMP END CRITICAL
                end if
                !$OMP END PARALLEL
               

                if(global_min_idx .eq. -1)then
                    rjct_cnt=rjct_cnt+1
                    
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
                    grad_fin%prev_erg,'               ',0.0,'             ',0.0,'        ',acpt_cnt,'          ',rjct_cnt
                    
                    temp%hjk(:,pick)=haml%hjk(:,pick); temp%hjk(pick,:)=haml%hjk(pick,:)
                    temp%ovrlp(:,pick)=haml%ovrlp(:,pick);temp%ovrlp(pick,:)=haml%ovrlp(pick,:)
                else 
                    t=b*(0.5**(global_min_idx-1))
                    acpt_cnt=acpt_cnt+1
                    chng_trk(acpt_cnt)=pick
                    lr_chng_trk(acpt_cnt)=t
                    erg_chng_trk(acpt_cnt)=global_min_erg
                    call grad_do_haml_transfer(global,haml,zstore(pick),dvecs,grad_fin%vars(pick,:),pick)

                    ! zstore(pick)=global_zom
                    call zombiewriter(zstore(pick),pick,0)
                    ! dvecs=global_dvecs
                    ! haml=global_ham
                    rjct_cnt=0
                    ! grad_fin%vars(pick,:)=global_min_erg%dx
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
                grad_fin%prev_erg,'               ',global_min_erg,'             ',t,'        ',acpt_cnt,'          ',rjct_cnt
                  
                    grad_fin%grad_avlb=0
                    grad_fin%grad_avlb(pick)=d_grad_flg
                    grad_fin%prev_erg=global_min_erg
                end if 
                flush(6)
              
            end do
          
          
            write(6,"(a,i0,a,f21.16,a,i0,a)") "Energy after epoc no. ",epoc_cnt,": ", &
            grad_fin%prev_erg, ". ", acpt_cnt, " Zombie state(s) altered."
      
            picker=scramble(ndet-1)
            if(acpt_cnt.gt.0)then
                call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,erg_chng_trk,lr_chng_trk,0)
                epoc_cnt=epoc_cnt+1
            else 
                ! if(d_grad_flg==1)then
                !     d_grad_flg=2
                ! else if(d_grad_flg==2)then
                !     d_grad_flg=1
                ! end if
            end if
           
            orb_cnt=orb_cnt-1
          
            ! if(acpt_cnt.eq.0)then
            if(rjct_cnt.ge.(ndet*4))then
                call orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,&
                epoc_cnt,alpha,b,picker,1,rjct_cnt,epoc_max)
                orb_cnt=orb_cnt+1
            else if((orb_cnt.le.0))then
                call orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,&
                epoc_cnt,alpha,b,picker,2,rjct_cnt,epoc_max)
                orb_cnt=100
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
       
        if (ierr/=0) then
            write(0,"(a,i0)") "Gradient descent deallocations . ierr had value ", ierr
            errorflag=1
            return
        end if 

        call dealloc_grad_do(temp)
        call dealloc_grad_do(thread)
        call dealloc_grad_do(global)

        if(epoc_cnt.lt.epoc_max)then
            call orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,&
            epoc_cnt,alpha,b,picker,epoc_max-epoc_cnt,rjct_cnt,epoc_max)
        end if
      
        
        ! call deallocham(temp_ham)
        ! call deallocham(global_ham)
        ! call deallocdv(temp_dvecs)
        ! call deallocdv(thread_d) 

        ! call dealloczf(global_zom)
        ! call dealloczf(thread_zom)
        ! call dealloczf(temp_zom)
       
        return

    end subroutine full_zs_gd


    subroutine zombie_alter(zstore,haml,elect,erg,dvecs)

        implicit none

        type(zombiest),dimension(:),intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(dual), dimension(:),intent(inout)::erg
        type(hamiltonian),intent(inout)::haml
        type(grad)::grad_fin
        integer,dimension(:),allocatable::picker
        integer::epoc_cnt,epoc_max,cnt,rjct_cnt
        real(dp)::alpha,b
        integer::ierr,j
        
        if (errorflag .ne. 0) return

        ierr=0
        call allocgrad(grad_fin,ndet,norb)
        grad_fin%prev_erg=erg(timesteps+1)
      
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
        else
            call epoc_writer(grad_fin%prev_erg,0,0,0.0d0,0)
        end if

      

        alpha=0.5 ! learning rate reduction
        b=3.2D2   !starting learning rate
       
        epoc_max=2000
        allocate(picker(ndet-1),stat=ierr)
      
        if(epoc_cnt.lt.epoc_max)then
            rjct_cnt=0
            picker=scramble(ndet-1)
            ! call orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,epoc_cnt,alpha,b,picker,1,rjct_cnt,epoc_max) 
            call full_zs_gd(zstore,elect,dvecs,haml,erg,epoc_cnt,alpha,b,grad_fin,epoc_max,picker) 
            ! call orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,epoc_cnt,alpha,b,picker,100,rjct_cnt,epoc_max) 
            ! call full_zs_gd(zstore,elect,dvecs,haml,erg,epoc_cnt,alpha,b,grad_fin,epoc_max,picker) 
            
        end if 

        !Brings phi values back within the normal 0-2pi range
        

        deallocate(picker,stat=ierr)
       
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


    subroutine grad_do_transfer(in,out,pick)

        implicit none
        type(grad_do),intent(in)::in
        type(grad_do),intent(inout)::out
        integer,intent(in)::pick

        if (errorflag .ne. 0) return

        out%hjk(:,pick)=in%hjk(:,pick); out%hjk(pick,:)=in%hjk(pick,:)
        out%ovrlp(:,pick)=in%ovrlp(:,pick); out%ovrlp(pick,:)=in%ovrlp(pick,:)
        ! out%inv=in%inv
        out%kinvh=in%kinvh

        out%diff_hjk_1=in%diff_hjk_1; out%diff_hjk_2=in%diff_hjk_2
        out%diff_ovrlp_1=in%diff_ovrlp_1; out%diff_ovrlp_2=in%diff_ovrlp_2

        out%zom=in%zom
        out%dvec=in%dvec
        out%erg=in%erg
        return 

    end subroutine grad_do_transfer

    subroutine grad_do_haml_transfer(in,haml,zstore,dvec,grads,pick)

        implicit none
        type(grad_do),intent(in)::in
        type(hamiltonian),intent(inout)::haml
        type(zombiest),intent(inout)::zstore
        type(dvector),intent(inout)::dvec
        real(wp),dimension(:),intent(inout)::grads
        integer,intent(in)::pick
        integer::j
        if (errorflag .ne. 0) return

        haml%hjk(:,pick)=in%hjk(:,pick); haml%hjk(pick,:)=in%hjk(pick,:)
        haml%ovrlp(:,pick)=in%ovrlp(:,pick); haml%ovrlp(pick,:)=in%ovrlp(pick,:)
        ! haml%inv=in%inv
        haml%kinvh=in%kinvh

        do j=1, ndet
            haml%diff_hjk(:,j,pick)=in%diff_hjk_1(:,j)
            haml%diff_ovrlp(:,j,pick)=in%diff_ovrlp_1(:,j)
            if(j.ne.pick)then
                haml%diff_hjk(:,pick,j)=in%diff_hjk_2(:,j)
                haml%diff_ovrlp(:,pick,j)=in%diff_ovrlp_2(:,j)
            end if 
        end do

        zstore=in%zom
        dvec=in%dvec
        grads=in%erg%dx
        return 

    end subroutine grad_do_haml_transfer


    ! An elemental subroutine that takes an 1D array and sets A(j)=j 
    subroutine set_array(A)
        implicit none
        integer,dimension(:),intent(inout)::A
        integer::j
        do j=1,size(A)
            A(j)=j
        end do
        return
    end subroutine set_array

END MODULE gradient_descent