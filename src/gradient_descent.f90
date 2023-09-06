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
    subroutine he_full_row(haml,zstore,zs_diff,elecs,size,an_cr,an2_cr2,diff_state)

        implicit none 

        type(hamiltonian), intent(inout)::haml
        type(zombiest),dimension(:),intent(in)::zstore
        type(dual2),dimension(0:),intent(in)::zs_diff
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::size,diff_state
        integer, allocatable,dimension(:)::IPIV1
        real(dp),allocatable,dimension(:)::WORK1
        
        integer::ierr


        if (errorflag .ne. 0) return
        ierr=0
       
        call haml_ovrlp_column(haml,zs_diff,zstore,ndet,an_cr%ham,an2_cr2%ham,elecs,diff_state)
        
        haml%inv=haml%ovrlp
       
        allocate(WORK1(size),IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV or WORK1 vector allocation . ierr had value ", ierr
            errorflag=1
        end if 

      
        Call dgetrf(size, size, haml%inv%x, size, IPIV1, ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if
        if (ierr==0) call dgetri(size,haml%inv%x,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if

        deallocate(WORK1,IPIV1)
        if(d_grad_flg==2)then 
            call kinvh_grad(haml%inv,haml%hjk,haml%ovrlp,haml%kinvh)
        else
            call DGEMM("N","N",size,size,size,1.d0,haml%inv%x,size,haml%hjk%x,size,0.d0,haml%kinvh%x,size)
        end if
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
            haml%hjk(pick,j)%dx=haml%diff_hjk(pick,j,:); haml%hjk(j,pick)%dx=haml%diff_hjk(pick,j,:)
            haml%ovrlp(pick,j)%dx=haml%diff_ovrlp(pick,j,:); haml%ovrlp(j,pick)%dx=haml%diff_ovrlp(pick,j,:)
        end do
       
        if(d_grad_flg==2)then 
            call kinvh_grad(haml%inv,haml%hjk,haml%ovrlp,haml%kinvh)
        end if
       
        if(grad_fin%grad_avlb(pick)==d_grad_flg)then 
            return
        else
            if(d_grad_flg==2)then
                call imaginary_time_prop2(dvec,erg,haml,ndet)
                grad_fin%vars(pick,:)=erg(1+timesteps)%dx
                grad_fin%grad_avlb(pick)=2
            else
                erg(timesteps+1)=ergcalc(haml%hjk,dvec%d)
                grad_fin%vars(pick,:)=erg(1+timesteps)%dx
                grad_fin%grad_avlb(pick)=1
            end if 

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


    subroutine orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,epoc_cnt,alpha,b,picker,maxloop,an_cr,an2_cr2,&
        rjct_cnt_in,epoc_max) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(dual),dimension(:), intent(inout)::erg
        type(oprts),intent(in)::an_cr,an2_cr2
        type(grad),intent(inout)::grad_fin
        real(dp),intent(in)::alpha,b
        integer,intent(inout)::epoc_cnt,rjct_cnt_in
        integer,intent(in)::maxloop,epoc_max
        integer,dimension(:),intent(inout)::picker

        type(dvector):: temp_dvecs,thread_d,global_dvecs       
        type(hamiltonian)::temp_ham,thread_ham,global_ham
        type(dual2),dimension(norb)::temp_zom_phi,thread_zom_phi,global_zom_phi
        type(dual2),dimension(0:2*norb)::temp_zom_val,thread_zom_val,global_zom_val
        type(dual)::global_min_fxtdk,min_fxtdk,fxtdk
       

        integer::rjct_cnt,acpt_cnt,pick,pickorb,rjct_cnt2,loops,lralt_zs,acpt_cnt_2,ierr
        real(dp)::t,erg_str
        integer::j,n,p,loop_max
        integer,dimension(:),allocatable::chng_trk,chng_trk2,pickerorb
        integer::global_min_fxtdk_idx,min_fxtdk_idx
        
        DOUBLE PRECISION, external::ZBQLU01

        if (errorflag .ne. 0) return
       
        ierr=0

        call allocdv(temp_dvecs,ndet)
        call allocdv(thread_d,ndet)
        call allocdv(global_dvecs,ndet)
        call allocham(thread_ham,ndet,norb)
        call allocham(temp_ham,ndet,norb)
        call allocham(global_ham,ndet,norb)

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
                call dual_2_dual2(zstore(pick)%phi(:),global_zom_phi,1)

                global_zom_val(1:norb)=sin(global_zom_phi(1:norb))
                global_zom_val(norb+1:2*norb)=cos(global_zom_phi(1:norb))
       
    
                do n=1,norb
                    rjct_cnt=0
                    pickorb=pickerorb(n)
                   
                    global_min_fxtdk=grad_fin%prev_erg
                    global_min_fxtdk_idx=-1
                    
                    !$OMP PARALLEL DEFAULT(NONE) SHARED(loop_max, b, alpha, zstore, grad_fin, haml, elect, ndet, &
                    !$OMP & an_cr, an2_cr2, timesteps,global_min_fxtdk,global_min_fxtdk_idx,global_zom_phi,global_zom_val,&
                    !$omp & global_ham,global_dvecs,pick,norb,pickorb) &
                    !$OMP & PRIVATE(lralt_zs, temp_zom_phi,temp_zom_val, temp_ham, erg, fxtdk, min_fxtdk, min_fxtdk_idx,& 
                    !$OMP & thread_zom_phi, thread_zom_val, thread_ham,thread_d,temp_dvecs,t)
                   
                    min_fxtdk = grad_fin%prev_erg !0
                    min_fxtdk_idx = -1

                    !$omp do !ordered schedule(static,1)
                    do lralt_zs=1,45
                        t=b*alpha**(lralt_zs-1)
                        temp_zom_phi=global_zom_phi
                        temp_zom_val=global_zom_val
                       
                        temp_zom_phi(pickorb)= temp_zom_phi(pickorb)-(t*grad_fin%vars(pick,pickorb))
                        temp_zom_val(pickorb)=sin(temp_zom_phi(pickorb))
                        temp_zom_val(pickorb+norb)=cos(temp_zom_phi(pickorb))
                  
                        temp_ham=haml
                        
                        call he_full_row(temp_ham,zstore,temp_zom_val,elect,ndet,an_cr,an2_cr2,pick)
                     
                        ! Imaginary time propagation for back tracing

                        erg=0
                        call imaginary_time_prop2(temp_dvecs,erg,temp_ham,ndet)
                        fxtdk=erg(timesteps+1)
            
                       
                        if((fxtdk .lt. min_fxtdk))then
                            min_fxtdk = fxtdk
                            min_fxtdk_idx = lralt_zs
                            thread_ham=temp_ham
                            thread_zom_phi = temp_zom_phi
                            thread_zom_val = temp_zom_val
                            thread_d = temp_dvecs
            
                            !$omp cancel do
                        end if
                       !$omp cancellation point do
                    end do
                    
                    !$omp end do
                    if(min_fxtdk_idx.ne.-1)then
                        !$OMP CRITICAL
                        ! Check the minimum fxtdk value across all threads
                        if ((min_fxtdk .lt. global_min_fxtdk)) then
                            ! Update the global minimum fxtdk and corresponding temp_zom and temp_ham
                            global_min_fxtdk = min_fxtdk
                            global_min_fxtdk_idx = min_fxtdk_idx
                            global_zom_phi = thread_zom_phi
                            global_zom_val = thread_zom_val
                            global_ham = thread_ham
                            global_dvecs = thread_d
                        end if
                        !$OMP END CRITICAL
                    end if
                    !$OMP END PARALLEL
                    
                    if(global_min_fxtdk_idx .ne. -1)then
                        t=b*0.5**(global_min_fxtdk_idx-1)
                        acpt_cnt=acpt_cnt+1
                        dvecs=global_dvecs
                        chng_trk2(acpt_cnt)=pickorb
                        rjct_cnt2=0
                        haml=global_ham
                        grad_fin%grad_avlb(pick)=d_grad_flg
                        grad_fin%prev_erg=global_min_fxtdk
                        grad_fin%vars(pick,:)=global_min_fxtdk%dx
                        rjct_cnt_in=0
                    end if 
                    
                    
                    write(6,'(1a)',advance='no') '|'
                    flush(6)
                 
                end do
               
                if(acpt_cnt.gt.0)then
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,*(i0:','))")'  ', pick,'          ', &
                    erg_str,'             ',grad_fin%prev_erg,'          ',chng_trk2(1:acpt_cnt) 
                  
                    zstore(pick)%phi%x=global_zom_phi%x
                    call dual2_2_dual(zstore(pick)%val,global_zom_val,1)
                    call zombiewriter(zstore(pick),pick,0)
                    acpt_cnt_2=acpt_cnt_2+1
                    chng_trk(acpt_cnt_2)=pick
                else 
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,i0)")'  ',pick,'          ', &
                    erg_str,'             ',grad_fin%prev_erg,'          ',0
                   
                    rjct_cnt2=rjct_cnt2+1
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
       
        call deallocham(global_ham)
        call deallocdv(thread_d)
        call deallocham(temp_ham)
        call deallocdv(temp_dvecs)

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


    subroutine full_zs_gd(zstore,elect,dvecs,haml,erg,epoc_cnt,alpha,b,an_cr,an2_cr2,grad_fin,epoc_max,picker) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(dual),dimension(:), intent(inout)::erg
        type(oprts),intent(in)::an_cr,an2_cr2
        type(grad),intent(inout)::grad_fin
        real(dp),intent(in)::b,alpha
        integer,intent(inout)::epoc_cnt
        integer,intent(in)::epoc_max
        integer,dimension(:),intent(inout)::picker

        type(dvector):: temp_dvecs,thread_d,global_dvecs       
        type(hamiltonian)::temp_ham,thread_ham,global_ham
        type(dual2),dimension(norb)::temp_zom_phi,thread_zom_phi,global_zom_phi
        type(dual2),dimension(0:2*norb)::temp_zom_val,thread_zom_val,global_zom_val
        type(dual)::global_min_fxtdk,min_fxtdk,fxtdk

        integer::lralt,rjct_cnt,rjct_cnt2,acpt_cnt,pick,lralt_temp,loop_max,orb_cnt,ierr
        real(wp)::newb,t
        integer::j
        integer,dimension(:),allocatable::chng_trk
        real(wp),dimension(:),allocatable::lr_chng_trk,erg_chng_trk

        DOUBLE PRECISION, external::ZBQLU01
       
        integer::global_min_fxtdk_idx,min_fxtdk_idx

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
        newb=b
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        loop_max=13!5

        if(epoc_cnt.eq.1)then
            orb_cnt=500
        else
            orb_cnt=100
        end if 
        
        call allocdv(temp_dvecs,ndet)
        call allocdv(thread_d,ndet)
        call allocdv(global_dvecs,ndet)
        call allocham(temp_ham,ndet,norb)
        call allocham(thread_ham,ndet,norb)
        call allocham(global_ham,ndet,norb)
       
       
        temp_ham=haml
        
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
                global_min_fxtdk=grad_fin%prev_erg
              
                global_min_fxtdk_idx=-1
                !$OMP PARALLEL DEFAULT(NONE) SHARED(loop_max, newb, alpha, zstore, grad_fin, haml, elect, ndet, &
                !$OMP & an_cr, an2_cr2, timesteps,global_min_fxtdk,global_min_fxtdk_idx,global_zom_phi,global_zom_val,global_ham,&
                !$omp & global_dvecs,pick,norb) &
                !$OMP & PRIVATE(lralt_temp, temp_zom_phi,temp_zom_val, temp_ham, erg, fxtdk, min_fxtdk, min_fxtdk_idx,&
                !$OMP & thread_zom_phi, thread_zom_val, thread_ham, thread_d,temp_dvecs,t)
                min_fxtdk = grad_fin%prev_erg !0
                min_fxtdk_idx = -1
                !$omp do !ordered schedule(static,1)
                do lralt_temp=1,45!24!loop_max
                    !$omp cancellation point do
                
                    t=newb*(alpha**(lralt_temp-1))
                 
                    call dual_2_dual2((zstore(pick)%phi(:)-(t*grad_fin%vars(pick,:))),temp_zom_phi,1)
                    
                
                    temp_zom_val(1:norb)=sin(temp_zom_phi)
                    temp_zom_val(norb+1:)=cos(temp_zom_phi)
                    temp_ham=haml
                
                    call he_full_row(temp_ham,zstore,temp_zom_val,elect,ndet,an_cr,an2_cr2,pick)

                
                    call imaginary_time_prop2(temp_dvecs,erg,temp_ham,ndet)

                    fxtdk=erg(timesteps+1)
                 
                    if((fxtdk .lt. min_fxtdk))then
                        min_fxtdk = fxtdk
                        min_fxtdk_idx = lralt_temp
                        thread_ham=temp_ham
                        thread_zom_phi = temp_zom_phi
                        thread_zom_val = temp_zom_val
                        thread_d = temp_dvecs
                       !$omp cancel do
                    end if
                    !$omp cancellation point do
                end do 
                !$omp end do
                if(min_fxtdk_idx.ne.-1)then
                !$OMP CRITICAL
                ! Check the minimum fxtdk value across all threads
                if ((min_fxtdk .lt. global_min_fxtdk)) then
                    ! Update the global minimum fxtdk and corresponding temp_zom and temp_ham
                    global_min_fxtdk = min_fxtdk
                    global_min_fxtdk_idx = min_fxtdk_idx
                    global_zom_phi = thread_zom_phi
                    global_zom_val = thread_zom_val
                    global_ham = thread_ham
                    global_dvecs = thread_d
                end if
                !$OMP END CRITICAL
                end if
                !$OMP END PARALLEL
               

                if(global_min_fxtdk_idx .eq. -1)then
                    rjct_cnt=rjct_cnt+1
                    
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
                    grad_fin%prev_erg,'               ',0.0,'             ',0.0,'        ',acpt_cnt,'          ',rjct_cnt
                   
                else 
                    t=newb*(0.5**(global_min_fxtdk_idx-1))
                    acpt_cnt=acpt_cnt+1
                    chng_trk(acpt_cnt)=pick
                    lr_chng_trk(acpt_cnt)=t
                    erg_chng_trk(acpt_cnt)=global_min_fxtdk
                    zstore(pick)%phi%x=global_zom_phi%x
                    call dual2_2_dual(zstore(pick)%val,global_zom_val,1)
                    call zombiewriter(zstore(pick),pick,0)
                    dvecs=global_dvecs
                    haml=global_ham
                    rjct_cnt=0
                    grad_fin%vars(pick,:)=global_min_fxtdk%dx
                    write(6,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
                grad_fin%prev_erg,'               ',global_min_fxtdk%x,'             ',t,'        ',acpt_cnt,'          ',rjct_cnt
                  
                    grad_fin%grad_avlb(pick)=d_grad_flg
                    grad_fin%prev_erg=global_min_fxtdk
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
                if(d_grad_flg==1)then
                    d_grad_flg=2
                else if(d_grad_flg==2)then
                    d_grad_flg=1
                end if
            end if
           
            orb_cnt=orb_cnt-1
          
            ! if(acpt_cnt.eq.0)then
            if(rjct_cnt.ge.(ndet*4))then
                call orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,&
                epoc_cnt,alpha,newb,picker,1,an_cr,an2_cr2,rjct_cnt,epoc_max)
                orb_cnt=orb_cnt+1
            else if((orb_cnt.le.0))then
                call orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,&
                epoc_cnt,alpha,newb,picker,2,an_cr,an2_cr2,rjct_cnt,epoc_max)
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

        
        call orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,&
        epoc_cnt,alpha,b,picker,epoc_max-epoc_cnt,an_cr,an2_cr2,rjct_cnt,epoc_max)
      
       
        call deallocham(temp_ham)
        call deallocham(global_ham)
        call deallocdv(temp_dvecs)
        call deallocdv(thread_d) 
       
        return

    end subroutine full_zs_gd


    subroutine zombie_alter(zstore,haml,elect,erg,dvecs,an_cr,an2_cr2)

        implicit none

        type(zombiest),dimension(:),intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(oprts),intent(in)::an_cr,an2_cr2
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

      

        alpha=0.5  ! learning rate reduction
        b=1.D1 !starting learning rate
       
        epoc_max=4 #2000
        allocate(picker(ndet-1),stat=ierr)
      
        if(epoc_cnt.lt.epoc_max)then
            rjct_cnt=0
            picker=scramble(ndet-1)
            call orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,epoc_cnt,alpha,b,picker,2,an_cr,an2_cr2,rjct_cnt,epoc_max) 
            call full_zs_gd(zstore,elect,dvecs,haml,erg,epoc_cnt,alpha,b,an_cr,an2_cr2,grad_fin,epoc_max,picker) 
            ! call orbital_gd(zstore,grad_fin,elect,dvecs,erg,haml,epoc_cnt,alpha,b,picker,100,an_cr,an2_cr2,rjct_cnt,epoc_max) 
            ! call full_zs_gd(zstore,elect,dvecs,haml,erg,epoc_cnt,alpha,b,an_cr,an2_cr2,grad_fin,epoc_max,picker) 
            
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


END MODULE gradient_descent