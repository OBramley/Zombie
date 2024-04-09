MODULE gradient_descent

    use mod_types
    use randgen
    use globvars
    use alarrays
    use ham
    use imgtp
    use outputs
    use infnan_mod
    use zom
    use operators 
    ! use neural_network

    implicit none 
   
    integer::epoc_cnt !epoc counter
    integer::rjct_cnt_global=0
    integer::pick !Chosen zombie state
    integer,dimension(:),allocatable::picker
    integer,dimension(:),allocatable::chng_trk
   
    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators

    ! Subroutine that controls and calcualtes all of the hamiltonian variables 
    subroutine he_full_row(temp,zstore,elecs,size,orb)

        implicit none 

        type(grad_do),intent(inout)::temp
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::size,orb
        integer, allocatable,dimension(:)::IPIV1
        real(dp),allocatable,dimension(:)::WORK1
        integer::ierr=0


        if (errorflag .ne. 0) return
        if(orb.eq.0)then
             call haml_ovrlp_column(temp,zstore,ndet,elecs,pick)
        else
            ! call haml_vals_2_orb_2(zstore,temp,elecs,pick,orb)
            call haml_ovrlp_column_orb(temp,zstore,ndet,elecs,pick,orb)    
        end if
        temp%inv=temp%ovrlp
       
        allocate(WORK1(size),IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in IPIV or WORK1 vector allocation . ierr had value ", ierr
            errorflag=1
        end if 

      
        Call dgetrf(size, size, temp%inv, size, IPIV1, ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)")"Error in DGETRF",ierr
        end if
        if (ierr==0) call dgetri(size,temp%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)")"Error in DGETRF",ierr
        end if

        deallocate(WORK1,IPIV1,stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in IPIV or WORK1 vector deallocation . ierr had value ", ierr
            errorflag=1
        end if
       
        call DGEMM("N","N",size,size,size,1.d0,temp%inv,size,temp%hjk,size,0.d0,temp%kinvh,size)
        
        return

    end subroutine he_full_row

    
    subroutine grad_calculate(haml,dvecs,zstore,grad_fin,orb)

        implicit none 

        type(grad),intent(inout)::grad_fin
        type(hamiltonian),intent(inout)::haml
        type(dvector),intent(inout)::dvecs
        type(zombiest),dimension(:),intent(in)::zstore
        real(wp)::ham_c_d
        real(wp),dimension(ndet)::temp
        integer,intent(in)::orb
        integer::j,k
       
        if (errorflag .ne. 0) return

        if(orb==0)then
            if(grad_fin%grad_avlb(1,pick)==0)then
                call DGEMV("N",ndet,ndet,1.d0,haml%hjk,ndet,dvecs%d,1,0.d0,temp,1)
                ham_c_d=(dvecs%d_o_d/(dvecs%norm*dvecs%norm*dvecs%norm))*dot_product(temp,dvecs%d_1)
                !$omp parallel private(j,k,temp) shared(zstore,pick,ndet,norb,dvecs,ham_c_d,grad_fin) 
                !$omp do collapse(2)
                do j=1,ndet
                    do k=1,norb
                        if(grad_fin%ovrlp_grad_avlb(k,j,pick).eq.0)then
                            grad_fin%ovrlp_grad(k,j,pick)=haml%ovrlp(j,pick)*&
                            (zstore(j)%val(k)*zstore(pick)%val(k+norb)-zstore(j)%val(k+norb)*zstore(pick)%val(k))/&
                            (zstore(j)%val(k)*zstore(pick)%val(k)+zstore(j)%val(k+norb)*zstore(pick)%val(k+norb))
                            grad_fin%ovrlp_grad_avlb(k,j,pick)=1
                        end if 
                    end do
                end do 
                !$omp end do
                !$omp critical
                !$omp end critical
                !$omp do
                do k=1,norb
                    temp=grad_fin%ovrlp_grad(k,:,pick)*dvecs%d_1     
                    temp(pick)= dot_product(grad_fin%ovrlp_grad(k,:,pick),dvecs%d_1)
                    grad_fin%vars(pick,k)=dot_product(temp,dvecs%d_1)*ham_c_d
                end do
                !$omp end do
                !$omp end parallel
                grad_fin%grad_avlb(:,pick)=1
            else 
                return
            end if 
        else
            if(grad_fin%grad_avlb(orb,pick)==0)then
                call DGEMV("N",ndet,ndet,1.d0,haml%hjk,ndet,dvecs%d,1,0.d0,temp,1)
                ham_c_d=(dvecs%d_o_d/(dvecs%norm*dvecs%norm*dvecs%norm))*dot_product(temp,dvecs%d_1)
                !$omp parallel do private(j) shared(zstore,grad_fin,haml,pick,ndet,norb,orb)
                do j=1,ndet
                    if(grad_fin%ovrlp_grad_avlb(orb,j,pick).eq.0)then
                        grad_fin%ovrlp_grad(orb,j,pick)=haml%ovrlp(j,pick)*&
                        (zstore(j)%val(orb)*zstore(pick)%val(orb+norb)-zstore(j)%val(orb+norb)*zstore(pick)%val(orb))/&
                        (zstore(j)%val(orb)*zstore(pick)%val(orb)+zstore(j)%val(orb+norb)*zstore(pick)%val(orb+norb))
                        grad_fin%ovrlp_grad_avlb(orb,j,pick)=1
                    end if
                end do  
                !$omp end parallel do
                temp=grad_fin%ovrlp_grad(orb,:,pick)*dvecs%d_1
                temp(pick)=dot_product(grad_fin%ovrlp_grad(orb,:,pick),dvecs%d_1)
                grad_fin%vars(pick,orb)=dot_product(temp,dvecs%d_1)*ham_c_d
                grad_fin%grad_avlb(orb,pick)=1
            else 
                return
            end if
        end if
    
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

    subroutine orbital_gd(zstore,grad_fin,elect,dvecs,haml,maxloop)!,neural_net) 

        implicit none 

        type(zombiest),dimension(:),allocatable,intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(grad),intent(inout)::grad_fin
        integer,intent(in)::maxloop
        type(grad_do)::temp,thread
        type(zombiest),dimension(:),allocatable::zstore_temp
        integer::rjct_cnt,acpt_cnt,pickorb,loops,lralt_zs,acpt_cnt_2,pow,pow_1
        integer::j,n,p,chng_chng,tracker,lralt_extra,extra_flag
        integer,dimension(:),allocatable::chng_trk2,pickerorb
        real(wp)::t,erg_str,num_av,num_av_temp
        integer::ierr=0
        ! type(neural_network_layer),dimension(:)::neural_net

        if (errorflag .ne. 0) return
     
        call alloc_grad_do(temp,ndet)
        call alloc_grad_do(thread,ndet)
        allocate(pickerorb(norb),stat=ierr)
        if(ierr==0) allocate(chng_trk2(norb),stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Gradient descent allocations . ierr had value ", ierr
            errorflag=1
            return
        end if 
       
        chng_trk2=0 !stores which orbitals in the ZS have changed 
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        acpt_cnt_2=0 !counts how many orbitals have been changed
        loops=0 !counts the number of loops through 
        lralt_zs=0 !For the learning rate
        lralt_extra=0 !remove higher learning rates
        tracker=0
        extra_flag=0
        p=70-norb
       
        if((epoc_cnt.gt.2).and.(blind_clone_num.gt.100))then
            chng_chng=blind_clone_num-100
        else
            chng_chng=blind_clone_num
        end if
       
        call haml_to_grad_do(haml,dvecs,temp)
        thread=temp
        num_av=0
        do j=1,ndet
            num_av=num_av+numf(zstore(j),zstore(j))
        end do 
        num_av=num_av/ndet
        
        ! Main GD Loop
        do while((rjct_cnt.lt.(norb*100)).and.(epoc_cnt.lt.epoc_max))
            loops=loops+1
            do p=1, ((norb-8)/2)
                write(stdout,'(1a)',advance='no') ' '
            end do 
            write(stdout,"(a)",advance='no') 'Progress'
            do p=1, ((norb-8)/2)
                write(stdout,'(1a)',advance='no') ' '
            end do 
            write(stdout,"(a)")'   | Zombie state | Previous Energy     | Energy after Gradient Descent steps   | Orbitals altered '
       
            chng_trk=0
            acpt_cnt_2=0  
            t=lr*(lr_alpha**lralt_zs)
            pow=(int(lralt_zs/2)+pow_1)
            do j=1,ndet-1
                write(stdout,'(i3)',advance='no') j
                erg_str=grad_fin%prev_erg
                pick=picker(j)
                chng_trk2=0
                acpt_cnt=0
                pickerorb=scramble_norb(norb)
                call haml_to_grad_do(haml,dvecs,thread)

                do n=1, norb
                    pickorb=n !pickerorb(n)
                    call grad_calculate(haml,dvecs,zstore,grad_fin,pickorb)
                    thread%zom=zstore(pick)
                    temp=thread

                    temp%zom%phi(pickorb) = thread%zom%phi(pickorb)-(t*grad_fin%vars(pick,pickorb))
                    call val_set(temp%zom,pickorb)
                    ! if(lr_loop_max.ge.10)then
                    if((abs(temp%zom%val(n)-zstore(pick)%val(n)).lt.1.0d-11).and.&
                    (abs(temp%zom%val(n+norb)-zstore(pick)%val(n+norb)).lt.1.0d-11))then
                        write(stdout,'(1a)',advance='no') '!'
                        cycle
                    end if
                    !     num_av_temp=num_av
                    !     num_av_temp=((num_av_temp*ndet)-numf(zstore(pick),zstore(pick))+numf(temp%zom,temp%zom))/ndet
                    !     if((abs(num_av_temp-nel).gt.abs(num_av-nel)+5.0d-3*(0.25)**(lralt_zs)))then
                    !         write(stdout,'(1a)',advance='no') '!'
                    !         cycle
                    !     end if
                    ! end if 
                    ! if((lralt_zs.gt.6).and.(abs(num_av_temp-nel).gt.abs(num_av-nel)+5.0d0**(-lralt_zs+2)))then
                    !     write(stdout,'(1a)',advance='no') '!'
                    !     cycle
                    ! end if

                    call he_full_row(temp,zstore,elect,ndet,pickorb)
                    call imaginary_time(temp,ndet)
                   
                    
                    if(grad_fin%prev_erg-temp%erg.ge.1.0d-14)then
                        ! num_av=num_av_temp
                        acpt_cnt=acpt_cnt+1
                        chng_trk2(acpt_cnt)=pickorb
                        rjct_cnt=0
                        rjct_cnt_global=0
                        call grad_do_haml_transfer(temp,haml,zstore(pick),dvecs)
                        thread=temp
                        grad_fin%grad_avlb=0
                        grad_fin%ovrlp_grad_avlb(:,:,pick)=0
                        grad_fin%ovrlp_grad_avlb(:,pick,:)=0
                        grad_fin%prev_erg=temp%erg
                    ! else
                    !     print*,'n',temp%zom%val(n)-zstore(pick)%val(n),temp%zom%val(n+norb)-zstore(pick)%val(n+norb)
                    end if
                
                    write(stdout,'(1a)',advance='no') '|'
                    flush(6)
                end do
                
                if(acpt_cnt.gt.0)then
                    write(stdout,"(a,i3,a,f21.16,a,f21.16,a,*(i0:','))")'  ', pick,'          ', &
                    erg_str,'             ',grad_fin%prev_erg,'          ',chng_trk2(1:acpt_cnt) 
                    
                    acpt_cnt_2=acpt_cnt_2+1
                    chng_trk(acpt_cnt_2)=pick
                else 
                    write(stdout,"(a,i3,a,f21.16,a,f21.16,a,i0)")'  ',pick,'          ', &
                    erg_str,'             ',grad_fin%prev_erg,'          ',0
                   
                    rjct_cnt=rjct_cnt+1
                    rjct_cnt_global=rjct_cnt_global+1   
                end if
            end do
            
        write(stdout,"(a,i0,a,f21.16,a,f10.5)") "Energy after epoch no. ",epoc_cnt,": ",grad_fin%prev_erg, "    Learning rate:",t
        ! write(stdout,"(a,i0,a,f21.16)") "Average number of electrons afer epoch no, ",epoc_cnt,": ",num_av
            if(acpt_cnt_2.gt.0)then
                do j=1,acpt_cnt_2
                    call zombiewriter(zstore(chng_trk(j)),chng_trk(j),zstore(chng_trk(j))%gram_num)
                end do
                call epoc_writer(grad_fin%prev_erg,epoc_cnt,t,chng_trk,0)
                epoc_cnt=epoc_cnt+1
            else
                loops=loops-1
            end if 
            
            if((acpt_cnt_2.lt.(0.15*ndet)).and.(extra_flag.eq.0).and.(lralt_zs.eq.lralt_extra).and.(tracker.gt.-1))then 
                lralt_extra=lralt_extra+1
                extra_flag=1
            end if 
           
            lralt_zs=lralt_zs+1
            chng_chng=chng_chng-1
            if(lralt_zs.gt.lr_loop_max)then
                if(lralt_extra==0)then 
                    lralt_extra=2
                else 
                    lralt_extra=0
                end if 
                lralt_zs=lralt_extra
                extra_flag=0
                if((acpt_cnt_2.lt.((ndet)/3)).or.((ndet.gt.5).and.(acpt_cnt_2.lt.3)).or.(tracker.lt.0))then
                    tracker=tracker+1
                end if
            end if

            picker=scramble(ndet-1)
           
            if(((tracker.ge.1).or.(chng_chng.le.0)))then
                if(ndet.lt.ndet_max)then
                    ! num_av=num_av*ndet
                    deallocate(picker,stat=ierr)
                    allocate(picker(ndet+ndet_increase-1),stat=ierr)
                    picker(ndet_increase+1:)=scramble(ndet-1)
                    do j=1,ndet_increase
                        picker(j)=ndet+j
                    end do
                    ndet=ndet+ndet_increase
                    call dealloc_grad_do(temp)
                    call alloc_grad_do(temp,ndet)
                    call dealloc_grad_do(thread)
                    call alloc_grad_do(thread,ndet)
                    call zstore_increase(zstore,elect,dvecs,haml,grad_fin,temp,thread,lralt_zs,extra_flag,lralt_extra,tracker)
                    if(blind_clone_num.gt.100)then
                        chng_chng=blind_clone_num-100
                    else
                        chng_chng=blind_clone_num
                    end if 
                    deallocate(chng_trk,stat=ierr)
                    allocate(chng_trk(ndet-1),stat=ierr)
                    do j=(ndet-ndet_increase+1),ndet
                        call zombiewriter(zstore(j),j,zstore(j)%gram_num)
                        ! num_av=num_av+numf(zstore(j),zstore(j))
                    end do
                    ! num_av=num_av/ndet
                    lralt_zs=0
                    if(blind_clone_num.gt.100)then
                        chng_chng=blind_clone_num-100
                    else
                        chng_chng=blind_clone_num
                    end if 
                    if(modulo(ndet,10).eq.0)then
                        chng_chng=50
                    end if
                    ! if(ndet==ndet_max)then
                        ! lr_loop_max=12
                        ! chng_chng=150
                    ! end if 
                    ! if(ndet.gt.ndet_max/2)then
                    !     pow_1=pow_1+1
                    ! end if
                else if(lr_loop_max.lt.10)then
                    ! if((tracker.ge.1))then
                        lr_loop_max=lr_loop_max+1
                        blind_clone_num=blind_clone_num+2
                    ! end if
                    tracker=0
                    extra_flag=1
                    ! lralt_extra=0
                    ! lralt_zs=0
                    if(blind_clone_num.gt.100)then
                        chng_chng=blind_clone_num-100
                    else
                        chng_chng=blind_clone_num
                    end if 
                    ! chng_chng=200
                  
                    ! num_av=0
                    ! do j=1,ndet
                    !     num_av=num_av+numf(zstore(j),zstore(j))
                    ! end do 
                    ! num_av=num_av/ndet
                else
                    ! pick=picker(1)
                    ! zstore(pick)%phi(:)=0
                    ! do j=1,ndet
                    !     if(j.ne.pick)then
                    !         zstore(pick)%phi(:)=zstore(pick)%phi(:)+zstore(j)%phi(:)
                    !     end if 
                    ! end do
                    ! zstore(pick)%phi(:)=(zstore(pick)%phi(:)/(ndet-1))
                    ! call val_set(zstore(pick))
                    ! temp%zom=zstore(pick)
                    ! call haml_ovrlp_column(temp,zstore,ndet,elect,pick)
                    ! haml%ovrlp=temp%ovrlp
                    ! haml%hjk=temp%hjk
                    ! call  hamgen_inv(haml,ndet)
                    ! call haml_to_grad_do(haml,dvecs,temp)
                    ! call imaginary_time_erg(temp,ndet)
                    ! grad_fin%prev_erg=temp%erg 
                    ! grad_fin%grad_avlb=0
                    ! grad_fin%ovrlp_grad_avlb=0
                    ! dvecs=temp%dvec

                    thread=temp
                    tracker=0
                    lralt_extra=0
                    lralt_zs=0
                    chng_chng=300
                    chng_chng=50
                    extra_flag=1
                end if  
            end if 
            
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
      
        deallocate(pickerorb,stat=ierr)
        if(ierr==0) deallocate(chng_trk2,stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Gradient descent deallocations . ierr had value ", ierr
            errorflag=1
            return
        end if 
        grad_fin%grad_avlb=0
        
        return

    end subroutine orbital_gd

    subroutine orbital_gd_gram(gramstore,zstore,grad_fin,elect,dvecs,haml,chng_trk2,temp,thread,acpt_cnt_2,lralt_zs,state)

        implicit none 
        type(gram),dimension(:)::gramstore
        type(zombiest),dimension(:),allocatable,intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(grad),intent(inout)::grad_fin
        integer,dimension(:),intent(inout)::chng_trk2
        type(grad_do),intent(inout)::temp,thread
        integer,intent(inout)::acpt_cnt_2,lralt_zs
        integer,intent(in)::state
        type(zombiest),dimension(:),allocatable::zstore_temp
        integer::rjct_cnt,acpt_cnt,pickorb
        integer::j,n,p
   
        real(wp)::t,erg_str
        integer::ierr=0
        
        if (errorflag .ne. 0) return

        chng_trk2=0 !stores which orbitals in the ZS have changed
        acpt_cnt_2=0 !counts how many orbitals have been changed
        chng_trk=0
        t=lr*(lr_alpha**lralt_zs)
        p=70-norb

        call haml_to_grad_do(haml,dvecs,temp)
        thread=temp
        
        do p=1, ((norb-8)/2)
            write(stdout,'(1a)',advance='no') ' '
        end do 
        write(stdout,"(a)",advance='no') 'Progress'
        do p=1, ((norb-8)/2)
            write(stdout,'(1a)',advance='no') ' '
        end do 
        write(stdout,"(a)")'   | Zombie state | Previous Energy     | Energy after Gradient Descent steps   | Orbitals altered '
    
        do j=1,ndet
            pick=picker(j)
            if((state.eq.1).and.(pick.eq.1))then 
                cycle
            end if 
            write(stdout,'(i3)',advance='no') j
            erg_str=grad_fin%prev_erg
            pick=picker(j)
            chng_trk2=0
            acpt_cnt=0
            call haml_to_grad_do(haml,dvecs,thread)
            if(state.gt.1)then
                thread%wf_ovrlp(1:state-1,1:ndet,1:ndet)=gramstore(state)%grads%ovrlp_grad(1:state-1,1:ndet,1:ndet)
            end if 
           
            do n=1, norb
                pickorb=n
                call grad_calculate(haml,dvecs,zstore,grad_fin,pickorb)
                thread%zom=zstore(pick)
                temp=thread
          
                temp%zom%phi(pickorb) = thread%zom%phi(pickorb)-(t*grad_fin%vars(pick,pickorb))
                call val_set(temp%zom,pickorb)
                call he_full_row(temp,zstore,elect,ndet,pickorb)

                if(state.eq.1)then
                    call imaginary_time(temp,ndet)
                else
                    call gram_ovrlp_var_temp(gramstore,temp%wf_ovrlp,temp%zom,state,pick)
                    call imaginary_time(temp,gramstore,ndet,state)
                end if
                if(((grad_fin%prev_erg-temp%erg.ge.1.0d-14).and.state.eq.1).or.&
                ((state.ne.1).and.(abs(grad_fin%prev_erg-temp%erg).le.1.0d-3)))then
                    acpt_cnt=acpt_cnt+1
                    chng_trk2(acpt_cnt)=pickorb
                    rjct_cnt=0
                    rjct_cnt_global=0
                    call grad_do_haml_transfer(temp,haml,zstore(pick),dvecs)
                    gramstore(state)%d_ovrlp_d=dot_product(gramstore(state)%dvecs%d,&
                        matmul(gramstore(state)%haml%ovrlp,gramstore(state)%dvecs%d))
                    thread=temp
                    if(state.gt.1)then
                        gramstore(state)%wf_ovrlp(1:state-1,1:ndet,1:ndet)=thread%wf_ovrlp(1:state-1,1:ndet,1:ndet)
                    end if
                    grad_fin%grad_avlb=0
                    grad_fin%ovrlp_grad_avlb(:,:,pick)=0
                    grad_fin%ovrlp_grad_avlb(:,pick,:)=0
                    grad_fin%prev_erg=temp%erg
                end if
            
                write(stdout,'(1a)',advance='no') '|'
                flush(6)
            end do
            
            if(acpt_cnt.gt.0)then
                write(stdout,"(a,i3,a,f21.16,a,f21.16,a,*(i0:','))")'  ', pick,'          ', &
                erg_str,'             ',grad_fin%prev_erg,'          ',chng_trk2(1:acpt_cnt) 
                
                acpt_cnt_2=acpt_cnt_2+1
                chng_trk(acpt_cnt_2)=pick
            else 
                write(stdout,"(a,i3,a,f21.16,a,f21.16,a,i0)")'  ',pick,'          ', &
                erg_str,'             ',grad_fin%prev_erg,'          ',0
                
                rjct_cnt=rjct_cnt+1
                rjct_cnt_global=rjct_cnt_global+1   
            end if
        end do
        write(stdout,"(a,i0)") "State no. ",state
        write(stdout,"(a,i0,a,f21.16,a,f10.5)") "Energy after epoch no. ",epoc_cnt,": ",grad_fin%prev_erg, "    Learning rate:",t
    
        return

    end subroutine orbital_gd_gram

    subroutine orbital_gd_gram_control(gramstore,elect,epoc_cnts)
        
        implicit none
        type(gram),dimension(:)::gramstore
        type(elecintrgl),intent(in)::elect
        integer,dimension(:),intent(inout)::epoc_cnts
        type(grad_do)::temp,thread
        integer,dimension(:),allocatable::chng_trk2,chng_chng,tracker,lralt_extra,extra_flag,lralt_zs
        integer::acpt_cnt_2,k,j,l
        integer::ierr=0


        if (errorflag .ne. 0) return
        call alloc_grad_do(temp,ndet)
        call alloc_grad_do(thread,ndet)
        allocate(temp%wf_ovrlp(gramnum,ndet,ndet),stat=ierr)
        allocate(thread%wf_ovrlp(gramnum,ndet,ndet),stat=ierr)
        allocate(picker(ndet),stat=ierr)
        if(ierr==0) allocate(chng_trk(ndet),stat=ierr)
        if(ierr==0) allocate(chng_trk2(norb),stat=ierr)
        if(ierr==0) allocate(tracker(gramnum+1),stat=ierr)
        if(ierr==0) allocate(chng_chng(gramnum+1),stat=ierr)
        if(ierr==0) allocate(lralt_extra(gramnum+1),stat=ierr)
        if(ierr==0) allocate(extra_flag(gramnum+1),stat=ierr)
        if(ierr==0) allocate(lralt_zs(gramnum+1),stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Gradient descent allocations . ierr had value ", ierr
            errorflag=1
            return
        end if
        do j=1, gramnum+1
            tracker(j)=0
            lralt_extra(j)=0
            extra_flag(j)=0
            if(epoc_cnts(j).gt.2)then
                chng_chng(j)=blind_clone_num-100
            else
                chng_chng(j)=blind_clone_num
            end if
        end do
        lralt_zs=0
        acpt_cnt_2=0
        picker=scramble_norb(ndet)

        rjct_cnt_global=0
        do while(rjct_cnt_global.lt.(norb*100))
            do j=1,gramnum+1
                do l=1, lr_loop_max
                epoc_cnt=epoc_cnts(j)
                call orbital_gd_gram(gramstore,gramstore(j)%zstore,gramstore(j)%grads,elect,gramstore(j)%dvecs,gramstore(j)%haml,&
                chng_trk2,temp,thread,acpt_cnt_2,lralt_zs(j),j)
                if(acpt_cnt_2.gt.0)then
                    do k=1,acpt_cnt_2
                        call zombiewriter(gramstore(j)%zstore(chng_trk(k)),chng_trk(k),j)
                    end do
                    call epoc_writer(gramstore(j)%grads%prev_erg,epoc_cnts(j),lr*(lr_alpha**lralt_zs(j)),chng_trk,0,j)
                    
                    epoc_cnts(j)=epoc_cnts(j)+1
                end if

                if((acpt_cnt_2.lt.(0.15*ndet)).and.(extra_flag(j).eq.0).and.&
                    (lralt_zs(j).eq.lralt_extra(j)).and.(tracker(j).gt.-1))then 
                    lralt_extra(j)=lralt_extra(j)+1
                    extra_flag(j)=1
                end if

                lralt_zs(j)=lralt_zs(j)+1
                chng_chng(j)=chng_chng(j)-1
                if(lralt_zs(j).gt.lr_loop_max)then
                    lralt_zs(j)=lralt_extra(j)
                    extra_flag(j)=0
                    if((acpt_cnt_2.lt.((ndet)/3)).or.((ndet.gt.5).and.(acpt_cnt_2.lt.3)).or.(tracker(j).lt.0))then
                        tracker(j)=tracker(j)+1
                    end if
                end if
            end do
            end do
           
            if((any(chng_chng.le.0)).or.(any(tracker.ge.1)))then
                if(ndet.lt.ndet_max)then
                    deallocate(picker,stat=ierr)
                    allocate(picker(ndet+ndet_increase),stat=ierr)
                    picker(ndet_increase+1:)=scramble_norb(ndet)
                    do j=1,ndet_increase
                        picker(j)=ndet+j
                    end do
                    ndet=ndet+ndet_increase
                    deallocate(temp%wf_ovrlp,stat=ierr)
                    deallocate(thread%wf_ovrlp,stat=ierr)
                    call dealloc_grad_do(temp)
                    call alloc_grad_do(temp,ndet)
                    call dealloc_grad_do(thread)
                    call alloc_grad_do(thread,ndet)
                    allocate(temp%wf_ovrlp(gramnum,ndet,ndet),stat=ierr)
                    allocate(thread%wf_ovrlp(gramnum,ndet,ndet),stat=ierr)
                    do j=1,gramnum+1
                        if(tracker(j).ge.1)then 
                            chng_chng(j)=blind_clone_num-100
                        else 
                            chng_chng(j)=chng_chng(j)+100
                        end if
                        call zstore_increase(gramstore(j)%zstore,elect,gramstore(j)%dvecs,gramstore(j)%haml,gramstore(j)%grads,&
                                temp,thread,lralt_zs(j),extra_flag(j),lralt_extra(j),tracker(j))
                        do k=(ndet-ndet_increase+1),ndet
                            gramstore(j)%zstore(k)%gram_num=j
                            call zombiewriter(gramstore(j)%zstore(k),k,j)
                        end do

                        if(j.eq.1)then
                            call imaginary_time(temp,ndet)
                            gramstore(1)%d_ovrlp_d=dot_product(gramstore(1)%dvecs%d,&
                                        matmul(gramstore(1)%haml%ovrlp,gramstore(1)%dvecs%d))
                        else 
                            call gram_ovrlp_fill(gramstore,j)
                            temp%wf_ovrlp(1:j-1,1:ndet,1:ndet)=gramstore(j)%grads%ovrlp_grad(1:j-1,1:ndet,1:ndet)
                            call imaginary_time(temp,gramstore,ndet,j)
                            gramstore(j)%d_ovrlp_d=dot_product(gramstore(j)%dvecs%d,&
                                        matmul(gramstore(j)%haml%ovrlp,gramstore(j)%dvecs%d))
                        end if 
                        gramstore(j)%grads%prev_erg=temp%erg
                       
                    end do 
                    deallocate(chng_trk,stat=ierr)
                    allocate(chng_trk(ndet),stat=ierr)
                else if(lr_loop_max.lt.10)then
                    if(any(tracker.ge.1))then
                        lr_loop_max=lr_loop_max+1
                    end if
                    tracker=0
                    lralt_extra=0
                    lralt_zs=0
                    chng_chng=200
                    extra_flag=1
                else 
                    tracker=0
                    lralt_extra=0
                    lralt_zs=0
                    chng_chng=300
                    extra_flag=1
                end if  
            end if 
        
        end do
        
        deallocate(temp%wf_ovrlp,stat=ierr)
        deallocate(thread%wf_ovrlp,stat=ierr)
        call dealloc_grad_do(temp)
        call dealloc_grad_do(thread)
        call dealloc_grad_do(temp)
        call dealloc_grad_do(thread)
        if(ierr==0) deallocate(picker,stat=ierr)
        if(ierr==0) deallocate(chng_trk,stat=ierr)
        if(ierr==0) deallocate(chng_trk2,stat=ierr)
        if(ierr==0) deallocate(tracker,stat=ierr)
        if(ierr==0) deallocate(chng_chng,stat=ierr)
        if(ierr==0) deallocate(lralt_extra,stat=ierr)
        if(ierr==0) deallocate(extra_flag,stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Gradient descent deallocations . ierr had value ", ierr
            errorflag=1
            return
        end if 
        



    end subroutine orbital_gd_gram_control

    subroutine zstore_increase(zstore,elect,dvecs,haml,grad_fin,temp,thread,lralt_zs,extra_flag,lralt_extra,tracker)
        implicit none
        type(zombiest),dimension(:),allocatable,intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(grad),intent(inout)::grad_fin
        type(grad_do),intent(inout)::temp,thread
        type(zombiest),dimension(:),allocatable::zstore_temp
        integer::lralt_zs,extra_flag,lralt_extra,tracker
        integer:: j
        tracker=-1
        extra_flag=1
        lralt_extra=0
      
        call alloczs(zstore_temp,ndet)
        do j=ndet-ndet_increase,ndet
            call biased_func(zstore_temp(j))
            call val_set(zstore_temp(j))
        end do
      
        zstore_temp(1:(ndet-ndet_increase))=zstore
        call dealloczs(zstore)
        call alloczs(zstore,ndet)
        zstore=zstore_temp
        call dealloczs(zstore_temp)
        do j=ndet-ndet_increase,ndet
            temp%zom=zstore(j)
            call haml_ovrlp_column(temp,zstore,ndet,elect,j)
        end do
        temp%ovrlp(1:ndet-ndet_increase,1:ndet-ndet_increase)=haml%ovrlp(1:ndet-ndet_increase,1:ndet-ndet_increase)
        temp%hjk(1:ndet-ndet_increase,1:ndet-ndet_increase)=haml%hjk(1:ndet-ndet_increase,1:ndet-ndet_increase)
        call deallocham(haml)
        call allocham(haml,ndet)
        haml%ovrlp=temp%ovrlp
        haml%hjk=temp%hjk
        call  hamgen_inv(haml,ndet)
        call deallocdv(dvecs)
        call allocdv(dvecs,ndet)
        call haml_to_grad_do(haml,dvecs,temp)
        call imaginary_time(temp,ndet)
        call deallocgrad(grad_fin)
        call allocgrad(grad_fin,ndet,norb)
        grad_fin%prev_erg=temp%erg 
        grad_fin%grad_avlb=0
        grad_fin%ovrlp_grad_avlb=0
        dvecs=temp%dvecs
        thread=temp
        lralt_zs=0
     
        return
    end subroutine zstore_increase


    subroutine full_zs_gd(zstore,elect,dvecs,haml,grad_fin) 

        implicit none 

        type(zombiest),dimension(:),allocatable,intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(grad),intent(inout)::grad_fin
        type(grad_do)::temp,thread
        integer::acpt_cnt,lralt_temp,orb_cnt,j
        real(wp)::t
        real(wp),dimension(:),allocatable::lr_chng_trk,erg_chng_trk
        integer::ierr=0

        if (errorflag .ne. 0) return
    
        if(ierr==0) allocate(lr_chng_trk(ndet-1),stat=ierr)
        if(ierr==0) allocate(erg_chng_trk(ndet-1),stat=ierr)
       
       
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Gradient descent allocations . ierr had value ", ierr
            errorflag=1
            return
        end if 
       
    
        acpt_cnt=0  !counts how many ZS have been changed
        lralt_temp=0
        if(epoc_cnt.eq.1)then
            orb_cnt=lr_loop_max+1
        else
            orb_cnt=lr_loop_max+1
        end if 

        call alloc_grad_do(temp,ndet)
        call alloc_grad_do(thread,ndet)
        call haml_to_grad_do(haml,dvecs,thread)
        temp=thread
       
        do while(rjct_cnt_global.lt.(ndet-1)*30)

            write(stdout,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent step &
            &   |       Learning rate      | Acceptance count | Rejection count'
            chng_trk=0
            lr_chng_trk=0
            erg_chng_trk=0
            call haml_to_grad_do(haml,dvecs,thread)
            t=lr*(lr_alpha**(lralt_temp))
            ! t=0.001d0
            do j=1,(ndet-1)
                
                pick=picker(j)
                call grad_calculate(haml,dvecs,zstore,grad_fin,0)
                temp=thread
                temp%zom=zstore(pick)
                temp%zom%phi=zstore(pick)%phi
                temp%zom%phi=zstore(pick)%phi-(t*grad_fin%vars(pick,:))
                call val_set(temp%zom)
           
                call he_full_row(temp,zstore,elect,ndet,0)
                call imaginary_time(temp,ndet)

                if(grad_fin%prev_erg-temp%erg.ge.1.0d-11)then
                    acpt_cnt=acpt_cnt+1
                    chng_trk(acpt_cnt)=pick
                    lr_chng_trk(acpt_cnt)=t
                    erg_chng_trk(acpt_cnt)=temp%erg
                    call grad_do_haml_transfer(temp,haml,zstore(pick),dvecs)
                    call zombiewriter(zstore(pick),pick,zstore(pick)%gram_num)
                    rjct_cnt_global=0
            
                    write(stdout,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
    grad_fin%prev_erg,'               ',temp%erg,'             ',t,'        ',acpt_cnt,'          ',rjct_cnt_global
        
                    grad_fin%grad_avlb=0
                    grad_fin%ovrlp_grad_avlb(:,:,pick)=0
                    grad_fin%ovrlp_grad_avlb(:,pick,:)=0
                    grad_fin%prev_erg=temp%erg
                    thread=temp
                else
                    rjct_cnt_global=rjct_cnt_global+1
                    write(stdout,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
                    grad_fin%prev_erg,'               ',0.0,'             ',0.0,'        ',acpt_cnt,'          ',rjct_cnt_global
                end if

                flush(6)
                
            end do
          
           
            write(stdout,"(a,i0,a,f21.16,a,i0,a)") "Energy after epoc no. ",epoc_cnt,": ", &
            grad_fin%prev_erg, ". ", acpt_cnt, " Zombie state(s) altered."
      
            picker=scramble(ndet-1)
            if(acpt_cnt.gt.0)then
                call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,erg_chng_trk,lr_chng_trk,0)
                epoc_cnt=epoc_cnt+1
            end if
            lralt_temp=lralt_temp+1
        
            orb_cnt=orb_cnt-1

            ! if((orb_cnt.le.0))then
            if((lralt_temp.gt.lr_loop_max))then
                lralt_temp=0
                if(acpt_cnt.eq.0)then
                    lralt_temp=lralt_temp+1
                
                    ! call orbital_gd(zstore,grad_fin,elect,dvecs,haml,loop_max+1)
                    call haml_to_grad_do(haml,dvecs,thread)
                    orb_cnt=lr_loop_max+1
                end if
            end if
 
            acpt_cnt=0
            if(epoc_cnt.gt.epoc_max)then 
                exit 
            end if
        end do

        if(ierr==0) deallocate(lr_chng_trk,stat=ierr)
        if(ierr==0) deallocate(erg_chng_trk,stat=ierr)
       
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Gradient descent deallocations . ierr had value ", ierr
            errorflag=1
            return
        end if 

        call dealloc_grad_do(temp)
        call dealloc_grad_do(thread)

        if(epoc_cnt.lt.epoc_max)then
            ! call orbital_gd(zstore,grad_fin,elect,dvecs,haml,(epoc_max-epoc_cnt))
        end if
      
        return

    end subroutine full_zs_gd


    subroutine zombie_alter(zstore,haml,elect,dvecs)!,neural_net)

        implicit none

        type(zombiest),dimension(:),allocatable,intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(grad)::grad_fin
        ! type(neural_network_layer),dimension(:)::neural_net
        integer::cnt,j
        integer::ierr=0
       
        if (errorflag .ne. 0) return

        call allocgrad(grad_fin,ndet,norb)
        grad_fin%prev_erg=ergcalc(haml%hjk,dvecs%d)
        grad_fin%grad_avlb=0
        grad_fin%ovrlp_grad_avlb=0
        epoc_cnt=1 !epoc counter
        if(rstrtflg.eq.'y')then 
            ierr=0
            open(unit=450,file='epoc.csv',status="old",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening epoc file to read in. ierr had value ", ierr
                errorflag=1
                return
            end if
            cnt=0
            do 
                read(450,*,iostat=ierr)
                if(ierr<0)then
                    exit
                else if(ierr>0)then
                    write(stderr,"(a,i0)") "Error in counting epocs. ierr had value ", ierr
                    errorflag=1
                    return
                else 
                    cnt=cnt+1
                end if
                
            end do
            close(450) 
            open(unit=450,file='epoc.csv',status="old",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening epoc file to read in. ierr had value ", ierr
                errorflag=1
                return
            end if
            do j=1,cnt-1
                read(450,*,iostat=ierr)
            end do 
            read(450,*,iostat=ierr)epoc_cnt
    
            close(450) 

            write(stderr,"(a,i0)") "Epoc read in as ", epoc_cnt
        
            call epoc_writer(grad_fin%prev_erg,epoc_cnt,0,0.0d0,1)
            epoc_cnt=epoc_cnt+1
        else
            call epoc_writer(grad_fin%prev_erg,0,0,0.0d0,0)
        end if

        allocate(picker(ndet-1),stat=ierr)
        if(ierr==0) allocate(chng_trk(ndet-1),stat=ierr)
        ! call omp_set_nested(.true.)
        if(epoc_cnt.lt.epoc_max)then
            picker=scramble(ndet-1)
            ! call full_zs_gd(zstore,elect,dvecs,haml,grad_fin) 
            ! call orbital_gd(zstore,grad_fin,elect,dvecs,haml,1000)
            ! call full_zs_gd(zstore,elect,dvecs,haml,grad_fin) 
            call orbital_gd(zstore,grad_fin,elect,dvecs,haml,epoc_max-epoc_cnt)!,neural_net)
            
        end if 

        !Brings phi values back within the normal 0-2pi range
        

        deallocate(picker,stat=ierr)
        deallocate(chng_trk,stat=ierr)
        return

    end subroutine zombie_alter

    ! Produces a random order for the ZS to be posisbly changed
    function scramble( number_of_values ) result(out)
        
        implicit none

        integer,intent(in)::number_of_values
        integer::out(number_of_values),array(number_of_values+1)
        integer::n,m,k,j,l,jtemp
        !real::r
        ! DOUBLE PRECISION, external::ZBQLU01

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
        ! DOUBLE PRECISION, external::ZBQLU01
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


    subroutine grad_do_haml_partial_transfer(in,haml,dvecs)

        implicit none
        type(grad_do),intent(inout)::in
        type(hamiltonian),intent(in)::haml
        type(dvector),intent(inout)::dvecs

        if (errorflag .ne. 0) return

        in%hjk(:,pick)= haml%hjk(:,pick); in%hjk(pick,:)=haml%hjk(pick,:)
        in%ovrlp(:,pick)= haml%ovrlp(:,pick); in%ovrlp(pick,:)=haml%ovrlp(pick,:)
        in%inv=haml%inv
        in%kinvh=haml%kinvh
        in%dvecs=dvecs
        return 

    end subroutine grad_do_haml_partial_transfer

    subroutine grad_do_haml_transfer(in,haml,zstore,dvecs)

        implicit none
        type(grad_do),intent(in)::in
        type(hamiltonian),intent(inout)::haml
        type(zombiest),intent(inout)::zstore
        type(dvector),intent(inout)::dvecs

        if (errorflag .ne. 0) return

        haml%hjk(:,pick)=in%hjk(:,pick); haml%hjk(pick,:)=in%hjk(pick,:)
        haml%ovrlp(:,pick)=in%ovrlp(:,pick); haml%ovrlp(pick,:)=in%ovrlp(pick,:)
        haml%inv=in%inv
        haml%kinvh=in%kinvh

        zstore=in%zom
        dvecs=in%dvecs
        return 

    end subroutine grad_do_haml_transfer

    subroutine haml_to_grad_do(haml,dvecs,out)

        implicit none 
        type(grad_do),intent(inout)::out
        type(hamiltonian),intent(in)::haml
        type(dvector),intent(in)::dvecs

        if (errorflag .ne. 0) return
        out%hjk=haml%hjk
        out%ovrlp=haml%ovrlp  
        out%inv=haml%inv
        out%kinvh=haml%kinvh
        out%dvecs%d=dvecs%d

        return 

    end subroutine 

    subroutine gram_ovrlp_var_temp(gramstore,wf_ovrlp,z1d,state,var)
        implicit none
        type(gram),dimension(:),intent(in)::gramstore
        real(wp),dimension(:,:,:),intent(inout)::wf_ovrlp
        type(zombiest),intent(in)::z1d
        integer,intent(in)::state,var
        integer::j,k

        if(errorflag.ne.0) return
        do j=1,state-1
            do k=1,ndet
                wf_ovrlp(j,k,var)=product(gramstore(state)%zstore(k)%val(1:norb)*&
                z1d%val(1:norb)+gramstore(state)%zstore(k)%val(1+norb:2*norb)*z1d%val(1+norb:2*norb))
              wf_ovrlp(j,var,k)=wf_ovrlp(j,k,var)
            end do
        end do

    end subroutine gram_ovrlp_var_temp

    ! subroutine nn_haml(zstore,neural_net,temp,pickorb,hnuc)
    !     implicit none
    !     type(neural_network_layer),dimension(:),intent(inout)::neural_net
    !     type(grad_do),intent(inout)::temp
    !     type(zombiest),dimension(:)::zstore
    !     real(dp),intent(in)::hnuc
    !     integer,intent(in)::pickorb
    !     real(wp),dimension((2*norb)+5)::input_features
    !     ! real(wp),dimension(5)::input_features
    !     integer::j
    !     integer, allocatable,dimension(:)::IPIV1
    !     real(dp),allocatable,dimension(:)::WORK1
    !     integer::ierr=0


    !     if (errorflag .ne. 0) return
        
    !     input_features(1)=temp%zom%val(pickorb+norb)
    !     input_features(2)=temp%zom%val(pickorb)
    !     input_features(4)=pickorb
    !     input_features(5)=0
    !     input_features((7+norb):)=temp%zom%val(1:norb)
    !     do j=1,ndet
    !         input_features(3)=temp%hjk(pick,j)
    !         if(j==pick)then
    !             input_features(5)=1
    !             temp%ovrlp(j,pick)=1.0d0
    !             input_features(6:(6+norb))=temp%zom%val(1:norb)
    !         else 
    !             input_features(5)=0
    !             input_features(6:(6+norb))=zstore(j)%val(1:norb)
    !             temp%ovrlp(j,pick)=product((zstore(j)%val(1:norb)*temp%zom%val(1:norb))+&
    !                                        (zstore(j)%val(1+norb:2*norb)*temp%zom%val(1+norb:2*norb)))
    !         end if
    !         temp%hjk(j,pick)=forward_result(neural_net,input_features)
    !     end do 
        
    !     temp%ovrlp(pick,:)=temp%ovrlp(:,pick)
    !     temp%hjk(:,pick)=temp%hjk(:,pick)+temp%ovrlp(:,pick)*hnuc
    !     temp%hjk(pick,:)= temp%hjk(:,pick)
    !     temp%inv=temp%ovrlp
       
    !     allocate(WORK1(ndet),IPIV1(ndet),stat=ierr)
    !     if (ierr/=0) then
    !         write(stderr,"(a,i0)") "Error in IPIV or WORK1 vector allocation . ierr had value ", ierr
    !         errorflag=1
    !     end if 

    !     Call dgetrf(ndet, ndet, temp%inv, ndet, IPIV1, ierr)
    !     if (ierr/=0) then
    !         write(stderr,"(a,i0)")"Error in DGETRF",ierr
    !     end if
    !     if (ierr==0) call dgetri(ndet,temp%inv,ndet,IPIV1,WORK1,ndet,ierr)
    !     if (ierr/=0) then
    !         write(stderr,"(a,i0)")"Error in DGETRF",ierr
    !     end if

    !     deallocate(WORK1,IPIV1,stat=ierr)
    !     if (ierr/=0) then
    !         write(stderr,"(a,i0)") "Error in IPIV or WORK1 vector deallocation . ierr had value ", ierr
    !         errorflag=1
    !     end if
       
    !     call DGEMM("N","N",ndet,ndet,ndet,1.d0,temp%inv,ndet,temp%hjk,ndet,0.d0,temp%kinvh,ndet)

    !     return

    ! end subroutine nn_haml


END MODULE gradient_descent