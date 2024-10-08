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
    real(wp),dimension(:),allocatable::orb_store
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

    subroutine orbital_gd(zstore,grad_fin,elect,dvecs,haml)

        implicit none 

        type(zombiest),dimension(:),allocatable,intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(grad),intent(inout)::grad_fin
        type(grad_do)::temp,thread
        integer::rjct_cnt,acpt_cnt,pickorb,loops,lralt_zs,acpt_cnt_2,lralt_extra2
        integer::j,n,p,chng_chng,tracker,lralt_extra,extra_flag,chng_chng2,reduc2
        integer,dimension(:),allocatable::chng_trk2,pickerorb
        real(wp)::t,erg_str,num_av,reduc,comp
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
        lralt_extra2=lr_loop_max !remove lower learning rates
        tracker=-1
        extra_flag=0
        p=70-norb
       
        reduc=reduc_in
        reduc2=0
       
        comp=grad_fin%prev_erg
        chng_chng=25  !blind_clone_num/4
        chng_chng2=blind_clone_num
        if(ndet.eq.ndet_max)then
            lr_loop_max=min_clone_lr
        end if
       
        call haml_to_grad_do(haml,dvecs,temp)
        if(gramflg.eq.'y')then
            ! if(epoc_cnt.eq.1)then
            !     gram_Store=gramnum 
            !     gramflg='n'
            ! end if 
            call imaginary_time(temp,ndet)
            grad_fin%prev_erg=temp%erg
            dvecs=temp%dvecs
        end if
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
                    if((abs(t*grad_fin%vars(pick,pickorb)).lt.reduc*1000).or.&
                        (abs(t*grad_fin%vars(pick,pickorb)).gt.3*pirl).or.&
                        (abs(t*grad_fin%vars(pick,pickorb)).lt.1.0d-12))then
                        write(stdout,'(1a)',advance='no') '!'
                        cycle
                    end if

                    thread%zom=zstore(pick)
                    temp=thread
                    temp%zom%phi(pickorb) = thread%zom%phi(pickorb)-(t*grad_fin%vars(pick,pickorb))
                    if(abs(temp%zom%phi(pickorb)).gt.2*pirl)then
                        write(stdout,'(1a)',advance='no') '!'
                        cycle
                    end if 
                    
                    call val_set(temp%zom,pickorb)
                    call he_full_row(temp,zstore,elect,ndet,pickorb)
                    call imaginary_time(temp,ndet)
     
                    !if((grad_fin%prev_erg-temp%erg.ge.1.0d-14))then
                    if((temp%erg.lt.grad_fin%prev_erg+reduc*(1**(-reduc2))))then
                    ! if((temp%erg.lt.grad_fin%prev_erg).or.((t.gt.0.1).and.(temp%erg.lt.grad_fin%prev_erg+reduc)))then
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
                    
                end if
            end do
            
           write(stdout,"(a,i0,a,f21.16,a,f10.5)") "Energy after epoch no. ",epoc_cnt,": ",grad_fin%prev_erg, "    Learning rate:",t
       
            if(acpt_cnt_2.gt.0)then
                do j=1,acpt_cnt_2
                    call zombiewriter(zstore(chng_trk(j)),chng_trk(j),zstore(chng_trk(j))%gram_num)
                end do
                call epoc_writer(grad_fin%prev_erg,epoc_cnt,t,chng_trk,0)
                epoc_cnt=epoc_cnt+1
                chng_chng2=chng_chng2-1
            else
                loops=loops-1
                rjct_cnt_global=rjct_cnt_global+1  
            end if 
            write(stdout,"(a,i0)") "Number of epochs until additional Zombie state, ",  chng_chng2
            print*,reduc,comp
            if((acpt_cnt_2.lt.(0.25*ndet)).and.(tracker.gt.-1))then
                if(lralt_zs.eq.lralt_extra)then
                    lralt_extra=lralt_extra+1
                ! else if(lralt_zs.eq.lralt_extra2)then
                !     lralt_extra2=lralt_extra2-1
                end if 
            end if 
        
            lralt_zs=lralt_zs+1
            chng_chng=chng_chng-1
            if(modulo(lralt_zs,2).eq.0)then
                reduc2=reduc2+1
            end if 

            if((chng_chng2.le.0).or.(rjct_cnt_global.gt.3*(lr_loop_max)))then
                if((abs(comp-grad_fin%prev_erg).lt.reduc*1000).or.(grad_fin%prev_erg.gt.comp))then
                    if((ndet.ge.ndet_max))then
                        reduc=reduc/10
                        ! if(reduc.lt.1.0d-13)then
                        !     reduc=0
                        ! end if 
                    end if
                else if(comp-grad_fin%prev_erg.gt.reduc*5000)then    
                    reduc=reduc*5
                    if(reduc.gt.1.0d-7)then
                        reduc=1.0d-7
                     end if
                end if
                if((ndet.lt.ndet_max))then
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
                    deallocate(chng_trk,stat=ierr)
                    allocate(chng_trk(ndet-1),stat=ierr)
                    do j=(ndet-ndet_increase+1),ndet
                        call zombiewriter(zstore(j),j,zstore(j)%gram_num)
                    end do
                    ! if(ndet.ge.ndet_max)then
                    !     blind_clone_num=blind_clone_num/2
                    ! end if
                else if(lr_loop_max.lt.min_clone_lr)then
                    lr_loop_max=lr_loop_max+1
                end if
                comp=grad_fin%prev_erg
                tracker=-1
                lralt_extra=0
                lralt_zs=0
                chng_chng=25  !blind_clone_num/4
                chng_chng2=blind_clone_num
                lralt_extra2=lr_loop_max
                reduc2=0
            else if((chng_chng.le.0))then
                lralt_zs=0
                reduc2=0
                lralt_extra=0
                chng_chng=25  !blind_clone_num/4
                lralt_extra2=lr_loop_max
                if(((comp-grad_fin%prev_erg).lt.reduc*1000).or.(grad_fin%prev_erg.gt.comp))then
                    reduc=reduc/10
                end if 
                comp=grad_fin%prev_erg
            end if 
           
            if((lralt_zs.gt.lralt_extra2))then
                picker=scramble(ndet-1)
                lralt_zs=lralt_extra
                extra_flag=0
                reduc2=0
                if((acpt_cnt_2.lt.((ndet)/3)).or.(tracker.lt.0))then
                    tracker=tracker+1
                end if    
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
        integer:: j,k
      
        tracker=-1
        extra_flag=1
        lralt_extra=0
        call alloczs(zstore_temp,ndet-ndet_increase)
        zstore_temp=zstore
        call dealloczs(zstore)
        call alloczs(zstore,ndet)
        zstore(1:(ndet-ndet_increase))=zstore_temp
        call dealloczs(zstore_temp)
        do j=(ndet-ndet_increase)+1,ndet
            ! if((ndet_increase.eq.1).and.(ndet.gt.3))then 
            !     zstore(j)%phi=orb_store
            ! else 
            if(zst.eq.'BB')then
                call biased_func(zstore(j))
            else
                do k=1,norb
                    zstore(j)%phi(k)=0.5*pirl*(ZBQLU01()) 
                end do
            end if 
            call val_set(zstore(j))
            zstore(j)%gram_num=zstore(1)%gram_num
        end do
        do j=(ndet-ndet_increase)+1,ndet
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

        if(gramflg.eq.'n')then
            grad_fin%prev_erg=ergcalc(haml%hjk,dvecs%d)
        end if
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
    
            ! close(450) 

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
            call orbital_gd(zstore,grad_fin,elect,dvecs,haml)

        end if 
      
        
       
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

        out=[(j,j=1,number_of_values)]
        array=[(j,j=1,number_of_values+1)]
        n=1; m=number_of_values
        do k=1,2
            do j=1,m+1
                l = n + FLOOR((m+1-n)*ZBQLU01())
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
        !real::r
        out=[(j,j=1,number_of_values)]
        
        n=1; m=number_of_values
        do k=1,2
            do j=1,m
                l = n + FLOOR((m-n)*ZBQLU01())
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

    ! subroutine gram_ovrlp_var_temp(gramstore,wf_ovrlp,z1d,state,var)
    !     implicit none
    !     type(gram),dimension(:),intent(in)::gramstore
    !     real(wp),dimension(:,:,:),intent(inout)::wf_ovrlp
    !     type(zombiest),intent(in)::z1d
    !     integer,intent(in)::state,var
    !     integer::j,k

    !     if(errorflag.ne.0) return
    !     do j=1,state-1
    !         do k=1,ndet
    !             wf_ovrlp(j,k,var)=product(gramstore(state)%zstore(k)%val(1:norb)*&
    !             z1d%val(1:norb)+gramstore(state)%zstore(k)%val(1+norb:2*norb)*z1d%val(1+norb:2*norb))
    !           wf_ovrlp(j,var,k)=wf_ovrlp(j,k,var)
    !         end do
    !     end do

    ! end subroutine gram_ovrlp_var_temp

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