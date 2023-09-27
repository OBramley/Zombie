MODULE gradient_descent

    use mod_types
    use globvars
    use alarrays
    use ham
    use imgtp
    use outputs
    use infnan_mod
    use zom 
   

    implicit none 
    real(wp)::alpha=0.4 ! learning rate reduction
    real(wp)::b=5.0D2   !starting learning rate
    integer::epoc_max=2000
    integer::epoc_cnt !epoc counter
    integer::loop_max=8 !max number of loops in gd
    integer::rjct_cnt_global=0
    integer::pick !Chosen zombie state
    integer,dimension(:),allocatable::picker
    integer,dimension(:),allocatable::chng_trk
   
    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators

    ! Subroutine that controls and calcualtes all of the hamiltonian variables 
    subroutine he_full_row(temp,zstore,elecs,size)

        implicit none 

        type(grad_do),intent(inout)::temp
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::size
        integer, allocatable,dimension(:)::IPIV1
        real(dp),allocatable,dimension(:)::WORK1
        integer::ierr=0


        if (errorflag .ne. 0) return
       
        call haml_ovrlp_column(temp,zstore,ndet,elecs,pick)
        ! call haml_ovrlp_column_gpu(temp_d,zstore_d,size,elecs_d,pick)
        ! temp=temp_d
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

    
    subroutine grad_calculate(haml,dvec,zstore,grad_fin,orb)

        implicit none 

        type(grad),intent(inout)::grad_fin
        type(hamiltonian),intent(inout)::haml
        type(dvector),intent(inout)::dvec
        type(zombiest),dimension(:),intent(in)::zstore
        real(wp)::ham_c_d
        real(wp),dimension(norb,ndet)::ovrlp_dx
        real(wp),dimension(ndet)::temp
        integer,intent(in)::orb
        integer::j,k
       
        if (errorflag .ne. 0) return

      
        if(orb==0)then
            if(grad_fin%grad_avlb(1,pick)==0)then
                call DGEMV("N",ndet,ndet,1.d0,haml%hjk,ndet,dvec%d,1,0.d0,temp,1)
                ham_c_d=(-dvec%d_o_d/(dvec%norm*dvec%norm*dvec%norm))*dot_product(temp,dvec%d_1)
                ovrlp_dx=1.0d0
                !$omp parallel private(j,k,temp) shared(zstore,ovrlp_dx,pick,ndet,norb,dvec,ham_c_d,grad_fin) 
                !$omp do collapse(2)
                do j=1,ndet
                    do k=1,norb
                        ovrlp_dx(k,j)=haml%ovrlp(j,pick)/(zstore(pick)%val(k)%x*zstore(j)%val(k)%x+&
                        zstore(pick)%val(k+norb)%x*zstore(j)%val(norb+k)%x)

                        ovrlp_dx(k,j)=ovrlp_dx(k,j)*(zstore(pick)%val(k)%dx(k)*zstore(j)%val(k)%x+&
                        zstore(pick)%val(k+norb)%dx(k)*zstore(j)%val(norb+k)%x)
                    end do
                end do 
                !$omp end do
                !$omp critical
                ovrlp_dx(:,pick)=0
                !$omp end critical
                !$omp do
                do k=1,norb
                    temp=ovrlp_dx(k,:)*dvec%d_1
                    temp(pick)=dot_product(ovrlp_dx(k,:),dvec%d_1)
                    grad_fin%vars(pick,k)=dot_product(temp,dvec%d_1)*ham_c_d
                end do
                !$omp end do
                !$omp end parallel
                grad_fin%grad_avlb(:,pick)=1
            else 
                return
            end if 
        else
            if(grad_fin%grad_avlb(orb,pick)==0)then
                call DGEMV("N",ndet,ndet,1.d0,haml%hjk,ndet,dvec%d,1,0.d0,temp,1)
                ham_c_d=(-dvec%d_o_d/(dvec%norm*dvec%norm*dvec%norm))*dot_product(temp,dvec%d_1)
                ovrlp_dx=1.0d0
                !!$omp parallel do reduction(*: ovrlp_dx)collapse(2) private(j,l) shared(zstore,pick,ndet,norb,orb)
                do j=1,ndet
                    ovrlp_dx(orb,j)=haml%ovrlp(j,pick)/(zstore(pick)%val(orb)%x*zstore(j)%val(orb)%x+&
                    zstore(pick)%val(orb+norb)%x*zstore(j)%val(norb+orb)%x)

                    ovrlp_dx(orb,j)=ovrlp_dx(orb,j)*(zstore(pick)%val(orb)%dx(orb)*zstore(j)%val(orb)%x+&
                    zstore(pick)%val(orb+norb)%dx(orb)*zstore(j)%val(norb+orb)%x)
                end do 
                !!$omp end parallel do
                ovrlp_dx(orb,pick)=0
                temp=ovrlp_dx(orb,:)*dvec%d_1
                temp(pick)=dot_product(ovrlp_dx(orb,:),dvec%d_1)
                grad_fin%vars(pick,orb)=dot_product(temp,dvec%d_1)*ham_c_d
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

    subroutine orbital_gd(zstore,grad_fin,elect,dvecs,haml,maxloop) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(grad),intent(inout)::grad_fin
        integer,intent(in)::maxloop
        type(grad_do)::temp,thread
        integer::rjct_cnt,acpt_cnt,pickorb,loops,lralt_zs,acpt_cnt_2
        real(dp)::t,erg_str
        integer::j,n,p
        integer,dimension(:),allocatable::chng_trk2,pickerorb
        integer::ierr=0
        
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
        acpt_cnt_2=0
        loops=0
        p=70-norb
     
    
        call haml_to_grad_do(haml,dvecs,temp)
        thread=temp

        do while(rjct_cnt.lt.(norb*100))
            loops=loops+1
           
            do p=1, ((norb-8)/2)
                write(stdout,'(1a)',advance='no') ' '
            end do 
            write(stdout,"(a)",advance='no') 'Progress'
            do p=1, ((norb-8)/2)
                write(stdout,'(1a)',advance='no') ' '
            end do 
            write(stdout,"(a)") ' | Zombie state | Previous Energy     | Energy after Gradient Descent steps   | Orbitals altered '
       
            chng_trk=0
            acpt_cnt_2=0  
        
            do j=1,ndet-1
               
                erg_str=grad_fin%prev_erg
                pick=picker(j)
                chng_trk2=0
                acpt_cnt=0
                pickerorb=scramble_norb(norb)
                call haml_to_grad_do(haml,dvecs,thread)

                do n=1,norb
                    pickorb=pickerorb(n)
                    call grad_calculate(haml,dvecs,zstore,grad_fin,pickorb)
                
                    thread%zom=zstore(pick)
            
                    do lralt_zs=1,loop_max
                        temp=thread
                        t=b*((alpha/(lralt_zs))**(lralt_zs-1))
                    
                        temp%zom%phi(pickorb)%x = thread%zom%phi(pickorb)%x-(t*grad_fin%vars(pick,pickorb))
                        call val_set(temp%zom,pickorb)
                        
                        call he_full_row(temp,zstore,elect,ndet)
                    
                        call imaginary_time_erg(temp,ndet)
                        
                        if((temp%erg .lt. grad_fin%prev_erg))then
                            acpt_cnt=acpt_cnt+1
                            chng_trk2(acpt_cnt)=pickorb
                            rjct_cnt=0
                            grad_fin%grad_avlb=0
                            grad_fin%prev_erg=thread%erg
                            rjct_cnt_global=0
                            call grad_do_haml_transfer(temp,haml,zstore(pick),dvecs)
                            thread=temp
                            grad_fin%grad_avlb=0
                            grad_fin%prev_erg=thread%erg
                            exit
                        end if
                    end do
                    
                    write(stdout,'(1a)',advance='no') '|'
                    flush(6)
                end do
               
                if(acpt_cnt.gt.0)then
                    write(stdout,"(a,i3,a,f21.16,a,f21.16,a,*(i0:','))")'  ', pick,'          ', &
                    erg_str,'             ',grad_fin%prev_erg,'          ',chng_trk2(1:acpt_cnt) 
                    
                    call zombiewriter(zstore(pick),pick,0)
                    acpt_cnt_2=acpt_cnt_2+1
                    chng_trk(acpt_cnt_2)=pick
                else 
                    write(stdout,"(a,i3,a,f21.16,a,f21.16,a,i0)")'  ',pick,'          ', &
                    erg_str,'             ',grad_fin%prev_erg,'          ',0
                   
                    rjct_cnt=rjct_cnt+1
                    rjct_cnt_global=rjct_cnt_global+1   
                end if
                
            end do
           
            write(stdout,"(a,i0,a,f21.16)") "Energy after epoch no. ",epoc_cnt,": ",grad_fin%prev_erg
          
            if(acpt_cnt_2.gt.0)then
                call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,0)
                epoc_cnt=epoc_cnt+1
            else
                loops=loops-1
            end if 

            picker=scramble(ndet-1)
           
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


    subroutine full_zs_gd(zstore,elect,dvecs,haml,grad_fin) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(grad),intent(inout)::grad_fin
        type(grad_do)::temp,thread
        integer::acpt_cnt,lralt_temp,orb_cnt,min_idx,j
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
       
        if(epoc_cnt.eq.1)then
            orb_cnt=2
        else
            orb_cnt=2
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
            do j=1,(ndet-1)
                
                pick=picker(j)
                
             
                call grad_calculate(haml,dvecs,zstore,grad_fin,0)
              
                do lralt_temp=1,loop_max
                    temp=thread
                    t=b*((alpha/(lralt_temp))**(lralt_temp-1))
                    ! t=b*(alpha**(lralt_temp-1))
                    temp%zom=zstore(pick)
                    temp%zom%phi=zstore(pick)%phi
                    temp%zom%phi%x=zstore(pick)%phi%x-(t*grad_fin%vars(pick,:))
                    call val_set(temp%zom)
                    call he_full_row(temp,zstore,elect,ndet)
                    call imaginary_time_erg(temp,ndet)
                    if(temp%erg .lt. grad_fin%prev_erg)then
                        t=b*(alpha**(min_idx-1))
                        acpt_cnt=acpt_cnt+1
                        chng_trk(acpt_cnt)=pick
                        lr_chng_trk(acpt_cnt)=t
                        erg_chng_trk(acpt_cnt)=thread%erg
                        !$acc update device(zstore(pick))
                        call grad_do_haml_transfer(thread,haml,zstore(pick),dvecs)
                        call zombiewriter(zstore(pick),pick,0)
                        rjct_cnt_global=0
                
                    write(stdout,"(a,i3,a,f21.16,a,f21.16,a,f21.16,a,i3,a,i3)") '       ', pick,'              ', &
        grad_fin%prev_erg,'               ',thread%erg,'             ',t,'        ',acpt_cnt,'          ',rjct_cnt_global
            
                        grad_fin%grad_avlb=0
                        grad_fin%prev_erg=thread%erg
                        call grad_do_haml_partial_transfer(temp,haml,dvecs)
                        thread=temp
                        min_idx=-1
                    end if
               
                end do
             
                if(min_idx .eq. -1)then
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
           
            orb_cnt=orb_cnt-1
          
            if(rjct_cnt_global.ge.(ndet+1))then
                call orbital_gd(zstore,grad_fin,elect,dvecs,haml,1)
                call haml_to_grad_do(haml,dvecs,thread)
                orb_cnt=orb_cnt+1
            else if((orb_cnt.le.0))then
                call orbital_gd(zstore,grad_fin,elect,dvecs,haml,2)
                call haml_to_grad_do(haml,dvecs,thread)
                orb_cnt=5
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
            call orbital_gd(zstore,grad_fin,elect,dvecs,haml,(epoc_max-epoc_cnt))
        end if
      
        return

    end subroutine full_zs_gd


    subroutine zombie_alter(zstore,haml,elect,dvecs)

        implicit none

        type(zombiest),dimension(:),intent(inout)::zstore
        type(elecintrgl),intent(in)::elect
        type(dvector),intent(inout)::dvecs
        type(hamiltonian),intent(inout)::haml
        type(grad)::grad_fin
        integer::cnt,j
        integer::ierr=0
       
        if (errorflag .ne. 0) return

        call allocgrad(grad_fin,ndet,norb)
        grad_fin%prev_erg=ergcalc(haml%hjk,dvecs%d)
      
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
        ! call omp_set_nested(.TRUE.)

        allocate(picker(ndet-1),stat=ierr)
        if(ierr==0) allocate(chng_trk(ndet-1),stat=ierr)
        ! call omp_set_nested(.true.)
        if(epoc_cnt.lt.epoc_max)then
            picker=scramble(ndet-1)
            ! call full_zs_gd(zstore,elect,dvecs,haml,grad_fin) 
            call orbital_gd(zstore,grad_fin,elect,dvecs,haml,1000)
            call full_zs_gd(zstore,elect,dvecs,haml,grad_fin) 
            !call orbital_gd(zstore,grad_fin,elect,dvecs,haml,100)
            
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


    subroutine grad_do_haml_partial_transfer(in,haml,dvec)

        implicit none
        type(grad_do),intent(inout)::in
        type(hamiltonian),intent(in)::haml
        type(dvector),intent(inout)::dvec

        if (errorflag .ne. 0) return

        in%hjk(:,pick)= haml%hjk(:,pick); in%hjk(pick,:)=haml%hjk(pick,:)
        in%ovrlp(:,pick)= haml%ovrlp(:,pick); in%ovrlp(pick,:)=haml%ovrlp(pick,:)
        in%inv=haml%inv
        in%kinvh=haml%kinvh
        in%dvec=dvec
        return 

    end subroutine grad_do_haml_partial_transfer

    subroutine grad_do_haml_transfer(in,haml,zstore,dvec)

        implicit none
        type(grad_do),intent(in)::in
        type(hamiltonian),intent(inout)::haml
        type(zombiest),intent(inout)::zstore
        type(dvector),intent(inout)::dvec

        if (errorflag .ne. 0) return

        haml%hjk(:,pick)=in%hjk(:,pick); haml%hjk(pick,:)=in%hjk(pick,:)
        haml%ovrlp(:,pick)=in%ovrlp(:,pick); haml%ovrlp(pick,:)=in%ovrlp(pick,:)
        haml%inv=in%inv
        haml%kinvh=in%kinvh

        zstore=in%zom
        dvec=in%dvec
        return 

    end subroutine grad_do_haml_transfer

    subroutine haml_to_grad_do(haml,dvec,out)

        implicit none 
        type(grad_do),intent(inout)::out
        type(hamiltonian),intent(in)::haml
        type(dvector),intent(in)::dvec

        if (errorflag .ne. 0) return
        out%hjk=haml%hjk
        out%ovrlp=haml%ovrlp  
        out%inv=haml%inv
        out%kinvh=haml%kinvh
        out%dvec%d=dvec%d

        return 

    end subroutine 



END MODULE gradient_descent