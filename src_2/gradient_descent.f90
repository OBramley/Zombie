MODULE gradient_descent

    use globvars
    use alarrays
    use ham
    use operators
    use grad_d
    use imgtp
    use outputs
    use infnan_mod

    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators
    subroutine he_full_row(haml,zstore,elecs,diff_state,size,occupancy_2an,occupancy_an_cr,occupancy_an)
        
        implicit none
        type(hamiltonian), intent(inout)::haml
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::size,diff_state
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an,occupancy_an_cr
        integer,dimension(:,:,:),intent(in)::occupancy_an
        complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
        complex(kind=8),allocatable, dimension(:,:,:)::z2l
        integer::j,k,l,m,ierr,equal
        complex(kind=8)::h1etot, h2etot
        real(kind=8),dimension(norb)::h1etot_diff_bra,h2etot_diff_bra
        real(kind=8),dimension(norb)::h1etot_diff_ket,h2etot_diff_ket
        integer, allocatable,dimension(:)::IPIV1
        complex(kind=8),allocatable,dimension(:)::WORK1
        
        if (errorflag .ne. 0) return
        ierr = 0

        allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        h1etot=cmplx(0.0,0.0)
        h2etot=cmplx(0.0,0.0)
 
    
        !$omp parallel shared(zstore,z1jk, z2l) private(j,k,l)
       
        !$omp do
        do l=1, norb
            z2l(l,1,:)=zstore(diff_state)%sin(:)
            z2l(l,2,:)=zstore(diff_state)%cos(:)
            z2l(l,2,l)=z2l(l,1,l)
            z2l(l,1,l)=cmplx(0.0,0.0)
        end do
        !$omp end do

        !$omp do
        do j=1, norb
            do k=1, norb
                z1jk(j,k,:,:)=z2l(j,:,:)
                z1jk(j,k,2,k)=z1jk(j,k,1,k)
                z1jk(j,k,1,k)=cmplx(0.0,0.0)
            end do
        end do
        !$omp end do
        !$omp end parallel
        z1jk=z1jk*occupancy_2an

        do m=1,size
      
            h1etot=cmplx(0.0,0.0)
            h2etot=cmplx(0.0,0.0)
         
            !$omp parallel shared(z2l)
            !$omp do
            do l=1, norb
                z2l(l,1,:)=zstore(m)%sin(:)
                z2l(l,2,:)=zstore(m)%cos(:)
                z2l(l,2,l)=z2l(l,1,l)
                z2l(l,1,l)=cmplx(0.0,0.0)
            end do
            !$omp end do
            !$omp end parallel
            
            equal=9
           
            if(m.eq.diff_state)then
                call one_elec_part(zstore(diff_state),z2l,h1etot,occupancy_an_cr,&
                                    elecs%h1ei,h1etot_diff_bra,h1etot_diff_ket,zstore(m),equal)
            else 
                call one_elec_part(zstore(diff_state),z2l,h1etot,occupancy_an_cr,&
                                    elecs%h1ei,h1etot_diff_bra,h1etot_diff_ket,zstore(m),equal)
            end if
           
            z2l=z2l*occupancy_an
            !$omp flush(z2l)
            if(m.eq.diff_state)then
                call two_elec_part(zstore(diff_state),z1jk,z2l,h2etot,occupancy_2an,occupancy_an,&
                                            elecs%h2ei,h2etot_diff_bra,h2etot_diff_ket,zstore(m),equal)
            else 
                call two_elec_part(zstore(diff_state),z1jk,z2l,h2etot,occupancy_2an,occupancy_an,&
                                            elecs%h2ei,h2etot_diff_bra,h2etot_diff_ket,zstore(m),equal)
            end if
            
            haml%ovrlp(diff_state,m)=overlap(zstore(diff_state),zstore(m))
            haml%ovrlp(m,diff_state)= haml%ovrlp(diff_state,m)
            haml%hjk(diff_state,m)=h1etot+h2etot+(elecs%hnuc*haml%ovrlp(diff_state,m))
            haml%hjk(m,diff_state)=haml%hjk(diff_state,m)
           
        end do
       
        deallocate(z1jk,stat=ierr)
        if(ierr==0) deallocate(z2l,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        haml%inv=haml%ovrlp
        allocate(IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
            errorflag=1
        end if 
        
        if (ierr==0) allocate(WORK1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
            errorflag=1
        end if   

        if (ierr==0) call ZGETRF(size,size,haml%inv,size,IPIV1,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if
        if (ierr==0) call ZGETRI(size,haml%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if

        if (ierr==0) deallocate(IPIV1,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV vector deallocation . ierr had value ", ierr
            errorflag=1
        end if

        if (ierr==0) deallocate(WORK1,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in WORK vector deallocation . ierr had value ", ierr
            errorflag=1
        end if

        !$omp parallel
        !$omp workshare
        haml%kinvh=matmul(haml%inv,haml%hjk)
        !$omp end workshare
        !$omp end parallel

        return

    end subroutine he_full_row

    subroutine gradient_row(haml,zstore,elecs,diff_state,size,occupancy_2an,occupancy_an_cr,occupancy_an)

        implicit none

        type(hamiltonian), intent(inout)::haml
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::size,diff_state
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an,occupancy_an_cr
        integer,dimension(:,:,:),intent(in)::occupancy_an
        complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
        complex(kind=8),allocatable, dimension(:,:,:)::z2l
        complex(kind=8),dimension(2,norb)::zomt
        real(kind=8),dimension(norb)::h1etot_diff_bra,h2etot_diff_bra
        real(kind=8),dimension(ndet,ndet,norb)::temp2
        real(kind=8),dimension(2,norb)::overlap_diff,h1etot_diff,h2etot_diff
        integer,dimension(2,norb)::occupancy
        integer::j,k,l,m,ierr, jspin
        


        if (errorflag .ne. 0) return
        ierr = 0

        allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        h1etot_diff_bra=0.0
        h2etot_diff_bra=0.0
       
        !$omp parallel shared(zstore,z1jk, z2l) private(j,k,l)
       
        !$omp do
        do l=1, norb
            z2l(l,1,:)=zstore(diff_state)%sin(:)
            z2l(l,2,:)=zstore(diff_state)%cos(:)
            z2l(l,2,l)=z2l(l,1,l)
            z2l(l,1,l)=cmplx(0.0,0.0)
        end do
        !$omp end do
    
        !$omp do
        do j=1, norb
            do k=1, norb
                z1jk(j,k,:,:)=z2l(j,:,:)
                z1jk(j,k,2,k)=z1jk(j,k,1,k)
                z1jk(j,k,1,k)=cmplx(0.0,0.0)
            end do
        end do
        !$omp end do
        !$omp end parallel
        z1jk=z1jk*occupancy_2an

        do m=1,size

            h1etot_diff_bra=0.0
            h2etot_diff_bra=0.0
         
            !$omp parallel shared(z2l)
            !$omp do
            do l=1, norb
                z2l(l,1,:)=zstore(m)%sin(:)
                z2l(l,2,:)=zstore(m)%cos(:)
                z2l(l,2,l)=z2l(l,1,l)
                z2l(l,1,l)=cmplx(0.0,0.0)
            end do
            !$omp end do
            !$omp end parallel

            h1etot_diff=0.0
            !$omp parallel private(j,k,zomt) shared(elecs,occupancy_an_cr,z2l,zstore,h1etot_diff)
            !$omp do schedule(dynamic) reduction(+:h1etot_diff)
            do j=1, norb
                do k=1, norb
                    zomt(:,:)=z2l(j,:,:)
                    zomt(1,k)=zomt(2,k)
                    zomt(2,k)=cmplx(0.0,0.0)
                    zomt=zomt*occupancy_an_cr(j,k,:,:)
                    if(diff_state.eq.m)then 
                        h1etot_diff = h1etot_diff+diff_overlap_cran(zstore(diff_state),zstore(m),&
                                        1,zomt,j,k,occupancy_an_cr(j,k,:,:))*elecs%h1ei(j,k)
                    else 
                        h1etot_diff = h1etot_diff+diff_overlap_cran(zstore(diff_state),zstore(m),&
                                        2,zomt,j,k,occupancy_an_cr(j,k,:,:))*elecs%h1ei(j,k)
                    end if
                end do
            end do
            !$omp end do
            !$omp end parallel
            h1etot_diff_bra=h1etot_diff(1,:)

            h2etot_diff=0.0
            do j=1, norb
                if(modulo(j,2)==0)then
                    jspin=2
                else
                    jspin=1
                end if
                do k=1, norb
                    if(j.eq.k)then
                        cycle
                    end if
                    do l=jspin, norb, 2
                        occupancy=occupancy_2an(j,k,:,:)*occupancy_an(l,:,:)
                        if(diff_state.eq.m)then 
                            h2etot_diff = h2etot_diff + z_an_z3_diff(z1jk(j,k,:,:),z2l(l,:,:),elecs%h2ei(j,k,l,:),1,&
                                occupancy,j,k,l,zstore(diff_state),zstore(m))
                        else
                            h2etot_diff = h2etot_diff + z_an_z3_diff(z1jk(j,k,:,:),z2l(l,:,:),elecs%h2ei(j,k,l,:),2,&
                            occupancy,j,k,l,zstore(diff_state),zstore(m))
                        end if
                    end do
                end do
            end do
            h2etot_diff_bra = h2etot_diff(1,:)*0.5

            haml%diff_hjk(diff_state,m,:)=h1etot_diff_bra+h2etot_diff_bra
            if(m.eq.diff_state)then
                haml%diff_ovrlp(diff_state,m,:) = 0
            else
                overlap_diff = diff_overlap(zstore(diff_state),zstore(m),2)
                haml%diff_ovrlp(diff_state,m,:) = overlap_diff(1,:)
            end if

        end do  

        deallocate(z1jk,stat=ierr)
        if(ierr==0) deallocate(z2l,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        j=diff_state 
        do k=1, ndet
            do l=1, ndet
                if(l.eq.j)then
                    temp2(k,j,:)=matmul(REAL(haml%inv(k,:)),haml%diff_ovrlp(j,:,:))
                else
                    temp2(k,l,:)=real(haml%inv(k,l))*haml%diff_ovrlp(j,l,:)
                end if
            end do
        end do
        do k=1, ndet
            do l=1, ndet
                haml%diff_invh(j,k,l,:)=matmul(transpose(temp2(k,:,:)),real(haml%kinvh(:,l)))*(-1)
            end do
        end do

        return

    end subroutine gradient_row

    subroutine zombie_alter(zstore,grad_fin,haml,elect,en,dvecs,chng_trk)

        implicit none

        type(zombiest),dimension(:),intent(inout)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml
        integer,dimension(:),intent(inout)::chng_trk
        type(zombiest),dimension(:),allocatable::temp_zom
        type(hamiltonian)::temp_ham
        integer,allocatable,dimension(:,:,:)::occupancy_an
        integer,allocatable,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an

        integer::lralt,rjct_cnt,epoc_cnt,next,acpt_cnt,acpt1,acpt2,pick,rpstk,pickorb,rjct_cnt2
        real(kind=8)::alpha,b,newb,t,fxtdk,l2_rglrstn
        integer(kind=8)::temp_int1,temp_int2
        integer,dimension(ndet-1)::picker,rsrtpass

        integer,dimension(ndet,norb)::lralt_zs
        integer,dimension(norb)::pickerorb
        real(kind=8),dimension(ndet,norb)::alpha_zs,newb_zs
        integer,dimension(norb)::chng_trk2
       
        
        logical::nanchk
        character(len=4)::ergerr
        integer::ierr,j,k,l,m,n
        
        
        if (errorflag .ne. 0) return

        allocate(occupancy_an(norb,2,norb),stat=ierr)
        if(ierr==0) allocate(occupancy_2an(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(occupancy_an_cr(norb,norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in occupancy vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        occupancy_an=1
        occupancy_2an=1
        occupancy_an_cr=1

        !$omp parallel private(j,k,l) shared(occupancy_an,occupancy_2an,occupancy_an_cr)
        !$omp do
        do j=1,norb
            occupancy_an(j,1,j)=0
            do l=j-1, 1, -1
                occupancy_an(j,1,l)=-1
            end do
        end do
        !$omp end do
        !$omp barrier
        !$omp do
        do j=1,norb
            do k=1, norb
                occupancy_2an(j,k,:,:)=occupancy_an(j,:,:)
                occupancy_2an(j,k,2,k)=occupancy_2an(j,k,1,k)
                occupancy_2an(j,k,1,k)=0
                occupancy_an_cr(j,k,:,:)=occupancy_an(j,:,:)
                occupancy_an_cr(j,k,1,k)=occupancy_an_cr(j,k,2,k)
                occupancy_an_cr(j,k,2,k)=0
                do l=k-1,1,-1
                    occupancy_2an(j,k,1,l)=occupancy_2an(j,k,1,l)*(-1)
                    occupancy_an_cr(j,k,1,l)=occupancy_an_cr(j,k,1,l)*(-1)
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    
        alpha=0.1  ! learning rate reduction
        lralt=0    ! power alpha is raised to  
        b=1.0D0 !starting learning rate
        newb=b
        epoc_cnt=0 !epoc counter
        chng_trk=0 !stores which if any ZS changed
        rjct_cnt=0 !tracks how many rejections 
        t=newb*(alpha**lralt) !learning rate
        acpt_cnt=0  !counts how many ZS have been changed
        l2_rglrstn=0.01 !L2 regularisation lambda paramter
        rsrtpass=0
        rpstk=0 !checks for computer error
        call alloczs(temp_zom,ndet)
        call allocham(temp_ham,ndet,norb)

        do j=1, ndet-1
            picker(j)=j+1
        end do
       
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

            call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,1)
            rsrtpass=1
            ierr=0
        end if
        

        do while(t.gt.(1.0d-10))
            t=newb*(alpha**lralt)
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent step &
                           &   | Learning rate | Acceptance count | Rejection count'
            do j=1,(ndet-1)
                pick=picker(j)
                
                do k=1,norb
                    nanchk=.false. 
                    if(is_nan(grad_fin%vars(pick,k)).eqv..true.)then
                        call gradient_row(haml,zstore,elect,pick,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                        dvecs(1)%d=cmplx(0.0,0.0)
                        en%erg=0
                        en%t=0
                        call imgtime_prop(dvecs,en,temp_ham,pick)
                        call final_grad(dvecs(1),haml,grad_fin,pick)
                        grad_fin%grad_avlb(pick)=1
                        exit 
                    end if 
                end do
                ! Setup temporary zombie state
                temp_zom=zstore
                temp_zom(pick)%phi(:)=zstore(pick)%phi(:)-(t*(grad_fin%vars(pick,:)))
                temp_zom(pick)%phi(:)=temp_zom(pick)%phi(:)+l2_rglrstn*((grad_fin%vars(pick,:))*grad_fin%vars(pick,:))
               
                do k=1,norb
                    if(is_nan(temp_zom(pick)%phi(k)).eqv..true.)then
                        temp_zom(pick)%phi=zstore(pick)%phi
                        exit
                    end if
                end do
                temp_zom(pick)%sin=sin(cmplx(temp_zom(pick)%phi,0.0d0,kind=8))
                temp_zom(pick)%cos=cos(cmplx(temp_zom(pick)%phi,0.0d0,kind=8))
                temp_ham=haml
                call he_full_row(temp_ham,temp_zom,elect,pick,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                ! Imaginary time propagation for back tracing
                dvecs(1)%d=cmplx(0.0,0.0)
                en%erg=0
                en%t=0
                call imgtime_prop(dvecs,en,temp_ham,0)
                fxtdk=real(en%erg(1,timesteps+1))
                temp_int2=int(fxtdk*(1.0d13),kind=8)
                temp_int1=int((grad_fin%prev_erg*1.0d13),kind=8)
                if(is_nan(fxtdk).eqv..true.)then 
                    nanchk=.true.
                    ergerr='NaN ' 
                else if(is_posinf(fxtdk).eqv..true.)then
                    nanchk=.true.
                    ergerr='+NaN'  
                else if(is_neginf(fxtdk).eqv..true.)then
                    nanchk=.true.
                    ergerr='-NaN'  
                end if
                
                write(6,"(a,i0,a,f21.16,a,f21.16,a,f12.10,a,i0,a,i0)") '       ', pick,'              ', &
                grad_fin%prev_erg,'               ',fxtdk,'            ',t,'          ',acpt_cnt,'                 ',rjct_cnt
                if(nanchk.eqv..false.)then
                    ! Check if energy is lower and accept or reject
                    if(fxtdk.lt.grad_fin%prev_erg)then
                        acpt_cnt=acpt_cnt+1
                        zstore=temp_zom
                        zstore(pick)%update_num=zstore(pick)%update_num+1
                        call zombiewriter(zstore(pick),pick,rsrtpass(pick))
                        rsrtpass(pick)=0
                        chng_trk(j)=pick
                        haml=temp_ham
                        rjct_cnt=0
                        grad_fin%grad_avlb=0
                        grad_fin%vars=0.0
                        haml%diff_hjk=0
                        haml%diff_invh=0
                        haml%diff_ovrlp=0
                        grad_fin%prev_erg=fxtdk
                        if(lralt.gt.0)then 
                            lralt=lralt-1
                        end if
                        t=newb*(alpha**lralt)
                    else 
                        rjct_cnt=rjct_cnt+1
                        if(rjct_cnt.eq.(ndet-1))then 
                            lralt=lralt+1
                            t=newb*(alpha**lralt)
                            rjct_cnt=0
                        end if
                    end if
                else 
                    write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                                               " for zombie state ", pick, ", on epoc ", epoc_cnt 
                end if

                ! if(temp_int1.eq.temp_int2)then 
                !     rpstk=rpstk+1 
                !     if(rpstk.eq.1)then 
                !         lralt=lralt+1 
                !     end if
                ! else 
                !     rpstk=0 
                ! end if

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
                    next=picker(j+1)
                end if
                
             
                if(grad_fin%grad_avlb(next).eq.1)then  
                    do k=1,norb
                        if(isnan(grad_fin%vars(next,k)).eqv..true.)then
                            grad_fin%grad_avlb(next)=0
                            exit 
                        end if 
                    end do
                end if
                !Set up gradients for next pass
                if(grad_fin%grad_avlb(next).eq.0)then 
                    call gradient_row(haml,zstore,elect,next,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                    do k=1,norb 
                        if(is_nan(grad_fin%vars(next,k)).eqv..true.)then
                           
                            nanchk=.true. 
                            grad_fin%vars=0 
                            grad_fin%grad_avlb=0
                            do l=1 ,10
                                call gradient_row(haml,zstore,elect,next,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                                do m=1,norb
                                    if(is_nan(grad_fin%vars(next,m)).eqv..true.)then 
                                        exit 
                                    end if 
                                    nanchk=.false. 
                                end do 
                                if(nanchk.eqv..false.)then 
                                    exit 
                                end if 
                            end do  
                            if(nanchk.eqv..true.)then 
                                errorflag=1 
                                write(0,"(a,i0)") "Error in gradient calculation for zombie state number ", next 
                                return 
                            else 
                                exit 
                            end if 
                        end if 
                    end do
                    grad_fin%grad_avlb(next)=1
                    dvecs(1)%d=cmplx(0.0,0.0)
                    en%erg=0
                    en%t=0
                    call imgtime_prop(dvecs,en,temp_ham,next)
                    call final_grad(dvecs(1),haml,grad_fin,next)
                end if

                
            
            end do

            call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,0)
            write(6,"(a,i0,a,f21.16,a,f12.10,a,i0,a)") "Energy after epoc no. ",epoc_cnt,": ", &
                grad_fin%prev_erg, ". The current learning rate is: ",t, ". ", acpt_cnt, " Zombie state(s) altered."

            rpstk=0 
            ! Increases learnign rate after each epoc if acceptance has occured
            ! do j=1,ndet-1
            !     if(chng_trk(j).ne.0)then
            !         ! rjct_cnt=0
            !         if(lralt.gt.0)then 
            !             lralt=lralt-1
            !         end if
            !         t=newb*(alpha**lralt)
            !         if(lralt.eq.0)then 
            !             exit 
            !         end if 
            !     end if 
            ! end do 
            if((acpt_cnt.gt.1).and.(newb.lt.b))then
                newb=(((newb+b)/2)+b)/2
                write(6,"(a,f12.10)") "b in learning rate equation t=b*(alpha^x) raised to ",newb 
                lralt=0    
                rjct_cnt=0
                acpt_cnt=0
            end if

            ! If only a single ZS is being altered the learning rate is lowered to allow others to be changed
            if(acpt_cnt.eq.1)then   
                do j=1,ndet-1
                    if(chng_trk(j).ne.0)then 
                        if(acpt1.eq.0)then 
                            acpt1=chng_trk(j)
                            acpt_cnt=0
                        else 
                            if(chng_trk(j).ne.acpt1)then 
                                acpt1=chng_trk(j)
                                acpt2=0
                                acpt_cnt=0 
                            else if(chng_trk(j).eq.acpt1)then
                                acpt2=acpt2+1
                                acpt_cnt=0
                                if(acpt2.gt.3)then
                                    if((t.gt.9.9d-5))then
                                        lralt=lralt+1
                                    end if
                                end if 
                            end if 
                            acpt_cnt=0
                        end if
                        exit
                    end if
                end do
            else if(acpt_cnt.gt.1)then
                acpt1=0
                acpt2=0
                acpt_cnt=0
            end if


            chng_trk=0

            t=newb*(alpha**lralt)
            !Redces b value in learning rate 
            if((t.lt.9.9d-7))then
                if(newb.gt.1.27d-1)then
                    newb=(((0.025+newb)/2)+newb)/2+0.025
                    rjct_cnt=0
                    lralt=0
                else if((newb.gt.1.1d-1))then
                    newb=0.1
                    rjct_cnt=0
                    lralt=0
                else 
                    if(newb.eq.0.1)then
                        if(alpha.le.0.6)then
                            alpha=alpha/2
                            rjct_cnt=0
                            lralt=0
                        else 
                            exit
                        end if
                    end if    
                end if
                write(6,"(a,f12.10,a,f12.10)") "b in learning rate equation t=b*(alpha^x) reduced to ",newb, "alpha is ", alpha 
            end if
            acpt_cnt=0
            if(epoc_cnt.gt.20000)then 
                exit 
            end if
        end do

        ! Now to change individual orbital values

        alpha_zs=0.1  ! learning rate reduction
        lralt_zs=0    ! power alpha is raised to 
        newb_zs=b
        chng_trk=0 !stores which if any ZS changed
        chng_trk2=0 !stores which orbitals in the ZS have changed 
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        l2_rglrstn=0.01 !L2 regularisation lambda paramter
        rpstk=0 !checks for computer error
        rjct_cnt2=0
        acpt_cnt=0
        do while(rjct_cnt2.lt.(norb*100))
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent steps &
                           &   | Orbitals altered '
            do j=1,(ndet-1)
                grad_fin%current_erg=grad_fin%prev_erg
                pick=picker(j)
                pickerorb=scramble_norb(norb)
                chng_trk2=0
                acpt_cnt=0
                do n=1,norb
                    rjct_cnt=0
                    pickorb=pickerorb(n)
                    t=newb_zs(pick,pickorb)*(alpha_zs(pick,pickorb)**lralt_zs(pick,pickorb))
                    do while(t.gt.(1.0d-10))
                        nanchk=.false.
                        if(is_nan(grad_fin%vars(pick,pickorb)).eqv..true.)then
                            call gradient_row(haml,zstore,elect,pick,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                            dvecs(1)%d=cmplx(0.0,0.0)
                            en%erg=0
                            en%t=0
                            call imgtime_prop(dvecs,en,temp_ham,pick)
                            call final_grad(dvecs(1),haml,grad_fin,pick)
                            grad_fin%grad_avlb(pick)=1 
                        end if 
                        
                        ! Setup temporary zombie state
                        temp_zom=zstore
                        temp_zom(pick)%phi(pickorb)=zstore(pick)%phi(pickorb)-(t*(grad_fin%vars(pick,pickorb)))
                        temp_zom(pick)%phi(pickorb)=temp_zom(pick)%phi(pickorb)+&
                                            l2_rglrstn*((grad_fin%vars(pick,pickorb))*grad_fin%vars(pick,pickorb))
                    
                        if(is_nan(temp_zom(pick)%phi(pickorb)).eqv..true.)then
                            t=(1.0d-11)
                        end if
                
                        temp_zom(pick)%sin(pickorb)=sin(cmplx(temp_zom(pick)%phi(pickorb),0.0d0,kind=8))
                        temp_zom(pick)%cos(pickorb)=cos(cmplx(temp_zom(pick)%phi(pickorb),0.0d0,kind=8))
                        temp_ham=haml
                        call he_full_row(temp_ham,temp_zom,elect,pick,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                        ! Imaginary time propagation for back tracing
                        dvecs(1)%d=cmplx(0.0,0.0)
                        en%erg=0
                        en%t=0
                        call imgtime_prop(dvecs,en,temp_ham,0)
                        fxtdk=real(en%erg(1,timesteps+1))
                        if(is_nan(fxtdk).eqv..true.)then 
                            ergerr='NaN ' 
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb 
                            t=(1.0d-11)
                            nanchk=.true.
                            rjct_cnt=1
                        else if(is_posinf(fxtdk).eqv..true.)then
                            ergerr='+NaN'  
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb  
                            t=(1.0d-11)
                            nanchk=.true.
                            rjct_cnt=1
                        else if(is_neginf(fxtdk).eqv..true.)then
                            ergerr='-NaN'  
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb 
                            t=(1.0d-11)
                            nanchk=.true.
                            rjct_cnt=1
                        end if
                        
                        if(nanchk.eqv..false.)then
                            ! Check if energy is lower and accept or reject
                            if(fxtdk.lt.grad_fin%prev_erg)then
                                t=(1.0d-11)
                                acpt_cnt=acpt_cnt+1
                                zstore=temp_zom
                                chng_trk(j)=pick
                                chng_trk2(acpt_cnt)=pickorb
                                rjct_cnt=0
                                rjct_cnt2=0
                                haml=temp_ham
                                grad_fin%grad_avlb=0
                                grad_fin%vars=0.0
                                haml%diff_hjk=0
                                haml%diff_invh=0
                                haml%diff_ovrlp=0
                                grad_fin%prev_erg=fxtdk
                                if(lralt_zs(pick,pickorb).gt.0)then 
                                    lralt_zs(pick,pickorb)=lralt_zs(pick,pickorb)-1
                                end if
                            else 
                                lralt_zs(pick,pickorb)=lralt_zs(pick,pickorb)+1
                                t=newb_zs(pick,pickorb)*(alpha_zs(pick,pickorb)**lralt_zs(pick,pickorb))
                                rjct_cnt=rjct_cnt+1
                            end if
                        end if
                    end do
               
                    if((rjct_cnt.eq.0).and.(n.ne.norb))then
                        call gradient_row(haml,zstore,elect,pick,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                        do k=1,norb 
                            if(is_nan(grad_fin%vars(pick,k)).eqv..true.)then
                                nanchk=.true. 
                                grad_fin%vars=0 
                                grad_fin%grad_avlb=0
                                do l=1 ,10
                                    call gradient_row(haml,zstore,elect,pick,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                                    do m=1,norb
                                        if(is_nan(grad_fin%vars(pick,m)).eqv..true.)then 
                                            exit 
                                        end if 
                                        nanchk=.false. 
                                    end do 
                                    if(nanchk.eqv..false.)then 
                                        exit 
                                    end if 
                                end do  
                                if(nanchk.eqv..true.)then 
                                    errorflag=1 
                                    write(0,"(a,i0)") "Error in gradient calculation for zombie state number ", pick 
                                    return 
                                else 
                                    exit 
                                end if 
                            end if 
                        end do
                        dvecs(1)%d=cmplx(0.0,0.0)
                        en%erg=0
                        en%t=0
                        call imgtime_prop(dvecs,en,temp_ham,pick)
                        call final_grad(dvecs(1),haml,grad_fin,next)
                        grad_fin%grad_avlb(pick)=1
                    else if((rjct_cnt.ne.0)) then
                        if(newb_zs(pick,pickorb).gt.1.27d-1)then
                            newb_zs(pick,pickorb)=(((0.025+newb_zs(pick,pickorb))/2)+newb_zs(pick,pickorb))/2+0.025
                            lralt_zs(pick,pickorb)=0
                        else if((newb_zs(pick,pickorb).gt.1.1d-1))then
                            newb_zs(pick,pickorb)=0.1
                            lralt_zs(pick,pickorb)=0
                        else 
                            if(newb_zs(pick,pickorb).eq.0.1)then
                                if(alpha_zs(pick,pickorb).le.0.6)then
                                    alpha_zs(pick,pickorb)=alpha_zs(pick,pickorb)/2
                                    lralt_zs(pick,pickorb)=0
                                else 
                                    exit
                                end if
                            end if    
                        end if
                    end if
                end do
                

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
                    next=picker(j+1)
                end if

                if(grad_fin%grad_avlb(next).eq.1)then  
                    do k=1,norb
                        if(isnan(grad_fin%vars(next,k)).eqv..true.)then
                            grad_fin%grad_avlb(next)=0
                            exit 
                        end if 
                    end do
                end if

                !Set up gradients for next pass
                if(grad_fin%grad_avlb(next).eq.0)then 
                    call gradient_row(haml,zstore,elect,next,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                    do k=1,norb 
                        if(is_nan(grad_fin%vars(next,k)).eqv..true.)then
                            nanchk=.true. 
                            grad_fin%vars=0 
                            grad_fin%grad_avlb=0
                            do l=1 ,10
                                call gradient_row(haml,zstore,elect,next,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                                do m=1,norb
                                    if(is_nan(grad_fin%vars(next,m)).eqv..true.)then 
                                        exit 
                                    end if 
                                    nanchk=.false. 
                                end do 
                                if(nanchk.eqv..false.)then 
                                    exit 
                                end if 
                            end do  
                            if(nanchk.eqv..true.)then 
                                errorflag=1 
                                write(0,"(a,i0)") "Error in gradient calculation for zombie state number ", next 
                                return 
                            else 
                                exit 
                            end if 
                        end if 
                    end do
                    grad_fin%grad_avlb(next)=1
                    dvecs(1)%d=cmplx(0.0,0.0)
                    en%erg=0
                    en%t=0
                    call imgtime_prop(dvecs,en,temp_ham,next)
                    call final_grad(dvecs(1),haml,grad_fin,next)
                end if
                if(acpt_cnt.gt.0)then
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,*(i0:','))") '       ', pick,'              ', &
                    grad_fin%current_erg,'               ',grad_fin%prev_erg,'            ',chng_trk2(1:acpt_cnt) 
                    rsrtpass(pick)=0
                    zstore(pick)%update_num=zstore(pick)%update_num+1
                    call zombiewriter(zstore(pick),pick,rsrtpass(pick))
                else 
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,i0)") '       ', pick,'              ', &
                    grad_fin%current_erg,'               ',grad_fin%prev_erg,'            ',0
                    rjct_cnt2=rjct_cnt2+1
                end if
            end do

            call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,0)
            write(6,"(a,i0,a,f21.16)") "Energy after epoch no. ",epoc_cnt,": ",grad_fin%prev_erg
            acpt_cnt=0
            chng_trk=0
            if(epoc_cnt.gt.30000)then 
                exit 
            end if
        end do

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
        deallocate(occupancy_an,stat=ierr)
        if(ierr==0) deallocate(occupancy_2an,stat=ierr)
        if(ierr==0) deallocate(occupancy_an_cr,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in occupancy vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

    end subroutine zombie_alter

    ! Produces a random order for the ZS to be posisbly changed
    function scramble( number_of_values ) result(out)
        
        implicit none

        integer,intent(in)    :: number_of_values
        integer,allocatable   :: out(:),array(:)
        integer::n,m,k,j,l,jtemp
        real::u

        out=[(j,j=1,number_of_values)]
        array=[(j,j=1,number_of_values+1)]
        n=1; m=number_of_values
        do k=1,2
            do j=1,m+1
                call random_number(u)
                l = n + FLOOR((m+1-n)*u)
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

        integer,intent(in)    :: number_of_values
        integer,allocatable   :: out(:)
        integer::n,m,k,j,l,jtemp
        real::u

        out=[(j,j=1,number_of_values)]
        n=1; m=number_of_values
        do k=1,2
            do j=1,m
                call random_number(u)
                l = n + FLOOR((m-n)*u)
                jtemp=out(l)
                out(l)=out(j)
                out(j)=jtemp
            end do
        end do

        return
    end function scramble_norb

END MODUlE gradient_descent