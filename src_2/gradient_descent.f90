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
        !$omp parallel do
        do k=1, ndet
            do l=1, ndet
                if(l.eq.diff_state)then
                    temp2(k,diff_state,:)=matmul(REAL(haml%inv(k,:)),haml%diff_ovrlp(diff_state,:,:))
                else
                    temp2(k,l,:)=real(haml%inv(k,l))*haml%diff_ovrlp(diff_state,l,:)
                end if
            end do
        end do
        !$omp parallel do
        do k=1, ndet
            do l=1, ndet
                haml%diff_invh(diff_state,k,l,:)=matmul(transpose(temp2(k,:,:)),real(haml%kinvh(:,l)))*(-1)
            end do
        end do

        return

    end subroutine gradient_row

    subroutine grad_calc(haml,zstore,elect,pick,occupancy_2an,occupancy_an_cr,occupancy_an,dvec,grad_fin,en,d_diff_flg)

        implicit none 

        type(zombiest),dimension(:),intent(in)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvec
        type(hamiltonian),intent(inout)::haml
        integer,dimension(:,:,:),intent(in)::occupancy_an
        integer,dimension(:,:,:,:),intent(in)::occupancy_an_cr,occupancy_2an
        integer,intent(in)::pick,d_diff_flg
        integer::k,l,m
        logical::nanchk
        type(energy),intent(inout)::en

        if (errorflag .ne. 0) return

        call gradient_row(haml,zstore,elect,pick,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
        nanchk=.false. 
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
        en%erg=0
        en%t=0

        if(d_diff_flg.eq.1)then
            call imgtime_prop(dvec,en,haml,pick)
        end if

        call final_grad(dvec(1),haml,grad_fin,pick,d_diff_flg)

       
        if(d_diff_flg.eq.0)then
            grad_fin%grad_avlb(pick)=1
        else if(d_diff_flg.eq.1)then 
            grad_fin%grad_avlb(pick)=2
        end if
        return

    end subroutine grad_calc

    subroutine orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
        occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker,maxloop) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore,temp_zom
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvecs,temp_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml,temp_ham
        integer,dimension(:),intent(inout)::chng_trk
        integer,intent(inout)::epoc_cnt
        integer,intent(in)::maxloop
        integer,dimension(:,:,:)::occupancy_an
        integer,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an
        real(kind=8),intent(in)::b,alphain
        integer,dimension(ndet-1),intent(inout)::picker
        integer::rjct_cnt,next,acpt_cnt,pick,pickorb,rjct_cnt2,loops,d_diff_flg
        real(kind=8)::t,fxtdk,l2_rglrstn
        integer,dimension(ndet-1)::rsrtpass
        integer,dimension(ndet,norb)::lralt_zs
        integer,dimension(norb)::pickerorb
        real(kind=8),dimension(ndet,norb)::alpha_zs,newb_zs
        integer,dimension(norb)::chng_trk2
        logical::nanchk
        character(len=4)::ergerr
        integer::j,k,l,n
        
        alpha_zs=alphain  ! learning rate reduction
        lralt_zs=0    ! power alpha is raised to 
        newb_zs=b
        chng_trk=0 !stores which if any ZS changed
        chng_trk2=0 !stores which orbitals in the ZS have changed 
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        l2_rglrstn=1! !L2 regularisation lambda paramter
        rjct_cnt2=0
        acpt_cnt=0
        loops=0
        d_diff_flg=1
        do while(rjct_cnt2.lt.(norb*100))
            loops=loops+1
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent steps &
                           &   | Orbitals altered '
            do j=1,(ndet-1)
                grad_fin%current_erg=grad_fin%prev_erg
                pick=picker(j)
                pickerorb=scramble_norb(norb)
                ! pickerorb=(/1,2,3,4,5,6,7,8,9,10/)
                chng_trk2=0
                acpt_cnt=0
                do n=1,norb
                    rjct_cnt=0
                    pickorb=pickerorb(n)
                    t=newb_zs(pick,pickorb)*(alpha_zs(pick,pickorb)**lralt_zs(pick,pickorb))
                    do while(t.gt.(1.0d-13))
                        nanchk=.false.
                        if(is_nan(grad_fin%vars(pick,pickorb)).eqv..true.)then
                            call grad_calc(haml,zstore,elect,pick,occupancy_2an,occupancy_an_cr,&
                            occupancy_an,dvecs,grad_fin,en,d_diff_flg) 
                        end if 
                        
                       
                        ! Setup temporary zombie state
                        temp_zom=zstore
                        temp_zom(pick)%phi(pickorb)=zstore(pick)%phi(pickorb)-(t*(grad_fin%vars(pick,pickorb)))
                        temp_zom(pick)%phi(pickorb)=temp_zom(pick)%phi(pickorb)+&
                                            l2_rglrstn*((grad_fin%vars(pick,pickorb))*grad_fin%vars(pick,pickorb))
                        if(is_nan(temp_zom(pick)%phi(pickorb)).eqv..true.)then
                            t=(1.0d-14)
                        end if
                        temp_zom(pick)%sin(pickorb)=sin(cmplx(temp_zom(pick)%phi(pickorb),0.0d0,kind=8))
                        temp_zom(pick)%cos(pickorb)=cos(cmplx(temp_zom(pick)%phi(pickorb),0.0d0,kind=8))
                        temp_ham=haml
                        call he_full_row(temp_ham,temp_zom,elect,pick,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                        ! Imaginary time propagation for back tracing
                        
                        en%erg=0
                        en%t=0
                        call imgtime_prop(temp_dvecs,en,temp_ham,0)
                        fxtdk=real(en%erg(1,timesteps+1))
                    
                        
                        if(is_nan(fxtdk).eqv..true.)then 
                            ergerr='NaN ' 
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb 
                            t=(1.0d-14)
                            nanchk=.true.
                            rjct_cnt=1
                        else if(is_posinf(fxtdk).eqv..true.)then
                            ergerr='+NaN'  
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb  
                            t=(1.0d-14)
                            nanchk=.true.
                            rjct_cnt=1
                        else if(is_neginf(fxtdk).eqv..true.)then
                            ergerr='-NaN'  
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb 
                            t=(1.0d-14)
                            nanchk=.true.
                            rjct_cnt=1
                        end if
                        
                        if(nanchk.eqv..false.)then
                            ! Check if energy is lower and accept or reject
                            if(fxtdk.lt.grad_fin%prev_erg)then
                                t=(1.0d-14)
                                acpt_cnt=acpt_cnt+1
                                zstore=temp_zom
                                dvecs=temp_dvecs
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
                                grad_fin%prev_mmntm(pick,:)=zstore(pick)%phi(:)
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
                        call grad_calc(haml,zstore,elect,pick,occupancy_2an,&
                        occupancy_an_cr,occupancy_an,dvecs,grad_fin,en,d_diff_flg) 
                    end if
                end do
            

                if(j.eq.(ndet-1))then 
                    epoc_cnt=epoc_cnt+1
                    picker=scramble(ndet-1)
                    next=picker(1)
                    lralt_zs=0
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
                        if(is_nan(grad_fin%vars(next,k)).eqv..true.)then
                            grad_fin%grad_avlb(next)=0
                            exit 
                        end if 
                    end do
                end if

                !Set up gradients for next pass
                if(grad_fin%grad_avlb(next).eq.0)then 
                    call grad_calc(haml,zstore,elect,next,occupancy_2an,occupancy_an_cr,occupancy_an,dvecs,grad_fin,en,d_diff_flg) 
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
            if(loops.ge.maxloop)then 
                exit 
            end if
        end do

    end subroutine orbital_gd


    subroutine full_zs_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
        occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore,temp_zom
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvecs
        type(dvector), dimension(:), allocatable:: temp_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml,temp_ham
        integer,dimension(:),intent(inout)::chng_trk
        integer,dimension(:,:,:),intent(in)::occupancy_an
        integer,dimension(:,:,:,:),intent(in)::occupancy_an_cr,occupancy_2an
        integer,intent(inout)::epoc_cnt
        real(kind=8),intent(in)::b,alphain
        integer,dimension(ndet-1),intent(inout)::picker
        integer::lralt,rjct_cnt,next,acpt_cnt,pick,orbitcnt,d_diff_flg,lralt_temp,mmntmflg,loop_max,loop_dwn
        real(kind=8)::newb,t,fxtdk,l2_rglrstn,alpha,mmnmtb
        real(kind=8),dimension(norb)::gradient_norm
        real(kind=8),dimension(ndet)::mmntm,mmntma
        integer,dimension(ndet-1)::rsrtpass
        !DOUBLE PRECISION, external::ZBQLU01
        real::r
        logical::nanchk
        character(len=4)::ergerr
        integer::j,k,l

       
        if (errorflag .ne. 0) return


        alpha=alphain  ! learning rate reduction
        lralt=0    ! power alpha is raised to  
        newb=b
        chng_trk=0 !stores which if any ZS changed
        rjct_cnt=0 !tracks how many rejections 
        t=newb*(alpha**lralt) !learning rate
        acpt_cnt=0  !counts how many ZS have been changed
        l2_rglrstn=1! !L2 regularisation lambda paramter
        orbitcnt=0
        mmntm=0
        mmntma=1
        mmntmflg=0
        loop_max=20
        loop_dwn=0 

        do while(rjct_cnt.lt.200)
            t=newb*(alpha**lralt)
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent step &
            &   | Learning rate | Acceptance count | Rejection count'
            
            do j=1,(ndet-1)
                
                pick=picker(j)
                lralt_temp=lralt
                do k=1,norb
                    if(is_nan(grad_fin%vars(pick,k)).eqv..true.)then
                        call grad_calc(haml,zstore,elect,pick,occupancy_2an,occupancy_an_cr,&
                        occupancy_an,dvecs,grad_fin,en,d_diff_flg)
                        nanchk=.false. 
                        exit 
                    end if 
                end do
                gradient_norm=((grad_fin%vars(pick,:))*grad_fin%vars(pick,:))
                do while(lralt_temp.lt.(loop_max))
                    t=newb*(alpha**lralt_temp)
                    nanchk=.false.
                    mmnmtb=(t*mmntm(pick))/mmntma(pick)
                    ! Setup temporary zombie state
                    temp_zom=zstore
                    temp_zom(pick)%phi(:)=zstore(pick)%phi(:)-(t*(grad_fin%vars(pick,:)))+&
                       mmnmtb*(zstore(pick)%phi(:)-grad_fin%prev_mmntm(pick,:))
                    temp_zom(pick)%phi(:)=temp_zom(pick)%phi(:)+l2_rglrstn*gradient_norm
                
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
                    en%erg=0
                    en%t=0
                    call imgtime_prop(temp_dvecs,en,temp_ham,0)
                
                    fxtdk=real(en%erg(1,timesteps+1))
                   
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
                            dvecs=temp_dvecs
                            rjct_cnt=0
                            grad_fin%grad_avlb=0
                            grad_fin%vars=0.0
                            haml%diff_hjk=0
                            haml%diff_invh=0
                            haml%diff_ovrlp=0
                            grad_fin%prev_erg=fxtdk
                            d_diff_flg=0
                            grad_fin%prev_mmntm(pick,:)=zstore(pick)%phi(:)
                            mmntma(pick)=t
                            mmntm(pick)=mmnmtb
                            loop_dwn=loop_dwn+1
                            if(orbitcnt.lt.0)then
                                orbitcnt=orbitcnt+1 
                            else
                                orbitcnt=0
                            end if
                            if((loop_max.gt.20).and.(loop_dwn.eq.5))then
                                loop_max=loop_max-1
                                loop_dwn=0
                            end if
                            ! if(lralt.eq.1)then 
                            !     lralt=0
                            ! else if(lralt.gt.1)then 
                            !     lralt=lralt-2
                            ! end if
                            write(6,"(a,i0,a,f21.16,a,f21.16,a,f12.10,a,i0,a,i0)") '       ', pick,'              ', &
                grad_fin%prev_erg,'               ',fxtdk,'            ',t,'          ',acpt_cnt,'                 ',rjct_cnt
                            Exit
                            ! t=newb*(alpha**lralt)
                        else 
                            lralt_temp=lralt_temp+1
                            ! rjct_cnt=rjct_cnt+1
                            ! if(rjct_cnt.eq.(ndet-1))then
                            !     lralt=lralt+1
                            !     t=newb*(alpha**lralt)
                            !     rjct_cnt=0
                            !     if(epoc_cnt.gt.100)then
                            !         orbitcnt=orbitcnt+1
                            !     end if
                            ! end if
                        end if
                    else 
                        write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                                                " for zombie state ", pick, ", on epoc ", epoc_cnt 
                    end if
                end do
                if(lralt_temp.ge.loop_max)then
                    rjct_cnt=rjct_cnt+1
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,f12.10,a,i0,a,i0)") '       ', pick,'              ', &
                grad_fin%prev_erg,'               ',fxtdk,'            ',0.0,'          ',acpt_cnt,'                 ',rjct_cnt
                    if(orbitcnt.ne.0)then
                        orbitcnt=orbitcnt+1
                    end if
                    if(loop_max.lt.35)then 
                        loop_max=loop_max+1
                    end if
                    loop_dwn=0 
                end if
                
               
                t=newb*(alpha**lralt)
                if(j.eq.(ndet-1))then 
                    epoc_cnt=epoc_cnt+1
                    if(acpt_cnt.gt.0)then 
                        picker=scramble(ndet-1)
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
        
                if(grad_fin%grad_avlb(next).eq.1)then 
                    d_diff_flg=1 
                else 
                    d_diff_flg=0
                    call random_number(r)
                    if(r.lt.0.3)then 
                    ! if(ZBQLU01(1).lt.0.3)then 
                        d_diff_flg=1
                    end if
                end if

                !Set up gradients for next pass
                if(grad_fin%grad_avlb(next).lt.2)then
                    call grad_calc(haml,zstore,elect,next,occupancy_2an,occupancy_an_cr,occupancy_an,dvecs,grad_fin,en,d_diff_flg) 
                end if

            end do

            call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,0)
            write(6,"(a,i0,a,f21.16,a,f12.10,a,i0,a)") "Energy after epoc no. ",epoc_cnt,": ", &
                grad_fin%prev_erg, ". The current learning rate is: ",t, ". ", acpt_cnt, " Zombie state(s) altered."

            
            ! If only a single ZS is being altered the learning rate is lowered to allow others to be changed
            if(acpt_cnt.eq.1)then
                if((orbitcnt.ge.0).and.(epoc_cnt.gt.100))then  
                    call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
                    occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker,1)
                    orbitcnt=-25
                    lralt=0
                end if   
            else if(acpt_cnt.eq.0)then
                if((epoc_cnt.gt.50))then  
                    call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
                    occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker,1)
                    orbitcnt=-25
                end if
            end if

            if((epoc_cnt.gt.500).and.(orbitcnt.ge.0).and.(modulo(epoc_cnt,100).eq.0))then
                call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
                occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker,1)
                orbitcnt=-25
                lralt=0
            end if

            chng_trk=0
            t=newb*(alpha**lralt)
            if(mmntmflg.eq.0)then 
                mmntm=0.4
                mmntmflg=1
            end if 
            acpt_cnt=0
            if(epoc_cnt.gt.30000)then 
                exit 
            end if
        end do

    end subroutine full_zs_gd

    subroutine zombie_alter(zstore,grad_fin,haml,elect,en,dvecs,chng_trk)

        implicit none

        type(zombiest),dimension(:),intent(inout)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvecs
        type(dvector), dimension(:), allocatable:: temp_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml
        integer,dimension(:),intent(inout)::chng_trk
        type(zombiest),dimension(:),allocatable::temp_zom
        type(hamiltonian)::temp_ham
        integer,allocatable,dimension(:,:,:)::occupancy_an
        integer,allocatable,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an
        integer::epoc_cnt,epoc_max
        real(kind=8)::alpha,b
        integer,dimension(ndet-1)::picker,rsrtpass

     
        integer::ierr,j,k,l
        
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

            call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,1)
            rsrtpass=1
            ierr=0
        end if

        alpha=0.5  ! learning rate reduction
        b=8.0D0 !starting learning rate
        
       
        chng_trk=0 !stores which if any ZS changed
        rsrtpass=0
        epoc_max=30000

        do j=1, ndet-1
            picker(j)=j+1
        end do

        call alloczs(temp_zom,ndet)
        call allocham(temp_ham,ndet,norb)
        call allocdv(temp_dvecs,1,ndet,norb)

        if(epoc_cnt.eq.0)then
            call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
            occupancy_an_cr,occupancy_2an,epoc_cnt,alpha,b,picker,1)
        end if 
     
        call full_zs_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,&
        chng_trk,temp_zom,occupancy_an,occupancy_an_cr,occupancy_2an,epoc_cnt,alpha,b,picker) 
        
        call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
        occupancy_an_cr,occupancy_2an,epoc_cnt,alpha,b,picker,(epoc_max-epoc_cnt)) 

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
        call deallocdv(temp_dvecs) 
        deallocate(occupancy_an,stat=ierr)
        if(ierr==0) deallocate(occupancy_2an,stat=ierr)
        if(ierr==0) deallocate(occupancy_an_cr,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in occupancy vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

    end subroutine zombie_alter


    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators
    subroutine he_full_row_gpu(ph1ei,ph2ei,phnuc,psin,pcos,occupancy_2an,&
        occupancy_an_cr,occupancy_an,phjk,povrlp,pinv,pkinvh,diff_state,size)

        implicit none

        real(kind=8),pointer,dimension(:,:),intent(in)::ph1ei
        real(kind=8), dimension(:,:,:,:),pointer,intent(in)::ph2ei
        real(kind=8),pointer,intent(in) :: phnuc
        complex(kind=8), dimension(:,:),intent(in)::psin
        complex(kind=8), dimension(:,:),intent(in)::pcos
        complex(kind=8), dimension(:,:), pointer,intent(inout)::phjk
        complex(kind=8), dimension(:,:), pointer,intent(inout)::povrlp
        complex(kind=8), dimension(:,:), pointer,intent(inout)::pinv
        complex(kind=8), dimension(:,:), pointer,intent(inout)::pkinvh
        integer,intent(in)::size,diff_state
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an,occupancy_an_cr
        integer,dimension(:,:,:),intent(in)::occupancy_an
        complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
        complex(kind=8),allocatable, dimension(:,:,:)::z2l
        integer::j,k,l,m,ierr,n,gmax,hmin,jspin
        complex(kind=8)::h1etot, h2etot
        complex(kind=8),dimension(2,norb)::zomt,vmult
        complex(kind=8),dimension(norb)::gg,hh
        integer, allocatable,dimension(:)::IPIV1
        complex(kind=8),allocatable,dimension(:)::WORK1
        complex(kind=8)::totc

        if (errorflag .ne. 0) return
        ierr = 0

        allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

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
      
        !!$omp target teams map(alloc:z1jk(norb,norb,2,norb),z2l(norb,2,norb),zomt(2,norb),&
        !!$omp vmult(2,norb),gg(norb),hh(norb)) &
        !!$omp & map(to:h1etot,h2etot,jspin,gmax,hmin,totc)&
        !!$omp & shared(ph1ei,z2l,z1jk,occupancy_2an,occupancy_an,occupancy_an_cr,ph2ei)&
        !!$omp & private(vmult,gg,hh,jspin,gmax,hmin,zomt) 

        !$omp target data map(alloc:z1jk(norb,norb,2,norb),z2l(norb,2,norb),zomt(2,norb),&
        !$omp vmult(2,norb),gg(norb),hh(norb)) &
        !$omp & map(to:h1etot,h2etot,jspin,gmax,hmin,totc)

        h1etot=cmplx(0.0,0.0)
        h2etot=cmplx(0.0,0.0)
        totc=(0.0,0.0)
        !$omp target teams 
        !$omp  distribute parallel do simd
        do l=1, norb
            z2l(l,1,:)=psin(diff_state,:) 
            z2l(l,2,:)=pcos(diff_state,:)
            z2l(l,2,l)=z2l(l,1,l)
            z2l(l,1,l)=cmplx(0.0,0.0)
        end do
       !$omp end  distribute parallel do simd

        !$omp  distribute parallel do simd collapse(2)
        do j=1, norb
            do k=1, norb
                z1jk(j,k,:,:)=z2l(j,:,:)
                z1jk(j,k,2,k)=z1jk(j,k,1,k)
                z1jk(j,k,1,k)=cmplx(0.0,0.0)
            end do
        end do
        !$omp end  distribute parallel do simd
        !$omp end target teams
        z1jk=z1jk*occupancy_2an

        do m=1,size
      
            h1etot=cmplx(0.0,0.0)
            h2etot=cmplx(0.0,0.0)
            totc=(0.0,0.0)
            !$omp target teams distribute parallel do simd
            do l=1, norb
                z2l(l,1,:)=psin(m,:) 
                z2l(l,2,:)=pcos(m,:)
                z2l(l,2,l)=z2l(l,1,l)
                z2l(l,1,l)=cmplx(0.0,0.0)
            end do
            !$omp end target teams distribute parallel do simd
            
            !$omp target teams distribute parallel do simd collapse(2) reduction(+:totc) &
            !$omp shared(ph1ei,occupancy_an_cr,z2l) private(zomt)
            do j=1, norb
                do k=1, norb
                    if(ph1ei(j,k).eq.0.0)then
                        cycle
                    end if
                    zomt(:,:)=z2l(j,:,:)
                    zomt(1,k)=zomt(2,k)
                    zomt(2,k)=cmplx(0.0,0.0)
                    zomt=zomt*occupancy_an_cr(j,k,:,:)
                    totc=totc+product((conjg(psin(diff_state,:))*zomt(1,:))+(conjg(pcos(diff_state,:))*zomt(2,:)))*ph1ei(j,k) 
                end do 
            end do
           !$omp end target teams distribute parallel do simd
            h1etot=totc
          
            z2l=z2l*occupancy_an
            totc=(0.0,0.0)
            !$omp flush
            !$omp target teams distribute parallel do reduction(+:h2etot,totc)&
            !$omp & private(vmult,gg,hh,jspin,gmax,hmin,zomt) shared(z2l,z1jk,ph2ei)
            do j=1, norb
                if(psin(diff_state,j)==(0.0,0.0))then
                    CYCLE
                end if
                if(modulo(j,2)==0)then
                    jspin=2
                else
                    jspin=1
                end if
                do k=1, norb
                    if(j.eq.k) cycle
                    if(occ_iszero(z1jk(j,k,:,:)).eqv..true.)then
                        CYCLE
                    end if
                    do l=jspin, norb, 2
                        if(psin(m,l)==(0.0,0.0))then
                            CYCLE
                        end if
    
                        vmult=conjg(z1jk(j,k,:,:))*(z2l(l,:,:))
        
                        gg(1:norb)=(0.0,0.0)
                        hh(1:norb)=(0.0,0.0)
                        gmax=norb
                        gg(1)=vmult(2,1)-vmult(1,1)
    
                        do n=2, norb
                            gg(n)=gg(n-1)*(vmult(2,n)-vmult(1,n))
                            if(gg(n)==(0.0,0.0))then
                                gmax=n
                                EXIT 
                            end if
                        end do
                        
                        hmin=0
                        hh(norb) = vmult(2,norb)+vmult(1,norb)
                        do n=(norb-1),1,-(1)
                            hh(n)=hh(n+1)*(vmult(2,n)+vmult(1,n))
                            if(hh(n)==(0.0,0.0))then
                                hmin=n
                                EXIT 
                            end if
                        end do
                        totc=(0.0,0.0)
                        if (gmax < hmin) then
                            h2etot=h2etot+totc
                            cycle
                        end if
    
                        if(ph2ei(j,k,l,1).ne.0) then
                            totc = totc+(conjg(z1jk(j,k,2,1))*z2l(l,1,1)*hh(2)*ph2ei(j,k,l,1))
                        end if
    
                        do n=2,norb-1
                            if(ph2ei(j,k,l,n).ne.0.0) then
                                totc = totc+ (gg(n-1)*conjg(z1jk(j,k,2,n))*z2l(l,1,n)*hh(n+1)*ph2ei(j,k,l,n))
                            end if
                        end do
    
                        if(ph2ei(j,k,l,norb).ne.0) then
                            totc = totc +(gg(norb-1)*conjg(z1jk(j,k,2,norb))*z2l(l,1,norb)*ph2ei(j,k,l,norb))
                        end if
                        h2etot=h2etot+totc
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            !$omp flush(h2etot)
            h2etot=h2etot*0.5
          
            povrlp(diff_state,m)=product(((conjg(psin(diff_state,:))*psin(m,:)))+((conjg(pcos(diff_state,:))*pcos(m,:))))
            povrlp(m,diff_state)= povrlp(diff_state,m)
            phjk(diff_state,m)=h1etot+h2etot+(phnuc*povrlp(diff_state,m))
            phjk(m,diff_state)=phjk(diff_state,m)
          

        end do
        !$omp end target data
        !$omp target update from(povrlp,phjk)
        pinv=povrlp
        
        if (ierr==0) call ZGETRF(size,size,pinv,size,IPIV1,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if
        if (ierr==0) call ZGETRI(size,pinv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if

        !$omp parallel
        !$omp workshare
        pkinvh=matmul(pinv,phjk)
        !$omp end workshare
        !$omp end parallel
        !$omp target update to(pinv, pkinvh)
       

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
        deallocate(z1jk,stat=ierr)
        if(ierr==0) deallocate(z2l,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        return

    end subroutine he_full_row_gpu

    subroutine gradient_row_gpu(psin,pcos,pphi,ph1ei,ph2ei,pkinvh,pinv,pdiff_hjk,&
        pdiff_ovrlp,pdiff_invh,occupancy_an,occupancy_an_cr,occupancy_2an,diff_state,size)


        implicit none

        real(kind=8),pointer,dimension(:,:),intent(in)::ph1ei
        real(kind=8), dimension(:,:,:,:),pointer,intent(in)::ph2ei
        complex(kind=8), dimension(:,:),intent(in)::psin
        complex(kind=8), dimension(:,:),intent(in)::pcos
        real(kind=8),dimension(:,:),intent(in)::pphi
        complex(kind=8), dimension(:,:), pointer,intent(inout)::pinv
        complex(kind=8), dimension(:,:), pointer,intent(inout)::pkinvh
        real(kind=8), dimension(:,:,:), pointer,intent(inout)::pdiff_hjk
        real(kind=8), dimension(:,:,:), pointer,intent(inout)::pdiff_ovrlp
        real(kind=8), dimension(:,:,:,:), pointer,intent(inout)::pdiff_invh
        integer,intent(in)::size,diff_state
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an,occupancy_an_cr
        integer,dimension(:,:,:),intent(in)::occupancy_an
        complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
        complex(kind=8),allocatable, dimension(:,:,:)::z2l
        complex(kind=8),dimension(2,norb)::zomt
        real(kind=8),dimension(norb)::h1etot_diff_bra,h2etot_diff_bra
        real(kind=8),dimension(ndet,ndet,norb)::temp2
        real(kind=8),dimension(2,norb)::vmult_dd,vmultr
        integer,dimension(2,norb)::occupancy
        integer::j,k,l,m,ierr, jspin
        real(kind=8),dimension(norb)::bra_prod,prod,temp1,gg_1,hh_1,temp
        integer::n,p,gmax1,hmin1,breakflag
        real(kind=8)::tot
        if (errorflag .ne. 0) return
        ierr = 0

        allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

      
        !$omp target teams map(alloc:z1jk(norb,norb,2,norb),z2l(norb,2,norb),zomt(2,norb),temp(norb),temp2(ndet,ndet,norb),&
        !$omp  bra_prod(norb),prod(norb),temp1(norb),vmultr(2,norb),vmult_dd(2,norb),occupancy(2,norb),&
        !$omp h1etot_diff_bra(norb),h2etot_diff_bra(norb)) & !map(to:j,k,l,diff_state,gg_1,hh_1)&
        !map(to:j,k,l,diff_state,gg_1,hh_1,gmax1,hmin1,jspin,breakflag,tot)&
        !$omp & shared(ph1ei,ph2ei,psin,pcos,pphi,occupancy_an_cr,z2l,z1jk,occupancy_2an,occupancy_an,diff_state)&
        !$omp & private(j,k,l,zomt,temp1,bra_prod,prod,gg_1,hh_1,gmax1,hmin1,jspin,occupancy,breakflag,vmultr,vmult_dd,temp,tot)
        !!$omp &num_teams(8) !call omp_set_num_teams(8)
       
        !$omp distribute parallel do simd
        do l=1, norb
            z2l(l,1,:)=psin(diff_state,:) 
            z2l(l,2,:)=pcos(diff_state,:)
            z2l(l,2,l)=z2l(l,1,l)
            z2l(l,1,l)=cmplx(0.0,0.0)
        end do
       !$omp end distribute parallel do simd

    
        !$omp distribute parallel do simd collapse(2)
        do j=1, norb
            do k=1, norb
                z1jk(j,k,:,:)=z2l(j,:,:)
                z1jk(j,k,2,k)=z1jk(j,k,1,k)
                z1jk(j,k,1,k)=cmplx(0.0,0.0)
            end do
        end do
        !$omp end distribute parallel do simd

        z1jk=z1jk*occupancy_2an

        do m=1,size

            h1etot_diff_bra=0.0
            h2etot_diff_bra=0.0
         
            !$omp distribute parallel do simd
            do l=1, norb
                z2l(l,1,:)=psin(m,:) 
                z2l(l,2,:)=pcos(m,:)
                z2l(l,2,l)=z2l(l,1,l)
                z2l(l,1,l)=cmplx(0.0,0.0)
            end do
            !$omp end distribute parallel do simd


            !$omp distribute parallel do simd collapse(2) reduction(+:h1etot_diff_bra)
            do j=1, norb
                do k=1, norb
                    if(ph1ei(j,k).eq.0.0)then
                        cycle
                    end if
                    zomt(:,:)=z2l(j,:,:)
                    zomt(1,k)=zomt(2,k)
                    zomt(2,k)=cmplx(0.0,0.0)
                    zomt=zomt*occupancy_an_cr(j,k,:,:)
                    prod=real(((conjg(psin(diff_state,:))*zomt(1,:)))+((conjg(pcos(diff_state,:)))*zomt(2,:)))
                        ! Differntial of overlap of the same ZS 0 unless an operator has acted on ZS
                    if(diff_state.eq.m)then
                        if(j.eq.k)then
                            if((real(pcos(m,j)).eq.0).and.(real(psin(m,j)).eq.1).or.(real(psin(m,j)).eq.0))then
                                h1etot_diff_bra(j)=h1etot_diff_bra(j)+0
                            else
                                bra_prod=prod           !dead amplitude is zero
                                bra_prod(j)=sin(2*pphi(diff_state,j))*occupancy_an_cr(j,k,1,j)
                                h1etot_diff_bra(j)= h1etot_diff_bra(j)+product(bra_prod)*ph1ei(j,k)    
                            end if 
                        else if(j.ne.k)then
                            bra_prod=prod
                            if((real(pcos(m,j)).eq.0).and.(real(psin(m,j)).eq.1)) then
                                bra_prod(j)=-1*occupancy_an_cr(j,k,2,j)
                            else if((real(psin(m,j)).eq.0).and.(real(pcos(m,j)).eq.1)) then
                                bra_prod(j)=occupancy_an_cr(j,k,2,j)
                            else
                                bra_prod(j)=cos(2*pphi(diff_state,j))*occupancy_an_cr(j,k,2,j)
                            end if
                            h1etot_diff_bra(j)= h1etot_diff_bra(j) + product(bra_prod)*ph1ei(j,k)      
                            bra_prod=prod
                            if((real(pcos(m,k)).eq.0).and.(real(psin(m,k)).eq.1)) then
                                bra_prod(k)=-1*occupancy_an_cr(j,k,1,k)
                            else if((real(psin(m,k)).eq.0).and.(real(pcos(m,k)).eq.1)) then
                                bra_prod(k)=occupancy_an_cr(j,k,1,k)
                            else
                                bra_prod(k)=cos(2*pphi(diff_state,k))*occupancy_an_cr(j,k,1,k)
                            end if               !dead amplitude is zero
                            h1etot_diff_bra(k)= h1etot_diff_bra(k) + product(bra_prod)*ph1ei(j,k)
                        end if
                    else
                        bra_prod=real(psin(diff_state,:)*psin(m,:)*occupancy_an_cr(j,k,1,:)-&
                        psin(diff_state,:)*pcos(m,:)*occupancy_an_cr(j,k,2,:))            
                        if(j.eq.k)then  !dead amplitude is zero
                            bra_prod(j)=real(pcos(diff_state,j)*psin(m,j)*occupancy_an_cr(j,k,1,j))
                        else                   
                            bra_prod(j)=-real(psin(diff_state,j)*psin(m,j))*occupancy_an_cr(j,k,2,j)
                            bra_prod(k)=real(pcos(diff_state,k)*pcos(m,k))*occupancy_an_cr(j,k,1,k) 
                        end if
                        do l=1,norb
                            temp1=prod
                            temp1(l)=bra_prod(l)
                            h1etot_diff_bra(l)= h1etot_diff_bra(l)+product(temp1)*ph1ei(j,k)   
                        end do 
                    end if
                end do
            end do
            !$omp end distribute parallel do simd

            !$omp distribute parallel do reduction(+:h2etot_diff_bra)
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
                        if(diff_state.eq.m)then !Differentiation when zombie states are the same
                            do n=1, norb    !Differentiating w.r.t to orbital n
                                breakflag=0
                                temp=0
                                vmultr=real(z1jk(j,k,:,:))*real(z2l(l,:,:))
                                vmult_dd=vmultr
                                if(j.eq.k)then
                                    temp(n)=0
                                    cycle
                                else if((j.eq.n).or.(k.eq.n))then
                                    if(l.eq.n)then    !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(1j)
                                        vmult_dd(1,n)=0
                                        vmult_dd(2,n)=2*real(psin(diff_state,n)*pcos(diff_state,n))*occupancy(2,n)
                                    else
                                        vmult_dd(1,n)=0        !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(0j)
                                        vmult_dd(2,n)=(1.0-2*((real(psin(diff_state,n)))**2))*occupancy(2,n)   
                                    end if
                                else if((j.ne.n).or.(k.ne.n))then
                                    if(l.eq.n)then !before diff alive:0 dead:a^(a)_(0j)*a^(a)_(1j)
                                        vmult_dd(1,n)=0
                                        vmult_dd(2,n)=(1.0-2*((real(psin(diff_state,n)))**2))*occupancy(2,n)
                                    else
                                        breakflag=n !before diff alive:a^(a)_(1j)*a^(a)_(1j) dead:a^(a)_(0j)*a^(a)_(0j)
                                    end if          !Unless an opeator acts at position j this evaluates to 0
                                end if
                            
                                gg_1(1:norb)=(0.0,0.0)
                                hh_1(1:norb)=(0.0,0.0)
                                gmax1=norb
                                gg_1(1)=vmult_dd(2,1)-vmult_dd(1,1)
                
                                do p=2, norb
                                    gg_1(p)=gg_1(p-1)*(vmult_dd(2,p)-vmult_dd(1,p))
                                    if(gg_1(p)==(0.0))then
                                        gmax1=p
                                        EXIT 
                                    end if
                                end do
                                
                                hmin1=0
                                hh_1(norb) = vmult_dd(2,norb)+vmult_dd(1,norb)
                                do p=(norb-1),1,-(1)
                                    hh_1(p)=hh_1(p+1)*(vmult_dd(2,p)+vmult_dd(1,p))
                                    if(hh_1(p)==(0.0))then
                                        hmin1=p
                                        EXIT 
                                    end if
                                end do
                
                                tot=0.0
                                if (gmax1 < hmin1) then
                                    temp(n)=tot
                                    cycle
                                end if
                
                                if((breakflag.eq.0).or.(breakflag.eq.1))then
                                    if(ph2ei(j,k,l,1).ne.0) then
                                        if(breakflag.eq.1)then
                                            temp(n)=(1.0-2*((real(psin(diff_state,1)))**2))*&
                                            occupancy(2,1)*hh_1(2)*ph2ei(j,k,l,1)
                                            cycle
                                        end if
                                        if(n.eq.1)then
                                            if((j.eq.1).or.(k.eq.1))then
                                                if(l.ne.1)then
                                                    tot = tot+&
                                                    (2*real(psin(diff_state,1)*pcos(diff_state,1))*occupancy(2,1)*&
                                                    hh_1(2)*ph2ei(j,k,l,1))
                                                end if
                                            else
                                                if(l.ne.1)then
                                                    tot = tot+&
                                                    ((1.0-2*((real(psin(diff_state,1)))**2))*occupancy(2,1)*&
                                                    hh_1(2)*ph2ei(j,k,l,1))
                                                end if 
                                            end if
                                        else
                                            tot = tot+&
                                            (REAL(conjg(z1jk(j,k,2,1))*z2l(l,1,1))*hh_1(2)*ph2ei(j,k,l,1))
                                        end if
                                    end if
                                else if((breakflag.lt.norb))then
                                    if(breakflag.ne.0)then
                                        temp(n)=(gg_1(n-1)*((1.0-2*((real(psin(diff_state,n)))**2))*&
                                        occupancy(2,n))*hh_1(n+1)*ph2ei(j,k,l,n))
                                        cycle
                                    end if
                                    do p=2,norb-1
                                        if(ph2ei(j,k,l,p).ne.0.0) then
                                            if(p.eq.n)then
                                                if((j.eq.p).or.(k.eq.p))then
                                                    if(l.ne.p)then
                                                        tot = tot+&
                                                        (gg_1(p-1)*(2*real(psin(diff_state,p)*pcos(diff_state,p))*&
                                                        occupancy(2,p))*hh_1(p+1)*ph2ei(j,k,l,p))
                                                    end if
                                                else
                                                    if(l.ne.p)then
                                                        tot = tot+&
                                                        (gg_1(p-1)*((1.0-2*((real(psin(diff_state,p)))**2))*&
                                                        occupancy(2,p))*hh_1(p+1)*ph2ei(j,k,l,p))
                                                    end if 
                                                end if
                                            else
                                                tot = tot+&
                                                (gg_1(p-1)*REAL(conjg(z1jk(j,k,2,p))*z2l(l,1,p))*hh_1(p+1)*ph2ei(j,k,l,p))
                                            end if
                                        end if
                                    end do
                                else
                                    if(breakflag.eq.norb)then
                                        temp(n)=(gg_1(norb-1)*((1.0-2*((real(psin(diff_state,norb)))**2))&
                                        *occupancy(2,norb))*ph2ei(j,k,l,norb))
                                        cycle
                                    end if
                                    if(ph2ei(j,k,l,norb).ne.0) then
                                        if(norb.eq.n)then
                                            if((j.eq.norb).or.(k.eq.norb))then
                                                if(l.ne.norb)then
                                                    tot= tot+&
                                                    (gg_1(norb-1)*(2*real(psin(diff_state,norb)*pcos(diff_state,norb))*&
                                                    occupancy(2,norb))*ph2ei(j,k,l,norb))
                                                end if
                                            else
                                                if(l.ne.norb)then
                                                    tot = tot+&
                                                     (gg_1(norb-1)*((1.0-2*((real(psin(diff_state,norb)))**2))*&
                                                     occupancy(2,norb))*ph2ei(j,k,l,norb))
                                                end if    
                                            end if
                                        else
                                            tot = tot+&
                                            (gg_1(norb-1)*REAL(conjg(z1jk(j,k,2,norb))*z2l(l,1,norb))*ph2ei(j,k,l,norb))
                                        end if
                                    end if
                                end if
                                temp(n)=tot
                            end do
                            h2etot_diff_bra(:)=h2etot_diff_bra(:)+temp(:)
                        else  !Differentiation when zombie states are not the same
                            temp=0
                            vmultr=real(z1jk(j,k,:,:))*real(z2l(l,:,:))
                            do n=1, norb
                                vmult_dd=vmultr
                                if(j.eq.k)then
                                    temp(n)=0
                                    cycle
                                else if((j.eq.n).or.(k.eq.n))then
                                    if(l.eq.n)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                        vmult_dd(1,n)=0
                                        vmult_dd(2,n)=real(pcos(diff_state,n)*psin(m,n))*occupancy(2,n)
                                    else !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(0j)
                                        vmult_dd(1,n)=0
                                        vmult_dd(2,n)=real(pcos(diff_state,n)*pcos(m,n))*occupancy(2,n)
                                    end if
                                else if((j.ne.n).or.(k.ne.n))then
                                    if(l.eq.n)then !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                        vmult_dd(1,n)=0
                                        vmult_dd(2,n)=-real(psin(diff_state,n)*psin(m,n))*occupancy(2,n)
                                    else    !before diff alive:a^(a)_(1j)*a^(b)_(1j) dead:a^(a)_(0j)*a^(b)_(0j)
                                        vmult_dd(1,n)=real(pcos(diff_state,n)*psin(m,n))*occupancy(1,n)
                                        vmult_dd(2,n)=-real(psin(diff_state,n)*pcos(m,n))*occupancy(2,n)
                                    end if
                                end if
                
                                gg_1(1:norb)=(0.0,0.0)
                                hh_1(1:norb)=(0.0,0.0)
                                gmax1=norb
                                gg_1(1)=(vmult_dd(2,1))-(vmult_dd(1,1))
                
                                do p=2, norb
                                    gg_1(p)=gg_1(p-1)*((vmult_dd(2,p)-vmult_dd(1,p)))
                                    if(gg_1(p)==(0.0,0.0))then
                                        gmax1=p
                                        EXIT 
                                    end if
                                end do
                
                                
                                hmin1=0
                                
                                hh_1(norb) = (vmult_dd(2,norb)+vmult_dd(1,norb))
                                
                                do p=(norb-1),1,-(1)
                                    hh_1(p)=hh_1(p+1)*(vmult_dd(2,p)+vmult_dd(1,p))
                                    if(hh_1(p)==(0.0,0.0))then
                                        hmin1=p
                                        EXIT 
                                    end if
                                end do
                
                                tot=0.0
                                if((gmax1 < hmin1))then
                                    temp(n)=tot
                                    cycle
                                end if
                
                                if(ph2ei(j,k,l,1).ne.0) then
                                    if(n.eq.1)then
                                        if((j.eq.1).or.(k.eq.1))then
                                            if(l.ne.1)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                                tot = tot+&
                                                ((REAL(pcos(diff_state,1)*psin(m,1))*occupancy(2,1))*hh_1(2)*ph2ei(j,k,l,1))
                                            end if
                                        else
                                            if(l.ne.1)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                                tot = tot+&
                                                ((REAL(-psin(diff_state,1)*psin(m,1))*occupancy(2,1))*hh_1(2)*ph2ei(j,k,l,1))
                                            end if 
                                        end if
                                    else
                                        tot= tot+&
                                        (REAL(conjg(z1jk(j,k,2,1))*z2l(l,1,1))*hh_1(2)*ph2ei(j,k,l,1))
                                    end if
                                end if
                                do p=2,norb-1
                                    if(ph2ei(j,k,l,p).ne.0.0) then
                                        if(p.eq.n)then
                                            if((j.eq.p).or.(k.eq.p))then
                                                if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                                    tot= tot+&
                                                    (gg_1(p-1)*(REAL(pcos(diff_state,p)*psin(m,p))*&
                                                    occupancy(2,p))*hh_1(p+1)*ph2ei(j,k,l,p))
                                                end if
                                            else
                                                if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                                    tot = tot+&
                                                    (gg_1(p-1)*(REAL(-psin(diff_state,p)*psin(m,p))*occupancy(2,p))*&
                                                    hh_1(p+1)*ph2ei(j,k,l,p))
                                                end if 
                                            end if
                                        else
                                            tot = tot+&
                                             (gg_1(p-1)*REAL(conjg(z1jk(j,k,2,p))*z2l(l,1,p))*hh_1(p+1)*ph2ei(j,k,l,p))
                                        end if
                                    end if
                                end do
                
                                if(ph2ei(j,k,l,norb).ne.0) then
                                    if(norb.eq.j)then
                                        if((j.eq.norb).or.(k.eq.norb))then
                                            if(l.ne.norb)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                                tot = tot+&
                                                (gg_1(norb-1)*(REAL(pcos(diff_state,norb)*psin(m,norb))*&
                                                occupancy(2,norb))*ph2ei(j,k,l,norb))
                                            end if
                                        else
                                            if(l.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                                tot =tot+&
                                                (gg_1(norb-1)*(REAL(-psin(diff_state,norb)*psin(m,norb))*&
                                                occupancy(2,norb))*ph2ei(j,k,l,norb))
                                            end if 
                                        end if
                                    else
                                        tot=tot+&
                                        (gg_1(norb-1)*REAL(conjg(z1jk(j,k,2,norb))*z2l(l,1,norb))*ph2ei(j,k,l,norb))
                                    end if
                                end if
                                temp(n)=tot
                            end do
                            h2etot_diff_bra=h2etot_diff_bra+temp
                        end if
                    end do
                end do
            end do
            !$omp end distribute parallel do 
            h2etot_diff_bra = h2etot_diff_bra*0.5

            temp =h1etot_diff_bra+h2etot_diff_bra
            pdiff_hjk(diff_state,m,:)=temp
            if(m.eq.diff_state)then
                pdiff_ovrlp(diff_state,m,:) = 0
            else
                prod=real((conjg(psin(diff_state,:))*psin(m,:))+(conjg(pcos(diff_state,:))*pcos(diff_state,:)))
                bra_prod=real(conjg(pcos(diff_state,:))*psin(m,:))-real(conjg(psin(diff_state,:))*pcos(m,:))
                !$omp parallel do
                do j=1,norb
                    temp1=prod
                    temp1(j)=bra_prod(j)
                    pdiff_ovrlp(diff_state,m,j)=product(temp1)
                end do
            end if
        end do  
        !$omp end target teams

        !$omp target update from(pinv,pkinvh)

        !$omp parallel do
        do k=1, ndet
            do l=1, ndet
                if(l.eq.diff_state)then
                    temp2(k,diff_state,:)=matmul(REAL(pinv(k,:)),pdiff_ovrlp(diff_state,:,:))
                else
                    temp2(k,l,:)=real(pinv(k,l))*pdiff_ovrlp(diff_state,l,:)
                end if
            end do
        end do
        !$omp parallel do
        do k=1, ndet
            do l=1, ndet
                pdiff_invh(diff_state,k,l,:)=matmul(transpose(temp2(k,:,:)),real(pkinvh(:,l)))*(-1)
            end do
        end do
        !$omp target update to(pdiff_invh)

        deallocate(z1jk,stat=ierr)
        if(ierr==0) deallocate(z2l,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        return

    end subroutine gradient_row_gpu

    subroutine grad_calc_gpu(psin,pcos,pphi,ph1ei,ph2ei,phjk,pkinvh,pinv,pdiff_hjk,pvars,pgrad_avlb,&
        pdiff_ovrlp,pdiff_invh,pick,occupancy_2an,occupancy_an_cr,occupancy_an,dvec,en,d_diff_flg,haml)

        implicit none 

        real(kind=8),pointer,dimension(:,:),intent(in)::ph1ei
        real(kind=8), dimension(:,:,:,:),pointer,intent(in)::ph2ei
        complex(kind=8), dimension(:,:),intent(in)::psin
        complex(kind=8), dimension(:,:),intent(in)::pcos
        real(kind=8),dimension(:,:),intent(in)::pphi
        complex(kind=8), dimension(:,:), pointer,intent(inout)::phjk
        complex(kind=8), dimension(:,:), pointer,intent(inout)::pinv
        complex(kind=8), dimension(:,:), pointer,intent(inout)::pkinvh
        real(kind=8), dimension(:,:,:), pointer,intent(inout)::pdiff_hjk
        real(kind=8), dimension(:,:,:), pointer,intent(inout)::pdiff_ovrlp
        real(kind=8), dimension(:,:,:,:), pointer,intent(inout)::pdiff_invh
        real(kind=8),dimension(:,:), pointer,intent(inout)::pvars
        integer,dimension(:),pointer,intent(inout)::pgrad_avlb
        type(hamiltonian),intent(inout)::haml
        type(dvector),dimension(:),intent(inout)::dvec
        type(energy),intent(inout)::en
        integer,dimension(:,:,:),intent(in)::occupancy_an
        integer,dimension(:,:,:,:),intent(in)::occupancy_an_cr,occupancy_2an
        integer,intent(in)::pick,d_diff_flg
        integer::k,l,m
        logical::nanchk
        
       
        if (errorflag .ne. 0) return

        call gradient_row_gpu(psin,pcos,pphi,ph1ei,ph2ei,pkinvh,pinv,pdiff_hjk,&
        pdiff_ovrlp,pdiff_invh,occupancy_an,occupancy_an_cr,occupancy_2an,pick,ndet)
        
        nanchk=.false.
        !$omp target update from(pvars)
        do k=1,norb
            if(is_nan(pvars(pick,k)).eqv..true.)then
                nanchk=.true. 
                pvars=0 
                pgrad_avlb=0
                do l=1 ,10
                    !!$omp declare target( gradient_row_gpu)
                    call gradient_row_gpu(psin,pcos,pphi,ph1ei,ph2ei,pkinvh,pinv,pdiff_hjk,&
                    pdiff_ovrlp,pdiff_invh,occupancy_an,occupancy_an_cr,occupancy_2an,pick,ndet)
                    !$omp target update from(pvars)
                    do m=1,norb
                        if(is_nan(pvars(pick,m)).eqv..true.)then 
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
        
        en%erg=0
        en%t=0
        !$omp target update if(d_diff_flg.eq.1) from(pdiff_hjk,pdiff_ovrlp,pdiff_invh)
        if(d_diff_flg.eq.1)then
            call imgtime_prop(dvec,en,haml,pick)
        end if

        call final_grad_gpu(pvars,phjk,pdiff_hjk,dvec(1)%d_diff,dvec(1)%d,pick,d_diff_flg)

       
       
        if(d_diff_flg.eq.0)then
            pgrad_avlb(pick)=1
        else if(d_diff_flg.eq.1)then 
            pgrad_avlb(pick)=2
        end if
       
        return

    end subroutine grad_calc_gpu

    subroutine orbital_gd_gpu(temp_ham,haml, psin,pcos,pphi,pvars,pprev_erg,pcurrent_erg,pgrad_avlb,pprev_mmntm,ph1ei,ph2ei,& 
        phnuc,dvecs,temp_dvecs,en,phjk,povrlp,pkinvh,pinv,pdiff_hjk,pdiff_ovrlp,pdiff_invh,phjkt,povrlpt,pkinvht,pinvt,chng_trk,&
        psint,pcost,pphit,occupancy_an,occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker,maxloop,zstore)
        
       
        implicit none 
        type(zombiest),dimension(:),intent(inout)::zstore
        real(kind=8),pointer,dimension(:,:),intent(in)::ph1ei
        real(kind=8), dimension(:,:,:,:),pointer,intent(in)::ph2ei
        real(kind=8),pointer,intent(in) :: phnuc
        complex(kind=8), dimension(ndet,norb),intent(inout)::psin,psint
        complex(kind=8), dimension(ndet,norb),intent(inout)::pcos,pcost
        real(kind=8),dimension(ndet,norb)::pphi, pphit
        complex(kind=8), dimension(:,:), pointer,intent(inout)::phjk,phjkt
        complex(kind=8), dimension(:,:), pointer,intent(inout)::povrlp,povrlpt
        complex(kind=8), dimension(:,:), pointer,intent(inout)::pinv,pinvt
        complex(kind=8), dimension(:,:), pointer,intent(inout)::pkinvh,pkinvht
        real(kind=8), dimension(:,:,:), pointer,intent(inout)::pdiff_hjk
        real(kind=8), dimension(:,:,:), pointer,intent(inout)::pdiff_ovrlp       
        real(kind=8), dimension(:,:,:,:), pointer,intent(inout)::pdiff_invh
        real(kind=8),dimension(:,:), pointer,intent(inout)::pvars
        real(kind=8),pointer,intent(inout):: pprev_erg
        real(kind=8),pointer,intent(inout):: pcurrent_erg
        integer,dimension(:),pointer,intent(inout)::pgrad_avlb
        real(kind=8),dimension(:,:),pointer,intent(inout)::pprev_mmntm
        type(dvector),dimension(:),intent(inout)::dvecs,temp_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::temp_ham,haml
        integer,dimension(:),intent(inout)::chng_trk
        integer,intent(inout)::epoc_cnt
        integer,intent(in)::maxloop
        integer,dimension(:,:,:)::occupancy_an
        integer,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an
        real(kind=8),intent(in)::b,alphain
        integer,dimension(ndet-1),intent(inout)::picker
        integer::rjct_cnt,next,acpt_cnt,pick,pickorb,rjct_cnt2,loops,d_diff_flg
        real(kind=8)::t,fxtdk,l2_rglrstn
        integer,dimension(ndet-1)::rsrtpass
        integer,dimension(ndet,norb)::lralt_zs
        integer,dimension(norb)::pickerorb
        real(kind=8),dimension(ndet,norb)::alpha_zs,newb_zs
        integer,dimension(norb)::chng_trk2
        logical::nanchk
        character(len=4)::ergerr
        integer::j,k,l,n
       
        alpha_zs=alphain  ! learning rate reduction
        lralt_zs=0    ! power alpha is raised to 
        newb_zs=b
        chng_trk=0 !stores which if any ZS changed
        chng_trk2=0 !stores which orbitals in the ZS have changed 
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        l2_rglrstn=1! !L2 regularisation lambda paramter
        rjct_cnt2=0
        loops=0
        d_diff_flg=1
        
        do while(rjct_cnt2.lt.(norb*100))
            loops=loops+1
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent steps &
                           &   | Orbitals altered '
            do j=1,(ndet-1)
                pcurrent_erg=pprev_erg
                pick=picker(j)
                pickerorb=scramble_norb(norb)
                ! pickerorb=(/1,2,3,4,5,6,7,8,9,10/)
                chng_trk2=0
                acpt_cnt=0
                do n=1,norb
                    rjct_cnt=0
                    pickorb=pickerorb(n)
                    t=newb_zs(pick,pickorb)*(alpha_zs(pick,pickorb)**lralt_zs(pick,pickorb))
                    do while(t.gt.(1.0d-13))
                        nanchk=.false.
                        if(is_nan(pvars(pick,pickorb)).eqv..true.)then
                            call grad_calc_gpu(psin,pcos,pphi,ph1ei,ph2ei,phjk,pkinvh,pinv,pdiff_hjk,pvars,pgrad_avlb,&
                            pdiff_ovrlp,pdiff_invh,pick,occupancy_2an,occupancy_an_cr,occupancy_an,dvecs,en,d_diff_flg,haml)
                        end if 
                       
                        !$omp target map(to:t,l2_rglrstn,pick,pickorb)
                       
                        ! Setup temporary zombie state
                        psint=psin
                        pcost=pcos
                        pphit=pphi
                      
                        pphit(pick,pickorb)=pphi(pick,pickorb)-(t*(pvars(pick,pickorb)))
                        pphit(pick,pickorb)=pphit(pick,pickorb)+l2_rglrstn*((pvars(pick,pickorb))*pvars(pick,pickorb))
                        if(is_nan(pphit(pick,pickorb)).eqv..true.)then
                            t=(1.0d-14)
                        end if
                        psint(pick,pickorb)=sin(cmplx(pphit(pick,pickorb),0.0d0,kind=8))
                        pcost(pick,pickorb)=cos(cmplx(pphit(pick,pickorb),0.0d0,kind=8))
                        phjkt=phjk
                        povrlpt=povrlp
                     
                       
                        !$omp end target
                        call he_full_row_gpu(ph1ei,ph2ei,phnuc,psint,pcost,occupancy_2an,&
                        occupancy_an_cr,occupancy_an,phjkt,povrlpt,pinvt,pkinvht,pick,ndet)
               
                        !$omp target update from(phjkt,povrlpt,pkinvht)
          
                        ! Imaginary time propagation for back tracing
                        
                        en%erg=0
                        en%t=0
                        call imgtime_prop(temp_dvecs,en,temp_ham,0)
                        fxtdk=real(en%erg(1,timesteps+1))
                    
                        
                        if(is_nan(fxtdk).eqv..true.)then 
                            ergerr='NaN ' 
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb 
                            t=(1.0d-14)
                            nanchk=.true.
                            rjct_cnt=1
                        else if(is_posinf(fxtdk).eqv..true.)then
                            ergerr='+NaN'  
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb  
                            t=(1.0d-14)
                            nanchk=.true.
                            rjct_cnt=1
                        else if(is_neginf(fxtdk).eqv..true.)then
                            ergerr='-NaN'  
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb 
                            t=(1.0d-14)
                            nanchk=.true.
                            rjct_cnt=1
                        end if
                        
                        if(nanchk.eqv..false.)then
                            ! Check if energy is lower and accept or reject
                            if(fxtdk.lt.pprev_erg)then
                                t=(1.0d-14)
                                acpt_cnt=acpt_cnt+1
                                !$omp target
                                psin=psint
                                pcos=pcost
                                pphi=pphit
                                phjk=phjkt
                                povrlp=povrlpt
                                pinv=pinvt
                                pkinvh=pkinvht
                                pvars=0.0
                                pdiff_hjk=0
                                pdiff_invh=0
                                pdiff_ovrlp=0
                                pprev_erg=fxtdk
                                pprev_mmntm(pick,:)=pphi(pick,:)
                                !$omp end target
                                pgrad_avlb=0
                                dvecs=temp_dvecs
                                chng_trk(j)=pick
                                chng_trk2(acpt_cnt)=pickorb
                                rjct_cnt=0
                                rjct_cnt2=0
                                !$omp target update from(psin,pcos,pphi,pprev_erg)
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
                        call grad_calc_gpu(psin,pcos,pphi,ph1ei,ph2ei,phjk,pkinvh,pinv,pdiff_hjk,pvars,pgrad_avlb,&
                        pdiff_ovrlp,pdiff_invh,pick,occupancy_2an,occupancy_an_cr,occupancy_an,dvecs,en,d_diff_flg,haml)
                    end if
                end do
                

                if(j.eq.(ndet-1))then 
                    epoc_cnt=epoc_cnt+1
                    picker=scramble(ndet-1)
                    next=picker(1)
                    lralt_zs=0
                    !Every 100 epoc brings phi values back within the normal 0-2pi range
                    if(modulo(epoc_cnt,100).eq.0)then 
                        do k=2,ndet 
                            do l=1,norb 
                                if((pphi(k,l).gt.2*pirl).or.(pphi(k,l).lt.0))then 
                                    pphi(k,l)=asin(real(psin(k,l)))
                                    if((pphi(k,l).gt.2*pirl))then
                                        pphi(k,l)=pphi(k,l)-2*pirl
                                    else if((pphi(k,l).lt.0))then
                                        pphi(k,l)=pphi(k,l)+2*pirl
                                    end if
                                end if
                            end do
                        end do
                    end if
                else 
                    next=picker(j+1)
                end if

                if(pgrad_avlb(next).eq.1)then  
                    do k=1,norb
                        if(is_nan(pvars(next,k)).eqv..true.)then
                            pgrad_avlb(next)=0
                            exit 
                        end if 
                    end do
                end if

                !Set up gradients for next pass
                if(pgrad_avlb(next).eq.0)then 
                    call grad_calc_gpu(psin,pcos,pphi,ph1ei,ph2ei,phjk,pkinvh,pinv,pdiff_hjk,pvars,pgrad_avlb,&
                    pdiff_ovrlp,pdiff_invh,next,occupancy_2an,occupancy_an_cr,occupancy_an,dvecs,en,d_diff_flg,haml)
                end if
                if(acpt_cnt.gt.0)then
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,*(i0:','))") '       ', pick,'              ', &
                    pcurrent_erg,'               ',pprev_erg,'            ',chng_trk2(1:acpt_cnt) 
                    rsrtpass(pick)=0
                    zstore(pick)%update_num=zstore(pick)%update_num+1
                    call zombiewriter(zstore(pick),pick,rsrtpass(pick))
                else 
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,i0)") '       ', pick,'              ', &
                    pcurrent_erg,'               ',pprev_erg,'            ',0
                    rjct_cnt2=rjct_cnt2+1
                end if
            end do
           
            call epoc_writer(pprev_erg,epoc_cnt,chng_trk,0)
            write(6,"(a,i0,a,f21.16)") "Energy after epoch no. ",epoc_cnt,": ",pprev_erg
            acpt_cnt=0
            chng_trk=0
            if(loops.ge.maxloop)then 
                exit 
            end if
        end do

    end subroutine orbital_gd_gpu


    subroutine full_zs_gd_gpu(temp_ham,haml,psin,pcos,pphi,pvars,pprev_erg,pcurrent_erg,pgrad_avlb,pprev_mmntm,ph1ei,ph2ei,& 
        phnuc,dvecs,temp_dvecs,en,phjk,povrlp,pkinvh,pinv,pdiff_hjk,pdiff_ovrlp,pdiff_invh,phjkt,povrlpt,pkinvht,pinvt,chng_trk,&
        psint,pcost,pphit,occupancy_an,occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker,zstore)

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore
        real(kind=8),pointer,dimension(:,:),intent(in)::ph1ei
        real(kind=8), dimension(:,:,:,:),pointer,intent(in)::ph2ei
        real(kind=8),pointer,intent(in) :: phnuc
        complex(kind=8), dimension(ndet,norb),intent(inout)::psin,psint
        complex(kind=8), dimension(ndet,norb),intent(inout)::pcos,pcost
        real(kind=8),dimension(ndet,norb)::pphi, pphit
        complex(kind=8), dimension(:,:), pointer,intent(inout)::phjk,phjkt
        complex(kind=8), dimension(:,:), pointer,intent(inout)::povrlp,povrlpt
        complex(kind=8), dimension(:,:), pointer,intent(inout)::pinv,pinvt
        complex(kind=8), dimension(:,:), pointer,intent(inout)::pkinvh,pkinvht
        real(kind=8), dimension(:,:,:), pointer,intent(inout)::pdiff_hjk
        real(kind=8), dimension(:,:,:), pointer,intent(inout)::pdiff_ovrlp
        real(kind=8), dimension(:,:,:,:), pointer,intent(inout)::pdiff_invh
        real(kind=8),dimension(:,:), pointer,intent(inout)::pvars
        real(kind=8),pointer,intent(inout):: pprev_erg
        real(kind=8),pointer,intent(inout):: pcurrent_erg
        integer,dimension(:),pointer,intent(inout)::pgrad_avlb
        real(kind=8),dimension(:,:),pointer,intent(inout)::pprev_mmntm
        type(dvector),dimension(:),intent(inout)::dvecs
        type(dvector), dimension(:),intent(inout) :: temp_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::temp_ham,haml
        integer,dimension(:),intent(inout)::chng_trk
        integer,dimension(:,:,:),intent(in)::occupancy_an
        integer,dimension(:,:,:,:),intent(in)::occupancy_an_cr,occupancy_2an
        integer,intent(inout)::epoc_cnt
        real(kind=8),intent(in)::b,alphain
        integer,dimension(ndet-1),intent(inout)::picker
        integer::lralt,rjct_cnt,next,acpt_cnt,pick,orbitcnt,d_diff_flg,lralt_temp,mmntmflg,loop_max,loop_dwn
        real(kind=8)::newb,t,fxtdk,l2_rglrstn,alpha,mmnmtb
        real(kind=8),dimension(norb)::gradient_norm
        real(kind=8),dimension(ndet)::mmntm,mmntma
        integer,dimension(ndet-1)::rsrtpass
        ! DOUBLE PRECISION, external::ZBQLU01
        real::r
        logical::nanchk
        character(len=4)::ergerr
        integer::j,k,l

        if (errorflag .ne. 0) return


        alpha=alphain  ! learning rate reduction
        lralt=0    ! power alpha is raised to  
        newb=b
        chng_trk=0 !stores which if any ZS changed
        rjct_cnt=0 !tracks how many rejections 
        t=newb*(alpha**lralt) !learning rate
        acpt_cnt=0  !counts how many ZS have been changed
        l2_rglrstn=1! !L2 regularisation lambda paramter
        orbitcnt=0
        mmntm=0
        mmntma=1
        mmntmflg=0
        loop_max=20
        loop_dwn=0 

        do while(rjct_cnt.lt.200)
            t=newb*(alpha**lralt)
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent step &
            &   | Learning rate | Acceptance count | Rejection count'
            
            do j=1,(ndet-1)
                
                pick=picker(j)
                lralt_temp=lralt
                do k=1,norb
                    if(is_nan(pvars(pick,k)).eqv..true.)then
                        call grad_calc_gpu(psin,pcos,pphi,ph1ei,ph2ei,phjk,pkinvh,pinv,pdiff_hjk,pvars,pgrad_avlb,&
                            pdiff_ovrlp,pdiff_invh,pick,occupancy_2an,occupancy_an_cr,occupancy_an,dvecs,en,d_diff_flg,haml)
                        nanchk=.false. 
                        exit 
                    end if 
                end do
                gradient_norm=((pvars(pick,:))*pvars(pick,:))
                do while(lralt_temp.lt.(loop_max))
                    t=newb*(alpha**lralt_temp)
                    nanchk=.false.


                    mmnmtb=(t*mmntm(pick))/mmntma(pick)
                    ! Setup temporary zombie state
                    !$omp target map(to:t,l2_rglrstn,pick,mmnmtb)
                    psint=psin
                    pcost=pcos
                    pphit=pphi
                    
                    pphit(pick,:)=pphi(pick,:)-(t*(pvars(pick,:)))+mmnmtb*(pphi(pick,:)-pprev_mmntm(pick,:))
                    pphit(pick,:)=pphit(pick,:)+l2_rglrstn*((pvars(pick,:))*pvars(pick,:))
                
                    do k=1,norb
                        if(is_nan(pphit(pick,k)).eqv..true.)then
                            pphit(pick,:)=pphi(pick,:)
                            exit
                        end if
                    end do
                    psint(pick,:)=sin(cmplx(pphit(pick,:),0.0d0,kind=8))
                    pcost(pick,:)=cos(cmplx(pphit(pick,:),0.0d0,kind=8))
                    phjkt=phjk
                    povrlpt=povrlp
                    pinvt=pinv
                    pkinvht=pkinvh
                   !$omp end target
                    call he_full_row_gpu(ph1ei,ph2ei,phnuc,psint,pcost,occupancy_2an,&
                        occupancy_an_cr,occupancy_an,phjkt,povrlpt,pinvt,pkinvht,pick,ndet)
                    !$omp target update from(phjkt,povrlpt,pkinvht)

                    ! Imaginary time propagation for back tracing
                    en%erg=0
                    en%t=0
                    call imgtime_prop(temp_dvecs,en,temp_ham,0)
                
                    fxtdk=real(en%erg(1,timesteps+1))
                   
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
            
                    if(nanchk.eqv..false.)then
                        ! Check if energy is lower and accept or reject
                        if(fxtdk.lt.pprev_erg)then
                            acpt_cnt=acpt_cnt+1
                            !$omp target
                            psin=psint
                            pcos=pcost
                            pphi=pphit
                            phjk=phjkt
                            povrlp=povrlpt
                            pinv=pinvt
                            pkinvh=pkinvht
                            pvars=0.0
                            pdiff_hjk=0
                            pdiff_invh=0
                            pdiff_ovrlp=0
                            pprev_erg=fxtdk
                            d_diff_flg=0
                            pprev_mmntm(pick,:)=pphi(pick,:)
                            rsrtpass(pick)=0
                            chng_trk(j)=pick
                            !$omp end target
                            pgrad_avlb=0
                            dvecs=temp_dvecs
                            rjct_cnt=0
                            mmntma(pick)=t
                            mmntm(pick)=mmnmtb
                            loop_dwn=loop_dwn+1
                            !$omp target update from(psin,pcos,pphi,pprev_erg)
                            zstore(pick)%update_num=zstore(pick)%update_num+1
                            call zombiewriter(zstore(pick),pick,rsrtpass(pick))
                            if(orbitcnt.lt.0)then
                                orbitcnt=orbitcnt+1 
                            else
                                orbitcnt=0
                            end if
                            if((loop_max.gt.20).and.(loop_dwn.eq.5))then
                                loop_max=loop_max-1
                                loop_dwn=0
                            end if
                       
                            write(6,"(a,i0,a,f21.16,a,f21.16,a,f12.10,a,i0,a,i0)") '       ', pick,'              ', &
                         pprev_erg,'               ',fxtdk,'            ',t,'          ',acpt_cnt,'                 ',rjct_cnt
                            Exit
                           
                        else 
                            lralt_temp=lralt_temp+1
                        end if
                    else 
                        write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                                                " for zombie state ", pick, ", on epoc ", epoc_cnt 
                    end if
                end do
                if(lralt_temp.ge.loop_max)then
                    rjct_cnt=rjct_cnt+1
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,f12.10,a,i0,a,i0)") '       ', pick,'              ', &
                pprev_erg,'               ',fxtdk,'            ',0.0,'          ',acpt_cnt,'                 ',rjct_cnt
                    if(orbitcnt.ne.0)then
                        orbitcnt=orbitcnt+1
                    end if
                    if(loop_max.lt.35)then 
                        loop_max=loop_max+1
                    end if
                    loop_dwn=0 
                end if
                
               
                t=newb*(alpha**lralt)
                if(j.eq.(ndet-1))then 
                    epoc_cnt=epoc_cnt+1
                    if(acpt_cnt.gt.0)then
                        picker=scramble(ndet-1)
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
        
                if(pgrad_avlb(next).eq.1)then 
                    d_diff_flg=1 
                else 
                    d_diff_flg=0
                    call random_number(r)
                    if(r.lt.0.3)then 
                    ! if(ZBQLU01(1).lt.0.3)then 
                        d_diff_flg=1
                    end if
                end if

                !Set up gradients for next pass
                if(pgrad_avlb(next).lt.2)then
                    call grad_calc_gpu(psin,pcos,pphi,ph1ei,ph2ei,phjk,pkinvh,pinv,pdiff_hjk,pvars,pgrad_avlb,&
                    pdiff_ovrlp,pdiff_invh,next,occupancy_2an,occupancy_an_cr,occupancy_an,dvecs,en,d_diff_flg,haml) 
                end if

            end do

            call epoc_writer(pprev_erg,epoc_cnt,chng_trk,0)
            write(6,"(a,i0,a,f21.16,a,f12.10,a,i0,a)") "Energy after epoc no. ",epoc_cnt,": ", &
                pprev_erg, ". The current max learning rate is: ",t, ". ", acpt_cnt, " Zombie state(s) altered."

            
            ! If only a single ZS is being altered the learning rate is lowered to allow others to be changed
            if(acpt_cnt.eq.1)then
                if((orbitcnt.ge.0).and.(epoc_cnt.gt.100))then  
                    call orbital_gd_gpu(temp_ham,haml,psin,pcos,pphi,pvars,pprev_erg,pcurrent_erg,pgrad_avlb,pprev_mmntm,& 
                    ph1ei,ph2ei,phnuc,dvecs,temp_dvecs,en,phjk,povrlp,pkinvh,pinv,pdiff_hjk,pdiff_ovrlp,pdiff_invh,&
                    phjkt,povrlpt,pkinvht,pinvt,chng_trk,psint,pcost,pphit,occupancy_an,occupancy_an_cr,occupancy_2an,&
                    epoc_cnt,alphain,b,picker,1,zstore)
                    orbitcnt=-(10*ndet)
                    lralt=0
                end if   
            else if(acpt_cnt.eq.0)then
                call orbital_gd_gpu(temp_ham,haml,psin,pcos,pphi,pvars,pprev_erg,pcurrent_erg,pgrad_avlb,pprev_mmntm,& 
                ph1ei,ph2ei,phnuc,dvecs,temp_dvecs,en,phjk,povrlp,pkinvh,pinv,pdiff_hjk,pdiff_ovrlp,pdiff_invh,&
                phjkt,povrlpt,pkinvht,pinvt,chng_trk,psint,pcost,pphit,occupancy_an,occupancy_an_cr,occupancy_2an,&
                epoc_cnt,alphain,b,picker,1,zstore)
                orbitcnt=-(10*ndet)
            end if

            if((epoc_cnt.gt.500).and.(orbitcnt.ge.0).and.(modulo(epoc_cnt,100).eq.0))then
                call orbital_gd_gpu(temp_ham,haml,psin,pcos,pphi,pvars,pprev_erg,pcurrent_erg,pgrad_avlb,pprev_mmntm,& 
                ph1ei,ph2ei,phnuc,dvecs,temp_dvecs,en,phjk,povrlp,pkinvh,pinv,pdiff_hjk,pdiff_ovrlp,pdiff_invh,&
                phjkt,povrlpt,pkinvht,pinvt,chng_trk,psint,pcost,pphit,occupancy_an,occupancy_an_cr,occupancy_2an,&
                epoc_cnt,alphain,b,picker,1,zstore)
                orbitcnt=-(10*ndet)
                lralt=0
            end if

            chng_trk=0
            t=newb*(alpha**lralt)
            if(mmntmflg.eq.0)then 
                mmntm=0.4
                mmntmflg=1
            end if 
            acpt_cnt=0
            if(epoc_cnt.gt.30000)then 
                exit 
            end if
        end do

    end subroutine full_zs_gd_gpu

    subroutine zombie_alter_gpu(zstore,grad_fin,haml,elect,en,dvecs,chng_trk)

        implicit none
       
        type(zombiest),dimension(:),intent(inout)::zstore
        type(grad),intent(inout),target::grad_fin
        type(elecintrgl),intent(in),target::elect
        type(dvector),dimension(:),intent(inout)::dvecs
        type(dvector), dimension(:), allocatable:: temp_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout),target::haml
        integer,dimension(:),intent(inout)::chng_trk
        type(hamiltonian),target::temp_ham
        integer,allocatable,dimension(:,:,:)::occupancy_an
        integer,allocatable,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an
        integer::epoc_cnt,epoc_max
        real(kind=8)::alpha,b
        integer,dimension(ndet-1)::picker,rsrtpass
       
        real(kind=8),pointer,dimension(:,:)::ph1ei
        real(kind=8), dimension(:,:,:,:), pointer::ph2ei
        real(kind=8),pointer :: phnuc

        complex(kind=8), dimension(ndet,norb)::psin,psint
        complex(kind=8), dimension(ndet,norb)::pcos,pcost
        real(kind=8),dimension(ndet,norb)::pphi, pphit

        complex(kind=8), dimension(:,:), pointer::phjk,phjkt
        complex(kind=8), dimension(:,:), pointer::povrlp,povrlpt
        complex(kind=8), dimension(:,:), pointer::pinv,pinvt
        complex(kind=8), dimension(:,:), pointer::pkinvh,pkinvht

        real(kind=8), dimension(:,:,:), pointer::pdiff_hjk
        real(kind=8), dimension(:,:,:), pointer::pdiff_ovrlp
        real(kind=8), dimension(:,:,:,:), pointer::pdiff_invh

        real(kind=8),dimension(:,:), pointer::pvars
        real(kind=8),pointer:: pprev_erg
        real(kind=8),pointer:: pcurrent_erg
        integer,dimension(:),pointer::pgrad_avlb
        real(kind=8),dimension(:,:),pointer::pprev_mmntm
        integer::ierr,j,k,l
        
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

            call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,1)
            rsrtpass=1
            ierr=0
        end if

        call allocham(temp_ham,ndet,norb)
        call allocdv(temp_dvecs,1,ndet,norb)
        temp_ham=haml 

        ph1ei=>elect%h1ei
        ph2ei=>elect%h2ei
        phnuc=>elect%hnuc 

        phjk=>haml%hjk
        povrlp=>haml%ovrlp
        pinv=>haml%inv
        pkinvh=>haml%kinvh
        pdiff_hjk=>haml%diff_hjk
        pdiff_ovrlp=>haml%diff_ovrlp
        pdiff_invh=>haml%diff_invh

        phjkt=>temp_ham%hjk
        povrlpt=>temp_ham%ovrlp
        pinvt=>temp_ham%inv
        pkinvht=>temp_ham%kinvh

        pvars=>grad_fin%vars
        pprev_erg=>grad_fin%prev_erg
        pcurrent_erg=>grad_fin%current_erg
        pgrad_avlb=>grad_fin%grad_avlb
        pprev_mmntm=>grad_fin%prev_mmntm
        

        do j=1,ndet
            psin(j,:)=zstore(j)%sin(:)
            pcos(j,:)=zstore(j)%cos(:)
            pphi(j,:)=zstore(j)%phi(:)
        end do

        psint=psin 
        pcost=pcos
        pphit=pphi

        alpha=0.5  ! learning rate reduction
        b=8.0D0 !starting learning rate
        
        chng_trk=0 !stores which if any ZS changed
        rsrtpass=0
        epoc_max=30000

        do j=1, ndet-1
            picker(j)=j+1
        end do
      
        !$omp target data map(to:ph1ei,ph2ei,phnuc,rsrtpass,occupancy_2an,occupancy_an_cr,occupancy_an,errorflag,&
        !$omp pprev_erg,pprev_mmntm,pdiff_hjk,pdiff_ovrlp,pdiff_invh,pkinvh,pinv,norb,ndet)&
        !$omp & map(alloc:psint(ndet,norb),pcost(ndet,norb),pphit(ndet,norb),phjkt(ndet,ndet),povrlpt(ndet,ndet),&
        !$omp pkinvht(ndet,ndet),pinvt(ndet,ndet)) map(tofrom:phjk,povrlp,psin,pcos,pphi,pvars) 
        
        if(epoc_cnt.eq.0)then
            call orbital_gd_gpu(temp_ham,haml,psin,pcos,pphi,pvars,pprev_erg,pcurrent_erg,pgrad_avlb,pprev_mmntm,& 
            ph1ei,ph2ei,phnuc,dvecs,temp_dvecs,en,phjk,povrlp,pkinvh,pinv,pdiff_hjk,pdiff_ovrlp,pdiff_invh,&
            phjkt,povrlpt,pkinvht,pinvt,chng_trk,psint,pcost,pphit,occupancy_an,occupancy_an_cr,occupancy_2an,&
            epoc_cnt,alpha,b,picker,1,zstore)
        end if 
        
        call full_zs_gd_gpu(temp_ham,haml,psin,pcos,pphi,pvars,pprev_erg,pcurrent_erg,pgrad_avlb,pprev_mmntm,ph1ei,ph2ei,& 
        phnuc,dvecs,temp_dvecs,en,phjk,povrlp,pkinvh,pinv,pdiff_hjk,pdiff_ovrlp,pdiff_invh,phjkt,povrlpt,pkinvht,pinvt,chng_trk,&
        psint,pcost,pphit,occupancy_an,occupancy_an_cr,occupancy_2an,epoc_cnt,alpha,b,picker,zstore)
            
        call orbital_gd_gpu(temp_ham,haml,psin,pcos,pphi,pvars,pprev_erg,pcurrent_erg,pgrad_avlb,pprev_mmntm,& 
        ph1ei,ph2ei,phnuc,dvecs,temp_dvecs,en,phjk,povrlp,pkinvh,pinv,pdiff_hjk,pdiff_ovrlp,pdiff_invh,&
        phjkt,povrlpt,pkinvht,pinvt,chng_trk,psint,pcost,pphit,occupancy_an,occupancy_an_cr,occupancy_2an,&
        epoc_cnt,alpha,b,picker,(epoc_max-epoc_cnt),zstore)
        
        !$omp end target data

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
        deallocate(occupancy_an,stat=ierr)
        if(ierr==0) deallocate(occupancy_2an,stat=ierr)
        if(ierr==0) deallocate(occupancy_an_cr,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in occupancy vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

    end subroutine zombie_alter_gpu


    ! Produces a random order for the ZS to be posisbly changed
    function scramble( number_of_values ) result(out)
        
        implicit none

        integer,intent(in)::number_of_values
        integer,allocatable::out(:),array(:)
        integer::n,m,k,j,l,jtemp
        real::r
        !DOUBLE PRECISION, external::ZBQLU01

        out=[(j,j=1,number_of_values)]
        array=[(j,j=1,number_of_values+1)]
        n=1; m=number_of_values
        do k=1,2
            do j=1,m+1
                call random_number(r)
                l = n + FLOOR((m+1-n)*r) !ZBQLU01(1))
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
        integer,allocatable::out(:)
        integer::n,m,k,j,l,jtemp
        !DOUBLE PRECISION, external::ZBQLU01
        real::r

        out=[(j,j=1,number_of_values)]
        n=1; m=number_of_values
        do k=1,2
            do j=1,m
                call random_number(r)
                l = n + FLOOR((m-n)*r) !ZBQLU01(1))
                jtemp=out(l)
                out(l)=out(j)
                out(j)=jtemp
            end do
        end do

        return
    end function scramble_norb
    
END MODUlE gradient_descent