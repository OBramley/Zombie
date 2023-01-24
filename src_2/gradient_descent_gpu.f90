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
    subroutine he_full_row_gpu(haml,zstore,elecs,diff_state,size,occupancy_2an,occupancy_an_cr,occupancy_an,&
        z2l_all,z1jk,z2l_inout)
        
        implicit none
        type(hamiltonian), intent(inout)::haml
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::size,diff_state
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an,occupancy_an_cr
        integer,dimension(:,:,:),intent(in)::occupancy_an
        complex(kind=8),dimension(:,:,:,:),intent(inout)::z2l_all
        complex(kind=8), dimension(:,:,:,:),intent(inout)::z1jk
        complex(kind=8), dimension(:,:,:),intent(inout)::z2l_inout
        ! complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
        complex(kind=8),allocatable, dimension(:,:,:)::z2l

        integer::j,k,l,m,ierr
        complex(kind=8)::h1etot, h2etot
        integer, allocatable,dimension(:)::IPIV1
        complex(kind=8),allocatable,dimension(:)::WORK1
        complex(kind=8),allocatable,dimension(:,:,:)::tot
        complex(kind=8),allocatable,dimension(:,:)::temp

        if (errorflag .ne. 0) return
        ierr = 0
        allocate(tot(norb,norb,norb),stat=ierr)
        allocate(temp(norb,norb),stat=ierr)
        ! allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        h1etot=cmplx(0.0,0.0)
        h2etot=cmplx(0.0,0.0)


        !$omp parallel shared(zstore, z1jk, z2l) private(j,k,l)
    
        !$omp do simd
        do l=1, norb
            z2l_inout(l,1,:)=zstore(diff_state)%sin(:)
            z2l_inout(l,2,:)=zstore(diff_state)%cos(:)
            z2l_inout(l,2,l)=z2l_inout(l,1,l)
            z2l_inout(l,1,l)=cmplx(0.0,0.0)
        end do


        !$omp do simd collapse(2)
        do j=1, norb
            do k=1, norb
                z1jk(j,k,:,:)=z2l_inout(j,:,:)
                z1jk(j,k,2,k)=z1jk(j,k,1,k)
                z1jk(j,k,1,k)=cmplx(0.0,0.0)
            end do
        end do

        !$omp end parallel
        z1jk=z1jk*occupancy_2an

        do m=1,size
    
            h1etot=cmplx(0.0,0.0)
            h2etot=cmplx(0.0,0.0)
            temp=cmplx(0.0,0.0)
            tot=cmplx(0.0,0.0)

            if(m.eq.diff_state)then
                z2l=z2l_inout
            else 
                z2l=z2l_all(m,:,:,:)
            end if 
    
            !!$omp parallel shared(z2l)
            !!$omp do simd
            ! do l=1, norb
            !     z2l(l,1,:)=zstore(m)%sin(:)
            !     z2l(l,2,:)=zstore(m)%cos(:)
            !     z2l(l,2,l)=z2l(l,1,l)
            !     z2l(l,1,l)=cmplx(0.0,0.0)
            ! end do
            !!$omp end parallel         

            call one_elec_intfc_gpu(zstore(diff_state)%sin,zstore(diff_state)%cos,z2l,elecs%h1ei,temp,occupancy_an_cr)
            
            z2l=z2l*occupancy_an
            !$omp flush(z2l)
        
            call two_elec_intfc_gpu(zstore(diff_state)%sin,zstore(m)%sin,z2l,z1jk,elecs%h2ei,tot)
                 
            do j=1, norb
                do k=1, norb
                    h1etot=h1etot+temp(j,k)
                    do l=1, norb
                        h2etot=h2etot+tot(j,k,l)
                    end do
                end do
            end do

            h2etot=h2etot*0.5

            if(is_nan(real(h1etot)).eqv..true.)then 
                print*,zstore(diff_state)%phi
                stop
            end if 
            if(is_nan(real(h2etot)).eqv..true.)then 
                print*,zstore(diff_state)%phi
                stop
            end if 
            !print*,h1etot,h2etot
            haml%ovrlp(diff_state,m)=overlap(zstore(diff_state),zstore(m))
            haml%ovrlp(m,diff_state)= haml%ovrlp(diff_state,m)
            haml%hjk(diff_state,m)=h1etot+h2etot+(elecs%hnuc*haml%ovrlp(diff_state,m))
            haml%hjk(m,diff_state)=haml%hjk(diff_state,m)
        
        end do
    
        ! deallocate(z1jk,stat=ierr)
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

    end subroutine he_full_row_gpu



    subroutine one_elec_intfc_gpu(zs1sin,zs1cos,z2l,h1ei,temp,occupancy)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs1cos
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        integer,dimension(:,:,:,:),intent(in)::occupancy
        real(kind=8), dimension(:,:), intent(in)::h1ei
        complex(kind=8),dimension(:,:),intent(inout)::temp
        complex(kind=8),allocatable,dimension(:,:)::zomt
        integer::j,k,len,ierr

       
        if (errorflag .ne. 0) return

        len=norb
       
        allocate(zomt(2,len),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in zomt  allocation . ierr had value ", ierr
            errorflag=1
        end if
        
  
        temp=cmplx(0.0,0.0)


        !$omp target teams distribute parallel do simd collapse(2) &
        !$omp & map(to:h1ei(:,:),occupancy(:,:,:,:),z2l(:,:,:),zs1sin(:),zs1cos(:),len) &
        !$omp & map(tofrom:temp(:,:)) map(alloc:zomt(2,len))&
        !$omp & private(j,k,zomt) shared(h1ei,occupancy,z2l,zs1sin,zs1cos,temp)
        do j=1, len
            do k=1, len
                if(h1ei(j,k).ne.(0.0)) then 
                    zomt=z2l(j,:,:)
                    zomt(1,k)=zomt(2,k)
                    zomt(2,k)=cmplx(0.0,0.0)
                    zomt=zomt*occupancy(j,k,:,:)
                    temp(j,k)=product((conjg(zs1sin)*zomt(1,:))+(conjg(zs1cos)*zomt(2,:)))*h1ei(j,k)
                end if
            end do
        end do
        !$omp end target teams distribute parallel do simd

        deallocate(zomt,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in zomt  deallocation . ierr had value ", ierr
            errorflag=1
        end if
      
        return

    end subroutine one_elec_intfc_gpu

    subroutine two_elec_intfc_gpu(zs1sin,zs2sin,z2l,z1jk,h2ei,tot)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs2sin
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:,:),intent(in)::z1jk
        real(kind=8), dimension(:,:,:,:),intent(in)::h2ei
        complex(kind=8),dimension(:,:,:),intent(inout)::tot
        logical,dimension(:,:,:),allocatable::iszero
        complex(kind=8),dimension(:,:),allocatable::vmult
        complex(kind=8),allocatable,dimension(:)::gg,hh
        integer::j,k,l,len,ierr,p,gmax,hmin,q
    
        if (errorflag .ne. 0) return
        ierr=0
        tot=cmplx(0.0,0.0)
        len=norb

        allocate(gg(len),stat=ierr)
        allocate(hh(len),stat=ierr)
        allocate(vmult(2,len),stat=ierr)
        allocate(iszero(norb,norb,(norb/2)),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in vmult allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        iszero=.true.
     
        do j=1, norb 
            if(zs1sin(j)==(0.0,0.0))then
                iszero(j,:,:)=.false.
                cycle
            end if
            do k=1, norb
                if(j.eq.k) then
                    iszero(j,k,:)=.false.
                    cycle
                end if
                if(occ_iszero(z1jk(j,k,:,:)).eqv..true.)then
                    iszero(j,k,:)=.false.
                    CYCLE 
                end if
                do q=1, norb/2
                    if(zs2sin(2*q-(modulo(j,2)) )==(0.0,0.0))then
                        iszero(j,k,q)=.false.
                    end if
                end do 
            end do 
        end do
    
        !$omp target teams distribute parallel do simd collapse(3) &
        !$omp & map(alloc:hh(len),gg(len),vmult(2,len)) map(tofrom:tot) &
        !$omp & map(to:h2ei,z1jk,z2l,gmax,hmin,len,iszero) &
        !$omp & private(gmax,hmin,gg,hh,p,vmult,j,k,l,q) shared(z1jk,z2l,tot,iszero)
        do j=1, len
            do k=1, len   
                do q=1, len/2
                    l=2*q-(modulo(j,2)) 
                    if(iszero(j,k,q).eqv..false.)then
                        tot(j,k,l)=(0.0,0.0)
                    else
                        vmult=conjg(z1jk(j,k,:,:))*(z2l(l,:,:))
                        gg(1:len)=(0.0,0.0)
                        hh(1:len)=(0.0,0.0)
                        gmax=len
                        gg(1)=vmult(2,1)-vmult(1,1)
    
                        do p=2, len
                            gg(p)=gg(p-1)*(vmult(2,p)-vmult(1,p))
                            if(gg(p)==(0.0,0.0))then
                                gmax=p
                                EXIT 
                            end if
                        end do
                        
                        hmin=0
                        hh(len) = vmult(2,len)+vmult(1,len)
                        do p=(len-1),1,-(1)
                            hh(p)=hh(p+1)*(vmult(2,p)+vmult(1,p))
                            if(hh(p)==(0.0,0.0))then
                                hmin=p
                                EXIT 
                            end if
                        end do
    
                        tot(j,k,l)=(0.0,0.0)
                        if (gmax < hmin) then
                            tot(j,k,l)=(0.0,0.0)
                            cycle
                        end if
    
                        if(h2ei(j,k,l,1).ne.0) then
                            tot(j,k,l) = tot(j,k,l)+(conjg(z1jk(j,k,2,1))*z2l(l,1,1)*hh(2)*h2ei(j,k,l,1))
                        end if
    
                        do p=2,len-1
                            if(h2ei(j,k,l,p).ne.0.0) then
                                tot(j,k,l) = tot(j,k,l)+ (gg(p-1)*conjg(z1jk(j,k,2,p))*z2l(l,1,p)*hh(p+1)*h2ei(j,k,l,p))
                            end if
                        end do
    
                        if(h2ei(j,k,l,len).ne.0) then
                            tot(j,k,l) = tot(j,k,l) +(gg(len-1)*conjg(z1jk(j,k,2,len))*z2l(l,1,len)*h2ei(j,k,l,len))
                        end if
                    end if
                end do
            
            end do
        end do
        !$omp end target teams distribute parallel do simd

        deallocate(gg,stat=ierr)
        deallocate(hh,stat=ierr)
        deallocate(vmult,stat=ierr)
        deallocate(iszero,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in vmult allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        return

    end subroutine two_elec_intfc_gpu

    subroutine gradient_row_gpu(haml,zstore,elecs,diff_state,size,occupancy_2an,occupancy_an_cr,occupancy_an,z2l_all,z1jk_all)

        implicit none

        type(hamiltonian), intent(inout)::haml
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::size,diff_state
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an,occupancy_an_cr
        integer,dimension(:,:,:),intent(in)::occupancy_an
        complex(kind=8),dimension(:,:,:,:,:),intent(inout)::z1jk_all
        complex(kind=8),dimension(:,:,:,:),intent(inout)::z2l_all
        ! complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
        complex(kind=8),allocatable, dimension(:,:,:)::z2l
        real(kind=8),allocatable,dimension(:)::h1etot_diff_bra,h2etot_diff_bra
        real(kind=8),allocatable,dimension(:,:,:)::temp2
        real(kind=8),allocatable,dimension(:,:,:)::h1etot_diff
        real(kind=8),allocatable,dimension(:,:,:,:)::h2etot_diff
        real(kind=8),allocatable,dimension(:)::overlap_diff
        integer::j,k,l,m,ierr, equal
        

        
        if (errorflag .ne. 0) return
        ierr = 0

        allocate(h1etot_diff_bra(norb),stat=ierr)
        if(ierr==0) allocate(h2etot_diff_bra(norb),stat=ierr)
        if(ierr==0) allocate(temp2(ndet,ndet,norb),stat=ierr)
        if(ierr==0) allocate(h1etot_diff(norb,norb,norb),stat=ierr)
        if(ierr==0) allocate(h2etot_diff(norb,norb,norb,norb),stat=ierr)
        if(ierr==0) allocate(overlap_diff(norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in gradient array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        ! allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        h1etot_diff_bra=0.0
        h2etot_diff_bra=0.0
    
        ! !$omp parallel shared(zstore,z1jk, z2l) private(j,k,l)
        ! !$omp do simd
        ! do l=1, norb
        !     z2l(l,1,:)=zstore(diff_state)%sin(:)
        !     z2l(l,2,:)=zstore(diff_state)%cos(:)
        !     z2l(l,2,l)=z2l(l,1,l)
        !     z2l(l,1,l)=cmplx(0.0,0.0)
        ! end do
   
        ! !$omp do simd collapse(2)
        ! do j=1, norb
        !     do k=1, norb
        !         z1jk(j,k,:,:)=z2l(j,:,:)
        !         z1jk(j,k,2,k)=z1jk(j,k,1,k)
        !         z1jk(j,k,1,k)=cmplx(0.0,0.0)
        !     end do
        ! end do

        ! !$omp end parallel
        ! z1jk=z1jk*occupancy_2an

        do m=1,size

            h1etot_diff_bra=0.0
            h2etot_diff_bra=0.0
            h1etot_diff=0.0
            h2etot_diff=0.0
            if(diff_state.eq.m)then 
                equal=1
            else
                equal=2
            end if
            
            ! !$omp parallel shared(z2l)
            ! !$omp do simd
            ! do l=1, norb
            !     z2l(l,1,:)=zstore(m)%sin(:)
            !     z2l(l,2,:)=zstore(m)%cos(:)
            !     z2l(l,2,l)=z2l(l,1,l)
            !     z2l(l,1,l)=cmplx(0.0,0.0)
            ! end do
            ! !$omp end parallel      

            call one_elec_intfc_grad_gpu(zstore(diff_state)%sin,zstore(diff_state)%cos,&
            zstore(m)%sin,zstore(m)%cos,z2l_all(m,:,:,:),occupancy_an_cr,elecs%h1ei,equal,h1etot_diff)
          
            z2l=z2l_all(m,:,:,:)*occupancy_an
            !$omp flush(z2l) 
            
            call two_elec_intfc_grad_gpu(zstore(diff_state)%sin,zstore(diff_state)%cos,&
            zstore(m)%sin,zstore(m)%cos,z2l_all(m,:,:,:),z1jk_all(diff_state,:,:,:,:),&
            elecs%h2ei,occupancy_2an,occupancy_an,equal,h2etot_diff)
            
            do j=1, norb
                do k=1,norb
                    h1etot_diff_bra=h1etot_diff_bra+h1etot_diff(j,k,:)
                    do l=1, norb
                        h2etot_diff_bra=h2etot_diff_bra+h2etot_diff(j,k,l,:)
                    end do 
                end do
            end do
                    
            haml%diff_hjk(diff_state,m,:)=h1etot_diff_bra+(h2etot_diff_bra*0.5)
            if(m.eq.diff_state)then
                haml%diff_ovrlp(diff_state,m,:) = 0
            else
                overlap_diff = diff_overlap_gpu(zstore(diff_state)%sin,zstore(diff_state)%cos,zstore(m)%sin,zstore(m)%cos,2) 
                haml%diff_ovrlp(diff_state,m,:) = overlap_diff
            end if

        end do  

        ! deallocate(z1jk,stat=ierr)
        if(ierr==0) deallocate(z2l,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        !$omp parallel shared(temp2,haml) private(k,l)
        !$omp do simd collapse(2)
        do k=1, ndet
            do l=1, ndet
                if(l.eq.diff_state)then
                    temp2(k,diff_state,:)=matmul(REAL(haml%inv(k,:)),haml%diff_ovrlp(diff_state,:,:))
                else
                    temp2(k,l,:)=real(haml%inv(k,l))*haml%diff_ovrlp(diff_state,l,:)
                end if
            end do
        end do
       
        !$omp do simd collapse(2)
        do k=1, ndet
            do l=1, ndet
                haml%diff_invh(diff_state,k,l,:)=matmul(transpose(temp2(k,:,:)),real(haml%kinvh(:,l)))*(-1)
            end do
        end do
        !$omp end parallel

        deallocate(h1etot_diff_bra,stat=ierr)
        if(ierr==0) deallocate(h2etot_diff_bra,stat=ierr)
        if(ierr==0) deallocate(overlap_diff,stat=ierr)
        if(ierr==0) deallocate(temp2,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in gradient array vector dedeallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        return

    end subroutine gradient_row_gpu


    subroutine one_elec_intfc_grad_gpu(zs1sin,zs1cos,zs2sin,zs2cos,z2l,occupancy,h1ei,equal,h1etot_diff)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        integer,dimension(:,:,:,:),intent(in)::occupancy
        real(kind=8), dimension(:,:), intent(in)::h1ei
        integer,intent(in)::equal
        real(kind=8),dimension(:,:,:),intent(inout)::h1etot_diff
        integer::j,k,l,len,ierr
        complex(kind=8),allocatable,dimension(:,:)::zomt
        real(kind=8),allocatable,dimension(:)::chng_prod,temp_prod,prod

       
        if (errorflag .ne. 0) return

        h1etot_diff=0.0

        len=norb
        allocate(zomt(2,len),stat=ierr)
        allocate(prod(len))
        allocate(chng_prod(len))
        allocate(temp_prod(len))
        h1etot_diff=0.0
      

        !$omp target teams distribute parallel do simd collapse(2) &
        !$omp & map(to:h1ei(:,:),occupancy(:,:,:,:),z2l(:,:,:),zs1sin(:),zs1cos(:),zs2sin(:),zs2cos(:),equal,len) &
        !$omp & map(tofrom:h1etot_diff(:,:,:)) map(alloc:prod(len),chng_prod(len),temp_prod(len),zomt(2,len)) &
        !$omp & private(j,k,l,prod,chng_prod,temp_prod,zomt) &
        !$omp & shared(h1ei,occupancy,z2l,zs1sin,zs1cos,zs2sin,zs2cos,equal,h1etot_diff)
        do j=1, len
            do k=1, len
                if(h1ei(j,k).ne.(0.0)) then
                    zomt=z2l(j,:,:)
                    zomt(1,k)=zomt(2,k)
                    zomt(2,k)=cmplx(0.0,0.0)
                    zomt=zomt*occupancy(j,k,:,:)
                    prod=real(((conjg(zs1sin)*zomt(1,:)))+((conjg(zs1cos)*zomt(2,:))))
                
                    ! Differntial of overlap of the same ZS 0 unless an operator has acted on ZS
                    if(equal.eq.1)then
                    
                        if(j.eq.k)then
                            if((real(zs2cos(j)).eq.0).and.(real(zs2sin(j)).eq.1).or.(real(zs2sin(j)).eq.0))then
                                h1etot_diff(j,k,j)=0.0
                            else
                                chng_prod=prod           !dead amplitude is zero
                                chng_prod(j)=real(2*zs1sin(j)*zs1cos(j)*occupancy(j,k,1,j))
                                h1etot_diff(j,k,j)=product(chng_prod)*h1ei(j,k)   
                            end if 
                        else if(j.ne.k)then
                            chng_prod=prod
                            if((real(zs2cos(j)).eq.0).and.(real(zs2sin(j)).eq.1)) then
                                chng_prod(j)=-1*occupancy(j,k,2,j)
                            else if((real(zs2sin(j)).eq.0).and.(real(zs2cos(j)).eq.1)) then
                                chng_prod(j)=occupancy(j,k,2,j)
                            else
                                chng_prod(j)=real(((zs1cos(j)**2)-(zs1sin(j)**2))*occupancy(j,k,2,j))
                            end if
                            h1etot_diff(j,k,j)=product(chng_prod)*h1ei(j,k)
                            chng_prod=prod
                            if((real(zs2cos(k)).eq.0).and.(real(zs2sin(k)).eq.1)) then
                                chng_prod(k)=-1*occupancy(j,k,1,k)
                            else if((real(zs2sin(k)).eq.0).and.(real(zs2cos(k)).eq.1)) then
                                chng_prod(k)=occupancy(j,k,1,k)
                            else
                                chng_prod(k)=real(((zs1cos(k)**2)-(zs1sin(k)**2))*occupancy(j,k,1,k))
                            end if               !dead amplitude is zero
                            h1etot_diff(j,k,k)=product(chng_prod)*h1ei(j,k)
                        end if
                    else if(equal.eq.2)then 
                        chng_prod=real(zs1cos*zs2sin*occupancy(j,k,1,:)-zs1sin*zs2cos*occupancy(j,k,2,:))            
                        if(j.eq.k)then  !dead amplitude is zero
                            chng_prod(j)=real(zs1cos(j)*zs2sin(j)*occupancy(j,k,1,j))
                        else                   
                            chng_prod(j)=-real(zs1sin(j)*zs2sin(j))*occupancy(j,k,2,j)
                            chng_prod(k)=real(zs1cos(k)*zs2cos(k))*occupancy(j,k,1,k) 
                        end if
                        do l=1,len
                            temp_prod=prod
                            temp_prod(l)=chng_prod(l)
                            h1etot_diff(j,k,l)=product(temp_prod)*h1ei(j,k)
                        end do
                    else if(equal.eq.3)then
                        chng_prod=real(zs1sin*zs2cos*occupancy(j,k,1,:)-zs1cos*zs2sin*occupancy(j,k,2,:))                  
                        if(j.eq.k)then  !dead amplitude is zero
                            chng_prod(j)=real(zs1sin(j)*zs2cos(j)*occupancy(j,k,1,j))
                        else                   
                            chng_prod(j)=real(zs1cos(j)*zs2cos(j))*occupancy(j,k,2,j)!alive amplitude is zero
                            chng_prod(k)=-real(zs1sin(k)*zs2sin(k))*occupancy(j,k,1,k) !dead amplitude is zero
                        end if
                     
                        do l=1,len
                            temp_prod=prod
                            temp_prod(l)=chng_prod(l)
                            h1etot_diff(j,k,l)=product(temp_prod)*h1ei(j,k)
                        end do
                      
                    end if
                end if
            end do
        end do
        !$omp end target teams distribute parallel do simd
       
     
        deallocate(zomt,stat=ierr)
        deallocate(prod,stat=ierr)
        deallocate(chng_prod,stat=ierr)
        deallocate(temp_prod,stat=ierr)
      
        return


    end subroutine one_elec_intfc_grad_gpu

    subroutine two_elec_intfc_grad_gpu(zs1sin,zs1cos,zs2sin,zs2cos,z2l,z1jk,h2ei,occupancy_2an,occupancy_an,equal,h2etot_diff_in)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:,:),intent(in)::z1jk
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an
        integer,dimension(:,:,:),intent(in)::occupancy_an
        real(kind=8), dimension(:,:,:,:), intent(in)::h2ei
        integer,intent(in)::equal
        real(kind=8),dimension(:,:,:,:),intent(inout)::h2etot_diff_in
        real(kind=8),allocatable,dimension(:,:,:,:)::h2etot_diff
        real(kind=8),allocatable,dimension(:)::gg_0,hh_0
        real(kind=8),dimension(:,:),allocatable::vmult_dd
        integer::j,k,l,len,ierr,q,p
        logical,dimension(:,:,:),allocatable::iszero_grad
        
        
        if (errorflag .ne. 0) return
        ierr=0
        len=norb
        allocate(h2etot_diff(len,len,len,len),stat=ierr)
        h2etot_diff=0.0
        allocate(gg_0(len),stat=ierr)
        allocate(hh_0(len),stat=ierr)
        allocate(vmult_dd(2,len),stat=ierr)
        allocate(iszero_grad(len,len,(len/2)),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in vmult_dd allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        

       

        iszero_grad=.true.
        do j=1, norb 
            if(zs1sin(j)==(0.0,0.0))then
                iszero_grad(j,:,:)=.false.
                cycle
            end if
            do k=1, norb
                if(j.eq.k) then
                    iszero_grad(j,k,:)=.false.
                    cycle
                end if
                do q=1, norb/2
                    if(zs2sin((2*q-(modulo(j,2))))==(0.0,0.0))then
                        iszero_grad(j,k,q)=.false.
                    end if
                end do 
            end do 
        end do
       
        !$omp target teams distribute parallel do simd collapse(4) &
        !$omp & map(alloc:hh_0(len),gg_0(len),vmult_dd(2,len))&
        !$omp & map(tofrom:h2etot_diff) &
        !$omp & map(to:h2ei,z1jk,z2l,zs1sin,zs1cos,zs2sin,zs2cos,len,iszero_grad,equal,occupancy_2an,occupancy_an) &
        !!$omp  parallel do  &
        !$omp & private(gg_0,hh_0,vmult_dd,j,k,l,q,p) &
        !$omp & shared(equal,z1jk,z2l,h2etot_diff,iszero_grad,zs1sin,zs1cos,zs2sin,zs2cos,occupancy_2an,occupancy_an)
        do j=1, len
            do k=1, len   
                do q=1, len/2
                    do p=1, len
                        l=2*q-(modulo(j,2)) 
                        if(iszero_grad(j,k,q).eqv..false.)then
                            h2etot_diff(j,k,l,p)=0.0
                        else
                            vmult_dd=real(conjg(z1jk(j,k,:,:))*(z2l(l,:,:)))
                            if(equal.eq.1)then
                                h2etot_diff(j,k,l,p)=z_an_z3_diff_gpu_1(len,z1jk(j,k,:,:),z2l(l,:,:),&
                                vmult_dd,h2ei(j,k,l,:),(occupancy_2an(j,k,:,:)*occupancy_an(l,:,:)),&
                                j,k,l,zs1sin,zs1cos,gg_0,hh_0)
                            else if(equal.eq.2)then
                                h2etot_diff(j,k,l,p)=z_an_z3_diff_gpu_2(len,z1jk(j,k,:,:),z2l(l,:,:),&
                                vmult_dd,h2ei(j,k,l,:),(occupancy_2an(j,k,:,:)*occupancy_an(l,:,:)),&
                                j,k,l,zs1sin,zs1cos,zs2sin,zs2cos,gg_0,hh_0)
                            else if(equal.eq.3)then
                                h2etot_diff(j,k,l,p)=z_an_z3_diff_gpu_3(len,z1jk(j,k,:,:),z2l(l,:,:),&
                                vmult_dd,h2ei(j,k,l,:),(occupancy_2an(j,k,:,:)*occupancy_an(l,:,:)),&
                                j,k,l,zs1sin,zs1cos,zs2sin,zs2cos,gg_0,hh_0)
                            end if 
                        end if 
                    end do
                end do
            end do
        end do 
        !!$omp end parallel do 
        !$omp end target teams distribute parallel do simd 
        
        h2etot_diff_in=h2etot_diff
        deallocate(h2etot_diff,stat=ierr)
        deallocate(gg_0,stat=ierr)
        deallocate(hh_0,stat=ierr)
        deallocate(vmult_dd,stat=ierr)
        deallocate(iszero_grad,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in vmult_dd deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        return

    end subroutine two_elec_intfc_grad_gpu

    subroutine grad_calc_gpu(haml,zstore,elect,pick,occupancy_2an,occupancy_an_cr,occupancy_an,&
        dvec,grad_fin,en,d_diff_flg,z2l_all,z1jk_all)

        implicit none 

        type(zombiest),dimension(:),intent(in)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvec
        type(hamiltonian),intent(inout)::haml
        integer,dimension(:,:,:),intent(in)::occupancy_an
        integer,dimension(:,:,:,:),intent(in)::occupancy_an_cr,occupancy_2an
        integer,intent(in)::pick,d_diff_flg
        complex(kind=8),dimension(:,:,:,:,:),intent(inout)::z1jk_all
        complex(kind=8),dimension(:,:,:,:),intent(inout)::z2l_all
        integer::k,l,m
        logical::nanchk
        type(energy),intent(inout)::en

        if (errorflag .ne. 0) return

        call gradient_row_gpu(haml,zstore,elect,pick,ndet,occupancy_2an,occupancy_an_cr,occupancy_an,z2l_all,z1jk_all)
        nanchk=.false. 
        do k=1,norb
            if(is_nan(grad_fin%vars(pick,k)).eqv..true.)then
                nanchk=.true. 
                grad_fin%vars=0 
                grad_fin%grad_avlb=0
                do l=1 ,10
                    call gradient_row_gpu(haml,zstore,elect,pick,ndet,occupancy_2an,occupancy_an_cr,occupancy_an,z2l_all,z1jk_all)
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

    end subroutine grad_calc_gpu

    subroutine orbital_gd_gpu(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
        occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker,maxloop,z2l_all,z1jk_all) 

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
        complex(kind=8),dimension(:,:,:,:,:),intent(inout)::z1jk_all
        complex(kind=8),dimension(:,:,:,:),intent(inout)::z2l_all
        complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
        complex(kind=8),allocatable, dimension(:,:,:)::z2l
        integer::rjct_cnt,next,acpt_cnt,pick,pickorb,rjct_cnt2,loops,d_diff_flg,lralt_zs
        real(kind=8)::t,fxtdk,l2_rglrstn
        integer,dimension(ndet-1)::rsrtpass
        integer,dimension(norb)::pickerorb
        integer,dimension(norb)::chng_trk2
        logical::nanchk
        character(len=4)::ergerr
        integer::j,k,l,n,ierr
        
    
        lralt_zs=0    ! power alpha is raised to 
        chng_trk=0 !stores which if any ZS changed
        chng_trk2=0 !stores which orbitals in the ZS have changed 
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        l2_rglrstn=1! !L2 regularisation lambda paramter
        rjct_cnt2=0
        loops=0
        d_diff_flg=1
        ierr=0

        allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        do while(rjct_cnt2.lt.(norb*100))
            loops=loops+1
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent steps &
                        &   | Orbitals altered '
            do j=1,(ndet-1)
               
                grad_fin%current_erg=grad_fin%prev_erg
                pick=picker(j)
                pickerorb=scramble_norb(norb)
                !pickerorb=(/1,2,3,4,5,6,7,8,9,10/)
                chng_trk2=0
                acpt_cnt=0
                do n=1,norb
                   
                    rjct_cnt=0
                    pickorb=pickerorb(n)
                    lralt_zs=0
                    t=b*(alphain**lralt_zs)
                    
                    do while(t.gt.(1.0d-13))
                    
                        nanchk=.false.
                        do k=1,norb
                            if(is_nan(grad_fin%vars(pick,k)).eqv..true.)then
                                call grad_calc_gpu(haml,zstore,elect,pick,occupancy_2an,occupancy_an_cr,&
                            occupancy_an,dvecs,grad_fin,en,d_diff_flg,z2l_all,z1jk_all) 
                                nanchk=.false. 
                                exit 
                            end if 
                        end do


                        if(is_nan(grad_fin%vars(pick,pickorb)).eqv..true.)then
                            call grad_calc_gpu(haml,zstore,elect,pick,occupancy_2an,occupancy_an_cr,&
                            occupancy_an,dvecs,grad_fin,en,d_diff_flg,z2l_all,z1jk_all) 
                        end if 
                        
                        ! Setup temporary zombie state
                        temp_zom=zstore
                        temp_zom(pick)%phi(pickorb)=zstore(pick)%phi(pickorb)-(t*(grad_fin%vars(pick,pickorb)))
                        temp_zom(pick)%phi(pickorb)=temp_zom(pick)%phi(pickorb)+&
                                            l2_rglrstn*((grad_fin%vars(pick,pickorb))*grad_fin%vars(pick,pickorb))
                        do k=1,norb
                            if(is_nan(temp_zom(pick)%phi(k)).eqv..true.)then
                                write(0,"(a,i0,a,i0,a)") "Error in Zombie state ", pick," in orbtial ", k, " which has value NaN"
                                t=(1.0d-14)
                                nanchk=.true.
                                rjct_cnt=1
                                fxtdk=0
                                exit 
                            end if 
                        end do
                  
                        if(nanchk.eqv..false.)then
                            temp_zom(pick)%sin(pickorb)=sin(cmplx(temp_zom(pick)%phi(pickorb),0.0d0,kind=8))
                            temp_zom(pick)%cos(pickorb)=cos(cmplx(temp_zom(pick)%phi(pickorb),0.0d0,kind=8))
                            temp_ham=haml
                            call he_full_row_gpu(temp_ham,temp_zom,elect,pick,ndet,occupancy_2an,&
                            occupancy_an_cr,occupancy_an,z2l_all,z1jk,z2l)
                            ! Imaginary time propagation for back tracing
                            
                            en%erg=0
                            en%t=0
                            call imgtime_prop(temp_dvecs,en,temp_ham,0)
                            fxtdk=real(en%erg(1,timesteps+1))
                        end if
                    
                
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
                                z1jk_all(pick,:,:,:,:)=z1jk
                                z2l_all(pick,:,:,:)=z2l
                            else 
                                lralt_zs=lralt_zs+1
                                t=b*(alphain**lralt_zs)
                                rjct_cnt=rjct_cnt+1
                            end if
                        end if
                    end do
                    
                    if((rjct_cnt.eq.0).and.(n.ne.norb))then
                        call grad_calc_gpu(haml,zstore,elect,pick,occupancy_2an,&
                        occupancy_an_cr,occupancy_an,dvecs,grad_fin,en,d_diff_flg,z2l_all,z1jk_all) 
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
                    call grad_calc_gpu(haml,zstore,elect,next,occupancy_2an,occupancy_an_cr,&
                    occupancy_an,dvecs,grad_fin,en,d_diff_flg,z2l_all,z1jk_all) 
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

        deallocate(z1jk,stat=ierr)
        if(ierr==0) deallocate(z2l,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        return

    end subroutine orbital_gd_gpu


    subroutine full_zs_gd_gpu(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
        occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker,z2l_all,z1jk_all) 

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
        complex(kind=8),dimension(:,:,:,:,:),intent(inout)::z1jk_all
        complex(kind=8),dimension(:,:,:,:),intent(inout)::z2l_all
        complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
        complex(kind=8),allocatable, dimension(:,:,:)::z2l
        integer::lralt,rjct_cnt,next,acpt_cnt,pick,orbitcnt,d_diff_flg,lralt_temp,mmntmflg,loop_max,loop_dwn
        real(kind=8)::newb,t,fxtdk,l2_rglrstn,alpha,mmnmtb
        real(kind=8),dimension(norb)::gradient_norm
        real(kind=8),dimension(ndet)::mmntm,mmntma
        integer,dimension(ndet-1)::rsrtpass
        DOUBLE PRECISION, external::ZBQLU01
        logical::nanchk
        character(len=4)::ergerr
        integer::j,k,l,ierr

    
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
        ierr=0

        allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        do while(rjct_cnt.lt.200)
            t=newb*(alpha**lralt)
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent step &
            &   | Learning rate | Acceptance count | Rejection count'
            
            do j=1,(ndet-1)
                
                pick=picker(j)
                lralt_temp=lralt
                do k=1,norb
                    if(is_nan(grad_fin%vars(pick,k)).eqv..true.)then
                        call grad_calc_gpu(haml,zstore,elect,pick,occupancy_2an,occupancy_an_cr,&
                        occupancy_an,dvecs,grad_fin,en,d_diff_flg,z2l_all,z1jk_all)
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
                            write(0,"(a,i0,a,i0,a)") "Error in Zombie state ", pick," in orbtial ", k, " which has value NaN"
                            t=(1.0d-14)
                            nanchk=.true.
                            rjct_cnt=1
                            fxtdk=0
                            exit
                        end if
                    end do
                    if(nanchk.eqv..false.)then
                        temp_zom(pick)%sin=sin(cmplx(temp_zom(pick)%phi,0.0d0,kind=8))
                        temp_zom(pick)%cos=cos(cmplx(temp_zom(pick)%phi,0.0d0,kind=8))
                        temp_ham=haml
                        call he_full_row_gpu(temp_ham,temp_zom,elect,pick,ndet,occupancy_2an,&
                        occupancy_an_cr,occupancy_an,z2l_all,z1jk,z2l)
                        
                        ! Imaginary time propagation for back tracing
                        en%erg=0
                        en%t=0
                        call imgtime_prop(temp_dvecs,en,temp_ham,0)
                        fxtdk=real(en%erg(1,timesteps+1))
                    end if

                    if(is_nan(fxtdk).eqv..true.)then 
                        nanchk=.true.
                        ergerr='NaN '
                        grad_fin%vars(pick,:)=0.0
                        grad_fin%grad_avlb(pick)=0
                    else if(is_posinf(fxtdk).eqv..true.)then
                        nanchk=.true.
                        ergerr='+NaN'
                        grad_fin%vars(pick,:)=0.0
                        grad_fin%grad_avlb(pick)=0  
                    else if(is_neginf(fxtdk).eqv..true.)then
                        nanchk=.true.
                        ergerr='-NaN'
                        grad_fin%vars(pick,:)=0.0
                        grad_fin%grad_avlb(pick)=0  
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
                            z1jk_all(pick,:,:,:,:)=z1jk
                            z2l_all(pick,:,:,:)=z2l
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
                grad_fin%prev_erg,'               ',fxtdk,'            ',t,'          ',acpt_cnt,'                 ',rjct_cnt
                            Exit
                        
                        else 
                            lralt_temp=lralt_temp+1
                        end if
                    else 
                        write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                                                " for zombie state ", pick, ", on epoc ", epoc_cnt 
                        Exit
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
                    if(ZBQLU01(1).lt.0.3)then 
                        d_diff_flg=1
                    end if
                end if

                !Set up gradients for next pass
                if(grad_fin%grad_avlb(next).lt.2)then
                    call grad_calc_gpu(haml,zstore,elect,next,occupancy_2an,occupancy_an_cr,&
                    occupancy_an,dvecs,grad_fin,en,d_diff_flg,z2l_all,z1jk_all) 
                end if

            end do

            call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,0)
            write(6,"(a,i0,a,f21.16,a,f12.10,a,i0,a)") "Energy after epoc no. ",epoc_cnt,": ", &
                grad_fin%prev_erg, ". The current learning rate is: ",t, ". ", acpt_cnt, " Zombie state(s) altered."

            
            ! If only a single ZS is being altered the learning rate is lowered to allow others to be changed
            if(acpt_cnt.eq.1)then
                if((orbitcnt.ge.0).and.(epoc_cnt.gt.100))then  
                    call orbital_gd_gpu(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
                    occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker,1,z2l_all,z1jk_all)
                    orbitcnt=-(10*ndet)
                    lralt=0
                end if   
            else if(acpt_cnt.eq.0)then
                if((epoc_cnt.gt.50))then  
                    call orbital_gd_gpu(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
                    occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker,1,z2l_all,z1jk_all)
                    orbitcnt=-(10*ndet)
                end if
            end if

            if((epoc_cnt.gt.500).and.(orbitcnt.ge.0).and.(modulo(epoc_cnt,100).eq.0))then
                call orbital_gd_gpu(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
                occupancy_an_cr,occupancy_2an,epoc_cnt,alphain,b,picker,1,z2l_all,z1jk_all)
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

        deallocate(z1jk,stat=ierr)
        if(ierr==0) deallocate(z2l,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        return

    end subroutine full_zs_gd_gpu

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
        complex(kind=8),allocatable,dimension(:,:,:,:,:)::z1jk_all
        complex(kind=8),allocatable,dimension(:,:,:,:)::z2l_all
        integer::epoc_cnt,epoc_max
        real(kind=8)::alpha,b
        integer,dimension(ndet-1)::picker,rsrtpass
    
        integer::ierr,j,k,l,p
        
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

        allocate(z1jk_all(ndet,norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l_all(ndet,norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        !$omp parallel do simd collapse(2)
        do p=1, ndet 
            do l=1, norb
                z2l_all(p,l,1,:)=zstore(p)%sin(:)
                z2l_all(p,l,2,:)=zstore(p)%cos(:)
                z2l_all(p,l,2,l)=z2l_all(p,l,1,l)
                z2l_all(p,l,1,l)=cmplx(0.0,0.0)
                do k=1, norb
                    z1jk_all(p,l,k,:,:)=z2l_all(p,l,:,:)
                    z1jk_all(p,l,k,2,k)=z1jk_all(p,l,k,1,k)
                    z1jk_all(p,l,k,1,k)=cmplx(0.0,0.0)
                end do
                z1jk_all(p,:,:,:,:)=z1jk_all(p,:,:,:,:)*occupancy_2an
            end do
        end do
        !$omp end parallel do simd 

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
        epoc_max=32000

        do j=1, ndet-1
            picker(j)=j+1
        end do

        call alloczs(temp_zom,ndet)
        call allocham(temp_ham,ndet,norb)
        call allocdv(temp_dvecs,1,ndet,norb)

        if(epoc_cnt.eq.0)then
            call orbital_gd_gpu(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
            occupancy_an_cr,occupancy_2an,epoc_cnt,alpha,b,picker,1,z2l_all,z1jk_all)
        end if 
        if(epoc_cnt.lt.epoc_max)then
            call full_zs_gd_gpu(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,&
            chng_trk,temp_zom,occupancy_an,occupancy_an_cr,occupancy_2an,epoc_cnt,alpha,b,picker,z2l_all,z1jk_all) 
            
            call orbital_gd_gpu(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,chng_trk,temp_zom,occupancy_an,&
            occupancy_an_cr,occupancy_2an,epoc_cnt,alpha,b,picker,(epoc_max-epoc_cnt),z2l_all,z1jk_all) 
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
                !call random_number(r)
                l = n + FLOOR((m+1-n)*ZBQLU01(1))
                !l = n + FLOOR((m+1-n)*r) !ZBQLU01(1))
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
                !call random_number(r)
                l = n + FLOOR((m-n)*ZBQLU01(1))
                !l = n + FLOOR((m-n)*r) 
                jtemp=out(l)
                out(l)=out(j)
                out(j)=jtemp

            end do
        end do

        return
    end function scramble_norb

END MODULE gradient_descent