MODULE ham

    use globvars
    use alarrays
    use operators
    use omp_lib

    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators
    subroutine he_row_gpu(ham,zstore,elecs,row,size,occupancy_2an,occupancy_an_cr,occupancy_an,passback)
        
        implicit none
        type(hamiltonian), intent(inout)::ham
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::row,size
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an,occupancy_an_cr
        integer,dimension(:,:,:),intent(in)::occupancy_an
        complex(kind=8),dimension(:,:,:),intent(inout)::passback
        complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
        complex(kind=8),allocatable, dimension(:,:,:)::z2l
        integer::j,k,l,m,ierr,equal
        complex(kind=8)::h1etot, h2etot
        real(kind=8),dimension(norb)::h1etot_diff_bra,h2etot_diff_bra
        real(kind=8),dimension(norb)::h1etot_diff_ket,h2etot_diff_ket
        real(kind=8),dimension(norb)::overlap_diff
        
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
        h1etot_diff_bra=0.0
        h1etot_diff_ket=0.0
        h2etot_diff_bra=0.0
        h2etot_diff_ket=0.0
    
   
        !$omp parallel shared(passback,zstore,z1jk, z2l) private(j,k,l)
        if(row.eq.1)then
            
            !$omp do simd
            do l=1, norb
                z2l(l,1,:)=zstore(1)%sin(:)
                z2l(l,2,:)=zstore(1)%cos(:)
                z2l(l,2,l)=z2l(l,1,l)
                z2l(l,1,l)=cmplx(0.0,0.0)
            end do
            !$omp end do simd
           
        else
            z2l=passback
          
        end if
        !$omp do simd collapse(2)
        do j=1, norb
            do k=1, norb
                z1jk(j,k,:,:)=z2l(j,:,:)
                z1jk(j,k,2,k)=z1jk(j,k,1,k)
                z1jk(j,k,1,k)=cmplx(0.0,0.0)
            end do
        end do
        !$omp end do simd
        !$omp end parallel
        z1jk=z1jk*occupancy_2an

        !!$omp flush(z1jk)
        !!$omp parallel private(j,k,l,z2l,h1etot,h2etot, &
        !!$omp h1etot_diff_bra,h2etot_diff_bra, h1etot_diff_ket,h2etot_diff_ket,&
        !!$omp  h1etot_diff,h2etot_diff,m) shared(z1jk,zstore,row,size,occupancy_2an,occupancy_an_cr,occupancy_an)
        !!$omp do schedule(dynamic)
        do m=row,size
      
            h1etot=cmplx(0.0,0.0)
            h2etot=cmplx(0.0,0.0)
            h1etot_diff_bra=0.0
            h1etot_diff_ket=0.0
            h2etot_diff_bra=0.0
            h2etot_diff_ket=0.0
    
          
            
            if(m.gt.row)then
               
                !$omp parallel shared(z2l)
                !$omp do simd
                do l=1, norb
                    z2l(l,1,:)=zstore(m)%sin(:)
                    z2l(l,2,:)=zstore(m)%cos(:)
                    z2l(l,2,l)=z2l(l,1,l)
                    z2l(l,1,l)=cmplx(0.0,0.0)
                end do
                !$omp end do simd
                !$omp end parallel
                if(m.eq.row+1)then 
                    passback=z2l
                end if
            end if
          
            h1etot=(0.0,0.0)
            h2etot=(0.0,0.0)
            h1etot_diff_bra=0.0
            h1etot_diff_ket=0.0
            h2etot_diff_bra=0.0
            h2etot_diff_ket=0.0
            
            equal=9
            if(2.eq.row)then 
                if(row.eq.m)then 
                    equal=1
                else if(row.ne.m)then 
                    equal=2
                end if
            else if(2.ne.row)then 
                if(m.eq.2)then 
                    equal=3
                end if
            else 
                equal = 9
            end if
        
            
            call one_elec_part_gpu(zstore(row)%sin,zstore(row)%cos,z2l,h1etot,occupancy_an_cr,&
                elecs%h1ei,h1etot_diff_bra,h1etot_diff_ket,zstore(m)%sin,zstore(m)%cos,equal)
            
           
            z2l=z2l*occupancy_an
            !!$omp flush(z2l)
         
            call two_elec_part_gpu(zstore(row)%sin,zstore(row)%cos,z1jk,z2l,h2etot,occupancy_2an,occupancy_an,&
                    elecs%h2ei,h2etot_diff_bra,h2etot_diff_ket,zstore(m)%sin,zstore(m)%cos,equal)
            
            ham%ovrlp(row,m)=overlap(zstore(row),zstore(m))
            ham%ovrlp(m,row)= ham%ovrlp(row,m)
            ham%hjk(row,m)=h1etot+h2etot+(elecs%hnuc*ham%ovrlp(row,m))
            ham%hjk(m,row)=ham%hjk(row,m)
            
            if(GDflg.eq.'y') then
                if(row.eq.2)then
                    ham%diff_hjk(row,m,:)=h1etot_diff_bra+h2etot_diff_bra
                    if(m.eq.row)then
                        ham%diff_ovrlp(row,m,:) = 0
                    else
                        overlap_diff = diff_overlap_gpu(zstore(row)%sin,zstore(row)%cos,zstore(m)%sin,zstore(m)%cos,equal)
                        ham%diff_ovrlp(row,m,:) = overlap_diff
                    end if
                else if(m.eq.2)then
                    ham%diff_hjk(m,row,:)=h1etot_diff_ket+h2etot_diff_ket
                    overlap_diff = diff_overlap_gpu(zstore(row)%sin,zstore(row)%cos,zstore(m)%sin,zstore(m)%cos,equal)
                    ham%diff_ovrlp(m,row,:) = overlap_diff
                end if
            end if
            ! write(6,"(a,i0,a,i0,a)") "Hamiltonian row/comun ",row,'/', m, " completed"
        end do
        !!$omp end do
        !!$omp end parallel
     

        deallocate(z1jk,stat=ierr)
        if(ierr==0) deallocate(z2l,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        
        return

    end subroutine he_row_gpu

    subroutine one_elec_part_gpu(zs1sin,zs1cos,z2l,h1etot,occupancy,h1ei,h1etot_diff_bra,h1etot_diff_ket,zs2sin,zs2cos,equal)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),intent(inout)::h1etot
        integer,dimension(:,:,:,:),intent(in)::occupancy
        real(kind=8),dimension(norb),intent(inout)::h1etot_diff_bra,h1etot_diff_ket
        real(kind=8), dimension(:,:), intent(in)::h1ei
        integer,intent(in)::equal
        real(kind=8),dimension(norb,norb,norb)::h1etot_diff
        complex(kind=8),dimension(norb,norb)::temp
        integer::j,k

        if (errorflag .ne. 0) return
        h1etot=(0.0,0.0)
        h1etot_diff_bra=0.0
        h1etot_diff_ket=0.0
        h1etot_diff=0.0
        temp=cmplx(0.0,0.0)
        !$omp target map(to:h1ei,occupancy,z2l,zs1sin,zs1cos,zs2sin,zs2cos,equal) map(tofrom:temp,h1etot_diff)
        !$omp teams distribute parallel do simd collapse(2) &
        !$omp & private(j,k) shared(h1ei,occupancy,z2l,zs1sin,temp)
        do j=1, norb
            do k=1, norb
                if(h1ei(j,k).ne.(0.0)) then 
                    temp(j,k)=one_elec_body_gpu(zs1sin,zs1cos,z2l(j,:,:),occupancy(j,k,:,:),h1ei(j,k),k)
                end if
            end do
        end do
        !$omp end teams distribute parallel do simd
        if(equal.lt.4)then 
            !$omp teams distribute parallel do simd collapse(2) &
            !$omp & private(j,k) shared(h1ei,occupancy,z2l,zs1sin,zs1cos,zs2sin,zs2cos,equal,h1etot_diff)
            do j=1, norb
                do k=1, norb
                    if(h1ei(j,k).ne.(0.0)) then 
                        h1etot_diff(j,k,:)=one_elec_body_grad_gpu(zs1sin,zs1cos,zs2sin,zs2cos,z2l(j,:,:),&
                        occupancy(j,k,:,:),h1ei(j,k),k,j,equal)
                    end if
                end do
            end do
            !$omp end teams distribute parallel do simd
        end if
        !$omp end target
       
        do j=1, norb
            do k=1, norb
                h1etot=h1etot+temp(j,k)
                if(equal.lt.4)then 
                    if(equal.eq.1)then
                        h1etot_diff_bra= h1etot_diff_bra+h1etot_diff(j,k,:)
                    else if(equal.eq.2)then
                        h1etot_diff_bra= h1etot_diff_bra+h1etot_diff(j,k,:)
                    else if(equal.eq.3)then
                        h1etot_diff_ket= h1etot_diff_ket+h1etot_diff(j,k,:)
                    end if
                end if
            end do
        end do

        if(equal.eq.1)then
            h1etot_diff_ket=h1etot_diff_bra
        end if 
     
        return

    end subroutine one_elec_part_gpu

    complex(kind=8) function one_elec_body_gpu(zs1sin,zs1cos,z2l,occupancy,h1ei,k)

        implicit none
        !$omp declare target
        complex(kind=8),dimension(:)::zs1sin,zs1cos
        integer,dimension(:,:),intent(in)::occupancy
        complex(kind=8),dimension(:,:),intent(in)::z2l
        real(kind=8), intent(in)::h1ei
        integer,intent(in)::k
        complex(kind=8),dimension(2,norb)::zomt

        zomt=z2l
        zomt(1,k)=zomt(2,k)
        zomt(2,k)=cmplx(0.0,0.0)
        zomt=zomt*occupancy
        one_elec_body_gpu=product((conjg(zs1sin)*zomt(1,:))+(conjg(zs1cos)*zomt(2,:)))*h1ei
     
        return

    end function one_elec_body_gpu

    function one_elec_body_grad_gpu(zs1sin,zs1cos,zs2sin,zs2cos,z2l,occupancy,h1ei,k,j,equal)

        implicit none
        !$omp declare target
        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        real(kind=8),dimension(norb)::one_elec_body_grad_gpu
        integer,dimension(:,:),intent(in)::occupancy
        complex(kind=8),dimension(:,:),intent(in)::z2l
        real(kind=8), intent(in)::h1ei
        integer,intent(in)::k,j,equal
        complex(kind=8),dimension(2,norb)::zomt

        zomt=z2l
        zomt(1,k)=zomt(2,k)
        zomt(2,k)=cmplx(0.0,0.0)
        zomt=zomt*occupancy
        one_elec_body_grad_gpu = diff_overlap_cran_gpu(zs1sin,zs1cos,zs2sin,zs2cos,equal,zomt,j,k,occupancy)*h1ei
        
        return

    end function one_elec_body_grad_gpu
        
    subroutine two_elec_part_gpu(zs1sin,zs1cos,z1jk,z2l,h2etot,occupancy_2an,occupancy_an,h2ei,&
        h2etot_diff_bra,h2etot_diff_ket,zs2sin,zs2cos,equal)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:,:),intent(in)::z1jk
        complex(kind=8),intent(inout)::h2etot
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an
        integer,dimension(:,:,:),intent(in)::occupancy_an
        real(kind=8),dimension(norb),intent(inout)::h2etot_diff_bra,h2etot_diff_ket
        real(kind=8), dimension(:,:,:,:), intent(in)::h2ei
        integer,intent(in)::equal
        real(kind=8),dimension(norb,norb,norb,norb)::h2etot_diff
        integer::j,k,l
        complex(kind=8),dimension(norb,norb,norb)::tot

        if (errorflag .ne. 0) return
        
        h2etot=(0.0,0.0)
        h2etot_diff_bra=0.0
        h2etot_diff_ket=0.0
        h2etot_diff=0.0
        tot=cmplx(0.0,0.0)
       

        !$omp parallel shared(z1jk,z2l,zs1sin,zs1cos,zs2sin,zs2cos,tot,occupancy_2an,occupancy_an,h2ei,equal,h2etot_diff) private(j)
        !$omp do
        do j=1, norb
            tot(j,:,:) = two_elec_part_body_gpu(zs1sin,zs2sin,z2l,z1jk(j,:,:,:),h2ei(j,:,:,:),j)
            if(equal.lt.4)then
                h2etot_diff(j,:,:,:) = two_elec_part_grad_gpu(zs1sin,zs1cos,zs2sin,zs2cos,z2l,z1jk(j,:,:,:),&
                h2ei(j,:,:,:),occupancy_2an(j,:,:,:),occupancy_an(:,:,:),j,equal)
            end if
        end do
        !$omp end  do
        !$omp end parallel
      
        do j=1, norb
            do k=1,norb
                do l=1, norb
                    h2etot=h2etot+tot(j,k,l)
                    if(equal.lt.4)then 
                        if(equal.eq.1)then
                            h2etot_diff_bra= h2etot_diff_bra+h2etot_diff(j,k,l,:)
                        else if(equal.eq.2)then
                            h2etot_diff_bra= h2etot_diff_bra+h2etot_diff(j,k,l,:)
                        else if(equal.eq.3)then
                            h2etot_diff_ket= h2etot_diff_ket+h2etot_diff(j,k,l,:)
                        end if
                    end if
                end do 
            end do
        end do

        if(equal.eq.1)then
            h2etot_diff_ket=h2etot_diff_bra
        end if 

        h2etot=h2etot*0.5
        h2etot_diff_bra = h2etot_diff_bra*0.5
        h2etot_diff_ket = h2etot_diff_ket*0.5
     
        return

    end subroutine two_elec_part_gpu

    function two_elec_part_body_gpu(zs1sin,zs2sin,z2l,z1jk,h2ei,j)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs2sin
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:),intent(in)::z1jk
        real(kind=8), dimension(:,:,:), intent(in)::h2ei
        integer,intent(in)::j
        complex(kind=8),dimension(norb,norb)::tot,two_elec_part_body_gpu
        integer::k,l,jspin

        tot=cmplx(0.0,0.0)

        if(zs1sin(j)==(0.0,0.0))then
            two_elec_part_body_gpu=(0.0,0.0)
            return
        end if
        if(modulo(j,2)==0)then
            jspin=2
        else
            jspin=1
        end if
        do k=1, norb
            if(j.eq.k) cycle
            if(occ_iszero(z1jk(k,:,:)).eqv..true.)then
                CYCLE
            end if
            !$omp target teams distribute parallel do simd 
            do l=jspin, norb, 2
                if(zs2sin(l)==(0.0,0.0))then
                    CYCLE
                end if
                tot(k,l) = z_an_z3_gpu(z1jk(k,:,:),z2l(l,:,:),h2ei(k,l,:))
            end do
            !$omp end target teams distribute parallel do simd 
        end do

        two_elec_part_body_gpu=tot

        return

    end function two_elec_part_body_gpu

    function two_elec_part_grad_gpu(zs1sin,zs1cos,zs2sin,zs2cos,z2l,z1jk,h2ei,occupancy_2an,occupancy_an,j,equal)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:),intent(in)::z1jk
        real(kind=8), dimension(:,:,:), intent(in)::h2ei
        integer,dimension(:,:,:),intent(in)::occupancy_2an
        integer,dimension(:,:,:),intent(in)::occupancy_an
        integer,intent(in)::j,equal
        real(kind=8),dimension(norb,norb,norb)::two_elec_part_grad_gpu
        integer,dimension(2,norb)::occupancy
        integer::k,l,jspin

        two_elec_part_grad_gpu=0.0

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
                occupancy=occupancy_2an(k,:,:)*occupancy_an(l,:,:)
                two_elec_part_grad_gpu(k,l,:) = z_an_z3_diff_gpu(z1jk(k,:,:),z2l(l,:,:),h2ei(k,l,:),&
                equal,occupancy,j,k,l,zs1sin,zs1cos,zs2sin,zs2cos)
            end do
        end do

        return

    end function two_elec_part_grad_gpu

    ! Top level routine to allocate hamiltonian and overlap matrices 
    subroutine hamgen(ham,zstore,elecs,size,verb)

        implicit none

        type(hamiltonian), intent(inout)::ham 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,allocatable,dimension(:,:,:)::occupancy_an
        integer,allocatable,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an
        complex(kind=8),allocatable, dimension(:,:,:)::passback
        integer, allocatable,dimension(:)::IPIV1
        integer,intent(in)::size,verb
        complex(kind=8),allocatable,dimension(:)::WORK1
        real(kind=8),dimension(ndet,ndet,norb)::temp2
        integer:: j,k,l,ierr

        ! real(kind=8),dimension(norb,norb)::ph1ei
        ! real(kind=8),dimension(norb,norb,norb,norb)::ph2ei
        ! real(kind=8) :: phnuc

        ! complex(kind=8), dimension(ndet,norb)::psin
        ! complex(kind=8), dimension(ndet,norb)::pcos
        ! real(kind=8),dimension(ndet,norb)::pphi

        ! complex(kind=8), dimension(ndet,ndet)::phjk
        ! complex(kind=8), dimension(ndet,ndet)::povrlp
        ! complex(kind=8), dimension(ndet,ndet)::pinv
        ! complex(kind=8), dimension(ndet,ndet)::pkinvh

        ! real(kind=8), dimension(ndet,ndet,norb)::pdiff_hjk
        ! real(kind=8), dimension(ndet,ndet,norb)::pdiff_ovrlp
        ! real(kind=8), dimension(ndet,ndet,ndet,norb)::pdiff_invh

        if (errorflag .ne. 0) return
        
        allocate(occupancy_an(norb,2,norb),stat=ierr)
        if(ierr==0) allocate(occupancy_2an(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(occupancy_an_cr(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(passback(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in occupancy vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        occupancy_an=1
        occupancy_2an=1
        occupancy_an_cr=1

        ! ph1ei=elecs%h1ei
        ! ph2ei=elecs%h2ei
        ! phnuc=elecs%hnuc 

        ! phjk=ham%hjk
        ! povrlp=ham%ovrlp
        ! pinv=ham%inv
        ! pkinvh=ham%kinvh
        ! if(GDflg.eq.'y')then
        !     pdiff_hjk=ham%diff_hjk
        !     pdiff_ovrlp=ham%diff_ovrlp
        !     pdiff_invh=ham%diff_invh
        ! end if
        

        ! do j=1,ndet
        !     psin(j,:)=zstore(j)%sin(:)
        !     pcos(j,:)=zstore(j)%cos(:)
        !     pphi(j,:)=zstore(j)%phi(:)
        ! end do
       
       
        !call omp_set_nested(.TRUE.)
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
        !$omp do simd collapse(2)
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
        !$omp end do simd
        !$omp end parallel
    
       
        do j=1, size
            call he_row_gpu(ham,zstore,elecs,j,size,occupancy_2an,occupancy_an_cr,occupancy_an,passback)
            if(verb.eq.1)then
                write(6,"(a,i0,a)") "Hamiltonian row ",j, " completed"
            end if  
        end do

        
     
       if(ierr==0) deallocate(passback,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in passback vector deallocation . ierr had value ", ierr
            errorflag=1
        end if

        ham%inv=ham%ovrlp
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

        if (ierr==0) call ZGETRF(size,size,ham%inv,size,IPIV1,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if
        if (ierr==0) call ZGETRI(size,ham%inv,size,IPIV1,WORK1,size,ierr)
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
        ham%kinvh=matmul(ham%inv,ham%hjk)
        !$omp end workshare
        !$omp end parallel
       
        if(GDflg.eq.'y')then
            do k=1, ndet
                do l=1, ndet
                    if(l.eq.2)then
                        temp2(k,2,:)=matmul(REAL(ham%inv(k,:)),ham%diff_ovrlp(2,:,:))
                    else
                        temp2(k,l,:)=real(ham%inv(k,l))*ham%diff_ovrlp(2,l,:)
                    end if
                end do
            end do
            do k=1, ndet
                do l=1, ndet
                    ham%diff_invh(2,k,l,:)=matmul(transpose(temp2(k,:,:)),real(ham%kinvh(:,l)))*(-1)
                end do
            end do
        end if
      
        return
        
    end subroutine hamgen

    function diff_overlap_cran_gpu(zs1sin,zs1cos,zs2sin,zs2cos,dtype,zomt,annihilate2,create2,occupancy)

        implicit none
        !$omp declare target
        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        real(kind=8),dimension(norb)::diff_overlap_cran_gpu
        complex(kind=8),dimension(:,:)::zomt
        integer,intent(in)::dtype,annihilate2,create2
        integer,dimension(:,:),intent(in)::occupancy
        real(kind=8),dimension(norb)::bra_prod,ket_prod,prod
        integer::j
        real(kind=8),dimension(norb)::temp,temp2

        if (errorflag .ne. 0) return

        prod=real(((conjg(zs1sin)*zomt(1,:)))+((conjg(zs1cos)*zomt(2,:))))
        diff_overlap_cran_gpu=0
        ! Differntial of overlap of the same ZS 0 unless an operator has acted on ZS
        if(dtype.eq.1)then
            if(annihilate2.eq.create2)then
                if((real(zs2cos(annihilate2)).eq.0).and.(real(zs2sin(annihilate2)).eq.1).or.(real(zs2sin(annihilate2)).eq.0))then
                    diff_overlap_cran_gpu(annihilate2)=0
                else
                    bra_prod=prod           !dead amplitude is zero
                    bra_prod(annihilate2)=real(2*zs1sin(annihilate2)*zs1cos(annihilate2)*occupancy(1,annihilate2))
                    diff_overlap_cran_gpu(annihilate2)=product(bra_prod)   
                end if 
            else if(annihilate2.ne.create2)then
                bra_prod=prod
                if((real(zs2cos(annihilate2)).eq.0).and.(real(zs2sin(annihilate2)).eq.1)) then
                    bra_prod(annihilate2)=-1*occupancy(2,annihilate2)
                else if((real(zs2sin(annihilate2)).eq.0).and.(real(zs2cos(annihilate2)).eq.1)) then
                    bra_prod(annihilate2)=occupancy(2,annihilate2)
                else
                    bra_prod(annihilate2)=real(((zs1cos(annihilate2)**2)-(zs1sin(annihilate2)**2))*occupancy(2,annihilate2))
                end if
                diff_overlap_cran_gpu(annihilate2)=product(bra_prod)   
                bra_prod=prod
                if((real(zs2cos(create2)).eq.0).and.(real(zs2sin(create2)).eq.1)) then
                    bra_prod(create2)=-1*occupancy(1,create2)
                else if((real(zs2sin(create2)).eq.0).and.(real(zs2cos(create2)).eq.1)) then
                    bra_prod(create2)=occupancy(1,create2)
                else
                    bra_prod(create2)=real(((zs1cos(create2)**2)-(zs1sin(create2)**2))*occupancy(1,create2))
                end if               !dead amplitude is zero
                diff_overlap_cran_gpu(create2)=product(bra_prod)  
            end if
        else if(dtype.eq.2)then 
            bra_prod=real(zs1cos*zs2sin*occupancy(1,:)-zs1sin*zs2cos*occupancy(2,:))            
            if(annihilate2.eq.create2)then  !dead amplitude is zero
                bra_prod(annihilate2)=real(zs1cos(annihilate2)*zs2sin(annihilate2)*occupancy(1,annihilate2))
            else                   
                bra_prod(annihilate2)=-real(zs1sin(annihilate2)*zs2sin(annihilate2))*occupancy(2,annihilate2)
                bra_prod(create2)=real(zs1cos(create2)*zs2cos(create2))*occupancy(1,create2) 
            end if
            !$omp parallel do simd shared(prod,bra_prod) private(temp,j)
            do j=1,norb
                temp=prod
                temp(j)=bra_prod(j)
                diff_overlap_cran_gpu(j)=product(temp)
            end do 
            !$omp end parallel do simd 
        else if(dtype.eq.3)then
            ket_prod=real(zs1sin*zs2cos*occupancy(1,:)-zs1cos*zs2sin*occupancy(2,:))                  
            if(annihilate2.eq.create2)then  !dead amplitude is zero
                ket_prod(annihilate2)=real(zs1sin(annihilate2)*zs2cos(annihilate2)*occupancy(1,annihilate2))
            else                   
                ket_prod(annihilate2)=real(zs1cos(annihilate2)*zs2cos(annihilate2))*occupancy(2,annihilate2)!alive amplitude is zero
                ket_prod(create2)=-real(zs1sin(create2)*zs2sin(create2))*occupancy(1,create2) !dead amplitude is zero
            end if
            !$omp parallel do simd  shared(prod,ket_prod) private(temp2,j)
            do j=1,norb
                temp2=prod
                temp2(j)=ket_prod(j)
                diff_overlap_cran_gpu(j)=product(temp2)
            end do
            !$omp end parallel do simd   
        end if
        
        RETURN


    end function diff_overlap_cran_gpu

     ! computes the vector of values formed by the derivative of the overlap with respect to each orbital. 
    ! Does not have capability to deal with states where creation and annihilation operators have acted
    function diff_overlap_gpu(zs1sin,zs1cos,zs2sin,zs2cos,dtype)
    
        implicit none
        !$omp declare target
        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        integer,intent(in)::dtype
        real(kind=8),dimension(norb)::diff_overlap_gpu
        real(kind=8),dimension(norb)::bra_prod,ket_prod,prod
        integer::j
        real(kind=8),dimension(norb)::temp,temp2

        if (errorflag .ne. 0) return
        prod=real((conjg(zs1sin)*zs2sin)+(conjg(zs1cos)*zs2cos))
        if(dtype.eq.2)then
            bra_prod=real(conjg(zs1cos)*zs2sin)-real(conjg(zs1sin)*zs2cos)
            !$omp parallel do simd  shared(prod,bra_prod) private(temp,j)
            do j=1,norb
                temp=prod
                temp(j)=bra_prod(j)
                diff_overlap_gpu(j)=product(temp)
            end do
            !$omp end parallel do simd  
        else if(dtype.eq.3)then
            ket_prod=real(conjg(zs1sin)*zs2cos)-real(conjg(zs1cos)*zs2sin)
            !$omp parallel do simd  shared(prod,ket_prod) private(temp2,j)
            do j=1,norb
                temp2=prod
                temp2(j)=ket_prod(j)
                diff_overlap_gpu(j)=product(temp2)
            end do
            !$omp end parallel do simd  
        end if

        return
        
    end function diff_overlap_gpu

    complex(kind=8) function z_an_z3_gpu(z1,z2,vec)

        implicit none
        !$omp declare target
        complex(kind=8),dimension(:,:),intent(in)::z1,z2
        complex(kind=8),dimension(:,:),allocatable::vmult
        real(kind=8),dimension(norb),intent(in)::vec
        complex(kind=8),dimension(norb)::gg,hh
        complex(kind=8)::tot
        integer::j,gmax,hmin,ierr

        if (errorflag .ne. 0) return
        ierr=0

      
        allocate(vmult(2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in vmult allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        
        vmult=conjg(z1)*(z2)
        
    
        gg(1:norb)=(0.0,0.0)
        hh(1:norb)=(0.0,0.0)
        gmax=norb
        gg(1)=vmult(2,1)-vmult(1,1)

        do j=2, norb
            gg(j)=gg(j-1)*(vmult(2,j)-vmult(1,j))
            if(gg(j)==(0.0,0.0))then
                gmax=j
                EXIT 
            end if
        end do
        
        hmin=0
        hh(norb) = vmult(2,norb)+vmult(1,norb)
        do j=(norb-1),1,-(1)
            hh(j)=hh(j+1)*(vmult(2,j)+vmult(1,j))
            if(hh(j)==(0.0,0.0))then
                hmin=j
                EXIT 
            end if
        end do


        tot=(0.0,0.0)
        if (gmax < hmin) then
            z_an_z3_gpu=tot
            return
        end if

        if(vec(1).ne.0) then
            tot = tot+(conjg(z1(2,1))*z2(1,1)*hh(2)*vec(1))
        end if

      
        do j=2,norb-1
            if(vec(j).ne.0.0) then
                tot = tot+ (gg(j-1)*conjg(z1(2,j))*z2(1,j)*hh(j+1)*vec(j))
            end if
        end do

        if(vec(norb).ne.0) then
            tot = tot +(gg(norb-1)*conjg(z1(2,norb))*z2(1,norb)*vec(norb))
        end if

        deallocate(vmult)
    
        z_an_z3_gpu=tot
    
        return

    end function z_an_z3_gpu

    function z_an_z3_diff_gpu(z1,z2,vec,dtype,occupancy,annihilate1,annihilate1_2,annihilate2,zs1sin,zs1cos,zs2sin,zs2cos)

        implicit none
        complex(kind=8),dimension(:,:),intent(in)::z1,z2
        real(kind=8),dimension(:,:),allocatable::vmult_dd,vmult_1d, vmult_2d,vmult
        real(kind=8),dimension(norb),intent(in)::vec
        integer,dimension(:,:),intent(in)::occupancy
        integer, intent(in)::dtype,annihilate1,annihilate1_2,annihilate2
        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        real(kind=8),dimension(norb)::z_an_z3_diff_gpu
        real(kind=8),dimension(norb)::temp,temp2,temp0
        real(kind=8),dimension(norb)::gg_1,hh_1,gg_2,hh_2,gg_0,hh_0
        real(kind=8)::tot1, tot2,tot0
        integer::j,k,gmax1,hmin1,gmax2,hmin2,ierr,breakflag,gmax0,hmin0

        if (errorflag .ne. 0) return
        ierr=0

        allocate(vmult(2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in vmult allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        vmult=real(conjg(z1)*z2)
        
        z_an_z3_diff_gpu=0
        if(dtype.eq.1) then !Differentiation when zombie states are the same
            allocate(vmult_dd(2,norb),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult_dd allocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
            temp0=0
            !$omp target teams distribute parallel do &
            !$omp & private(vmult_dd,gmax0,hmin0,breakflag,tot0,j,k) &
            !$omp shared(annihilate1,annihilate1_2,annihilate2,temp0,zs1sin,zs1cos,zs2sin,zs2cos,vmult,occupancy,vec,z1,z2)
            do j=1, norb    !Differentiating w.r.t to orbital j
                breakflag=0
                vmult_dd=vmult
                if(annihilate1.eq.annihilate1_2)then
                    temp0(j)=0
                    cycle
                else if((annihilate1.eq.j).or.(annihilate1_2.eq.j))then
                    if(annihilate2.eq.j)then    !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(1j)
                        vmult_dd(1,j)=0
                        vmult_dd(2,j)=2*real(zs1sin(j)*zs1cos(j))*occupancy(2,j) !(sin2x=2sinxcosx)
                    else
                        vmult_dd(1,j)=0        !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(0j)
                        vmult_dd(2,j)=(1.0-2*((real(zs1sin(j)))**2))*occupancy(2,j)    !(cos2x = cos^x -sin^2x)
                    end if
                else if((annihilate1.ne.j).or.(annihilate1_2.ne.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(0j)*a^(a)_(1j)
                        vmult_dd(1,j)=0
                        vmult_dd(2,j)=(1.0-2*((real(zs1sin(j)))**2))*occupancy(2,j)
                    else
                        breakflag=j !before diff alive:a^(a)_(1j)*a^(a)_(1j) dead:a^(a)_(0j)*a^(a)_(0j)
                    end if          !Unless an opeator acts at position j this evaluates to 0
                end if
            
                gg_0(1:norb)=(0.0,0.0)
                hh_0(1:norb)=(0.0,0.0)
                gmax0=norb
                gg_0(1)=vmult_dd(2,1)-vmult_dd(1,1)

                do k=2, norb
                    gg_0(k)=gg_0(k-1)*(vmult_dd(2,k)-vmult_dd(1,k))
                    if(gg_0(k)==(0.0))then
                        gmax0=k
                        EXIT 
                    end if
                end do
                
                hmin0=0
                hh_0(norb) = vmult_dd(2,norb)+vmult_dd(1,norb)
                do k=(norb-1),1,-(1)
                    hh_0(k)=hh_0(k+1)*(vmult_dd(2,k)+vmult_dd(1,k))
                    if(hh_0(k)==(0.0))then
                        hmin0=k
                        EXIT 
                    end if
                end do

                tot0=(0.0)
                if (gmax0 < hmin0) then
                    temp0(j)=tot0
                    cycle
                end if

                if((breakflag.eq.0).or.(breakflag.eq.1))then
                    if(vec(1).ne.0) then
                        if(breakflag.eq.1)then
                            temp0(j)=(1.0-2*((real(zs1sin(1)))**2))*occupancy(2,1)*hh_0(2)*vec(1)
                            cycle
                        end if
                        if(j.eq.1)then
                            if((annihilate1.eq.1).or.(annihilate1_2.eq.1))then
                                if(annihilate2.ne.1)then
                                    tot0 = tot0+(2*real(zs1sin(1)*zs1cos(1))*occupancy(2,1)*hh_0(2)*vec(1))
                                end if
                            else
                                if(annihilate2.ne.1)then
                                    tot0 = tot0+((1.0-2*((real(zs1sin(1)))**2))*occupancy(2,1)*hh_0(2)*vec(1))
                                end if 
                            end if
                        else
                            tot0 = tot0+(REAL(conjg(z1(2,1))*z2(1,1))*hh_0(2)*vec(1))
                        end if
                    end if
                else if((breakflag.lt.norb))then
                    if(breakflag.ne.0)then
                        temp0(j)=(gg_0(j-1)*((1.0-2*((real(zs1sin(j)))**2))*occupancy(2,j))*hh_0(j+1)*vec(j))
                        cycle
                    end if
                    do k=2,norb-1
                        if(vec(k).ne.0.0) then
                            if(k.eq.j)then
                                if((annihilate1.eq.k).or.(annihilate1_2.eq.k))then
                                    if(annihilate2.ne.k)then
                                        tot0 = tot0+(gg_1(k-1)*(2*real(zs1sin(k)*zs1cos(k))*occupancy(2,k))*hh_0(k+1)*vec(k))
                                    end if
                                else
                                    if(annihilate2.ne.k)then
                                        tot0 = tot0+(gg_1(k-1)*((1.0-2*((real(zs1sin(k)))**2))*occupancy(2,k))*hh_0(k+1)*vec(k))
                                    end if 
                                end if
                            else
                                tot0 = tot0+ (gg_0(k-1)*REAL(conjg(z1(2,k))*z2(1,k))*hh_0(k+1)*vec(k))
                            end if
                        end if
                    end do
                else
                    if(breakflag.eq.norb)then
                        temp0(j)=(gg_0(norb-1)*((1.0-2*((real(zs1sin(norb)))**2))*occupancy(2,norb))*vec(norb))
                        cycle
                    end if
                    if(vec(norb).ne.0) then
                        if(norb.eq.j)then
                            if((annihilate1.eq.norb).or.(annihilate1_2.eq.norb))then
                                if(annihilate2.ne.norb)then
                                    tot0 = tot0 +(gg_0(norb-1)*(2*real(zs1sin(norb)*zs1cos(norb))*occupancy(2,norb))*vec(norb))
                                end if
                            else
                                if(annihilate2.ne.norb)then
                                    tot0 = tot0+ (gg_0(norb-1)*((1.0-2*((real(zs1sin(norb)))**2))*occupancy(2,norb))*vec(norb))
                                end if    
                            end if
                        else
                            tot0 = tot0 +(gg_0(norb-1)*REAL(conjg(z1(2,norb))*z2(1,norb))*vec(norb))
                        end if
                    end if
                end if

                temp0(j)=tot0
            end do
            !$omp end target teams distribute parallel do

            z_an_z3_diff_gpu=temp0(:)
            z_an_z3_diff_gpu=temp0(:)
            deallocate(vmult,stat=ierr)
            if(ierr==0) deallocate(vmult_dd,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult deallocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
        
        else if(dtype.eq.2)then !Differentiation when zombie states are not the same
            allocate(vmult_1d(2,norb),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult_1d/vmult_2d allocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
            temp=0
            !$omp target teams distribute parallel do &
            !$omp & private(vmult_1d,gmax1,hmin1,tot1,j,k) shared(annihilate1,&
            !$omp annihilate1_2,annihilate2,temp,zs1sin,zs1cos,zs2sin,zs2cos,vmult,occupancy,vec,z1,z2,norb)
            do j=1, norb
                vmult_1d=vmult
                if(annihilate1.eq.annihilate1_2)then
                    temp(j)=0
                    cycle
                else if((annihilate1.eq.j).or.(annihilate1_2.eq.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                        vmult_1d(1,j)=0
                        vmult_1d(2,j)=real(zs1cos(j)*zs2sin(j))*occupancy(2,j)
                    else !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(0j)
                        vmult_1d(1,j)=0
                        vmult_1d(2,j)=real(zs1cos(j)*zs2cos(j))*occupancy(2,j)
                    end if
                else if((annihilate1.ne.j).or.(annihilate1_2.ne.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                        vmult_1d(1,j)=0
                        vmult_1d(2,j)=-real(zs1sin(j)*zs2sin(j))*occupancy(2,j)
                    else    !before diff alive:a^(a)_(1j)*a^(b)_(1j) dead:a^(a)_(0j)*a^(b)_(0j)
                        vmult_1d(1,j)=real(zs1cos(j)*zs2sin(j))*occupancy(1,j)
                        vmult_1d(2,j)=-real(zs1sin(j)*zs2cos(j))*occupancy(2,j)
                    end if
                end if

                gg_1(1:norb)=(0.0,0.0)
                hh_1(1:norb)=(0.0,0.0)
                gmax1=norb
                gg_1(1)=(vmult_1d(2,1))-(vmult_1d(1,1))

                do k=2, norb
                    gg_1(k)=gg_1(k-1)*((vmult_1d(2,k)-vmult_1d(1,k)))
                    if(gg_1(k)==(0.0,0.0))then
                        gmax1=k
                        EXIT 
                    end if
                end do
                
                hmin1=0
               
                hh_1(norb) = (vmult_1d(2,norb)+vmult_1d(1,norb))
              
                do k=(norb-1),1,-(1)
                    hh_1(k)=hh_1(k+1)*(vmult_1d(2,k)+vmult_1d(1,k))
                    if(hh_1(k)==(0.0,0.0))then
                        hmin1=k
                        EXIT 
                    end if
                end do

                tot1=(0.0)
               
                if((gmax1 < hmin1))then
                    temp(j)=tot1
                    cycle
                end if

                if(vec(1).ne.0) then
                    if(j.eq.1)then
                        if((annihilate1.eq.1).or.(annihilate1_2.eq.1))then
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot1 = tot1+((REAL(zs1cos(1)*zs2sin(1))*occupancy(2,1))*hh_1(2)*vec(1))
                            end if
                        else
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot1 = tot1+((REAL(-zs1sin(1)*zs2sin(1))*occupancy(2,1))*hh_1(2)*vec(1))
                            end if 
                        end if
                    else
                        tot1 = tot1+(REAL(conjg(z1(2,1))*z2(1,1))*hh_1(2)*vec(1))
                    end if
                end if
                do k=2,norb-1
                    if(vec(k).ne.0.0) then
                        if(k.eq.j)then
                            if((annihilate1.eq.k).or.(annihilate1_2.eq.k))then
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                    tot1 = tot1+(gg_1(k-1)*(REAL(zs1cos(k)*zs2sin(k))*occupancy(2,k))*hh_1(k+1)*vec(k))
                                end if
                            else
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                    tot1 = tot1+(gg_1(k-1)*(REAL(-zs1sin(k)*zs2sin(k))*occupancy(2,k))*hh_1(k+1)*vec(k))
                                end if 
                            end if
                        else
                            tot1 = tot1+ (gg_1(k-1)*REAL(conjg(z1(2,k))*z2(1,k))*hh_1(k+1)*vec(k))
                        end if
                    end if
                end do

                if(vec(norb).ne.0) then
                    if(norb.eq.j)then
                        if((annihilate1.eq.norb).or.(annihilate1_2.eq.norb))then
                            if(annihilate2.ne.norb)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot1 = tot1+(gg_1(norb-1)*(REAL(zs1cos(norb)*zs2sin(norb))*occupancy(2,norb))*vec(norb))
                            end if
                        else
                            if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot1 = tot1+(gg_1(norb-1)*(REAL(-zs1sin(norb)*zs2sin(norb))*occupancy(2,norb))*vec(norb))
                            end if 
                        end if
                    else
                        tot1 = tot1 +(gg_1(norb-1)*REAL(conjg(z1(2,norb))*z2(1,norb))*vec(norb))
                    end if
                end if

                temp(j)=tot1
            end do
            !$omp end target teams distribute parallel do
            z_an_z3_diff_gpu=temp(:)
            deallocate(vmult,stat=ierr)
            if(ierr==0) deallocate(vmult_1d,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult deallocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
        else if(dtype.eq.3)then
            allocate(vmult_2d(2,norb),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult_1d/vmult_2d allocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
            temp2=0
            !$omp target teams distribute parallel do &
            !$omp & private(vmult_2d,gmax2,hmin2,tot2,j,k) shared(annihilate1,&
            !$omp annihilate1_2,annihilate2,temp2,zs1sin,zs1cos,zs2sin,zs2cos,vmult,occupancy,vec,z1,z2)
            do j=1, norb
                vmult_2d=vmult
                if(annihilate1.eq.annihilate1_2)then
                    temp(j)=0
                    cycle
                else if((annihilate1.eq.j).or.(annihilate1_2.eq.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                        vmult_2d(1,j)=0
                        vmult_2d(2,j)=real(zs1sin(j)*zs2cos(j))*occupancy(2,j)
                    else !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(0j)
                        vmult_2d(1,j)=0
                        vmult_2d(2,j)=-real(zs1sin(j)*zs2sin(j))*occupancy(2,j)
                    end if
                else if((annihilate1.ne.j).or.(annihilate1_2.ne.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                        vmult_2d(1,j)=0
                        vmult_2d(2,j)=real(zs1cos(j)*zs2cos(j))*occupancy(2,j)
                    else    !before diff alive:a^(a)_(1j)*a^(b)_(1j) dead:a^(a)_(0j)*a^(b)_(0j)
                        vmult_2d(1,j)=real(zs1sin(j)*zs2cos(j))*occupancy(1,j)
                        vmult_2d(2,j)=-real(zs1cos(j)*zs2sin(j))*occupancy(2,j)
                    end if
                end if

                gg_2(1:norb)=(0.0,0.0)
                hh_2(1:norb)=(0.0,0.0)
                gmax2=norb
                gg_2(1)=(vmult_2d(2,1))-(vmult_2d(1,1))

                do k=2, norb
                    gg_2(k)=gg_2(k-1)*((vmult_2d(2,k)-vmult_2d(1,k)))
                    if(gg_2(k)==(0.0,0.0))then
                        gmax2=k
                        EXIT 
                    end if
                end do
                

                hmin2=0
                hh_2(norb) = (vmult_2d(2,norb)+vmult_2d(1,norb))
        
                do k=(norb-1),1,-(1)
                    hh_2(k)=hh_2(k+1)*(vmult_2d(2,k)+vmult_2d(1,k))
                    if(hh_2(k)==(0.0,0.0))then
                        hmin2=k
                        EXIT 
                    end if
                end do

   
                tot2=(0.0)
                if((gmax2 < hmin2))then
                    temp2(j)=tot2
                    cycle
                end if

                if(vec(1).ne.0) then
                    if(j.eq.1)then
                        if((annihilate1.eq.1).or.(annihilate1_2.eq.1))then
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot2 = tot2 +((REAL(zs2sin(1)*zs1cos(1))*occupancy(2,1))*hh_2(2)*vec(1))
                            end if
                        else
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot2 = tot2 +((REAL(zs2cos(1)*zs1cos(1))*occupancy(2,1))*hh_2(2)*vec(1))
                            end if 
                        end if
                    else
                        tot2 = tot2+(REAL(conjg(z1(2,1))*z2(1,1))*hh_2(2)*vec(1))
                    end if
                end if
                do k=2,norb-1
                    if(vec(k).ne.0.0) then
                        if(k.eq.j)then
                            if((annihilate1.eq.k).or.(annihilate1_2.eq.k))then
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                    tot2 = tot2 +(gg_2(k-1)*(REAL(zs2sin(k)*zs1cos(k))*occupancy(2,k))*hh_2(k+1)*vec(k))
                                end if
                            else
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                    tot2 = tot2 +(gg_2(k-1)*(REAL(zs2cos(k)*zs1cos(k))*occupancy(2,k))*hh_2(k+1)*vec(k))
                                end if 
                            end if
                        else
                            tot2 = tot2+ (gg_2(k-1)*REAL(conjg(z1(2,k))*z2(1,k))*hh_2(k+1)*vec(k))
                        end if
                    end if
                end do

                if(vec(norb).ne.0) then
                    if(norb.eq.j)then
                        if((annihilate1.eq.norb).or.(annihilate1_2.eq.norb))then
                            if(annihilate2.ne.norb)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot2 = tot2 +(gg_2(norb-1)*(REAL(zs2sin(norb)*zs1cos(norb))*occupancy(2,norb))*vec(norb))
                            end if
                        else
                            if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot2 = tot2 +(gg_2(norb-1)*(REAL(zs2cos(norb)*zs1cos(norb))*occupancy(2,norb))*vec(norb))
                            end if 
                        end if
                    else
                        tot2 = tot2 +(gg_2(norb-1)*REAL(conjg(z1(2,norb))*z2(1,norb))*vec(norb))
                    end if
                end if

                temp2(j)=tot2
             
            end do
            !$omp end target teams distribute parallel do
     
            z_an_z3_diff_gpu=temp2(:)
            deallocate(vmult,stat=ierr)
            if(ierr==0) deallocate(vmult_2d,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult deallocation . ierr had value ", ierr
                errorflag=1
                return
            end if 

        end if

        return

    end function z_an_z3_diff_gpu

END MODULE ham