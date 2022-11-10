MODULE ham 

    use globvars
    use alarrays
    use operators
    use omp_lib
    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators
    subroutine he_row(ham,zstore,elecs,row,size,occupancy_2an,occupancy_an_cr,occupancy_an,passback)
        
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
        real(kind=8),dimension(2,norb)::overlap_diff
        
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
            
            !$omp do
            do l=1, norb
                z2l(l,1,:)=zstore(1)%sin(:)
                z2l(l,2,:)=zstore(1)%cos(:)
                z2l(l,2,l)=z2l(l,1,l)
                z2l(l,1,l)=cmplx(0.0,0.0)
            end do
            !$omp end do
           
        else
            z2l=passback
          
        end if
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
                !$omp do
                do l=1, norb
                    z2l(l,1,:)=zstore(m)%sin(:)
                    z2l(l,2,:)=zstore(m)%cos(:)
                    z2l(l,2,l)=z2l(l,1,l)
                    z2l(l,1,l)=cmplx(0.0,0.0)
                end do
                !$omp end do
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
        
            if(m.eq.row)then
                call one_elec_part(zstore(row),z2l,h1etot,occupancy_an_cr,&
                                    elecs%h1ei,h1etot_diff_bra,h1etot_diff_ket,zstore(m),equal)
            else 
                call one_elec_part(zstore(row),z2l,h1etot,occupancy_an_cr,&
                                    elecs%h1ei,h1etot_diff_bra,h1etot_diff_ket,zstore(m),equal)
            end if
            
            z2l=z2l*occupancy_an
            !$omp flush(z2l)
            if(m.eq.row)then
                call two_elec_part(zstore(row),z1jk,z2l,h2etot,occupancy_2an,occupancy_an,&
                                            elecs%h2ei,h2etot_diff_bra,h2etot_diff_ket,zstore(m),equal)
            else 
                call two_elec_part(zstore(row),z1jk,z2l,h2etot,occupancy_2an,occupancy_an,&
                                            elecs%h2ei,h2etot_diff_bra,h2etot_diff_ket,zstore(m),equal)
            end if
          
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
                        overlap_diff = diff_overlap(zstore(row),zstore(m),equal)
                        ham%diff_ovrlp(row,m,:) = overlap_diff(1,:)
                    end if
                else if(m.eq.2)then
                    ham%diff_hjk(m,row,:)=h1etot_diff_ket+h2etot_diff_ket
                    overlap_diff = diff_overlap(zstore(row),zstore(m),equal)
                    ham%diff_ovrlp(m,row,:) = overlap_diff(2,:)
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

    end subroutine he_row

    subroutine one_elec_part(zs1,z2l,h1etot,occupancy,h1ei,h1etot_diff_bra,h1etot_diff_ket,zs2,equal)

        implicit none

        type(zombiest),intent(in)::zs1,zs2
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),intent(inout)::h1etot
        integer,dimension(:,:,:,:),intent(in)::occupancy
        real(kind=8),dimension(norb),intent(inout)::h1etot_diff_bra,h1etot_diff_ket
        real(kind=8), dimension(:,:), intent(in)::h1ei
        integer,intent(in)::equal
        complex(kind=8),dimension(2,norb)::zomt
        real(kind=8),dimension(2,norb)::h1etot_diff
        complex(kind=8)::temp
        integer::j,k

        if (errorflag .ne. 0) return
        h1etot=(0.0,0.0)
        h1etot_diff_bra=0.0
        h1etot_diff_ket=0.0
        h1etot_diff=0.0
        temp=cmplx(0.0,0.0)
        !$omp parallel private(j,k,zomt) shared(h1ei,occupancy,z2l,zs1,zs2,equal,temp,h1etot_diff)
        !$omp do schedule(dynamic) reduction(+:temp,h1etot_diff)
        do j=1, norb
            do k=1, norb
                zomt(:,:)=z2l(j,:,:)
                zomt(1,k)=zomt(2,k)
                zomt(2,k)=cmplx(0.0,0.0)
                zomt=zomt*occupancy(j,k,:,:)
                temp=temp+product((conjg(zs1%sin)*zomt(1,:))+(conjg(zs1%cos)*zomt(2,:)))*h1ei(j,k)
                if(GDflg.eq.'y')then
                    if(equal.lt.4)then 
                        h1etot_diff = h1etot_diff+ diff_overlap_cran(zs1,zs2,equal,zomt,j,k,occupancy(j,k,:,:))*h1ei(j,k)
                    end if
                end if
            end do
        end do
        !$omp end do
        !$omp end parallel
       
        h1etot=temp
        h1etot_diff_bra=  h1etot_diff(1,:)
        h1etot_diff_ket=  h1etot_diff(2,:)
    
        return

    end subroutine one_elec_part
        
    subroutine two_elec_part(zs1,z1jk,z2l,h2etot,occupancy_2an,occupancy_an,h2ei,h2etot_diff_bra,h2etot_diff_ket,zs2,equal)

        implicit none

        type(zombiest),intent(in)::zs1,zs2
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:,:),intent(in)::z1jk
        complex(kind=8),intent(inout)::h2etot
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an
        integer,dimension(:,:,:),intent(in)::occupancy_an
        real(kind=8),dimension(norb),intent(inout)::h2etot_diff_bra,h2etot_diff_ket
        real(kind=8), dimension(:,:,:,:), intent(in)::h2ei
        integer,intent(in)::equal
        integer,dimension(2,norb)::occupancy
        real(kind=8),dimension(2,norb)::h2etot_diff
        integer::j,k,l,jspin
        complex(kind=8)::tot

        if (errorflag .ne. 0) return
        
        h2etot=(0.0,0.0)
        h2etot_diff_bra=0.0
        h2etot_diff_ket=0.0
        h2etot_diff=0.0
        tot=cmplx(0.0,0.0)
        !$omp parallel shared(z1jk,z2l,zs1,zs2,tot,occupancy_2an,&
        !$omp           occupancy_an,h2ei,equal,h2etot_diff) private(j,k,l,jspin,occupancy)
        !$omp do schedule(dynamic) reduction(+:tot)
        do j=1, norb
            if(zs1%sin(j)==(0.0,0.0))then
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
                    if(zs2%sin(l)==(0.0,0.0))then
                        CYCLE
                    end if
                    tot = tot+z_an_z3(z1jk(j,k,:,:),z2l(l,:,:),h2ei(j,k,l,:))
                end do
            end do
        end do
        !$omp end do
      
       !$omp end parallel
        if(GDflg.eq.'y')then
            if(equal.lt.4)then 
            !!$omp do schedule(dynamic) reduction(+:h2etot_diff)
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
                            h2etot_diff = h2etot_diff + z_an_z3_diff(z1jk(j,k,:,:),z2l(l,:,:),h2ei(j,k,l,:),equal,&
                                occupancy,j,k,l,zs1,zs2)
                        end do
                    end do
                   
                end do
            end if
            !!$omp end do
        end if
        !!$omp end parallel
      
        h2etot=tot*0.5
        h2etot_diff_bra = h2etot_diff(1,:)*0.5
        h2etot_diff_ket = h2etot_diff(2,:)*0.5
   
        return

    end subroutine two_elec_part
    
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
       
       
        !$ call omp_set_nested(.TRUE.)
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
    
       
        do j=1, size
            call he_row(ham,zstore,elecs,j,size,occupancy_2an,occupancy_an_cr,occupancy_an,passback)
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

    
    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators
    subroutine he_row_gpu(ph1ei,ph2ei,phnuc,psin,pcos,pphi,occupancy_2an,&
        occupancy_an_cr,occupancy_an,passback,phjk,povrlp,pdiff_hjk,pdiff_ovrlp,row,size)
        
        implicit none
        
        real(kind=8),pointer,dimension(:,:),intent(in)::ph1ei
        real(kind=8), dimension(:,:,:,:),pointer,intent(in)::ph2ei
        real(kind=8),pointer,intent(in) :: phnuc
        complex(kind=8), dimension(:,:),intent(in)::psin
        complex(kind=8), dimension(:,:),intent(in)::pcos
        real(kind=8),dimension(:,:),intent(in)::pphi
        complex(kind=8), dimension(:,:), pointer,intent(inout)::phjk
        complex(kind=8), dimension(:,:), pointer,intent(inout)::povrlp
        real(kind=8), dimension(:,:,:), pointer,intent(inout)::pdiff_hjk
        real(kind=8), dimension(:,:,:), pointer,intent(inout)::pdiff_ovrlp
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
        real(kind=8),dimension(norb)::bra_prod,ket_prod,prod,temp
        
        !!$omp declare target 
        if (errorflag .ne. 0) return
        ierr = 0

        
        allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if
  
        
        !$omp target teams map(alloc:z1jk(norb,norb,2,norb),z2l(norb,2,norb),prod(norb),bra_prod(norb),temp(norb),ket_prod(norb),&
        !$omp h1etot_diff_bra(norb),h2etot_diff_bra(norb),h1etot_diff_ket(norb),h2etot_diff_ket(norb)) &
        !$omp map(to:h1etot,h2etot,equal,row,m)
        h1etot=cmplx(0.0,0.0)
        h2etot=cmplx(0.0,0.0)
        h1etot_diff_bra(:)=0.0
        h1etot_diff_ket(:)=0.0
        h2etot_diff_bra(:)=0.0
        h2etot_diff_ket(:)=0.0
    
    
        if(row.eq.1)then
            !$omp distribute parallel do simd
            do l=1, norb
                z2l(l,1,:)=psin(1,:) 
                z2l(l,2,:)=pcos(1,:)
                z2l(l,2,l)=z2l(l,1,l)
                z2l(l,1,l)=cmplx(0.0,0.0)
            end do
            !$omp end distribute parallel do simd
        else
            z2l=passback
        end if
        
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
  

        do m=row,size
        
            h1etot=cmplx(0.0,0.0)
            h2etot=cmplx(0.0,0.0)
            h1etot_diff_bra(:)=0.0
            h1etot_diff_ket(:)=0.0
            h2etot_diff_bra(:)=0.0
            h2etot_diff_ket(:)=0.0
    
            if(m.gt.row)then
                !$omp distribute parallel do simd
                do l=1, norb
                    z2l(l,1,:)=psin(m,:) 
                    z2l(l,2,:)=pcos(m,:)
                    z2l(l,2,l)=z2l(l,1,l)
                    z2l(l,1,l)=cmplx(0.0,0.0)
                end do
                !$omp end distribute parallel do simd
                if(m.eq.row+1)then 
                    passback=z2l
                end if
            end if
            
            
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
          
         
            call one_elec_part_gpu(psin(row,:),psin(m,:),pcos(row,:),pcos(m,:),pphi(row,:),&
                z2l,h1etot,occupancy_an_cr,ph1ei,h1etot_diff_bra,h1etot_diff_ket,equal) 

            z2l=z2l*occupancy_an
          
            call two_elec_part_gpu(psin(row,:),psin(m,:),pcos(row,:),pcos(m,:),z1jk,z2l,h2etot,occupancy_2an,&
            occupancy_an,ph2ei,h2etot_diff_bra,h2etot_diff_ket,equal)
            
            povrlp(row,m)=product(((conjg(psin(row,:))*psin(m,:)))+((conjg(pcos(row,:))*pcos(m,:))))
            povrlp(m,row)= povrlp(row,m)
            phjk(row,m)=h1etot+h2etot+(phnuc*povrlp(row,m))
            phjk(m,row)=phjk(row,m)
       
           

            if(GDflg.eq.'y') then
                if(row.eq.2)then
                    temp=h1etot_diff_bra(:)+h2etot_diff_bra(:)
                    pdiff_hjk(row,m,:)= temp !h1etot_diff_bra(:)+h2etot_diff_bra(:)
                    if(m.eq.row)then
                        pdiff_ovrlp(row,m,:) = 0
                    else
                        prod=real((conjg(psin(row,:))*psin(m,:))+(conjg(pcos(row,:))*pcos(m,:)))
                        bra_prod=real(conjg(pcos(row,:))*psin(m,:))-real(conjg(psin(row,:))*pcos(m,:))
                        !$omp parallel do
                        do j=1,norb
                            temp=prod
                            temp(j)=bra_prod(j)
                            pdiff_ovrlp(row,m,j)=product(temp)
                        end do
                 
                    end if
                else if(m.eq.2)then
                    temp=h1etot_diff_ket(:)+h2etot_diff_ket(:)
                    pdiff_hjk(m,row,:)=temp!h1etot_diff_ket(:)+h2etot_diff_ket(:)
                    prod=real((conjg(psin(row,:))*psin(m,:))+(conjg(pcos(row,:))*pcos(m,:)))
                    ket_prod=real(conjg(psin(row,:))*pcos(m,:))-real(conjg(pcos(row,:))*psin(m,:))
                    !$omp parallel do
                    do j=1,norb
                        temp=prod
                        temp(j)=ket_prod(j)
                        pdiff_ovrlp(m,row,j)=product(temp)
                    end do
                    
                end if
             
            end if
          
        end do
        !$omp end target teams
       
        deallocate(z1jk,stat=ierr)
        if(ierr==0) deallocate(z2l,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

       
        return

    end subroutine he_row_gpu

    subroutine one_elec_part_gpu(psin1,psin2,pcos1,pcos2,pphi1,z2l,h1etot,occupancy,&
        ph1ei,h1etot_diff_bra,h1etot_diff_ket,equal)
        !!$omp declare target 
        implicit none

        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),intent(inout)::h1etot
        integer,dimension(:,:,:,:),intent(in)::occupancy
        real(kind=8),dimension(:),intent(inout)::h1etot_diff_bra,h1etot_diff_ket
        real(kind=8),pointer,dimension(:,:),intent(in)::ph1ei
        integer,intent(in)::equal
        complex(kind=8), dimension(:),intent(in)::psin1,psin2
        complex(kind=8), dimension(:),intent(in)::pcos1,pcos2
        real(kind=8),dimension(:),intent(in)::pphi1
        complex(kind=8),dimension(2,norb)::zomt
        real(kind=8),dimension(norb)::temp1,temp2,bra_prod,ket_prod,prod
        integer::j,k,l

        
       
        !!$omp declare target
        if (errorflag .ne. 0) return
        h1etot=(0.0,0.0)
        h1etot_diff_bra=0.0
        h1etot_diff_ket=0.0

        
        
        !$omp target teams distribute parallel do simd collapse(2) reduction(+:h1etot,h1etot_diff_bra,h1etot_diff_ket)&
        !$omp & map(alloc:zomt(2,norb),bra_prod(norb),ket_prod(norb),prod(norb),temp1(norb),temp2(norb)) &
        !$omp & private(zomt,j,k,l,temp1,temp2,bra_prod,ket_prod,prod) &
        !$omp & shared(ph1ei,psin1,psin2,pcos1,pcos2,pphi1,occupancy,equal,z2l)
        do j=1, norb
            do k=1, norb
                if(ph1ei(j,k).eq.0.0)then
                    cycle
                end if
                zomt(:,:)=z2l(j,:,:)
                zomt(1,k)=zomt(2,k)
                zomt(2,k)=cmplx(0.0,0.0)
                zomt=zomt*occupancy(j,k,:,:)
                h1etot=h1etot+product((conjg(psin1)*zomt(1,:))+(conjg(pcos1)*zomt(2,:)))*ph1ei(j,k) 
                if(GDflg.eq.'y')then
                    if(equal.lt.4)then
                        prod=real(((conjg(psin1)*zomt(1,:)))+((conjg(pcos1)*zomt(2,:))))
                        ! Differntial of overlap of the same ZS 0 unless an operator has acted on ZS
                        if(equal.eq.1)then
                            if(j.eq.k)then
                                if((real(pcos2(j)).eq.0).and.(real(psin2(j)).eq.1).or.(real(psin2(j)).eq.0))then
                                    h1etot_diff_bra(j)=h1etot_diff_bra(j)+0
                                    h1etot_diff_ket(j)=h1etot_diff_ket(j)+0
                                else
                                    bra_prod=prod           !dead amplitude is zero
                                    bra_prod(j)=sin(2*pphi1(j))*occupancy(j,k,1,j)
                                    h1etot_diff_bra(j)= h1etot_diff_bra(j)+product(bra_prod)*ph1ei(j,k)
                                    h1etot_diff_ket(j)= h1etot_diff_ket(j)+product(bra_prod)*ph1ei(j,k)      
                                end if 
                            else if(j.ne.k)then
                                bra_prod=prod
                                if((real(pcos2(j)).eq.0).and.(real(psin2(j)).eq.1)) then
                                    bra_prod(j)=-1*occupancy(j,k,2,j)
                                else if((real(psin2(j)).eq.0).and.(real(pcos2(j)).eq.1)) then
                                    bra_prod(j)=occupancy(j,k,2,j)
                                else
                                    bra_prod(j)=cos(2*pphi1(j))*occupancy(j,k,2,j)
                                end if
                                h1etot_diff_bra(j)= h1etot_diff_bra(j) + product(bra_prod)*ph1ei(j,k) 
                                h1etot_diff_ket(j)= h1etot_diff_ket(j) + product(bra_prod)*ph1ei(j,k)        
                                bra_prod=prod
                                if((real(pcos2(k)).eq.0).and.(real(psin2(k)).eq.1)) then
                                    bra_prod(k)=-1*occupancy(j,k,1,k)
                                else if((real(psin2(k)).eq.0).and.(real(pcos2(k)).eq.1)) then
                                    bra_prod(k)=occupancy(j,k,1,k)
                                else
                                    bra_prod(k)=cos(2*pphi1(k))*occupancy(j,k,1,k)
                                end if               !dead amplitude is zero
                                h1etot_diff_bra(k)= h1etot_diff_bra(k) + product(bra_prod)*ph1ei(j,k)
                                h1etot_diff_ket(k)= h1etot_diff_ket(k) + product(bra_prod)*ph1ei(j,k)          
                            end if
                        else if(equal.eq.2)then 
                            bra_prod=real(pcos1*psin2*occupancy(j,k,1,:)-psin1*pcos2*occupancy(j,k,2,:))            
                            if(j.eq.k)then  !dead amplitude is zero
                                bra_prod(j)=real(pcos1(j)*psin2(j)*occupancy(j,k,1,j))
                            else                   
                                bra_prod(j)=-real(psin1(j)*psin2(j))*occupancy(j,k,2,j)
                                bra_prod(k)=real(pcos1(k)*pcos2(k))*occupancy(j,k,1,k) 
                            end if
                            do l=1,norb
                                temp1=prod
                                temp1(l)=bra_prod(l)
                                h1etot_diff_bra(l)= h1etot_diff_bra(l)+product(temp1)*ph1ei(j,k)   
                            end do 
                        else if(equal.eq.3)then
                            ket_prod=real(psin1*pcos2*occupancy(j,k,1,:)-pcos1*psin2*occupancy(j,k,2,:))                  
                            if(j.eq.k)then  !dead amplitude is zero
                                ket_prod(j)=real(psin1(j)*pcos2(j)*occupancy(j,k,1,j))
                            else                   
                                ket_prod(j)=real(pcos1(j)*pcos2(j))*occupancy(j,k,2,j)!alive amplitude is zero
                                ket_prod(k)=-real(psin1(k)*psin2(k))*occupancy(j,k,1,k) !dead amplitude is zero
                            end if
                            do l=1,norb
                                temp2=prod
                                temp2(l)=ket_prod(l)
                                h1etot_diff_ket(l)= h1etot_diff_ket(l) + product(temp2)*ph1ei(j,k)   
                            end do   
                        end if  
                    end if
                end if
            end do
        end do
        !$omp end target teams distribute parallel do simd

    
        return

    end subroutine one_elec_part_gpu
        
    subroutine two_elec_part_gpu(psin1,psin2,pcos1,pcos2,z1jk,z2l,h2etot,occupancy_2an,&
        occupancy_an,ph2ei,h2etot_diff_bra,h2etot_diff_ket,equal)
        !!$omp declare target 
        implicit none

        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:,:),intent(in)::z1jk
        complex(kind=8),intent(inout)::h2etot
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an
        integer,dimension(:,:,:),intent(in)::occupancy_an
        real(kind=8),dimension(norb),intent(inout)::h2etot_diff_bra,h2etot_diff_ket
        real(kind=8), pointer,dimension(:,:,:,:), intent(in)::ph2ei
        integer,intent(in)::equal
        complex(kind=8), dimension(:),intent(in)::psin1,psin2
        complex(kind=8), dimension(:),intent(in)::pcos1,pcos2
        integer,dimension(2,norb)::occupancy
        integer::j,k,l,jspin
        real(kind=8)::tot
        complex(kind=8)::totc
        complex(kind=8),dimension(2,norb)::vmult
        real(kind=8),dimension(2,norb)::vmult_dd,vmultr
        complex(kind=8),dimension(norb)::gg,hh
        real(kind=8),dimension(norb)::gg_1,hh_1,gg_2,hh_2
        integer::n,p,gmax,hmin,gmax1,hmin1,gmax2,hmin2,breakflag

        real(kind=8),dimension(norb)::temp
        if (errorflag .ne. 0) return

        
        h2etot=(0.0,0.0)
        h2etot_diff_bra=0.0
        h2etot_diff_ket=0.0
     
        !$omp target teams distribute parallel do reduction(+:h2etot) map(to:gmax,hmin,jspin,tot,totc) &
        !$omp & map(alloc:vmult(2,norb),gg(norb),hh(norb),occupancy(2,norb)) &
        !$omp & private(j,k,l,jspin,occupancy,gg,hh,gmax,hmin,vmult,tot,toc,temp) shared(z2l,z1jk,occupancy_2an,occupancy_an,ph2ei)
        do j=1, norb
            if(psin1(j)==(0.0,0.0))then
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
                    if(psin2(l)==(0.0,0.0))then
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
        h2etot=h2etot*0.5
       
        if(GDflg.eq.'y')then
            if(equal.lt.4)then 
                h2etot_diff_ket=0.0
                h2etot_diff_bra=0.0
                !$omp target teams distribute parallel do reduction(+:h2etot_diff_bra,h2etot_diff_ket) &
                !$omp map(to:gg_1,hh_1,gg_2,hh_2,gmax1,hmin1,gmax2,hmin2,jspin,breakflag) &
                !$omp & map(alloc:vmultr(2,norb),vmult_dd(2,norb),occupancy(2,norb),temp(norb)) &
                !$omp & private(gg_1,hh_1,gg_2,hh_2,gmax1,gmax2,hmin1,hmin2,jspin,occupancy,breakflag,vmultr,vmult_dd,temp)&
                !$omp & shared(ph2ei,z1jk,z2l)
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
                             !j,k,l=annihilate1,annihilate1_2,annihilate2
                            if(equal.eq.1) then !Differentiation when zombie states are the same
                                temp=0
                                vmultr=real(z1jk(j,k,:,:))*real(z2l(l,:,:))
                                do n=1, norb    !Differentiating w.r.t to orbital n
                                    breakflag=0
                                    vmult_dd=vmultr
                                    if(j.eq.k)then
                                        temp(n)=0
                                        cycle
                                    else if((j.eq.n).or.(k.eq.n))then
                                        if(l.eq.n)then    !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(1j)
                                            vmult_dd(1,n)=0
                                            vmult_dd(2,n)=2*real(psin1(n)*pcos1(n))*occupancy(2,n) !(sin2x=2sinxcosx)
                                        else
                                            vmult_dd(1,n)=0        !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(0j)
                                            vmult_dd(2,n)=(1.0-2*((real(psin1(n)))**2))*occupancy(2,n)    !(cos2x = cos^x -sin^2x)
                                        end if
                                    else if((j.ne.n).or.(k.ne.n))then
                                        if(l.eq.n)then !before diff alive:0 dead:a^(a)_(0j)*a^(a)_(1j)
                                            vmult_dd(1,n)=0
                                            vmult_dd(2,n)=(1.0-2*((real(psin1(n)))**2))*occupancy(2,n)
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
                                                temp(n)=(1.0-2*((real(psin1(1)))**2))*&
                                                occupancy(2,1)*hh_1(2)*ph2ei(j,k,l,1)
                                                cycle
                                            end if
                                            if(n.eq.1)then
                                                if((j.eq.1).or.(k.eq.1))then
                                                    if(l.ne.1)then
                                                        tot = tot+&
                                                        (2*real(psin1(1)*pcos1(1))*occupancy(2,1)*hh_1(2)*ph2ei(j,k,l,1))
                                                    end if
                                                else
                                                    if(l.ne.1)then
                                                        tot = tot+&
                                                        ((1.0-2*((real(psin1(1)))**2))*occupancy(2,1)*hh_1(2)*ph2ei(j,k,l,1))
                                                    end if 
                                                end if
                                            else
                                                tot = tot+&
                                                (REAL(conjg(z1jk(j,k,2,1))*z2l(l,1,1))*hh_1(2)*ph2ei(j,k,l,1))
                                            end if
                                        end if
                                    else if((breakflag.lt.norb))then
                                        if(breakflag.ne.0)then
                                            temp(n)=(gg_1(n-1)*((1.0-2*((real(psin1(n)))**2))*&
                                            occupancy(2,n))*hh_1(n+1)*ph2ei(j,k,l,n))
                                            cycle
                                        end if
                                        do p=2,norb-1
                                            if(ph2ei(j,k,l,p).ne.0.0) then
                                                if(p.eq.n)then
                                                    if((j.eq.p).or.(k.eq.p))then
                                                        if(l.ne.p)then
                                                            tot = tot+&
                                                            (gg_1(p-1)*(2*real(psin1(p)*pcos1(p))*&
                                                            occupancy(2,p))*hh_1(p+1)*ph2ei(j,k,l,p))
                                                        end if
                                                    else
                                                        if(l.ne.p)then
                                                            tot = tot+&
                                                            (gg_1(p-1)*((1.0-2*((real(psin1(p)))**2))*&
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
                                            temp(n)=(gg_1(norb-1)*((1.0-2*((real(psin1(norb)))**2))&
                                            *occupancy(2,norb))*ph2ei(j,k,l,norb))
                                            cycle
                                        end if
                                        if(ph2ei(j,k,l,norb).ne.0) then
                                            if(norb.eq.n)then
                                                if((j.eq.norb).or.(k.eq.norb))then
                                                    if(l.ne.norb)then
                                                        tot = tot+&
                                                        (gg_1(norb-1)*(2*real(psin1(norb)*pcos1(norb))*&
                                                        occupancy(2,norb))*ph2ei(j,k,l,norb))
                                                    end if
                                                else
                                                    if(l.ne.norb)then
                                                        tot = tot+&
                                                         (gg_1(norb-1)*((1.0-2*((real(psin1(norb)))**2))*&
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
                                h2etot_diff_ket=h2etot_diff_bra
                               
                                
                        
                            else if(equal.eq.2)then !Differentiation when zombie states are not the same
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
                                            vmult_dd(2,n)=real(pcos1(n)*psin2(n))*occupancy(2,n)
                                        else !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(0j)
                                            vmult_dd(1,n)=0
                                            vmult_dd(2,n)=real(pcos1(n)*pcos2(n))*occupancy(2,n)
                                        end if
                                    else if((j.ne.n).or.(k.ne.n))then
                                        if(l.eq.n)then !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                            vmult_dd(1,n)=0
                                            vmult_dd(2,n)=-real(psin1(n)*psin2(n))*occupancy(2,n)
                                        else    !before diff alive:a^(a)_(1j)*a^(b)_(1j) dead:a^(a)_(0j)*a^(b)_(0j)
                                            vmult_dd(1,n)=real(pcos1(n)*psin2(n))*occupancy(1,n)
                                            vmult_dd(2,n)=-real(psin1(n)*pcos2(n))*occupancy(2,n)
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
                                                    ((REAL(pcos1(1)*psin2(1))*occupancy(2,1))*hh_1(2)*ph2ei(j,k,l,1))
                                                end if
                                            else
                                                if(l.ne.1)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                                    tot = tot+&
                                                    ((REAL(-psin1(1)*psin2(1))*occupancy(2,1))*hh_1(2)*ph2ei(j,k,l,1))
                                                end if 
                                            end if
                                        else
                                            tot = tot+&
                                            (REAL(conjg(z1jk(j,k,2,1))*z2l(l,1,1))*hh_1(2)*ph2ei(j,k,l,1))
                                        end if
                                    end if
                                    do p=2,norb-1
                                        if(ph2ei(j,k,l,p).ne.0.0) then
                                            if(p.eq.n)then
                                                if((j.eq.p).or.(k.eq.p))then
                                                    if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                                        tot = tot+&
                                                        (gg_1(p-1)*(REAL(pcos1(p)*psin2(p))*&
                                                        occupancy(2,p))*hh_1(p+1)*ph2ei(j,k,l,p))
                                                    end if
                                                else
                                                    if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                                        tot = tot+&
                                                        (gg_1(p-1)*(REAL(-psin1(p)*psin2(p))*occupancy(2,p))*&
                                                        hh_1(p+1)*ph2ei(j,k,l,p))
                                                    end if 
                                                end if
                                            else
                                                tot= tot+&
                                                 (gg_1(p-1)*REAL(conjg(z1jk(j,k,2,p))*z2l(l,1,p))*hh_1(p+1)*ph2ei(j,k,l,p))
                                            end if
                                        end if
                                    end do
                    
                                    if(ph2ei(j,k,l,norb).ne.0) then
                                        if(norb.eq.n)then
                                            if((j.eq.norb).or.(k.eq.norb))then
                                                if(l.ne.norb)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                                    tot =tot+&
                                                    (gg_1(norb-1)*(REAL(pcos1(norb)*psin2(norb))*&
                                                    occupancy(2,norb))*ph2ei(j,k,l,norb))
                                                end if
                                            else
                                                if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                                    tot = tot+&
                                                    (gg_1(norb-1)*(REAL(-psin1(norb)*psin2(norb))*&
                                                    occupancy(2,norb))*ph2ei(j,k,l,norb))
                                                end if 
                                            end if
                                        else
                                            tot = tot+&
                                            (gg_1(norb-1)*REAL(conjg(z1jk(j,k,2,norb))*z2l(l,1,norb))*ph2ei(j,k,l,norb))
                                        end if
                                    end if
                                    temp(n)=tot
                                end do
                                h2etot_diff_bra=h2etot_diff_bra+temp
                                cycle
                            else if(equal.eq.3)then
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
                                            vmult_dd(2,n)=real(psin1(n)*pcos2(n))*occupancy(2,n)
                                        else !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(0j)
                                            vmult_dd(1,n)=0
                                            vmult_dd(2,n)=-real(psin1(n)*psin2(n))*occupancy(2,n)
                                        end if
                                    else if((j.ne.n).or.(k.ne.n))then
                                        if(l.eq.n)then !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                            vmult_dd(1,n)=0
                                            vmult_dd(2,n)=real(pcos1(n)*pcos2(n))*occupancy(2,n)
                                        else    !before diff alive:a^(a)_(1j)*a^(b)_(1j) dead:a^(a)_(0j)*a^(b)_(0j)
                                            vmult_dd(1,n)=real(psin1(n)*pcos2(n))*occupancy(1,n)
                                            vmult_dd(2,n)=-real(pcos1(n)*psin2(n))*occupancy(2,n)
                                        end if
                                    end if
                    
                                    gg_2(1:norb)=(0.0,0.0)
                                    hh_2(1:norb)=(0.0,0.0)
                                    gmax2=norb
                                    gg_2(1)=(vmult_dd(2,1))-(vmult_dd(1,1))
                    
                                    do p=2, norb
                                        gg_2(p)=gg_2(p-1)*((vmult_dd(2,p)-vmult_dd(1,p)))
                                        if(gg_2(p)==(0.0,0.0))then
                                            gmax2=p
                                            EXIT 
                                        end if
                                    end do
                                    
                    
                                    hmin2=0
                                    hh_2(norb) = (vmult_dd(2,norb)+vmult_dd(1,norb))
                            
                                    do p=(norb-1),1,-(1)
                                        hh_2(p)=hh_2(p+1)*(vmult_dd(2,p)+vmult_dd(1,p))
                                        if(hh_2(p)==(0.0,0.0))then
                                            hmin2=p
                                            EXIT 
                                        end if
                                    end do
                                    tot=0.0
                                    if((gmax2 < hmin2))then
                                        temp(n)=tot
                                        cycle
                                    end if
                    
                                    if(ph2ei(j,k,l,1).ne.0) then
                                        if(n.eq.1)then
                                            if((j.eq.1).or.(k.eq.1))then
                                                if(l.ne.1)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                                    tot = tot+&
                                                    ((REAL(psin2(1)*pcos1(1))*occupancy(2,1))*hh_2(2)*ph2ei(j,k,l,1))
                                                end if
                                            else
                                                if(l.ne.1)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                                    tot =tot+&
                                                    ((REAL(pcos2(1)*pcos1(1))*occupancy(2,1))*hh_2(2)*ph2ei(j,k,l,1))
                                                end if 
                                            end if
                                        else
                                            tot = tot+&
                                            (REAL(conjg(z1jk(j,k,2,1))*z2l(l,1,1))*hh_2(2)*ph2ei(j,k,l,1))
                                        end if
                                    end if
                                    do p=2,norb-1
                                        if(ph2ei(j,k,l,p).ne.0.0) then
                                            if(p.eq.n)then
                                                if((j.eq.p).or.(k.eq.p))then
                                                    if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                                        tot=tot+&
                                                        (gg_2(p-1)*(REAL(psin2(p)*pcos1(p))*&
                                                        occupancy(2,p))*hh_2(p+1)*ph2ei(j,k,l,p))
                                                    end if
                                                else
                                                    if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                                        tot = tot+&
                                                        (gg_2(p-1)*(REAL(pcos2(p)*pcos1(p))*&
                                                        occupancy(2,p))*hh_2(p+1)*ph2ei(j,k,l,p))
                                                    end if 
                                                end if
                                            else
                                                tot = tot+&
                                                (gg_2(p-1)*REAL(conjg(z1jk(j,k,2,p))*z2l(l,1,p))*hh_2(p+1)*ph2ei(j,k,l,p))
                                            end if
                                        end if
                                    end do
                    
                                    if(ph2ei(j,k,l,norb).ne.0) then
                                        if(norb.eq.n)then
                                            if((j.eq.norb).or.(k.eq.norb))then
                                                if(l.ne.norb)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                                    tot= tot+&
                                                    (gg_2(norb-1)*(REAL(psin2(norb)*pcos1(norb))*&
                                                    occupancy(2,norb))*ph2ei(j,k,l,norb))
                                                end if
                                            else
                                                if(l.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                                    tot = tot+&
                                                    (gg_2(norb-1)*(REAL(pcos2(norb)*pcos1(norb))*&
                                                    occupancy(2,norb))*ph2ei(j,k,l,norb))
                                                end if 
                                            end if
                                        else
                                            tot = tot+&
                                            (gg_2(norb-1)*REAL(conjg(z1jk(j,k,2,norb))*z2l(l,1,norb))*ph2ei(j,k,l,norb))
                                        end if
                                    end if

                                    temp(n)=tot

                                end do
                                h2etot_diff_ket=h2etot_diff_ket+temp
                                cycle
                            end if
                           
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do 
                h2etot_diff_bra = h2etot_diff_bra*0.5
                h2etot_diff_ket = h2etot_diff_ket*0.5
            end if
        end if

    
        
        return

    end subroutine two_elec_part_gpu
    
    ! Top level routine to allocate hamiltonian and overlap matrices 
    subroutine hamgen_gpu(ham,zstore,elecs,size,verb)

        implicit none

        
        type(hamiltonian), intent(inout),target::ham 
        type(zombiest),dimension(:),intent(in),target::zstore
        type(elecintrgl),intent(in),target::elecs
        integer,allocatable,dimension(:,:,:)::occupancy_an
        integer,allocatable,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an
        complex(kind=8),allocatable, dimension(:,:,:)::passback
        integer, allocatable,dimension(:)::IPIV1
        integer,intent(in)::size,verb
        complex(kind=8),allocatable,dimension(:)::WORK1
        real(kind=8),dimension(ndet,ndet,norb)::temp2
        
        real(kind=8),pointer,dimension(:,:)::ph1ei
        real(kind=8), dimension(:,:,:,:), pointer::ph2ei
        real(kind=8),pointer :: phnuc

        complex(kind=8), dimension(ndet,norb)::psin
        complex(kind=8), dimension(ndet,norb)::pcos
        real(kind=8),dimension(ndet,norb)::pphi

        complex(kind=8), dimension(:,:), pointer::phjk
        complex(kind=8), dimension(:,:), pointer::povrlp
        complex(kind=8), dimension(:,:), pointer::pinv
        complex(kind=8), dimension(:,:), pointer::pkinvh

        real(kind=8), dimension(:,:,:), pointer::pdiff_hjk
        real(kind=8), dimension(:,:,:), pointer::pdiff_ovrlp
        real(kind=8), dimension(:,:,:,:), pointer::pdiff_invh
        integer:: j,k,l,ierr

       
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
       
        
        occupancy_an=1
        occupancy_2an=1
        occupancy_an_cr=1
        !$omp parallel
        !$omp do
        do j=1,norb
            occupancy_an(j,1,j)=0
            do l=j-1, 1, -1
                occupancy_an(j,1,l)=-1
            end do
        end do
        !$omp end do
        !$omp  do
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
        
    
        ph1ei=>elecs%h1ei
        ph2ei=>elecs%h2ei
        phnuc=>elecs%hnuc 

        phjk=>ham%hjk
        povrlp=>ham%ovrlp
        pinv=>ham%inv
        pkinvh=>ham%kinvh
        pdiff_hjk=>ham%diff_hjk
        pdiff_ovrlp=>ham%diff_ovrlp
        pdiff_invh=>ham%diff_invh
        

        do j=1,ndet
            psin(j,:)=zstore(j)%sin(:)
            pcos(j,:)=zstore(j)%cos(:)
            pphi(j,:)=zstore(j)%phi(:)
        end do
        

        !$omp target data map(to:ph1ei,ph2ei,phnuc,psin,pcos,pphi,occupancy_2an,&
        !$omp occupancy_an_cr,occupancy_an,passback,errorflag,ierr,GDflg) &
        !$omp & map(tofrom:phjk,povrlp) if(GDflg.eq.'y') map(tofrom:pdiff_hjk,pdiff_ovrlp)
        do j=1, size
            call he_row_gpu(ph1ei,ph2ei,phnuc,psin,pcos,pphi,occupancy_2an,&
            & occupancy_an_cr,occupancy_an,passback,phjk,povrlp,pdiff_hjk,pdiff_ovrlp,j,size)
            if(verb.eq.1)then
                write(6,"(a,i0,a)") "Hamiltonian row ",j, " completed"
            end if  
        end do
        !$omp end target data
        
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
        

        deallocate(occupancy_an,stat=ierr)
        if(ierr==0) deallocate(passback,stat=ierr)
        if(ierr==0) deallocate(occupancy_2an,stat=ierr)
        if(ierr==0) deallocate(occupancy_an_cr,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in passback vector deallocation . ierr had value ", ierr
            errorflag=1
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
        
        return
       

    end subroutine hamgen_gpu
    
    
END MODULE ham