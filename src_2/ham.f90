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
        ! print*,'h1ei complete'
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
                    ! print*,j
                end do
            end if
            !!$omp end do
        end if
        !!$omp end parallel
      
        h2etot=tot*0.5
        h2etot_diff_bra = h2etot_diff(1,:)*0.5
        h2etot_diff_ket = h2etot_diff(2,:)*0.5
        ! print*,'h2ei complete'
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
       
        !!$omp parallel private(j) shared(occupancy_an,occupancy_2an,occupancy_an_cr,ham,zstore,elecs)
        !!$omp do ordered 
        !!$omp declare mapper(elcmap : type(elecintrgl)::elecs) map(elecs%h1ei,elecs%h2ei,elecs%hnuc)
        !!$omp declare mapper(hamap : type(hamiltonian)::ham) map(ham%hjk,ham%inv,ham%inv,ham%ovrlp)
        !!$omp declare mapper(zmap : type(zombiest)::zstore(1:ndet)) map(zs%alive,zs%dead,zs%phi,zs%cos,zs%sin,zs%update_num)
        ! !$omp target data map(alloc:ham%hjk(1:ndet,1:ndet))
        !!$omp target data map(tofrom:ham%hjk(1:ndet,1:ndet))
        !$omp target data map(to:occupancy_2an,occupancy_an_cr,occupancy_an,passback)
        do j=1, size
            call he_row(ham,zstore,elecs,j,size,occupancy_2an,occupancy_an_cr,occupancy_an,passback)
            if(verb.eq.1)then
                write(6,"(a,i0,a)") "Hamiltonian row ",j, " completed"
            end if  
        end do
        !$omp end target data
        !!$omp end target data
        ! print*, ham%diff_hjk(2,1,:)
        ! print*, ham%diff_ovrlp(2,1,:)
        !$ call omp_set_nested(.FALSE.)
      
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
                    if(l.eq.j)then
                        temp2(k,j,:)=matmul(REAL(ham%inv(k,:)),ham%diff_ovrlp(j,:,:))
                    else
                        temp2(k,l,:)=real(ham%inv(k,l))*ham%diff_ovrlp(j,l,:)
                    end if
                end do
            end do
            do k=1, ndet
                do l=1, ndet
                    ham%diff_invh(j,k,l,:)=matmul(transpose(temp2(k,:,:)),real(ham%kinvh(:,l)))*(-1)
                end do
            end do
        end if
        return
        
    end subroutine hamgen



END MODULE ham