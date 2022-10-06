MODULE ham 

    use globvars
    use alarrays
    use operators
    use omp_lib
    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators
    subroutine he_row(ham,zstore,elecs,row,size,occupancy_2an,occupancy_an_cr,occupancy_an)
        
        implicit none
        type(hamiltonian), intent(inout)::ham
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::row,size
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an,occupancy_an_cr
        integer,dimension(:,:,:),intent(in)::occupancy_an
        complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
        complex(kind=8),allocatable, dimension(:,:,:)::z2l
        integer::j,k,l,m,ierr
        complex(kind=8)::h1etot, h2etot,temp
        real(kind=8),dimension(norb)::h1etot_diff_bra,h2etot_diff_bra
        real(kind=8),dimension(norb)::h1etot_diff_ket,h2etot_diff_ket
        real(kind=8),dimension(2,norb)::h1etot_diff,h2etot_diff
        
        if (errorflag .ne. 0) return
        ierr = 0

    
        allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        h1etot=(0.0,0.0)
        h2etot=(0.0,0.0)
        temp=(0.0,0.0)
        h1etot_diff_bra=0.0
        h1etot_diff_ket=0.0
        h2etot_diff_bra=0.0
        h2etot_diff_ket=0.0
        h1etot_diff=0.0
        h2etot_diff=0.0
        
        !$ call omp_set_nested(.TRUE.)
        !!$omp  parallel private(j,k,l,z2l,h1etot,h2etot,temp, &
        !!$omp h1etot_diff_bra,h2etot_diff_bra, h1etot_diff_ket,h2etot_diff_ket,&
        !!$omp  h1etot_diff,h2etot_diff,m) shared(z1jk,zstore,row,size,occupancy_2an,occupancy_an_cr,occupancy_an)
        !$omp parallel
        !$omp do
        do j=1, norb
            do k=1, norb
                z1jk(j,k,1,:)=zstore(row)%sin(:)
                z1jk(j,k,2,:)=zstore(row)%cos(:)
                z1jk(j,k,2,j)=z1jk(j,k,1,j)
                z1jk(j,k,1,j)=(0.0,0.0)
                z1jk(j,k,2,k)=z1jk(j,k,1,k)
                z1jk(j,k,1,k)=(0.0,0.0)
            end do
        end do
        !$omp end do
        !$omp end parallel
        z1jk=z1jk*occupancy_2an
       
        !$omp parallel private(j,k,l,z2l,h1etot,h2etot, &
        !$omp h1etot_diff_bra,h2etot_diff_bra, h1etot_diff_ket,h2etot_diff_ket,&
        !$omp  h1etot_diff,h2etot_diff,m) shared(z1jk,zstore,row,size,occupancy_2an,occupancy_an_cr,occupancy_an)
        !$omp do
        do m=row,size
            h1etot=(0.0,0.0)
            h2etot=(0.0,0.0)
            temp=(0.0,0.0)
            h1etot_diff_bra=0.0
            h1etot_diff_ket=0.0
            h2etot_diff_bra=0.0
            h2etot_diff_ket=0.0
            h1etot_diff=0.0
            h2etot_diff=0.0
            print*,row,m, h1etot, h2etot
            !$omp parallel shared(z2l)
            !$omp do
            do l=1, norb
                z2l(l,1,:)=zstore(m)%sin(:)
                z2l(l,2,:)=zstore(m)%cos(:)
                z2l(l,2,l)=z2l(l,1,l)
                z2l(l,1,l)=(0.0,0.0)
            end do
            !$omp end do
            !$omp end parallel

            !!$omp parallel private(z2l) shared(row,m,z1jk,zstore,size,occupancy_2an,occupancy_an_cr,&
            !!$omp       occupancy_an,h1etot_diff_bra,h2etot_diff_bra,h1etot_diff_ket,h2etot_diff_ket)
            h1etot=(0.0,0.0)
            h2etot=(0.0,0.0)
            h1etot_diff_bra=0.0
            h1etot_diff_ket=0.0
            h2etot_diff_bra=0.0
            h2etot_diff_ket=0.0
            h1etot_diff=0.0
            h2etot_diff=0.0
            !!$omp sections
            !!$omp section 
            if(m.eq.row)then
                print*,"one elec call", row,m
                call one_elec_part(zstore(row),z2l,h1etot,occupancy_an_cr,&
                                    elecs%h1ei,h1etot_diff_bra,h1etot_diff_ket,zstore(m),1)
            else 
                print*,"one elec call", row,m
                call one_elec_part(zstore(row),z2l,h1etot,occupancy_an_cr,&
                                    elecs%h1ei,h1etot_diff_bra,h1etot_diff_ket,zstore(m),9)
            end if
            print*,row,m, h1etot, h2etot
            !!$omp section
            z2l=z2l*occupancy_an
            !$omp flush(z2l)
            if(m.eq.row)then
                print*,"two elec call", row,m
                call two_elec_part(zstore(row),z1jk,z2l,h2etot,occupancy_2an,occupancy_an,&
                                            elecs%h2ei,h2etot_diff_bra,h2etot_diff_ket,zstore(m),1)
            else 
                print*,"two elec call", row,m
                call two_elec_part(zstore(row),z1jk,z2l,h2etot,occupancy_2an,occupancy_an,&
                                            elecs%h2ei,h2etot_diff_bra,h2etot_diff_ket,zstore(m),9)
            end if
            !!$omp end sections
            !!$omp end parallel

            ham%ovrlp(row,m)=overlap(zstore(row),zstore(m))
            ham%ovrlp(m,row)= ham%ovrlp(row,m)
            ham%hjk(row,m)=h1etot+h2etot+(elecs%hnuc*ham%ovrlp(row,m))
            ham%hjk(m,row)=ham%hjk(row,m)
            print*,row,m, h1etot, h2etot
            if(GDflg.eq.'y') then
                ham%diff_hjk_bra(row,m,:)=h1etot_diff_bra+h2etot_diff_bra
                ham%diff_hjk_ket(m,row,:)=h1etot_diff_ket+h2etot_diff_ket
                if(m.eq.row)then
                    h1etot_diff = diff_overlap(zstore(row),zstore(m),0)
                    ham%diff_ovrlp_bra(row,m,:) = h1etot_diff(1,:)
                    ham%diff_ovrlp_ket(row,m,:) = h1etot_diff(1,:)
                else 
                    h1etot_diff = diff_overlap(zstore(row),zstore(m),1)
                    ham%diff_ovrlp_bra(row,m,:) =h1etot_diff(1,:)
                    ham%diff_ovrlp_ket(m,row,:) = h1etot_diff(2,:)
                end if
            end if
        end do
        !$omp end do
        !$omp end parallel
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
        !$omp parallel private(j,k,zomt,temp,h1etot_diff) shared(h1etot,h1ei,occupancy,z2l,&
        !$omp   h1etot_diff_bra,h1etot_diff_ket,zs1,zs2,equal)
        !$omp do
        do j=1, norb
            do k=1, norb
                zomt(:,:)=z2l(j,:,:)
                zomt(1,k)=zomt(2,k)
                zomt(2,k)=(0.0,0.0)
                zomt=zomt*occupancy(j,k,:,:)
                temp=product((conjg(zs1%sin)*zomt(1,:))+(conjg(zs1%cos)*zomt(2,:)))*h1ei(j,k)
                !$omp critical
                h1etot=h1etot+temp
                !$omp end critical
                if(GDflg.eq.'y')then
                    if(equal.eq.1)then
                        h1etot_diff = diff_overlap_cran(zs1,zs2,0,zomt,j,k,occupancy(j,k,:,:))*h1ei(j,k) 
                    else 
                        h1etot_diff = diff_overlap_cran(zs1,zs2,1,zomt,j,k,occupancy(j,k,:,:))*h1ei(j,k) 
                    end if
                    !$omp critical
                    h1etot_diff_bra= h1etot_diff_bra + h1etot_diff(1,:)
                    h1etot_diff_ket= h1etot_diff_ket + h1etot_diff(2,:)
                    !$omp end critical
                end if
            end do
        end do
        !$omp end do
        !$omp end parallel
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
        ! print*,'h2etotstart',equal, h2etot
        ! !$omp parallel private(j,k,l,h2etot_diff,tot,occupancy,jspin) shared(h2ei,occupancy_2an,occupancy_an,&
        ! !$omp   z2l,h2etot_diff_bra,h2etot_diff_ket,zs1,zs2,equal,h2etot,z1jk) 
        h2etot=(0.0,0.0)
        h2etot_diff_bra=0.0
        h2etot_diff_ket=0.0
        h2etot_diff=0.0
        ! !$omp flush
        ! $omp do reduction(+:h2etot)
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
                if(occ_iszero(z1jk(j,k,:,:)).eqv..true.)then
                    CYCLE
                end if
                !!$omp do ordered
                do l=jspin, norb, 2
                    if(zs2%sin(l)==(0.0,0.0))then
                        CYCLE
                    end if
                    tot = z_an_z3(z1jk(j,k,:,:),z2l(l,:,:),h2ei(j,k,l,:))
                    !!$omp critical
                    h2etot=h2etot+tot
                    !!$omp end critical
                    if(GDflg.eq.'y')then
                        occupancy=occupancy_2an(j,k,:,:)*occupancy_an(l,:,:)
                        if(equal.eq.1)then
                            h2etot_diff = z_an_z3_diff(z1jk(j,k,:,:),z2l(l,:,:),h2ei(j,k,l,:),0,&
                            occupancy,j,k,l,zs1,zs2)
                        else 
                            h2etot_diff = z_an_z3_diff(z1jk(j,k,:,:),z2l(l,:,:),h2ei(j,k,l,:),1,&
                            occupancy,j,k,l,zs1,zs2) 
                        end if
                        !!$omp critical
                        h2etot_diff_bra = h2etot_diff_bra+h2etot_diff(1,:)
                        h2etot_diff_ket = h2etot_diff_ket+ h2etot_diff(2,:)
                        !!$omp end critical
                    end if
                end do
                !!$omp end do
            end do
        end do
        ! !$omp end do 
        ! !$omp end parallel
        h2etot=h2etot*0.5
        h2etot_diff_bra=h2etot_diff_bra*0.5
        h2etot_diff_ket=h2etot_diff_ket*0.5
        !print*,'h2etotfinish',equal, h2etot

    end subroutine two_elec_part
    
    ! Top level routine to allocate hamiltonian and overlap matrices 
    subroutine hamgen(ham,zstore,elecs,size)

        implicit none

        type(hamiltonian), intent(inout)::ham 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,allocatable,dimension(:,:,:)::occupancy_an
        integer,allocatable,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an
        integer, allocatable,dimension(:)::IPIV1
        complex(kind=8),allocatable,dimension(:)::WORK1
        real(kind=8),dimension(ndet,ndet,norb)::temp2
        real(kind=8),dimension(ndet,norb)::matrix
        integer:: j,k,l,size,ierr,thread_1,thread_2

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
        ! call omp_set_nested(.TRUE.)
        thread_1=0
        thread_2=0
        !$  if(omp_get_num_threads().le.2)then 
        !$      thread_1=omp_get_num_threads()
        !$      thread_2=1
        !$  else if(omp_get_num_threads().le.4)then
        !$      thread_1=2
        !$      thread_2=1
        !$  else
        !$      thread_1=2
        !$      thread_2=2
        !$  end if
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
        
        !$omp parallel private(j) shared(occupancy_an,occupancy_2an,occupancy_an_cr,ham,zstore,elecs)
        !$omp do 
        do j=1, size
            call he_row(ham,zstore,elecs,j,size,occupancy_2an,occupancy_an_cr,occupancy_an)  
            write(6,"(a,i0,a)") "Hamiltonian row ",j, " completed"
        end do
        !$omp end do
        !$omp end parallel
        !$ call omp_set_nested(.FALSE.)
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
            write(0,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
            errorflag=1
        end if

        if (ierr==0) deallocate(WORK1,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
            errorflag=1
        end if
        !$omp parallel
        !$omp workshare
        ham%kinvh=matmul(ham%inv,ham%hjk)
        !$omp end workshare
        !$omp end parallel

        ! To find the derrivative of the inverse of a matrix we use the identity
        ! d(omega^-1)= omega^-1d(omega)omega^-1. Since this result is only used when being multiplied with the hamiltonian 
        ! we will actually calculate omega^-1d(omega)omega^-1*H 
        ! j decides which zombie state we differentiate by
        ! omega^-1d(omega) is first calcualted. creating an ndet*ndet*norb array 
        ! this tensor is then 
        if(GDflg.eq.'y')then
            do j=1, ndet
                ! do l=1, ndet
                    do k=1, ndet
                        if(j.eq.k)then
                            ! vector=REAL(ham%ovrlp(k,:))
                            ! matrix=ham%diff_ovrlp_ket(:,j,:)
                            temp2(k,j,:)=matmul(REAL(ham%ovrlp(k,:)),ham%diff_ovrlp_ket(:,j,:))   !matmul(vector,matrix)
                        else
                            do l=1, ndet
                                temp2(l,k,:)=REAL(ham%ovrlp(l,j))*ham%diff_ovrlp_bra(j,k,:)
                            end do   
                        end if
                    end do
                do k=1, ndet
                    ! temp(1:norb)=0
                    matrix=temp2(k,:,:)
                    do l=1, ndet
                        ham%diff_invh(j,k,l,:)=matmul(matrix,REAL(ham%kinvh(1:ndet,l)))
                    end do
                end do
            end do
        end if
        return
        
    end subroutine hamgen



END MODULE ham