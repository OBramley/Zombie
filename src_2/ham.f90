MODULE ham 

    use globvars
    use alarrays
    use operators

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
        ! real(kind=8),dimension(:)
        ! type(zombiest),allocatable, dimension(:,:)::z1jk
        ! type(zombiest),allocatable, dimension(:)::z2l
        complex(kind=8),allocatable,dimension(:,:)::zomt
        ! type(zombiest)::zomt
        integer::j,k,l,m,jspin,ierr
        complex(kind=8)::h1etot, h2etot,temp
        real(kind=8),dimension(norb)::h1etot_diff_bra,h2etot_diff_bra
        real(kind=8),dimension(norb)::h1etot_diff_ket,h2etot_diff_ket
        real(kind=8),dimension(2,norb)::h1etot_diff,h2etot_diff
        if (errorflag .ne. 0) return
        ierr = 0

        ! call alloczf(zomt)
        ! call alloczs2d(z1jk,norb)
        ! call alloczs(z2l,norb)
        allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if(ierr==0) allocate(zomt(2,norb),stat=ierr)
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
        !!$omp parallel private(j,k,l,jspin,zomt,z2l,h1etot,h2etot,temp, &
        !!$omp h1etot_diff_bra,h2etot_diff_bra, h1etot_diff_ket,h2etot_diff_ket,&
        !!$omp  h1etot_diff,h2etot_diff,m) shared(z1jk,zstore,row,size,occupancy_2an,occupancy_an_cr,occupancy_an)
        !!$omp do
        do j=1, norb
            do k=1, norb
                z1jk(j,k,1,:)=zstore(row)%sin(:)
                z1jk(j,k,2,:)=zstore(row)%cos(:)
                z1jk(j,k,2,j)=z1jk(j,k,1,j)
                z1jk(j,k,1,j)=(0.0,0.0)
                z1jk(j,k,2,k)=z1jk(j,k,1,k)
                z1jk(j,k,1,k)=(0.0,0.0)
                ! zomt=zstore(row)
                ! call an(zomt,j)
                ! call an(zomt,k)
                ! z1jk(j,k)=zomt
            end do
        end do
        z1jk=z1jk*occupancy_2an
        !!$omp end do
        !!$omp flush(z1jk)
        !!$omp barrier
        !!$omp do
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
            
            !!$omp do
            do l=1, norb
                z2l(l,1,:)=zstore(m)%sin(:)
                z2l(l,2,:)=zstore(m)%cos(:)
                z2l(l,2,l)=z2l(l,1,l)
                z2l(l,1,l)=(0.0,0.0)
                ! zomt=zstore(m)
                ! call an(zomt,l)
                ! z2l(l)=zomt
            end do
            z2l=z2l*occupancy_an
            
            do j=1, norb
                do k=1, norb
                    l=l+1
                    temp=(0.0,0.0)
                    zomt(:,:)=z2l(j,:,:)
                    zomt(1,k)=zomt(2,k)
                    zomt(2,k)=(0.0,0.0)
                    zomt=zomt*occupancy_an_cr(j,k,:,:)
                    temp=product((conjg(zstore(row)%sin)*zomt(1,:))+(conjg(zstore(row)%cos)*zomt(2,:)))*elecs%h1ei(j,k)
                    print*,temp
                    ! call cr(zomt,k)
                    ! temp = (overlap(zstore(row),zomt)*elecs%h1ei(j,k))
                    h1etot=h1etot+temp
                    if(GDflg.eq.'y')then
                        if(m.eq.row)then
                            h1etot_diff = diff_overlap_cran(zstore(row),zstore(m),0,&
                                            zomt,j,k,occupancy_an_cr(j,k,:,:))*elecs%h1ei(j,k) 
                        else 
                            h1etot_diff = diff_overlap_cran(zstore(row),zstore(m),1,&
                                            zomt,j,k,occupancy_an_cr(j,k,:,:))*elecs%h1ei(j,k) 
                        end if
                        h1etot_diff_bra= h1etot_diff_bra + h1etot_diff(1,:)
                        h1etot_diff_ket= h1etot_diff_ket + h1etot_diff(2,:)
                    end if
                end do
            end do

            print*, row, m
            print*, h1etot
            stop
            do j=1, norb
                if(zstore(row)%sin(j)==(0.0,0.0))then
                    CYCLE
                end if

                if(modulo(j,2)==0)then
                    jspin=2
                else
                    jspin=1
                end if

                do k=1, norb
                    if(occ_iszero(occupancy_2an(j,k,:,:)).eqv..true.)then
                        CYCLE
                    end if
                    do l=jspin, norb, 2
                        if(zstore(m)%sin(l)==(0.0,0.0))then
                            CYCLE
                        end if
                        h2etot = h2etot + z_an_z3(z1jk(j,k,:,:),z2l(l,:,:),elecs%h2ei(j,k,l,:))
                        if(GDflg.eq.'y')then
                            if(m.eq.row)then
                                h2etot_diff = z_an_z3_diff(z1jk(j,k,:,:),z2l(l,:,:),elecs%h2ei(j,k,l,:),0,&
                                (occupancy_2an(j,k,:,:)*occupancy_an(l,:,:)),j,k,l,zstore(row),zstore(m))
                            else 
                                h2etot_diff = z_an_z3_diff(z1jk(j,k,:,:),z2l(l,:,:),elecs%h2ei(j,k,l,:),1,&
                                (occupancy_2an(j,k,:,:)*occupancy_an(l,:,:)),j,k,l,zstore(row),zstore(m)) 
                            end if
                            h2etot_diff_bra = h2etot_diff_bra+h2etot_diff(1,:)
                            h2etot_diff_ket = h2etot_diff_ket+ h2etot_diff(2,:)
                        end if
                    end do
                end do
            end do
            h2etot=h2etot*0.5
            h2etot_diff_bra=h2etot_diff_bra*0.5
            h2etot_diff_ket=h2etot_diff_ket*0.5
            
            print*, row,m
            ham%ovrlp(row,m)=overlap(zstore(row),zstore(m))
            ham%ovrlp(m,row)= ham%ovrlp(row,m)
            ham%hjk(row,m)=h1etot+h2etot+(elecs%hnuc*ham%ovrlp(row,m))
            ham%hjk(m,row)=ham%hjk(row,m)

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
                ! ham%diff_ovrlp_bra(m,row,:)=ham%diff_ovrlp_bra(row,m,:)
                ! ham%diff_ovrlp_ket(m,row,:)=ham%diff_ovrlp_ket(row,m,:)
                ! ham%diff_hjk_bra(m,row,:)=ham%diff_hjk_bra(row,m,:)
                ! ham%diff_hjk_ket(row,m,:)=
                ! ham%diff_hjk_ket(row,m,:)

            end if
        end do
        ! !$omp end do
       !! $omp end parallel

        deallocate(z1jk,stat=ierr)
        if(ierr==0) deallocate(z2l,stat=ierr)
        if(ierr==0) deallocate(zomt,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        ! call dealloczs2d(z1jk)
        ! call dealloczs(z2l)
        ! call dealloczf(zomt)

        return

    end subroutine he_row


        
    
    ! Function to generate 1 electron hamiltonian element part
    ! complex(kind=8) function h1et(z1,z2,elecs)

    !     implicit none

    !     type(zombiest),intent(in)::z1,z2
    !     type(elecintrgl),intent(in)::elecs
    !     type(zombiest)::zomt
    !     integer::j,k,ierr
    !     complex(kind=8)::tot
        
    !     if (errorflag .ne. 0) return
    !     ierr = 0

    !     call alloczf(zomt)
    !     tot=(0.0,0.0)
    !     !$omp parallel private(j,k,zomt) shared(tot,elecs,z1,z2)
    !     !$omp do reduction(+:tot) 
    !     do j=1, norb
    !         do k=1, norb
    !             zomt=z2
        
    !             call an(zomt,j)
    !             call cr(zomt,k)
    !             tot = tot + (overlap(z1,zomt)*elecs%h1ei(j,k))
    !         end do
    !     end do
    !     !$omp end do
    !     !$omp end parallel
    !     call dealloczf(zomt)
        
    !     h1et=tot
    !     ! print*,h1et
    !     return

    ! end function h1et

    ! Function to generate 2 electron hamiltonian element part
    ! complex(kind=8) function h2et(z1,z2,elecs)

    !     implicit none

    !     type(zombiest),intent(in)::z1,z2
    !     type(elecintrgl),intent(in)::elecs
    !     type(zombiest),allocatable, dimension(:,:)::z1jk
    !     type(zombiest),allocatable, dimension(:)::z2l
    !     type(zombiest)::zomt
    !     integer::j,k,l,jspin, ierr
    !     complex(kind=8)::tot

    !     if (errorflag .ne. 0) return
    !     ierr = 0

      
    !     call alloczs2d(z1jk,norb)
    !     call alloczs(z2l,norb)
    !     call alloczf(zomt)
    !     tot=(0.0,0.0)

    !     !$omp parallel shared(z1jk,z1,z2,tot) private(j,k,l,jspin,zomt)
    !     !$omp do
    !     do j=1, norb
    !         do k=1, norb
    !             zomt=z1
    !             call an(zomt,j)
    !             call an(zomt,k)
    !             z1jk(j,k)=zomt
        
    !         end do
    !     end do
    !     !$omp end do NOWAIT
    !     !$omp do
    !     do l=1, norb
    !         zomt=z2
    
    !         call an(zomt,l)
    !         z2l(l)=zomt
            
    !     end do
    !     !$omp end do
        
    !     !$omp do reduction(+:tot)
    !     do j=1, norb
    !         if(z1%alive(j)==(0.0,0.0))then
    !             CYCLE
    !         end if

    !         if(modulo(j,2)==0)then
    !             jspin=2
    !         else
    !             jspin=1
    !         end if

    !         do k=1, norb
    !             if(iszero(z1jk(j,k)).eqv..true.)then
    !                 CYCLE
    !             end if
    !             do l=jspin, norb, 2
    !                 if(z2%alive(l)==(0.0,0.0))then
    !                     CYCLE
    !                 end if
    !                 tot = tot + z_an_z3(z1jk(j,k),z2l(l),elecs%h2ei(j,k,l,:))
    !             end do
    !         end do
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    !     call dealloczs2d(z1jk)
    !     call dealloczs(z2l)
    !     call dealloczf(zomt)

    !     h2et=tot*0.5
    !     ! print*,h2et
    !     return

    ! end function h2et

    
    ! Function to generate indivdual hamiltonian elements
    ! Old version not currently used but useful to keep
    ! complex(kind=8) function hamval(z1,z2,elecs,ovrl)

    !     implicit none

    !     type(zombiest),intent(in)::z1,z2
    !     type(elecintrgl),intent(in)::elecs
    !     complex(kind=8), intent(in)::ovrl

    !     hamval= h1et(z1,z2,elecs)+h2et(z1,z2,elecs)+((elecs%hnuc)*ovrl)
    !     return
    
    ! end function hamval



    ! Top level routine to allocate hamiltonian and overlap matrices 
    subroutine hamgen(ham,zstore,elecs,size)

        implicit none

        type(hamiltonian), intent(inout)::ham 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,allocatable,dimension(:,:,:)::occupancy_an
        integer,allocatable,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an
        ! real(kind=8),dimension(:)::shift_an
        ! real(kind=8)(:,:)::shift_2an, shift_an_cr
        integer, allocatable,dimension(:)::IPIV1
        complex(kind=8),allocatable,dimension(:)::WORK1
        ! real(kind=8),dimension(norb)::temp
        real(kind=8),dimension(ndet,ndet,norb)::temp2
        ! real(kind=8),dimension(ndet)::vector
        real(kind=8),dimension(ndet,norb)::matrix
        integer:: j,k,l,size,ierr

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
      
        !$omp parallel private(j) shared(occupancy_an,occupancy_2an,occupancy_an_cr,ham,zstore,elecs)
        !$omp do 
        do j=1,norb
            occupancy_an(j,2,j)=occupancy_an(j,1,j)
            occupancy_an(j,1,j)=0
            occupancy_an(j,1,1:(j-1))=occupancy_an(j,1,1:(j-1))*(-1)
        end do
        !$omp end do
        !$omp barrier
        !$omp do
        do j=1,norb
            do k=1,norb
                occupancy_2an(j,k,:,:)=occupancy_an(j,:,:)
                occupancy_2an(j,k,2,k)=occupancy_2an(j,k,1,k)
                occupancy_2an(j,k,1,k)=0
                occupancy_2an(j,k,1,1:(k-1))=occupancy_2an(j,k,1,1:(k-1))*(-1)
                occupancy_an_cr(j,k,:,:)=occupancy_an(j,:,:)
                occupancy_an_cr(j,k,1,k)=occupancy_an_cr(j,k,2,k)
                occupancy_an_cr(j,k,2,k)=0
                occupancy_an_cr(j,k,1,1:(k-1))=occupancy_an_cr(j,k,1,1:(k-1))*(-1)
            end do
        end do
        !$omp end do
        !$omp do
        do j=1, 1
            call he_row(ham,zstore,elecs,j,size,occupancy_2an,occupancy_an_cr,occupancy_an)  
            write(6,"(a,i0,a)") "Hamiltonian row ",j, " completed"
        end do
        !$omp end do
        !$omp end parallel
        ham%inv=ham%ovrlp
        allocate(IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        
        allocate(WORK1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if   

        call ZGETRF(size,size,ham%inv,size,IPIV1,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if
        call ZGETRI(size,ham%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if

        deallocate(IPIV1,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if

        deallocate(WORK1,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
            errorflag=1
            return
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
            print*,ham%diff_ovrlp_bra(2,2,:)
            print*,ham%diff_ovrlp_bra(2,4,:) 
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