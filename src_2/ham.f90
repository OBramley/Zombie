MODULE ham 

    use globvars
    use alarrays
    use operators

    contains
    
    ! Function to generate 1 electron hamiltonian element part
    complex(kind=8) function h1et(z1,z2,elecs)

        implicit none

        type(zombiest),intent(in)::z1,z2
        type(elecintrgl),intent(in)::elecs
        type(zombiest)::zomt
        integer::j,k,ierr
        complex(kind=8)::tot
        
        if (errorflag .ne. 0) return
        ierr = 0

        call alloczf(zomt)
        tot=(0.0,0.0)
        !$omp parallel private(j,k,zomt) shared(tot,elecs,z1,z2)
        !$omp do reduction(+:tot) 
        do j=1, norb
            do k=1, norb
                zomt%alive(1:norb)=z2%alive(1:norb)
                zomt%dead(1:norb)=z2%dead(1:norb)
                call an(zomt,j)
                call cr(zomt,k)
                tot = tot + (overlap(z1,zomt)*elecs%h1ei(j,k))
            end do
        end do
        !$omp end do
        !$omp end parallel
        call dealloczf(zomt)
        
        h1et=tot
        return

    end function h1et

    ! Function to generate 2 electron hamiltonian element part
    complex(kind=8) function h2et(z1,z2,elecs)

        implicit none

        type(zombiest),intent(in)::z1,z2
        type(elecintrgl),intent(in)::elecs
        type(zombiest),allocatable, dimension(:,:)::z1jk
        type(zombiest),allocatable, dimension(:)::z2l
        type(zombiest)::zomt
        integer::j,k,l,jspin, ierr
        complex(kind=8)::tot

        if (errorflag .ne. 0) return
        ierr = 0

      
        call alloczs2d(z1jk,norb)
        call alloczs(z2l,norb)
        call alloczf(zomt)
        tot=(0.0,0.0)

        !$omp parallel shared(z1jk,z1,z2,tot) private(j,k,l,jspin,zomt)
        !$omp do
        do j=1, norb
            do k=1, norb
                zomt%alive(1:norb)=z1%alive(1:norb)
                zomt%dead(1:norb)=z1%dead(1:norb)
                call an(zomt,j)
                call an(zomt,k)
                do l=1, norb
                    z1jk(j,k)%alive(l)=zomt%alive(l)
                    z1jk(j,k)%dead(l)=zomt%dead(l)
                end do
            end do
        end do
        !$omp end do NOWAIT
        !$omp do
        do l=1, norb
            zomt%alive(1:norb)=z2%alive(1:norb)
            zomt%dead(1:norb)=z2%dead(1:norb)
            call an(zomt,l)
            z2l(l)%alive(1:norb)=zomt%alive(1:norb)
            z2l(l)%dead(1:norb)=zomt%dead(1:norb)
        end do
        !$omp end do
        
        !$omp do reduction(+:tot)
        do j=1, norb
            if(z1%alive(j)==(0.0,0.0))then
                CYCLE
            end if

            if(modulo(j,2)==0)then
                jspin=2
            else
                jspin=1
            end if

            do k=1, norb
                if(iszero(z1jk(j,k)).eqv..true.)then
                    CYCLE
                end if
                do l=jspin, norb, 2
                    if(z2%alive(l)==(0.0,0.0))then
                        CYCLE
                    end if
                    tot = tot + z_an_z3(z1jk(j,k),z2l(l),elecs%h2ei(j,k,l,:))
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel

        call dealloczs2d(z1jk)
        call dealloczs(z2l)
        call dealloczf(zomt)

        h2et=tot*0.5
        return

    end function h2et

    
    ! Function to generate indivdual hamiltonian elements
    complex(kind=8) function hamval(z1,z2,elecs,ovrl)

        implicit none

        type(zombiest),intent(in)::z1,z2
        type(elecintrgl),intent(in)::elecs
        complex(kind=8), intent(in)::ovrl

        hamval= h1et(z1,z2,elecs)+h2et(z1,z2,elecs)+((elecs%hnuc)*ovrl)
        return
    
    end function hamval



    ! Top level routine to allocate hamiltonian and overlap matrices 
    subroutine hamgen(ham,zstore,elecs,size)

        implicit none

        type(hamiltonian), intent(inout)::ham 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer, allocatable,dimension(:)::IPIV1
        complex(kind=8),allocatable,dimension(:)::WORK1
        integer:: j,k,size,ierr

        !$omp parallel shared(ham,zstore) private(j,k)
        !$omp do 
        do j=1, size
            do k=j,size
                ham%ovrlp(j,k)=overlap(zstore(j),zstore(k))
                ham%ovrlp(k,j)= ham%ovrlp(j,k)
                ham%inv(j,k)=ham%ovrlp(j,k)
                ham%inv(k,j)=ham%ovrlp(j,k)
                ham%hjk(j,k)= hamval(zstore(j),zstore(k),elecs,ham%ovrlp(j,k))
                ham%hjk(k,j)=ham%hjk(j,k)
            end do
            write(6,"(a,i0,a)") "Hamiltonian row ",j, " completed"
        end do
        !$omp end do
        !$omp end parallel
        
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
        return
        
    end subroutine hamgen



END MODULE ham