MODULE ham 

    use globvars
    use alarrays
    use operators

    contains

    subroutine he_row(ham,zstore,elecs,row,size)
        
        implicit none
        type(hamiltonian), intent(inout)::ham
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::row,size
        type(zombiest),allocatable, dimension(:,:)::z1jk
        type(zombiest),allocatable, dimension(:)::z2l
        type(zombiest)::zomt
        integer::j,k,l,m,jspin,ierr
        complex(kind=8)::h1etot, h2etot,temp
        real(kind=8),dimension(norb)::h1etot_diff,h2etot_diff,diff_temp 

        if (errorflag .ne. 0) return
        ierr = 0

        call alloczf(zomt)
        call alloczs2d(z1jk,norb)
        call alloczs(z2l,norb)
        h1etot=(0.0,0.0)
        h2etot=(0.0,0.0)
        temp=(0.0,0.0)
        ! z1=zstore(row)
        !$omp parallel shared(z1jk,zstore) private(j,k,l,jspin,zomt,z2l,h1etot,h2etot,h1etot_diff,h2etot_diff)
        !$omp do
        do j=1, norb
            do k=1, norb
                zomt=zstore(row)
                call an(zomt,j)
                call an(zomt,k)
                z1jk(j,k)=zomt
            end do
        end do
        !$omp end do
        !$omp do
        do m=row,size
            h1etot=(0.0,0.0)
            h2etot=(0.0,0.0)
            temp=(0.0,0.0)
            
            do l=1, norb
                zomt=zstore(m)
                call an(zomt,l)
                z2l(l)=zomt
            end do
            
            do j=1, norb
                do k=1, norb
                    temp=(0.0,0.0)
                    zomt=z2l(j)
                    call cr(zomt,k)
                    temp = (overlap(zstore(row),zomt)*elecs%h1ei(j,k))
                    h1etot=h1etot+temp
                end do
            end do
  
            do j=1, norb
                if(zstore(row)%alive(j)==(0.0,0.0))then
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
                        if(zstore(m)%alive(l)==(0.0,0.0))then
                            CYCLE
                        end if
                        h2etot = h2etot + z_an_z3(z1jk(j,k),z2l(l),elecs%h2ei(j,k,l,:))
                    end do
                end do
            end do
            h2etot=h2etot*0.5
            

            ham%ovrlp(row,m)=overlap(zstore(row),zstore(m))
            ham%ovrlp(m,row)= ham%ovrlp(row,m)
            ham%hjk(row,m)=h1etot+h2etot+(elecs%hnuc*ham%ovrlp(row,m))
            ham%hjk(m,row)=ham%hjk(row,m)
        end do
        !$omp end do
        !$omp end parallel

        call dealloczs2d(z1jk)
        call dealloczs(z2l)
        call dealloczf(zomt)

        return

    end subroutine he_row


        
    
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
                zomt=z2
                ! zomt%alive(1:norb)=z2%alive(1:norb)
                ! zomt%dead(1:norb)=z2%dead(1:norb)
                call an(zomt,j)
                call cr(zomt,k)
                tot = tot + (overlap(z1,zomt)*elecs%h1ei(j,k))
            end do
        end do
        !$omp end do
        !$omp end parallel
        call dealloczf(zomt)
        
        h1et=tot
        ! print*,h1et
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
                zomt=z1
                ! zomt%alive(1:norb)=z1%alive(1:norb)
                ! zomt%dead(1:norb)=z1%dead(1:norb)
                call an(zomt,j)
                call an(zomt,k)
                z1jk(j,k)=zomt
                ! do l=1, norb
                !     z1jk(j,k)%alive(l)=zomt%alive(l)
                !     z1jk(j,k)%dead(l)=zomt%dead(l)
                ! end do
            end do
        end do
        !$omp end do NOWAIT
        !$omp do
        do l=1, norb
            zomt=z2
            ! zomt%alive(1:norb)=z2%alive(1:norb)
            ! zomt%dead(1:norb)=z2%dead(1:norb)
            call an(zomt,l)
            z2l(l)=zomt
            ! z2l(l)%alive(1:norb)=zomt%alive(1:norb)
            ! z2l(l)%dead(1:norb)=zomt%dead(1:norb)
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
        ! print*,h2et
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
        integer:: j,size,ierr


        do j=1, size
                call he_row(ham,zstore,elecs,j,size)
            write(6,"(a,i0,a)") "Hamiltonian row ",j, " completed"
        end do

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
        return
        
    end subroutine hamgen



END MODULE ham