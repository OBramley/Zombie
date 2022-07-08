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
         
        do j=1, norb
            do k=1, norb
                zomt%alive(1:norb)=z2%alive(1:norb)
                zomt%dead(1:norb)=z2%dead(1:norb)
                call an(zomt,j)
                call cr(zomt,k)
                tot = tot + (overlap(z1,zomt)*elecs%h1ei(j,k))
            end do
        end do
        h1et=tot
        return

    end function h1et

    ! Function to generate 1 electron hamiltonian element part
    complex(kind=8) function h2et(z1,z2,elecs)

        implicit none

        type(zombiest),intent(in)::z1,z2
        type(elecintrgl),intent(in)::elecs
        type(zombiest),allocatable, dimension(:,:)::z1jk
        type(zombiest),allocatable, dimension(:)::z2l
        type(zombiest)::zomt
        ! complex(kind=8),dimension(norb,norb,norb,2)::z1jk
        ! complex(kind=8),dimension(norb,norb,2)
        integer::j,k,l,jspin, ierr
        complex(kind=8)::tot

        if (errorflag .ne. 0) return
        ierr = 0


        call alloczs2d(z1jk,norb)
        call alloczs(z2l,norb)
        call alloczf(zomt)

        do j=1, norb
            do k=1, norb
                zomt%alive(1:norb)=z1%alive(1:norb)
                zomt%dead(1:norb)=z1%dead(1:norb)
                call an(zomt,j)
                call an(zomt,k)
                z1jk(j,k)%alive(1:norb)=zomt%alive(1:norb)
                z1jk(j,k)%dead(1:norb)=zomt%dead(1:norb)
            end do
        end do

        do l=1, norb
            zomt%alive(1:norb)=z2%alive(1:norb)
            zomt%dead(1:norb)=z2%dead(1:norb)
            call an(zomt,l)
            z2l(l)%alive(1:norb)=zomt%alive(1:norb)
            z2l(l)%dead(1:norb)=zomt%dead(1:norb)
        end do

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
                    tot = tot + z_an_z3(z1jk(j,k),z2l(l),elecs%h2ei(j,k,l,1:norb))
                end do
            end do
        end do

        call dealloczs2d(z1jk)
        call dealloczs(z2l)
        call dealloczf(zomt)

        h2et=tot
        return






    end function h2et

    
    ! Function to generate indivdual hamiltonian elements
    complex(kind=8) function hamval(z1,z2,elecs,ovrl)

        implicit none

        type(zombiest),intent(in)::z1,z2
        type(elecintrgl),intent(in)::elecs
        complex(kind=8), intent(in)::ovrl

        hamval= h1et(z1,z2,elecs)+h2et(z1,z2,elecs)+(elecs%hnuc*ovrl)

        return
    
        end function hamval



    ! Top level routine to allocate hamiltonian and overlap matrices 
    subroutine hamgen(ham,zstore,elecs,size)

        implicit none

        type(hamiltonian), intent(inout)::ham 
        type(zombiest),dimension(:),intent(inout)::zstore
        type(elecintrgl),intent(inout)::elecs
        integer:: j,k,size

        do j=1, size
            do k=j,size
                ham%ovrlp(j,k)=overlap(zstore(j),zstore(k))
                ham%ovrlp(k,j)= ham%ovrlp(j,k)
                ham%hjk(j,k)= hamval(zstore(j),zstore(k),elecs,ham%ovrlp(j,k))
                ham%hjk(k,j)=ham%hjk(j,k)
            end do
        end do

    end subroutine hamgen



END MODULE ham