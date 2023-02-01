MODULE ham

    use globvars
    use ham_2
    contains

    ! DO I NEED THE DERIVATIVE OF THE OVERLAP WHEN CALCUALTING HAMILTONIAN DERIVATIVE!!!
    !sums all values in a a 1D array
    real(kind=8) function arr_sum_self(to_sum)

        implicit none
        real(kind=8),dimension(:),intent(in)::to_sum
        integer::j

        arr_sum_self=0.0

        do j=1,size(to_sum)
            arr_sum_self=arr_sum_self+to_sum(j)
        end do

        return

    end function arr_sum_self

    !multiplies all values in an array
    real(kind=8) function arr_mult_self(to_sum)

        implicit none
        real(kind=8),dimension(:),intent(in)::to_sum
        integer::j

        arr_mult_self=1.0

        do j=1,size(to_sum)
            arr_mult_self=arr_mult_self*to_sum(j)
        end do

        return

    end function arr_mult_self

    !multilies two arrays elementally
    subroutine arr_mult(a,b,c)

        implicit none
        real(kind=8),dimension(:),intent(in)::a,b
        real(kind=8),dimension(:),intent(out)::c
        integer::j

        do j=1, size(a)
            c(j)=a(j)*b(j)
        end do 

    end subroutine arr_mult

    !adds elements of two arrays together
    subroutine arr_sum(a,b,c)

        implicit none
        real(kind=8),dimension(:),intent(in)::a,b
        real(kind=8),dimension(:),intent(out)::c
        integer::j

        do j=1, size(a)
            c(j)=a(j)+b(j)
        end do 

    end subroutine arr_sum

    subroutine one_electron_vals(alive1,dead1,alive2,dead2,el,h1etot) 
 
        implicit none
  
        real(kind=8), dimension(:),intent(in)::alive1,dead1,alive2,dead2,el
        real(kind=8),dimension(:),intent(out)::h1etot 
        integer::j
 
       
        do j=1,size(h1etot)
            h1etot(j)=h1e(alive1,dead1,alive2,dead2,el(j),j)
        end do  
        
        return

    end subroutine one_electron_vals

    subroutine two_electron_vals(alive1,dead1,alive2,dead2,el,h2etot) 

        implicit none
    
        real(kind=8), dimension(:),intent(in)::alive1,dead1,alive2,dead2,el
        real(kind=8),dimension(:),intent(out)::h2etot 
        integer::j
    
        
        do j=1,size(h2etot)
            h2etot(j)=h2e(alive1,dead1,alive2,dead2,el(j),j)
        end do  
        
        return
        
    end subroutine two_electron_vals

    subroutine ham_make(ham,ovrlp,alive,dead,el1,el2,hnuc) 

        implicit none
        
        real(kind=8),dimension(:,:),intent(inout)::ham 
        real(kind=8),dimension(:,:),intent(in)::ovrlp
        real(kind=8),dimension(:,:),intent(in)::alive,dead
        real(kind=8),dimension(:),intent(in)::el1,el2
        real(kind=8),allocatable,dimension(:)::h1etot,h2etot
        real(kind=8),intent(in)::hnuc
        real(kind=8)::h1,h2
        integer::j,k,ierr
    
        if (errorflag .ne. 0) return 
    
        ierr=0
        allocate(h1etot(size(el1)),stat=ierr)
        allocate(h2etot(size(el2)),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation. ierr had value ,", ierr 
            errorflag=1
            return
        end if
    
        do j=1,ndet
            do k=j,ndet
                call one_electron_vals(alive(:,j),dead(:,j),alive(:,k),dead(:,k),el1,h1etot)
                call two_electron_vals(alive(:,j),dead(:,j),alive(:,k),dead(:,k),el2,h2etot)
                h1=arr_sum_self(h1etot)
                h2=arr_sum_self(h2etot)
                ham(j,k)=h1+0.5*h2+hnuc*ovrlp(j,k); ham(k,j)=ham(j,k)
            end do
        end do 
    
    
    end subroutine ham_make

    subroutine ovrlp_make(ovrlp,alive,dead) 

        implicit none
        
        real(kind=8),dimension(:,:),intent(inout)::ovrlp
        real(kind=8),dimension(:,:),intent(in)::alive,dead
        real(kind=8),allocatable,dimension(:)::a1,d1,s1
        integer::j,k,ierr
    
        if (errorflag .ne. 0) return 
    
        ierr=0
        allocate(a1(norb),stat=ierr)
        allocate(d1(norb),stat=ierr)
        allocate(d1(norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation. ierr had value ,", ierr 
            errorflag=1
            return
        end if
    
        do j=1,ndet
            do k=j,ndet
                call arr_mult(alive(:,j),alive(:,k),a1)
                call arr_mult(dead(:,j),dead(:,k),d1)
                call arr_sum(a1,d1,s1)
                ovrlp(j,k)=arr_mult_self(s1); ovrlp(k,j)=ovrlp(j,k)
                
            end do
        end do 
    
    end subroutine ovrlp_make

    subroutine hamgen(ham,zstore,elecs,size,verb)

        implicit none 

        type(hamiltonian), intent(inout)::ham 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::size,verb

        real(kind=8),allocatable,dimension(:,:)::ovrlp,temp_ham
        real(kind=8),allocatable,dimension(:,:)::alive,dead
        integer, allocatable,dimension(:)::IPIV1
        real(kind=8),allocatable,dimension(:)::WORK1
        
        integer::j,ierr


        if (errorflag .ne. 0) return

        allocate(alive(size,norb),dead(size,norb),stat=ierr)
        if(ierr==0) allocate(ovrlp(size,size),temp_ham(size,size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in dead and alive vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        do j=1,size
            alive(:,j)=zstore(j)%sin
            dead(:,j)=zstore(j)%cos
        end do

        call ovrlp_make(ovrlp,alive,dead) 
        call ham_make(temp_ham,ovrlp,alive,dead,elecs%h1ei,elecs%h2ei,elecs%hnuc)
        ham%ovrlp=ovrlp
        ham%hjk=temp_ham
        ham%inv=ovrlp
        allocate(WORK1(size),IPIV1(size))
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV or WORK1 vector allocation . ierr had value ", ierr
            errorflag=1
        end if 

        Call dgetrf(size, size, ham%inv, size, IPIV1, ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if
        if (ierr==0) call dgetri(size,ham%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if

        deallocate(WORK1,IPIV1)

        call DGEMM("N","N",size,size,size,1.d0,ham%inv,size,ham%hjk,size,0.d0,ham%kinvh,size)

        deallocate(alive,dead,ovrlp,temp_ham,stat=ierr)

    end subroutine hamgen



    

END MODULE ham