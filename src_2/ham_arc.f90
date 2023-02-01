MODULE ham 

    use globvars
    use ham_2
    contains


    subroutine one_electron_vals(alive1,dead1,alive2,dead2,el,h1etot,n) 
        !$omp declare simd(one_electron_vals) uniform(n)
        implicit none
  
        real(kind=8), dimension(:),intent(in)::alive1,dead1,alive2,dead2,el
        real(kind=8),intent(out)::h1etot
        integer::j,n
 
        h1etot=0.0
        !$omp parallel do simd private(j) shared(alive1,alive2,dead1,dead2,el) reduction(+:h1etot)
        do j=1,n
            h1etot=h1etot+h1e(alive1,dead1,alive2,dead2,el(j),j)
        end do  
        !$omp end parallel do simd
        return

    end subroutine one_electron_vals

  
    subroutine two_electron_vals(alive1,dead1,alive2,dead2,el,h2etot,n) 
        !$omp declare simd(two_electron_vals) uniform(n)
        implicit none
    
        real(kind=8), dimension(:),intent(in)::alive1,dead1,alive2,dead2,el
        real(kind=8),intent(out)::h2etot
        integer::j,n
    
        h2etot=0.0
        !$omp parallel do simd private(j) shared(alive1,alive2,dead1,dead2,el) reduction(+:h2etot)
        do j=1,n
            h2etot=h2etot+h2e(alive1,dead1,alive2,dead2,el(j),j)
        end do  
        !$omp end parallel do simd
        return
        
    end subroutine two_electron_vals

    subroutine ham_make(ham,ovrlp,alive,dead,el1,el2,hnuc) 

        implicit none
        
        real(kind=8),dimension(:,:),intent(inout)::ham 
        real(kind=8),dimension(:,:),intent(in)::ovrlp
        real(kind=8),dimension(:,:),intent(in)::alive,dead
        real(kind=8),dimension(:),intent(in)::el1,el2
        real(kind=8),intent(in)::hnuc
        real(kind=8)::h1etot,h2etot
        integer::j,k,h1,h2
    
        if (errorflag .ne. 0) return 
        h1=size(el1)
        h2=size(el2)
        !$omp parallel do simd collapse(2) &
        !$omp & private(j,k,h1etot,h2etot) &
        !$omp & shared(alive,dead,el1,el2,ovrlp,hnuc) 
        do j=1,ndet
            do k=j,ndet
                call one_electron_vals(alive(:,j),dead(:,j),alive(:,k),dead(:,k),el1,h1etot,h1)
                call two_electron_vals(alive(:,j),dead(:,j),alive(:,k),dead(:,k),el2,h2etot,h2)
                ham(j,k)=h1etot+(0.5*h2etot)+(hnuc*ovrlp(j,k)); ham(k,j)=ham(j,k)
            end do
        end do 
        !$omp end parallel do simd
    
    end subroutine ham_make

    subroutine ovrlp_make(ovrlp,alive,dead) 

        implicit none
        
        real(kind=8),dimension(:,:),intent(inout)::ovrlp
        real(kind=8),dimension(:,:),intent(in)::alive,dead
        integer::j,k
    
        if (errorflag .ne. 0) return

        !$omp parallel do simd collapse(2) &
        !$omp & private(j,k) &
        !$omp & shared(alive,dead,ovrlp) 
        do j=1,ndet
            do k=j,ndet
                ovrlp(j,k)=product((alive(:,j)*alive(:,k))+(dead(:,j)*dead(:,k))); ovrlp(k,j)=ovrlp(j,k)
            end do
        end do 
        !$omp end parallel do simd

    end subroutine ovrlp_make


    subroutine hamgen(ham,zstore,elecs,size,verb)

        implicit none 

        type(hamiltonian), intent(inout)::ham 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::size,verb
        real(kind=8),allocatable,dimension(:,:)::alive,dead
        integer, allocatable,dimension(:)::IPIV1
        real(kind=8),allocatable,dimension(:)::WORK1
        
        integer::j,ierr


        if (errorflag .ne. 0) return

        allocate(alive(size,norb),dead(size,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in dead and alive vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        !$omp parallel do 
        do j=1,size
            alive(:,j)=zstore(j)%sin
            dead(:,j)=zstore(j)%cos
        end do
        !$omp end parallel do

        call ovrlp_make(ham%ovrlp,alive,dead) 
        call ham_make(ham%hjk,ham%ovrlp,alive,dead,elecs%h1ei,elecs%h2ei,elecs%hnuc)
        ham%inv=ham%ovrlp
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

        deallocate(alive,dead,stat=ierr)

    end subroutine hamgen


END MODULE ham

