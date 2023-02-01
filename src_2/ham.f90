MODULE ham 

    use globvars
    contains

    ! calculates indvidual hamiltonian elements taking in two Zombie states and a set of 
    ! creation and annihilation operations
    real(kind=8) function ham_vals(z1d,z2d,ops,el,len)

        implicit none 
        real(kind=8),dimension(0:),intent(in)::z1d,z2d
        real(kind=8),dimension(:),intent(in)::el
        integer,intent(in)::len
        type(oprts),intent(in)::ops
        real(kind=8)::ov
        integer::j,k

        
        ham_vals=0.0
        !!$omp parallel private(j,k,ov) shared(ops,z1d,z2d)
        !!$omp do simd  reduction(+:ham_vals) 
        do j=1,len
            ov=1.0
            !!$omp do simd reduction(*:ov)
            do k=1, norb
                ov=ov*((z1d(k)*z2d(ops%alive(k,j))*ops%neg_alive(k,j))+(z1d(k+norb)*z2d(ops%dead(k,j))*ops%neg_dead(k,j))) 
            end do
            !!$omp end do simd
            ham_vals=ham_vals+(ov*el(j))
        end do
        !!$omp end do simd
        !!$omp end parallel 
        
        return 
      
    end function ham_vals

    ! Calcualates a column of a hamiltonian Start specifies the row the column
    ! is started to be calcualted 
    function ham_column(hcol,z1d,zstore,an_cr,an2_cr2,elecs,start)

        implicit none
        real(kind=8),dimension(:),intent(inout)::hcol 
        type(zombiest),dimension(:),intent(in)::zstore
        real(kind=8),dimension(0:),intent(in)::z1d
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::start
        real(kind=8),dimension(ndet-(start-1))::ham_column
        real(kind=8)::h1etot,h2etot
        integer::j
        
        ham_column=hcol
        !!$omp parallel 
        !!$omp single
        do j=1,(ndet-(start-1))
            !!$omp task firstprivate(h1etot,j) shared(ham_column,zstore,hcol,an_cr,an2_cr2,elecs,z1d,start)
            h1etot = ham_vals(z1d,zstore(start+j-1)%val,an_cr,elecs%h1ei,elecs%h1_num)
            !!$omp atomic
            ham_column(j)=ham_column(j)+h1etot
            !!$omp end atomic
            !!$omp end task
            !!$omp task firstprivate(h2etot,j) shared(ham_column,zstore,hcol,an_cr,an2_cr2,elecs,z1d)
            h2etot = ham_vals(z1d,zstore(start+j-1)%val,an2_cr2,elecs%h2ei,elecs%h2_num)
            !!$omp atomic
            ham_column(j)=ham_column(j)+(0.5*h2etot)
            !!$omp end atomic
            !!$omp end task
            
        end do 
        !!$omp end single
        !!$omp end parallel  
        return

    end function ham_column

    ! Hamiltonian calcualtion - calcualtes the whole hamiltonian 
    subroutine ham_make(ham,zstore,elecs,an_cr,an2_cr2,verb) 

        implicit none
        
        real(kind=8),dimension(:,:),intent(inout)::ham 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::verb
        integer::j,ierr
    
        if (errorflag .ne. 0) return 
        ierr=0
    
        !!$omp parallel do &
        !!$omp & private(j) &
        !!$omp & shared(elecs,zstore,an_cr,an2_cr2,ham) 
        do j=1,ndet
            ham(j:,j) = ham_column(ham(j:,j),zstore(j)%val,zstore,an_cr,an2_cr2,elecs,j)
            if(verb.eq.1)then
                write(6,"(a,i0,a)") "Hamiltonian column ",j, " completed"
            end if 
        end do
        !!$omp end parallel do

        do j=1,ndet
            ham(j,:)=ham(:,j)
        end do
        
    end subroutine ham_make

    ! calculates individual overlaps where no creation and annihilation operations are needed
    real(kind=8) function overlap_1(z1d,z2d)

        implicit none
        real(kind=8),dimension(0:)::z1d,z2d
        integer::j
    
    
        overlap_1=1.0
        !!$omp parallel do simd reduction(*:overlap_1)
        do j=1,norb
            overlap_1=overlap_1*((z1d(j)*z2d(j))+(z1d(j+norb)*z2d(norb+j)))
        end do
        !!$omp end parallel do simd
        

        return 
    end function overlap_1

    ! function to calcualte an entire column of the overlap 
    function ovrlp_column(z1d,zstore,row)

        implicit none
        
        type(zombiest),dimension(:),intent(in)::zstore
        real(kind=8),dimension(0:),intent(in)::z1d
        integer,intent(in)::row
        real(kind=8),dimension(ndet)::ovrlp_column
        integer::j
        
        ovrlp_column=0.0
        !!$omp parallel do &
        !!$omp & shared(z1d,zstore,ovrlp_column) &
        !!$omp & private(j)
        do j=1,ndet
            if(j.ne.row)then 
                ovrlp_column(j)=overlap_1(z1d,zstore(j)%val)
            else
                ovrlp_column(j)=1.0
            end if 
        end do 
        !!$omp end parallel do
        return

    end function ovrlp_column

    !subroutine calcualates whole overlap matrix
    subroutine ovrlp_make(ovrlp,zstore)

        implicit none
        
        real(kind=8),dimension(:,:),intent(inout)::ovrlp
        type(zombiest),dimension(:),intent(in)::zstore
        integer::j,k
    
        if (errorflag .ne. 0) return

        !!$omp parallel do &
        !!$omp & private(j,k) &
        !!$omp & shared(zstore,ovrlp)
        do j=1,ndet
            do k=j,ndet
                if(k.ne.j)then 
                    ovrlp(j,k)=overlap_1(zstore(j)%val,zstore(k)%val); ovrlp(k,j)=ovrlp(j,k)
                else
                    ovrlp(j,k)=1.0
                end if 
            end do
        end do 
        !!$omp end parallel do 

    end subroutine ovrlp_make

    ! Subroutine that controls and calcualtes all of the hamiltonian variables 
    subroutine hamgen(ham,zstore,elecs,size,an_cr,an2_cr2,verb)

        implicit none 

        type(hamiltonian), intent(inout)::ham 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::size,verb
        ! real(kind=8),allocatable,dimension(:,:)::alive,dead
        integer, allocatable,dimension(:)::IPIV1
        real(kind=8),allocatable,dimension(:)::WORK1
        
        integer::ierr


        if (errorflag .ne. 0) return
        ierr=0
        ! allocate(alive(size,norb),dead(size,norb),stat=ierr)
        ! if (ierr/=0) then
        !     write(0,"(a,i0)") "Error in dead and alive vector allocation . ierr had value ", ierr
        !     errorflag=1
        !     return
        ! end if 
        ! !$omp parallel do 
        ! do j=1,size
        !     alive(:,j)=zstore(j)%sin
        !     dead(:,j)=zstore(j)%cos
        ! end do
        ! !$omp end parallel do

        call ovrlp_make(ham%ovrlp,zstore) !alive,dead)
       
        ham%hjk=ham%ovrlp*elecs%hnuc 
        call ham_make(ham%hjk,zstore,elecs,an_cr,an2_cr2,verb)
        ham%inv=ham%ovrlp
        allocate(WORK1(size),IPIV1(size),stat=ierr)
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

        ! deallocate(alive,dead,stat=ierr)

    end subroutine hamgen


END MODULE ham
