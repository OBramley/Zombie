MODULE ham 

    use mod_types
    use globvars
    use alarrays
    use omp_lib
    use dnad
    contains

    !Level 0 Hamiltonian Routine
    ! Subroutine that controls and calcualtes all of the hamiltonian variables 
    subroutine hamgen(haml,zstore,elecs,size,verb)

        implicit none 

        type(hamiltonian), intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::size,verb
        integer, allocatable,dimension(:)::IPIV1
        real(kind=8),allocatable,dimension(:)::WORK1
        integer::ierr

      
        

        if (errorflag .ne. 0) return
        ierr=0
        ! integer(kind=8)::beginning,rate,end
        ! call system_clock(beginning, rate)
        ! call system_clock(end)
        ! print *, "elapsed time: ", real(end - beginning) / real(rate)

        call haml_ovrlp_comb(haml,zstore,elecs,size,verb)
        
        haml%inv=haml%ovrlp
        allocate(WORK1(size),IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV or WORK1 vector allocation . ierr had value ", ierr
            errorflag=1
        end if 

        Call dgetrf(size, size, haml%inv, size, IPIV1, ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF ",ierr
        end if
        if (ierr==0) call dgetri(size,haml%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRI ",ierr
        end if

        deallocate(WORK1,IPIV1)

        call DGEMM("N","N",size,size,size,1.d0,haml%inv,size,haml%hjk,size,0.d0,haml%kinvh,size)

        return 

    end subroutine hamgen

    !##############################################################################################################################
    
    
    !Level 1 Routines to make the Hamiltonian and Overlap matrices

    subroutine haml_ovrlp_comb(haml,zstore,elecs,size,verb)
        implicit none

        type(hamiltonian), intent(inout)::haml 
        ! real(kind=8),dimension(:,:),intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(dual2),dimension(0:2*norb)::z1d
        integer,intent(in)::verb,size
        integer::j,k,ierr
        type(dual2)::ovlptot,hamtot
    
        if (errorflag .ne. 0) return 
        ierr=0

        !$omp parallel do simd&
        !$omp & private(j,k,ovlptot,hamtot,z1d) &
        !$omp & shared(elecs,zstore,haml) 
        do j=1,size
            z1d =typ2_2_typ1(zstore(j)%val)

            hamtot=haml_vals_2(z1d,z1d,elecs)+(elecs%hnuc)
            haml%hjk(j,j)=hamtot%x
            haml%diff_hjk(:,j,j)=hamtot%dx(1:norb)

            haml%ovrlp(j,j)=1.0d0
            haml%diff_ovrlp(:,j,j)=0.0d0 
           
            do k=j+1,size
                ovlptot=overlap_2(z1d,zstore(k)%val)
                haml%ovrlp(j,k)=ovlptot%x; haml%ovrlp(k,j)=haml%ovrlp(j,k)
                haml%diff_ovrlp(:,k,j)=ovlptot%dx(1:norb)
                haml%diff_ovrlp(:,j,k)=ovlptot%dx(1+norb:2*norb)

                hamtot=haml_vals_2(z1d,zstore(k)%val,elecs)+(ovlptot*elecs%hnuc)
                haml%hjk(j,k)=hamtot%x; haml%hjk(k,j)=haml%hjk(j,k)
                haml%diff_hjk(:,k,j)=hamtot%dx(1:norb)
                haml%diff_hjk(:,j,k)=hamtot%dx(1+norb:2*norb)
            end do 
            if(verb.eq.1)then
                write(6,"(a,i0,a)") "hamliltonian column ",j, " completed"
            end if 
        end do
        !$omp end parallel do simd

    end subroutine haml_ovrlp_comb

    

    !##############################################################################################################################

    !Level 2 routines to make an Overlap and Hamiltonian matrix column

    ! Calcualates a column of a hamiltonian Start specifies the row the column
    ! is started to be calcualted 

    subroutine haml_ovrlp_column(temp,zstore,size,elecs,row)

        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        type(grad_do),intent(inout)::temp
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::row,size
        type(dual2),dimension(0:2*norb)::z1d
        type(dual2)::ovlptot,hamtot
        integer::j

        if (errorflag .ne. 0) return
        z1d = typ2_2_typ1(temp%zom%val)
        !$omp parallel do  default(none) &
        !$omp & private(j,ovlptot,hamtot) &
        !$omp & shared(elecs,zstore,temp,z1d,row,norb) 
        do j=1,size
            if (j.ne.row) then
                ovlptot=overlap_2(z1d,zstore(j)%val)
                temp%ovrlp(row,j)=ovlptot%x; temp%ovrlp(row,j)%dx=ovlptot%dx(1:norb)
                temp%ovrlp(j,row)=temp%ovrlp(row,j)

                temp%diff_ovrlp_1(:,j)=ovlptot%dx(1:norb)
                temp%diff_ovrlp_2(:,j)=ovlptot%dx(1+norb:2*norb)

             
                hamtot=haml_vals_2(z1d,zstore(j)%val,elecs)+(ovlptot*elecs%hnuc)
                temp%hjk(row,j)%x=hamtot%x; temp%hjk(row,j)%dx=hamtot%dx(1:norb)
                temp%hjk(j,row)=temp%hjk(row,j)

                temp%diff_hjk_1(:,j)=hamtot%dx(1:norb)
                temp%diff_hjk_2(:,j)=hamtot%dx(1+norb:2*norb)
            else 
                temp%ovrlp(row,row)=1.0d0
                temp%ovrlp(row,row)%dx=0.0d0
                temp%diff_ovrlp_1(:,j)=0.0d0
             
                hamtot=haml_vals_2(z1d,z1d,elecs)+(elecs%hnuc)
                temp%hjk(row,row)%x=hamtot%x; temp%hjk(row,row)%dx=hamtot%dx(1:norb)

                temp%diff_hjk_1(:,j)=hamtot%dx(1:norb)
            end if 
        end do
        !$omp end parallel do 

        return

    end subroutine haml_ovrlp_column

    !##############################################################################################################################

    !Level 3 routine to make individual value in overlap and hamiltonian matrix

    ! calculates indvidual hamiltonian elements taking in two Zombie states and a set of 
    ! creation and annihilation operations
    function haml_vals(z1d,z2d,elecs) result(ham_tot)

        implicit none 
        type(dual2),dimension(0:),intent(in)::z1d,z2d
        type(elecintrgl),intent(in)::elecs
        type(dual2)::ham_tot
        type(dual2)::ov
        integer::j,k
        
        if (errorflag .ne. 0) return

       
        ham_tot=0.0d0
        !$omp simd
        do j=1,elecs%num
            ov=elecs%integrals(j)
            do k=1, norb
                ov=ov*((z1d(k)*z2d(elecs%alive(k,j))*elecs%neg_a(k,j))+(z1d((k+norb))*z2d((elecs%dead(k,j)))*elecs%neg_d(k,j)))
            end do
            ham_tot=ham_tot+ov
        end do
        !$omp end simd
        return 
      
    end function haml_vals

    ! calculates individual overlaps where no creation and annihilation operations are needed
    function overlap_1(z1d,z2d) result(ovrlp_tot)

        implicit none
        type(dual2),dimension(0:)::z1d,z2d
        type(dual2)::ovrlp_tot
        integer::j

        if (errorflag .ne. 0) return

        ovrlp_tot=1.0d0
      
        !$omp simd
        do j=1,norb
            ovrlp_tot=ovrlp_tot*((z1d(j)*z2d(j))+(z1d(j+norb)*z2d(norb+j)))
        end do
        !$omp end simd
      
        return 
    end function overlap_1

    !##############################################################################################################################
    
    
   
    function overlap_2(z1d,z2d) result(ovrlp_tot)

        implicit none
        type(dual2),dimension(0:)::z1d,z2d
        type(dual2)::ovrlp_tot
        integer::j,k

        if (errorflag .ne. 0) return

       
        ovrlp_tot=1.0d0
       
        do j=1,norb
            do k=1,j
                ovrlp_tot%dx(k)=((z1d(j)%x*z2d(j)%dx(k)+z1d(j)%dx(k)*z2d(j)%x)+&
                                (z1d(j+norb)%x*z2d(norb+j)%dx(k)+z1d(j+norb)%dx(k)*z2d(norb+j)%x))*ovrlp_tot%x + &
                                (z1d(j)%x*z2d(j)%x+z1d(j+norb)%x*z2d(norb+j)%x)*ovrlp_tot%dx(k)
                ovrlp_tot%dx(k+norb)=((z1d(j)%x*z2d(j)%dx(k+norb)+z1d(j)%dx(k+norb)*z2d(j)%x)+&
                                (z1d(j+norb)%x*z2d(norb+j)%dx(k+norb)+z1d(j+norb)%dx(k+norb)*z2d(norb+j)%x))*ovrlp_tot%x + &
                                (z1d(j)%x*z2d(j)%x+z1d(j+norb)%x*z2d(norb+j)%x)*ovrlp_tot%dx(k+norb)
            end do
            ovrlp_tot%x=ovrlp_tot%x*(z1d(j)%x*z2d(j)%x+z1d(j+norb)%x*z2d(norb+j)%x)
        end do
     
      
        return 
    end function overlap_2

    function overlap_3(val,z1d,z2d,alive,dead,aneg,dneg,j) result(ovrlp_tot)

        implicit none
        type(dual2),dimension(0:)::z1d,z2d
        type(dual2)::val
        integer(int16)::alive,dead
        integer(int8)::aneg,dneg
        type(dual2)::ovrlp_tot
        integer::j
        integer::k

        if (errorflag .ne. 0) return

        do k=1,norb
            ovrlp_tot%dx(k)=((z1d(j)%x*aneg*z2d(alive)%dx(k)+z1d(j)%dx(k)*aneg*z2d(alive)%x)+&
                            (z1d(j+norb)%x*dneg*z2d(dead)%dx(k)+z1d(j+norb)%dx(k)*dneg*z2d(dead)%x))*val%x + &
                            (z1d(j)%x*aneg*z2d(alive)%x+z1d(j+norb)%x*dneg*z2d(dead)%x)*val%dx(k)
            ovrlp_tot%dx(k+norb)=((z1d(j)%x*aneg*z2d(alive)%dx(k+norb)+z1d(j)%dx(+norb)*aneg*z2d(alive)%x)+&
                            (z1d(j+norb)%x*dneg*z2d(dead)%dx(k+norb)+z1d(j+norb)%dx(k+norb)*dneg*z2d(dead)%x))*val%x + &
                            (z1d(j)%x*aneg*z2d(alive)%x+z1d(j+norb)%x*dneg*z2d(dead)%x)*val%dx(k+norb)
        end do
        ovrlp_tot%x=val%x*(z1d(j)%x*aneg*z2d(alive)%x+z1d(j+norb)%x*dneg*z2d(dead)%x)
        
        return 
    end function overlap_3


    function haml_vals_2(z1d,z2d,elecs) result(ham_tot)

        implicit none 
        type(dual2),dimension(0:),intent(in)::z1d,z2d
        type(elecintrgl),intent(in)::elecs
        type(dual2)::ham_tot
        type(dual2)::ov
        integer::j,k
        
        if (errorflag .ne. 0) return

       
        ham_tot=0.0d0
       
        do j=1,elecs%num
            ov=elecs%integrals(j)
            do k=1, norb
                ov=overlap_3(ov,z1d,z2d,elecs%alive(k,j),elecs%dead(k,j),elecs%neg_a(k,j),elecs%neg_d(k,j),k)
                
            end do
            ham_tot=ham_tot+ov
            ! ham_tot%x=ham_tot%x+ov%x
            ! do l=1,2*norb
            !     ham_tot%dx(l)=ham_tot%dx(l)+ov%dx(l)
            ! end do
        end do
       
        return 
      
    end function haml_vals_2


!     attributes(global) subroutine overlap_2_kernel(z1d, z2d, ovrlp_tot)
!         type(dual2), device :: z1d(:), z2d(:)
!         type(dual2), device :: ovrlp_tot
!         integer :: j, k

!         integer :: tid

!         tid = threadIdx%x + (blockIdx%x - 1) * blockDim%x

!         if (tid <= norb) then
!         type(dual2) temp

!         temp%x = one
!         temp%dx = 0.0d0

!         do j = 1, norb
!             do k = 1, 2 * norb
!             temp%dx(k) = ((z1d(j)%x * z2d(j)%dx(k) + z1d(j)%dx(k) * z2d(j)%x) + &
!                             (z1d(j+norb)%x * z2d(norb+j)%dx(k) + z1d(j+norb)%dx(k) * z2d(norb+j)%x)) * temp%x + &
!                         (z1d(j)%x * z2d(j)%x + z1d(j+norb)%x * z2d(norb+j)%x) * temp%dx(k)
!             end do
!             temp%x = temp%x * (z1d(j)%x * z2d(j)%x + z1d(j+norb)%x * z2d(norb+j)%x)
!         end do

!         ovrlp_tot = temp
!         end if

!     end subroutine overlap_2_kernel

!   subroutine overlap_2(z1d, z2d, ovrlp_tot)
!     type(dual2), dimension(0:norb-1), intent(in) :: z1d, z2d
!     type(dual2), intent(out) :: ovrlp_tot
!     type(dual2), device :: d_z1d(norb+2*norb), d_z2d(norb+2*norb), d_ovrlp_tot
!     integer :: numBlocks, blockSize

!     ! Copy z1d and z2d from host to device
!     d_z1d(1:norb+2*norb) = z1d
!     d_z2d(1:norb+2*norb) = z2d

!     numBlocks = (norb + 255) / 256  ! Adjust block size as needed
!     blockSize = 256                 ! Adjust block size as needed

!     d_ovrlp_tot%x = one
!     d_ovrlp_tot%dx = 0.0d0

!     call overlap_2_kernel<<<numBlocks, blockSize>>>(d_z1d, d_z2d, d_ovrlp_tot)

!     ! Copy ovrlp_tot from device to host
!     ovrlp_tot = d_ovrlp_tot

!   end subroutine overlap_2



END MODULE ham
