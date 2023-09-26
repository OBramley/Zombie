MODULE ham 

    use mod_types
    use globvars
    use alarrays
    use omp_lib
  
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
        real(wp),allocatable,dimension(:)::WORK1
        integer::ierr=0

        if (errorflag .ne. 0) return
       
        call haml_ovrlp_comb(haml,zstore,elecs,size,verb)
        
        haml%inv=haml%ovrlp
        allocate(WORK1(size),IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in IPIV or WORK1 vector allocation . ierr had value ", ierr
            errorflag=1
        end if 

        Call dgetrf(size, size, haml%inv, size, IPIV1, ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)")"Error in DGETRF ",ierr
        end if
        if (ierr==0) call dgetri(size,haml%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)")"Error in DGETRI ",ierr
        end if

        deallocate(WORK1,IPIV1)

        call DGEMM("N","N",size,size,size,1.d0,haml%inv,size,haml%hjk,size,0.d0,haml%kinvh,size)

        return 

    end subroutine hamgen

    !##############################################################################################################################
    
   
    !##############################################################################################################################

    
    !Level 1 Routines to make the Hamiltonian and Overlap matrices

    subroutine haml_ovrlp_comb(haml,zstore,elecs,size,verb)
        implicit none

        type(hamiltonian), intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::verb,size
        integer::j,k
        real(wp)::ovlptot,hamtot

        if (errorflag .ne. 0) return 
       

        !$omp parallel do &
        !$omp & private(j,k,ovlptot,hamtot) &
        !$omp & shared(elecs,zstore,haml) 
        do j=1,size
           
            hamtot=haml_vals(zstore(j)%val%x,zstore(j)%val%x,elecs)+(elecs%hnuc)
            haml%hjk(j,j)=hamtot
            haml%ovrlp(j,j)=1.0d0
          
            do k=j+1,size
                ovlptot=overlap_1(zstore(j)%val%x,zstore(k)%val%x)
                haml%ovrlp(j,k)=ovlptot; haml%ovrlp(k,j)=haml%ovrlp(j,k)
                hamtot=haml_vals(zstore(j)%val%x,zstore(k)%val%x,elecs)+(ovlptot*elecs%hnuc)
                haml%hjk(j,k)=hamtot; haml%hjk(k,j)=haml%hjk(j,k)
              
            end do 
            if(verb.eq.1)then
                write(stdout,"(a,i0,a)") "hamliltonian column ",j, " completed"
            end if 
        end do
        !$omp end parallel do

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
       
        real(wp)::ovlptot,hamtot
        integer::j

        if (errorflag .ne. 0) return
      
        !$omp parallel do  default(none) &
        !$omp & private(j,ovlptot,hamtot) &
        !$omp & shared(elecs,zstore,temp,row,norb)
        !!$acc data present(zstore,elecs,norb,ndet) copy(temp) create(ovlptot,hamtot)
        !!$acc parallel loop
        do j=1,ndet
            if (j.ne.row) then
                ovlptot=overlap_1(temp%zom%val%x,zstore(j)%val%x)
                temp%ovrlp(j,row)=ovlptot
                ! temp%ovrlp(j,row)=temp%ovrlp(row,j)

                hamtot=haml_vals(temp%zom%val%x,zstore(j)%val%x,elecs)+(ovlptot*elecs%hnuc)
                temp%hjk(j,row)=hamtot
                ! temp%hjk(j,row)=temp%hjk(row,j)

            else 
                temp%ovrlp(row,row)=1.0d0
                hamtot=haml_vals(temp%zom%val%x,temp%zom%val%x,elecs)+(elecs%hnuc)
                temp%hjk(row,row)=hamtot
            end if 
        end do
        !$omp end parallel do 
        !!$acc end parallel loop
        !!$acc end data
        temp%hjk(row,:)=temp%hjk(:,row) 
        temp%ovrlp(row,:)=temp%ovrlp(:,row)
        
        return

    end subroutine haml_ovrlp_column

    !##############################################################################################################################

    !Level 3 routine to make individual value in overlap and hamiltonian matrix

   

    !##############################################################################################################################
    
    function overlap_1(z1d,z2d) result(ovrlp_tot)

        implicit none
        real(wp),dimension(0:)::z1d,z2d
        real(wp)::ovrlp_tot
        integer::j

        ovrlp_tot=1.0d0
        !$omp simd
        do j=1,norb
            ovrlp_tot=ovrlp_tot*(z1d(j)*z2d(j)+z1d(j+norb)*z2d(norb+j))
        end do
        !$omp end simd
      
        return 
    end function overlap_1

    function haml_vals(z1d,z2d,elecs) result(ham_tot)

        implicit none 
        real(wp),dimension(0:),intent(in)::z1d,z2d
        type(elecintrgl),intent(in)::elecs
        real(wp)::ham_tot
        integer::j
        
        ham_tot=0.0d0
        !$omp simd
        do j=1,elecs%num
            ham_tot=ham_tot+haml_singular(z1d,z2d,elecs%alive(:,j),elecs%dead(:,j),&
            elecs%neg_a(:,j),elecs%neg_d(:,j),elecs%integrals(j))
        end do
        !$omp end simd
        return 
      
    end function haml_vals

    real(wp) function haml_singular(z1d,z2d,alive,dead,neg_a,neg_d,ov)

        implicit none 
        real(wp),dimension(0:),intent(in)::z1d,z2d
        integer(int16),dimension(:),intent(in)::alive,dead
        integer(int8),dimension(:),intent(in)::neg_a,neg_d
        real(wp),value::ov
        integer::k
        
        !$omp simd
        do k=1, norb
            ov=ov*(z1d(k)*neg_a(k)*z2d(alive(k))+z1d(k+norb)*neg_d(k)*z2d(dead(k)))
        end do
        !$omp end simd
        haml_singular=ov
        return 
      
    end function haml_singular


END MODULE ham
