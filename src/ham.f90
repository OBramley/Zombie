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
            write(0,"(a,i0)")"Error in DGETRF ",ierr
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

        !$omp parallel do &
        !$omp & private(j,k,ovlptot,hamtot,z1d) &
        !$omp & shared(elecs,zstore,haml) 
        do j=1,size
            z1d =typ2_2_typ1(zstore(j)%val)

            hamtot=haml_vals(z1d,z1d,elecs)+(elecs%hnuc)
            haml%hjk(j,j)=hamtot%x
            haml%diff_hjk(:,j,j)=hamtot%dx(1:norb)

            haml%ovrlp(j,j)=1.0d0
            haml%diff_ovrlp(:,j,j)=0.0d0 
           
            do k=j+1,size
                ovlptot=overlap_1(z1d,zstore(k)%val)
                haml%ovrlp(j,k)=ovlptot%x; haml%ovrlp(k,j)=haml%ovrlp(j,k)
                haml%diff_ovrlp(:,k,j)=ovlptot%dx(1:norb)
                haml%diff_ovrlp(:,j,k)=ovlptot%dx(1+norb:2*norb)

                hamtot=haml_vals(z1d,zstore(k)%val,elecs)+(ovlptot*elecs%hnuc)
                haml%hjk(j,k)=hamtot%x; haml%hjk(k,j)=haml%hjk(j,k)
                haml%diff_hjk(:,k,j)=hamtot%dx(1:norb)
                haml%diff_hjk(:,j,k)=hamtot%dx(1+norb:2*norb)
            end do 
            if(verb.eq.1)then
                write(6,"(a,i0,a)") "hamliltonian column ",j, " completed"
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
        ! type(zombiest),intent(in)::z1
        type(elecintrgl),intent(in)::elecs
        ! type(hamiltonian),intent(inout)::haml
        integer,intent(in)::row,size
        type(dual2),dimension(0:2*norb)::z1d
        type(dual2)::ovlptot,hamtot
        integer::j

        if (errorflag .ne. 0) return
        z1d = typ2_2_typ1(temp%zom%val)
        !$omp parallel do &
        !$omp & private(j,ovlptot,hamtot) &
        !$omp & shared(elecs,zstore,haml,z1d) 
        do j=1,size
            if (j.ne.row) then
                ovlptot=overlap_1(z1d,zstore(j)%val)
                temp%ovrlp(row,j)=ovlptot%x; temp%ovrlp(row,j)%dx=ovlptot%dx(1:norb)
                temp%ovrlp(j,row)=temp%ovrlp(row,j)

                temp%diff_ovrlp_1(:,j)=ovlptot%dx(1:norb)
                temp%diff_ovrlp_2(:,j)=ovlptot%dx(1+norb:2*norb)

             
                hamtot=haml_vals(z1d,zstore(j)%val,elecs)+(ovlptot*elecs%hnuc)
                temp%hjk(row,j)%x=hamtot%x; temp%hjk(row,j)%dx=hamtot%dx(1:norb)
                temp%hjk(j,row)=temp%hjk(row,j)

                temp%diff_hjk_1(:,j)=hamtot%dx(1:norb)
                temp%diff_hjk_2(:,j)=hamtot%dx(1+norb:2*norb)
            else 
                temp%ovrlp(row,row)=1.0d0
                temp%ovrlp(row,row)%dx=0.0d0
                temp%diff_ovrlp_1(:,j)=0.0d0
             
                hamtot=haml_vals(z1d,z1d,elecs)+(elecs%hnuc)
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
     
        do j=1,elecs%num
            ov=elecs%integrals(j)
            do k=1, norb
                ov=ov*((z1d(k)*z2d(elecs%alive(k,j))*elecs%neg_a(k,j))+&
                (z1d((k+norb))*z2d((elecs%dead(k,j)))*elecs%neg_d(k,j))) 
            end do
            ham_tot=ham_tot+ov
        end do

     
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
      
       
        do j=1,norb
            ovrlp=ovrlp*sparse_mult_ovrlp(z1d(j),z1d(j+norb),z2d(j),z2d(j+norb),j)
            ! ovrlp_tot=ovrlp_tot*((z1d(j)*z2d(j))+(z1d(j+norb)*z2d(norb+j)))
        end do
   
      
        

        return 
    end function overlap_1

    !##############################################################################################################################
    
    function sparse_mult_ovrlp(z1d_a,z1d_d,z2d_a,z2d_d,j)result(mult)
        implicit none
        type(dual2),value::z1d_a,z1d_d,z2d_a,z2d_d
        type(dual2)::mult
        integer::j

       
        mult%x          = (z1d_a%x * z2d_a%x) + (z1d_d%x*z2d_d%x)
        mult%dx(j)      = (z1d_a%dx(j) * z2d_a%x) + (z1d_d%dx(j)*z2d_d%x)
        mult%dx(j+norb) = (z1d_a%x  * z2d_a%dx(j+norb)) + (z1d_d%x*z2d_d%dx(j+norb))
       
    end function sparse_mult_ovrlp

END MODULE ham
