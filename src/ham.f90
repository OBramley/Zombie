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

        Call dgetrf(size, size, haml%inv%x, size, IPIV1, ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF ",ierr
        end if
        if (ierr==0) call dgetri(size,haml%inv%x,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF ",ierr
        end if

        deallocate(WORK1,IPIV1)

        call DGEMM("N","N",size,size,size,1.d0,haml%inv%x,size,haml%hjk%x,size,0.d0,haml%kinvh%x,size)

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
        type(dual2)::htot,ovlptot,hamtot
    
        if (errorflag .ne. 0) return 
        ierr=0

        !$omp parallel do &
        !$omp & private(j,k,h1etot,h2etot,ovlptot,hamtot,z1d) &
        !$omp & shared(elecs,zstore,haml) 
        do j=1,size
            ovlptot=0.0d0; hamtot=0.0d0
            z1d =typ2_2_typ1(zstore(j)%val)
            ovlptot=1.0d0
            htot=haml_vals(z1d,z1d,elecs)
            hamtot=htot+(ovlptot*elecs%hnuc)

            call val_filler(haml%hjk,haml%ovrlp,haml%diff_hjk,haml%diff_ovrlp,hamtot,ovlptot,j,j)
            do k=j+1,size
                ovlptot=0.0d0; hamtot=0.0d0
                ovlptot=overlap_1(z1d,zstore(k)%val)
                htot=haml_vals(z1d,zstore(k)%val,elecs)
                hamtot=htot+(ovlptot*elecs%hnuc)
                call val_filler(haml%hjk,haml%ovrlp,haml%diff_hjk,haml%diff_ovrlp,hamtot,ovlptot,j,k)
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

    subroutine haml_ovrlp_column(haml,z1,zstore,size,elecs,row)

        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        type(zombiest),intent(in)::z1
        type(elecintrgl),intent(in)::elecs
        type(hamiltonian),intent(inout)::haml
        integer,intent(in)::row,size
        type(dual2),dimension(0:2*norb)::z1d
        type(dual2)::htot,ovlptot,hamtot
        integer::j

        if (errorflag .ne. 0) return
        z1d = typ2_2_typ1(z1%val)
        !$omp parallel do &
        !$omp & private(j,h1etot,h2etot,ovlptot,hamtot) &
        !$omp & shared(elecs,zstore,haml,z1d) 
        do j=1,size
            ovlptot=0.0d0; hamtot=0.0d0
            if (j.ne.row) then
                ovlptot=overlap_1(z1d,zstore(j)%val)
                htot=haml_vals(z1d,zstore(j)%val,elecs)
                hamtot=htot+(ovlptot*elecs%hnuc)
                call val_filler(haml%hjk,haml%ovrlp,haml%diff_hjk,haml%diff_ovrlp,hamtot,ovlptot,row,j)
            else 
                ovlptot=1.0d0
                htot=haml_vals(z1d,z1d,elecs)
                hamtot=htot+(ovlptot*elecs%hnuc)
                call val_filler(haml%hjk,haml%ovrlp,haml%diff_hjk,haml%diff_ovrlp,hamtot,ovlptot,row,j)
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
                ov=ov*((z1d(k)*z2d(elecs%ali_dead(k,j))*elecs%negs(k,j))+&
                (z1d((k+norb))*z2d((elecs%ali_dead(k+norb,j)))*elecs%negs(k+norb,j))) 
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
            ovrlp_tot=ovrlp_tot*((z1d(j)*z2d(j))+(z1d(j+norb)*z2d(norb+j)))
        end do
   
      
        

        return 
    end function overlap_1

    !##############################################################################################################################
    
    subroutine val_filler(hjk,ovrlp,diff_hjk,diff_ovrlp,hamtot,ovlptot,j,k)

        implicit none 
        type(dual), dimension(:,:),intent(inout)::hjk
        type(dual), dimension(:,:),intent(inout)::ovrlp
        real(wp), dimension(:,:,:),intent(inout)::diff_hjk
        real(wp), dimension(:,:,:),intent(inout)::diff_ovrlp
        type(dual2),intent(in)::hamtot,ovlptot
        integer,intent(in)::j,k

        if (errorflag .ne. 0) return

        hjk(j,k)%x=hamtot%x; hjk(j,k)%dx=hamtot%dx(1:norb)
        hjk(k,j)=hjk(j,k)
        
        ovrlp(j,k)=ovlptot%x; ovrlp(j,k)%dx=ovlptot%dx(1:norb)
        ovrlp(k,j)=ovrlp(j,k)
        
        diff_hjk(j,k,:)=hamtot%dx(1:norb)
        diff_ovrlp(k,j,:)=ovlptot%dx(1+norb:2*norb)
        if(j.ne.k)then 
            diff_hjk(k,j,:)=hamtot%dx(1+norb:2*norb)
            diff_ovrlp(j,k,:)=ovlptot%dx(1:norb)
        end if
        
       
        return
        
    end subroutine val_filler

END MODULE ham
