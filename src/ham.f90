MODULE ham 

    use mod_types
    use globvars
    use alarrays
    use omp_lib
    use dnad
    contains

    !Level 0 Hamiltonian Routine
    ! Subroutine that controls and calcualtes all of the hamiltonian variables 
    subroutine hamgen(haml,zstore,elecs,size,an_cr,an2_cr2,verb)

        implicit none 

        type(hamiltonian), intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
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

        call haml_ovrlp_comb(haml,zstore,elecs,size,an_cr%ham,an2_cr2%ham,verb)
        
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

    subroutine haml_ovrlp_comb(haml,zstore,elecs,size,an_cr,an2_cr2,verb)
        implicit none

        type(hamiltonian), intent(inout)::haml 
        ! real(kind=8),dimension(:,:),intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts_2),intent(in)::an_cr,an2_cr2
        type(dual2),dimension(0:2*norb)::z1d,z2d
        integer,intent(in)::verb,size
        integer::j,k,ierr
        type(dual2)::h1etot,h2etot,ovlptot,hamtot
    
        if (errorflag .ne. 0) return 
        ierr=0

        !$omp parallel do &
        !$omp & private(j,k,h1etot,h2etot,ovlptot,hamtot,z1d,z2d) &
        !$omp & shared(elecs,zstore,an_cr,an2_cr2,haml) 
        do j=1,size
            h1etot=0.0d0; h2etot=0.0d0; ovlptot=0.0d0; hamtot=0.0d0
            call dual_2_dual2(zstore(j)%val,z1d,1)
            ovlptot=1.0d0
            h1etot = haml_vals(z1d,z1d,an_cr,elecs%h1ei,elecs%h1_num)
            h2etot = haml_vals(z1d,z1d,an2_cr2,elecs%h2ei,elecs%h2_num)
            hamtot=h1etot+(0.5d0*h2etot)+(ovlptot*elecs%hnuc)

            call val_filler(haml%hjk,haml%ovrlp,haml%diff_hjk,haml%diff_ovrlp,hamtot,ovlptot,j,j)
            do k=j+1,size
                call dual_2_dual2(zstore(k)%val,z2d,2)
                h1etot=0.0d0; h2etot=0.0d0; ovlptot=0.0d0; hamtot=0.0d0
                ovlptot=overlap_1(z1d,z2d)
                h1etot = haml_vals(z1d,z2d,an_cr,elecs%h1ei,elecs%h1_num)
                h2etot = haml_vals(z1d,z2d,an2_cr2,elecs%h2ei,elecs%h2_num)
                hamtot=h1etot+(0.5d0*h2etot)+(ovlptot*elecs%hnuc)
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

    subroutine haml_ovrlp_column(haml,z1d,zstore,size,an_cr,an2_cr2,elecs,row)

        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        type(dual2),dimension(0:),intent(in)::z1d
        type(elecintrgl),intent(in)::elecs
        type(oprts_2),intent(in)::an_cr,an2_cr2
        type(hamiltonian),intent(inout)::haml
        type(dual2),dimension(0:2*norb)::z2d
        integer,intent(in)::row,size
        type(dual2)::h1etot,h2etot,ovlptot,hamtot
        integer::j

        if (errorflag .ne. 0) return

        !$omp parallel do &
        !$omp & private(j,h1etot,h2etot,ovlptot,hamtot,z2d) &
        !$omp & shared(elecs,zstore,an_cr,an2_cr2,haml,z1d) 
        do j=1,size
            h1etot=0.0d0; h2etot=0.0d0; ovlptot=0.0d0; hamtot=0.0d0
            if (j.ne.row) then
                call dual_2_dual2(zstore(j)%val,z2d,2)
                ovlptot=overlap_1(z1d,z2d)
                h1etot = haml_vals(z1d,z2d,an_cr,elecs%h1ei,elecs%h1_num)
                h2etot = haml_vals(z1d,z2d,an2_cr2,elecs%h2ei,elecs%h2_num)
                hamtot=h1etot+(0.5d0*h2etot)+(ovlptot*elecs%hnuc)
                call val_filler(haml%hjk,haml%ovrlp,haml%diff_hjk,haml%diff_ovrlp,hamtot,ovlptot,row,j)
            else 
                ovlptot=1.0d0
                h1etot = haml_vals(z1d,z1d,an_cr,elecs%h1ei,elecs%h1_num)
                h2etot = haml_vals(z1d,z1d,an2_cr2,elecs%h2ei,elecs%h2_num)
                hamtot=h1etot+(0.5d0*h2etot)+(ovlptot*elecs%hnuc)
                call val_filler(haml%hjk,haml%ovrlp,haml%diff_hjk,haml%diff_ovrlp,hamtot,ovlptot,row,j)
            end if 
        end do
        !$omp end parallel do 

        return

    end subroutine haml_ovrlp_column

    function haml_column(z1d,zstore,size,an_cr,an2_cr2,elecs,row) result(ham_tot)


        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        type(dual2),dimension(0:),intent(in)::z1d
        type(elecintrgl),intent(in)::elecs
        type(oprts_2),intent(in)::an_cr,an2_cr2
        integer,intent(in)::row,size
        type(dual2),dimension(size)::ham_tot
        type(dual2),dimension(0:2*norb)::z2d
        type(dual2)::h1etot,h2etot
        integer::j

        if (errorflag .ne. 0) return
        ham_tot=0.d0
        !$omp parallel do private(h1etot,h2etot,j,z2d) shared(an_cr,an2_cr2,elecs,z1d,row,ham_tot)
        do j=1,size
            h1etot=0.0d0; h2etot=0.0d0
            if (j.ne.row) then
                call dual_2_dual2(zstore(j)%val,z2d,2) 
                h1etot = haml_vals(z1d,z2d,an_cr,elecs%h1ei,elecs%h1_num)
                h2etot = haml_vals(z1d,z2d,an2_cr2,elecs%h2ei,elecs%h2_num)
                ham_tot(j)=ham_tot(j)+(0.5d0*h2etot)+h1etot
            else 
                h1etot = haml_vals(z1d,z1d,an_cr,elecs%h1ei,elecs%h1_num)
                h2etot = haml_vals(z1d,z1d,an2_cr2,elecs%h2ei,elecs%h2_num)
                ham_tot(j)=ham_tot(j)+(0.5d0*h2etot)+h1etot
            end if 
        end do 
        !$omp end parallel do
        
        return

    end function haml_column

    ! function to calcualte an entire column of the overlap 
    function ovrlp_column(z1d,size,zstore,row) result(ovrlp_tot)

        implicit none
        
        type(zombiest),dimension(:),intent(in)::zstore
        type(dual2),dimension(0:),intent(in)::z1d
        integer,intent(in)::row,size
        type(dual2),dimension(size)::ovrlp_tot
        type(dual2),dimension(0:2*norb)::z2d
       
        integer::j
        
        if (errorflag .ne. 0) return
       
        ovrlp_tot=0.0d0
  
        do j=1,size
            if(j.ne.row)then
                call dual_2_dual2(zstore(j)%val,z2d,2) 
                ovrlp_tot(j)=overlap_1(z1d,z2d)
            else
                ovrlp_tot(j)=1.0d0
            end if 
        end do 

        return

    end function ovrlp_column

    !##############################################################################################################################

    !Level 3 routine to make individual value in overlap and hamiltonian matrix

    ! calculates indvidual hamiltonian elements taking in two Zombie states and a set of 
    ! creation and annihilation operations
    function haml_vals(z1d,z2d,ops,el,len) result(ham_tot)

        implicit none 
        type(dual2),dimension(0:),intent(in)::z1d,z2d
        real(wp),dimension(:),intent(in)::el
        integer,intent(in)::len
        type(oprts_2),intent(in)::ops
        type(dual2)::ham_tot
        type(dual2)::ov
        integer::j,k
        
        if (errorflag .ne. 0) return

       
        ham_tot=0.0d0
     
       
   
        do j=1,len
            ov=1.0d0
            do k=1, norb
                ov=ov*((z1d(k)*z2d(ops%alive(k,j))*ops%neg_alive(k,j))+(z1d((k+norb))*z2d((ops%dead(k,j)))*ops%neg_dead(k,j))) 
            end do
            ham_tot=ham_tot+(ov*el(j))
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
