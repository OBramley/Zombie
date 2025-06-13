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

    subroutine hamgen_inv(haml,size)

        implicit none 

        type(hamiltonian), intent(inout)::haml
        integer,intent(in)::size
        integer, allocatable,dimension(:)::IPIV1
        real(wp),allocatable,dimension(:)::WORK1
        integer::ierr=0

        if (errorflag .ne. 0) return
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

    end subroutine hamgen_inv

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
        ! INTEGER(int64) :: c1, c2, cr, cm
        ! REAL(real64) :: rate

        if (errorflag .ne. 0) return 
      
        ! call omp_set_nested(.TRUE.)
        ! CALL SYSTEM_CLOCK(count_rate=cr)
        ! CALL SYSTEM_CLOCK(count_max=cm)
        ! rate = REAL(cr)
        ! CALL SYSTEM_CLOCK(c1)
        do j=1,size
            !$omp  parallel do & 
            !$omp & private(k) &
            !$omp & shared(j,elecs,size,zstore,haml)
            do k=j,size
                call haml_vals_2(zstore(j)%val,zstore(k)%val,haml%ovrlp(j,k),haml%hjk(j,k),elecs,(k-j))
            end do
            !$OMP end parallel do
            if(verb.eq.1)then
                write(stdout,"(a,i0,a)") "hamliltonian column ",j, " completed"
            end if 
        end do
       
        !$omp  parallel do
        do j=1,size
            do k=j,size
                haml%hjk(k,j)=haml%hjk(j,k)
                haml%ovrlp(k,j)=haml%ovrlp(j,k)
            end do
        end do 
        !$omp end parallel do
        ! CALL SYSTEM_CLOCK(c2)
        ! PRINT *, "Elapsed time (system_clock): ", real((c2 - c1),real64) / rate
     
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
        integer::j

        if (errorflag .ne. 0) return
    
        !$omp parallel do  default(none) &
        !$omp & private(j) &
        !$omp & shared(elecs,zstore,temp,row,norb,size)
        do j=1,size
            if (j.ne.row) then
                call haml_vals_2(temp%zom%val,zstore(j)%val,temp%ovrlp(j,row),temp%hjk(j,row),elecs,abs(j-row))
            else
                call haml_vals_2(temp%zom%val,temp%zom%val,temp%ovrlp(row,row),temp%hjk(row,row),elecs,abs(j-row))
            end if 
        end do
        !$omp end parallel do 
      
        temp%hjk(row,:)=temp%hjk(:,row) 
        temp%ovrlp(row,:)=temp%ovrlp(:,row)
        
        return

    end subroutine haml_ovrlp_column

    !Level 2 routines to make an Overlap and Hamiltonian matrix column

    ! Calcualates a column of a hamiltonian Start specifies the row the column
    ! is started to be calcualted 

    subroutine haml_ovrlp_column_orb(temp,zstore,size,elecs,row,orb)

        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        type(grad_do),intent(inout)::temp
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::row,orb,size
        integer::j

        if (errorflag .ne. 0) return
      
        !$omp parallel do  default(none) &
        !$omp & private(j) &
        !$omp & shared(elecs,zstore,temp,row,norb,orb,size)
        !!$omp target teams distribute parallel do
        do j=1,size
            if (j.ne.row) then
                temp%ovrlp(j,row)=temp%ovrlp(j,row)/(&
                (zstore(j)%val(orb)*zstore(row)%val(orb))+(zstore(j)%val(orb+norb)*zstore(row)%val(orb+norb)))
                call haml_vals_2_orb(temp%zom%val,zstore(j)%val,temp%ovrlp(j,row),temp%hjk(j,row),elecs,abs(j-row),orb)
            else
                call haml_vals_2_orb(temp%zom%val,temp%zom%val,temp%ovrlp(row,row),temp%hjk(row,row),elecs,abs(j-row),orb)
            end if 
        end do
        !$omp end parallel do 
        temp%hjk(row,:)=temp%hjk(:,row) 
        temp%ovrlp(row,:)=temp%ovrlp(:,row)
        
        return

    end subroutine haml_ovrlp_column_orb

    !##############################################################################################################################

    !Level 3 routine to make individual value in overlap and hamiltonian matrix


    !##############################################################################################################################
    
    subroutine haml_vals_2(z1d,z2d,ovrlp,ham_tot,elecs,sm)
        implicit none 
        real(wp),dimension(0:),intent(in)::z1d,z2d
        real(wp),intent(inout)::ovrlp,ham_tot
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::sm
        real(wp),dimension(4,norb)::perts
        real(wp),dimension(norb)::div,bth
        real(wp)::ov,aa,dd,ad,da
        integer::j,k,l
        real(wp),dimension(elecs%num)::ovrlp_vec
        
        ovrlp=1.0d0; ov=1.0d0
        if(sm==0)then
            !$omp simd
            do j=1,norb
                aa=z1d(j)*z2d(j)
                dd=z1d(j+norb)*z2d(norb+j)
                ad=z1d(j)*z2d(norb+j)
                da=z1d(j+norb)*z2d(j)
                div(j)=aa+dd
                perts(1,j)=da/div(j)
                perts(2,j)=ad/div(j) 
                perts(3,j)=aa/div(j)
                perts(4,j)=(-aa+dd)/div(j)
            end do
            !$omp end simd
            ham_tot=elecs%hnuc
        else 
            do j=1,norb
                aa=z1d(j)*z2d(j)
                dd=z1d(j+norb)*z2d(norb+j)
                ad=z1d(j)*z2d(norb+j)
                da=z1d(j+norb)*z2d(j)
                div(j)=aa+dd
                if(div(j).eq.0.0d0)then
                    ovrlp=0.0d0
                    ham_tot=0.0d0
                    return
                end if
                perts(1,j)=da/div(j)
                perts(2,j)=ad/div(j) 
                perts(3,j)=aa/div(j)
                perts(4,j)=(-aa+dd)/div(j)
                ovrlp=ovrlp*div(j)
            end do
           
            ham_tot=ovrlp*elecs%hnuc
        end if 
        
        ovrlp_vec=ovrlp
        !$omp simd
        do k=1,norb
            do l=1,elecs%orbital_choice2(0,k)
                ov=ovrlp_vec(elecs%orbital_choice2(k,(l*2)-1))*&
                perts(elecs%orbital_choice(k,elecs%orbital_choice2(k,(l*2)-1)),elecs%orbital_choice3(k))
                do j=elecs%orbital_choice2(k,(l*2)-1),elecs%orbital_choice2(k,l*2)
                    ovrlp_vec(j)=ov
                end do
            end do
           
        end do
        !$omp end simd
        
        !$omp simd
        do j=1,elecs%num
            ham_tot=ham_tot+(ovrlp_vec(j)*elecs%integrals(j))
        end do
        !$omp end simd
       

    
        return 
      
    end subroutine haml_vals_2

    subroutine haml_vals_2_orb(z1d,z2d,ovrlp,ham_tot,elecs,sm,orb)
        implicit none 
        real(wp),dimension(0:),intent(in)::z1d,z2d
        real(wp),intent(inout)::ovrlp,ham_tot
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::sm,orb
        real(wp),dimension(4,norb)::perts
        real(wp),dimension(norb)::div,bth
        real(wp)::ov,aa,dd,ad,da
        integer::j,k,l
        real(wp),dimension(elecs%num)::ovrlp_vec
        
        ov=1.0d0
        if(sm==0)then
            ovrlp=1.0d0
            !$omp simd
            do j=1,norb
                aa=z1d(j)*z2d(j)
                dd=z1d(j+norb)*z2d(norb+j)
                ad=z1d(j)*z2d(norb+j)
                da=z1d(j+norb)*z2d(j)
                perts(1,j)=da
                perts(2,j)=ad 
                perts(3,j)=aa
                perts(4,j)=(-aa+dd)
            end do
            !$omp end simd
            ham_tot=elecs%hnuc
        else 
            do j=1,norb
                aa=z1d(j)*z2d(j)
                dd=z1d(j+norb)*z2d(norb+j)
                ad=z1d(j)*z2d(norb+j)
                da=z1d(j+norb)*z2d(j)
                div(j)=aa+dd
                if(div(j).eq.0.0d0)then
                    ovrlp=0.0d0
                    ham_tot=0.0d0
                    return
                end if
                perts(1,j)=da/div(j)
                perts(2,j)=ad/div(j) 
                perts(3,j)=aa/div(j)
                perts(4,j)=(-aa+dd)/div(j)
            end do
            ovrlp=ovrlp*div(orb)
            ham_tot=ovrlp*elecs%hnuc
        end if 
        
        ovrlp_vec=ovrlp

        !$omp simd
        do k=1,norb
            do l=1,elecs%orbital_choice2(0,k)
                ov=ovrlp_vec(elecs%orbital_choice2(k,(l*2)-1))*&
                perts(elecs%orbital_choice(k,elecs%orbital_choice2(k,(l*2)-1)),elecs%orbital_choice3(k))
                do j=elecs%orbital_choice2(k,(l*2)-1),elecs%orbital_choice2(k,l*2)
                    ovrlp_vec(j)=ov
                end do
            end do
        end do
        !$omp end simd
       
        do j=1,elecs%num
            ham_tot=ham_tot+(ovrlp_vec(j)*elecs%integrals(j))
        end do
       
       
        return 
      
    end subroutine haml_vals_2_orb

    subroutine haml_vals_3_orb(z1d_old,z2d_old,z_new,ovrlp,ham_tot,elecs,sm,orb)
        implicit none 
        type(zombiest),intent(in)::z1d_old,z2d_old,z_new
        real(wp),intent(inout)::ovrlp,ham_tot
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::sm,orb
        real(wp),dimension(0:4)::perts
        real(wp)::aa,dd,ad,da,div
        integer::j,k,l
        
        div=ham_tot/elecs%num
        if(sm==0)then
            ! ovrlp=1.0d0
            ! perts(0)=(z_new%val(orb)*z_new%val(orb)+z_new%val(orb+norb)*z_new%val(orb+norb))/&
            ! (z1d_old%val(orb)*z1d_old%val(orb)+z1d_old%val(orb+norb)*z1d_old%val(orb+norb))
            ! perts(1)=(z_new%val(orb+norb)*z_new%val(orb))/(z1d_old%val(orb+norb)*z1d_old%val(orb))
            ! perts(2)=(z_new%val(orb)*z_new%val(orb+norb))/(z1d_old%val(orb)*z1d_old%val(orb+norb)) 
            ! perts(3)=(z_new%val(orb)*z_new%val(orb))/(z1d_old%val(orb)*z1d_old%val(orb))
            ! perts(4)=(-z_new%val(orb)*z_new%val(orb)+z_new%val(orb+norb)*z_new%val(orb+norb))/&
            !             (-z1d_old%val(orb)*z1d_old%val(orb)+z1d_old%val(orb+norb)*z1d_old%val(orb+norb))

            perts(0)=(z1d_old%val(orb)*z1d_old%val(orb)+z1d_old%val(orb+norb)*z1d_old%val(orb+norb))
            if(perts(0).ne.0.0)then
                perts(0)=(z_new%val(orb)*z_new%val(orb)+z_new%val(orb+norb)*z_new%val(orb+norb))/perts(0)
            end if
            perts(1)=(z1d_old%val(orb+norb)*z1d_old%val(orb))
            if(perts(1).ne.0.0)then
                perts(1)=(z_new%val(orb+norb)*z_new%val(orb))/perts(1)
            end if
            perts(2)=(z1d_old%val(orb)*z1d_old%val(orb+norb)) 
            if(perts(2).ne.0.0)then
                perts(2)=(z_new%val(orb)*z_new%val(orb+norb))/perts(2)
            end if 
            perts(3)=(z1d_old%val(orb)*z1d_old%val(orb))
            if(perts(3).ne.0.0)then
                perts(3)=(z_new%val(orb)*z_new%val(orb))/perts(3)
            end if
            perts(4)= (-z1d_old%val(orb)*z1d_old%val(orb)+z1d_old%val(orb+norb)*z1d_old%val(orb+norb))
            if(perts(4).ne.0.0)then
                perts(4)= (-z_new%val(orb)*z_new%val(orb)+z_new%val(orb+norb)*z_new%val(orb+norb))/perts(4)
            end if 
        else
            perts(0)=(z1d_old%val(orb)*z2d_old%val(orb)+z1d_old%val(orb+norb)*z2d_old%val(orb+norb))
            if(perts(0).ne.0.0)then
                perts(0)=(z_new%val(orb)*z2d_old%val(orb)+z_new%val(orb+norb)*z2d_old%val(orb+norb))/perts(0)
            end if
            perts(1)=(z1d_old%val(orb+norb)*z2d_old%val(orb))
            if(perts(1).ne.0.0)then
              perts(1)=(z_new%val(orb+norb)*z2d_old%val(orb))/perts(1)
            end if
            perts(2)=(z1d_old%val(orb)*z2d_old%val(orb+norb)) 
            if(perts(2).ne.0.0)then
                perts(2)=(z_new%val(orb)*z2d_old%val(orb+norb))/perts(2)
            end if 
            perts(3)=(z1d_old%val(orb)*z2d_old%val(orb))
            if(perts(3).ne.0.0)then
                perts(3)=(z_new%val(orb)*z2d_old%val(orb))/perts(3)
            end if
            perts(4)= (-z1d_old%val(orb)*z2d_old%val(orb)+z1d_old%val(orb+norb)*z2d_old%val(orb+norb))
            if(perts(4).ne.0.0)then
                perts(4)= (-z_new%val(orb)*z2d_old%val(orb)+z_new%val(orb+norb)*z2d_old%val(orb+norb))/perts(4)
            end if 
            
            ovrlp=ovrlp*perts(0)
            ham_tot=ovrlp*elecs%hnuc
        end if
       
        
        do j=1,elecs%num
            ham_tot=ham_tot+(div*perts(elecs%orbital_choice(orb,j)))
        end do
       
        return 
      
    end subroutine haml_vals_3_orb


  

END MODULE ham
