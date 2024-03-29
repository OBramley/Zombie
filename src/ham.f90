MODULE ham 

    use globvars
    use alarrays
    use omp_lib
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

        Call dgetrf(size, size, haml%inv, size, IPIV1, ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if
        if (ierr==0) call dgetri(size,haml%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if

        deallocate(WORK1,IPIV1)

        call DGEMM("N","N",size,size,size,1.d0,haml%inv,size,haml%hjk,size,0.d0,haml%kinvh,size)

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
        integer,intent(in)::verb,size
        integer::j,k,ierr
        real(kind=8)::h1etot,h2etot
    
        if (errorflag .ne. 0) return 
        ierr=0

        !$omp parallel do &
        !$omp & private(j,k,h1etot,h2etot) &
        !$omp & shared(elecs,zstore,an_cr,an2_cr2,haml) 
        do j=1,size
            haml%ovrlp(j,j)=1.d0
            h1etot = haml_vals(zstore(j)%val,zstore(j)%val,an_cr,elecs%h1ei,elecs%h1_num)
            h2etot = haml_vals(zstore(j)%val,zstore(j)%val,an2_cr2,elecs%h2ei,elecs%h2_num)
            haml%hjk(j,j)=h1etot+(0.5*h2etot)+(haml%ovrlp(j,j)*elecs%hnuc)
            do k=j+1,size
                haml%ovrlp(j,k)=overlap_1(zstore(j)%val,zstore(k)%val);haml%ovrlp(k,j)=haml%ovrlp(j,k)
                h1etot = haml_vals(zstore(j)%val,zstore(k)%val,an_cr,elecs%h1ei,elecs%h1_num)
                h2etot = haml_vals(zstore(j)%val,zstore(k)%val,an2_cr2,elecs%h2ei,elecs%h2_num)
                haml%hjk(j,k)=h1etot+(0.5*h2etot)+(haml%ovrlp(j,k)*elecs%hnuc);haml%hjk(k,j)=haml%hjk(j,k)
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
        real(kind=8),dimension(0:),intent(in)::z1d
        type(elecintrgl),intent(in)::elecs
        type(oprts_2),intent(in)::an_cr,an2_cr2
        type(hamiltonian),intent(inout)::haml
        integer,intent(in)::row,size
        real(kind=8)::h1etot,h2etot
        integer::j

        if (errorflag .ne. 0) return

        !$omp parallel do private(h1etot,h2etot,j) shared(haml,zstore,an_cr,an2_cr2,elecs,z1d,row)
        do j=1,size
            if (j.ne.row) then
                haml%ovrlp(j,row)=overlap_1(z1d,zstore(j)%val)
                haml%ovrlp(row,j)=haml%ovrlp(j,row)
                h1etot = haml_vals(z1d,zstore(j)%val,an_cr,elecs%h1ei,elecs%h1_num)
                h2etot = haml_vals(z1d,zstore(j)%val,an2_cr2,elecs%h2ei,elecs%h2_num)
                haml%hjk(j,row)=h1etot+(0.5*h2etot)+(haml%ovrlp(j,row)*elecs%hnuc)
                haml%hjk(row,j)=haml%hjk(j,row)
            else 
                haml%ovrlp(row,row)=1.d0
                h1etot = haml_vals(z1d,z1d,an_cr,elecs%h1ei,elecs%h1_num)
                h2etot = haml_vals(z1d,z1d,an2_cr2,elecs%h2ei,elecs%h2_num)
                haml%hjk(row,row)=h1etot+(0.5*h2etot)+elecs%hnuc
            end if 
        end do 
        !$omp end parallel do

        return

    end subroutine haml_ovrlp_column

    function haml_column(z1d,zstore,size,an_cr,an2_cr2,elecs,row)

        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        real(kind=8),dimension(0:),intent(in)::z1d
        type(elecintrgl),intent(in)::elecs
        type(oprts_2),intent(in)::an_cr,an2_cr2
        real(kind=8),dimension(ndet)::haml_column
        integer,intent(in)::row,size
        real(kind=8)::h1etot,h2etot
        integer::j

        if (errorflag .ne. 0) return
        haml_column=0.d0
       !$omp parallel do private(h1etot,h2etot,j) shared(zstore,an_cr,an2_cr2,elecs,z1d,row)
        do j=1,size
            if (j.ne.row) then 
                h1etot = haml_vals(z1d,zstore(j)%val,an_cr,elecs%h1ei,elecs%h1_num)
                haml_column(j)=haml_column(j)+h1etot
                h2etot = haml_vals(z1d,zstore(j)%val,an2_cr2,elecs%h2ei,elecs%h2_num)
                haml_column(j)=haml_column(j)+(0.5*h2etot)
            else 
                h1etot = haml_vals(z1d,z1d,an_cr,elecs%h1ei,elecs%h1_num)
                haml_column(j)=haml_column(j)+h1etot
                h2etot = haml_vals(z1d,z1d,an2_cr2,elecs%h2ei,elecs%h2_num)
                haml_column(j)=haml_column(j)+(0.5*h2etot)
            end if 
        end do 
        !$omp end parallel do
        
        return

    end function haml_column

    ! function to calcualte an entire column of the overlap 
    function ovrlp_column(z1d,size,zstore,row)

        implicit none
        
        type(zombiest),dimension(:),intent(in)::zstore
        real(kind=8),dimension(0:),intent(in)::z1d
        integer,intent(in)::row,size
        real(kind=8),dimension(ndet)::ovrlp_column
        integer::j
        
        if (errorflag .ne. 0) return
       
        ovrlp_column=1.0
        !!$omp parallel do & 
        !!$omp & shared(z1d,zstore,ovrlp_column,row) &
        !!$omp & private(j)
        do j=1,size
            if(j.ne.row)then 
                ovrlp_column(j)=overlap_1(z1d,zstore(j)%val)
            else
                ovrlp_column(j)=1.0
            end if 
        end do 
        !!$omp end parallel do
        return

    end function ovrlp_column

    !##############################################################################################################################

    !Level 3 routine to make individual value in overlap and hamiltonian matrix

    ! calculates indvidual hamiltonian elements taking in two Zombie states and a set of 
    ! creation and annihilation operations
    real(kind=8) function haml_vals(z1d,z2d,ops,el,len)

        implicit none 
        real(kind=8),dimension(0:),intent(in)::z1d,z2d
        real(kind=8),dimension(:),intent(in)::el
        integer,intent(in)::len
        type(oprts_2),intent(in)::ops
        real(kind=8)::ov
        integer::j,k
        
        if (errorflag .ne. 0) return

        haml_vals=0.0
        
     
       
        !!$omp parallel do reduction(+:haml_vals) private(k,ov)
        do j=1,len
            ov=1.0
            do k=1, norb
                ov=ov*((z1d(k)*z2d(ops%alive(k,j))*ops%neg_alive(k,j))+(z1d(k+norb)*z2d(ops%dead(k,j))*ops%neg_dead(k,j))) 
            end do
            haml_vals=haml_vals+(ov*el(j))
        end do
        !!$omp end parallel do 
     
        
        return 
      
    end function haml_vals

    ! calculates individual overlaps where no creation and annihilation operations are needed
    real(kind=8) function overlap_1(z1d,z2d)

        implicit none
        real(kind=8),dimension(0:)::z1d,z2d
        integer::j

        if (errorflag .ne. 0) return

        overlap_1=1.0
      
        !!$omp parallel do reduction(*:overlap_1)
        do j=1,norb
            overlap_1=overlap_1*((z1d(j)*z2d(j))+(z1d(j+norb)*z2d(norb+j)))
        end do
        !!$omp end parallel do 
      
        

        return 
    end function overlap_1

    !##############################################################################################################################
    

   
    ! ! hamliltonian calcualtion - calcualtes the whole hamliltonian 
    ! subroutine haml_make(haml,zstore,elecs,an_cr,an2_cr2,verb) 

    !     implicit none
        
    !     real(kind=8),dimension(:,:),intent(inout)::haml 
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     type(elecintrgl),intent(in)::elecs
    !     type(oprts_2),intent(in)::an_cr,an2_cr2
    !     integer,intent(in)::verb
    !     integer::j,ierr
      
    
    !     if (errorflag .ne. 0) return 
    !     ierr=0

    !     call omp_set_nested(.true.)
    
    !     !$omp parallel do schedule(guided) num_threads(6) &
    !     !$omp & private(j) &
    !     !$omp & shared(elecs,zstore,an_cr,an2_cr2,haml) 

    !     do j=ndet,1,-1
    !         call haml_column(haml(j:,j),zstore(j)%val,zstore,an_cr,an2_cr2,elecs,j) 
    !         if(verb.eq.1)then
    !             write(6,"(a,i0,a)") "hamliltonian column ",j, " completed"
    !         end if 
    !     end do
      
    !     !$omp end parallel do

    !     do j=1,ndet
    !         haml(j,:)=haml(:,j)
    !     end do

    !     return
        
    ! end subroutine haml_make

    ! !subroutine calcualates whole overlap matrix
    ! subroutine ovrlp_make(ovrlp,zstore)

    !     implicit none
        
    !     real(kind=8),dimension(:,:),intent(inout)::ovrlp
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     integer::j,k
        
    !     if (errorflag .ne. 0) return
    !     ovrlp=1.0d0
    !     !!$omp parallel do private(j,k) shared(ovrlp,zstore)
    !     do j=1,ndet
    !         do k=j+1,ndet
    !             ovrlp(j,k)=overlap_1(zstore(j)%val,zstore(k)%val); ovrlp(k,j)=ovrlp(j,k)
    !         end do
    !     end do  
    !     !!$omp end parallel do 
  

    !     return

    ! end subroutine ovrlp_make
    !Level 2 routines to calcualte overlap and hamiltonian gradient and hessian matrices

    ! ! subroutine that finds the gradient of the overlap w.r.t a specified state
    ! subroutine ovrlp_make_grad(zstore,state,ovrlp_grad,cmplt,loops)

    !     implicit none 
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     real(kind=8),dimension(:,:),intent(inout)::ovrlp_grad
    !     integer,intent(in)::state,loops
    !     real(kind=8),dimension(0:2*norb)::z1d
    !     integer,dimension(:),intent(in)::cmplt
    !     integer::j 

    !     if (errorflag .ne. 0) return
       
    !     do j=1, norb
    !         z1d(0:2*norb)=zstore(state)%val(0:2*norb)
    !         z1d(j)=zstore(state)%cos(j)
    !         z1d(j+norb)=(-1)*zstore(state)%sin(j)
    !         ovrlp_grad(j,:)=ovrlp_column_grad(z1d,zstore,state,cmplt,loops)
    !     end do
       

    !     return 

    ! end subroutine ovrlp_make_grad

    !  ! subroutine that finds the gradient of the overlap w.r.t a specified state
    ! subroutine ovrlp_make_grad_one_elec(zstore,state,ovrlp_grad,orb,cmplt,loops)

    !     implicit none 
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     real(kind=8),dimension(:,:),intent(inout)::ovrlp_grad
    !     integer,intent(in)::state,orb,loops
    !     real(kind=8),dimension(0:2*norb)::z1d
    !     integer,dimension(:),intent(in)::cmplt
        

    !     if (errorflag .ne. 0) return

    !     z1d(0:2*norb)=zstore(state)%val(0:2*norb)
    !     z1d(orb)=zstore(state)%cos(orb)
    !     z1d(orb+norb)=(-1)*zstore(state)%sin(orb)
    !     ovrlp_grad(orb,:)=ovrlp_column_grad(z1d,zstore,state,cmplt,loops)
      
    
    !     return 

    ! end subroutine ovrlp_make_grad_one_elec

    ! !  ! Hamiltonian calcualtion - calcualtes the gradient of ther hamliltonian w.r.t one zombie state 
    ! subroutine haml_grad(haml_diff,zstore,elecs,an_cr,an2_cr2,state,cmplt,loops) 

    !     implicit none
        
    !     real(kind=8),dimension(:,:),intent(inout)::haml_diff 
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     type(elecintrgl),intent(in)::elecs
    !     type(oprts),intent(in)::an_cr,an2_cr2
    !     integer,intent(in)::state,loops
    !     integer,dimension(:),intent(in)::cmplt
    !     real(kind=8),dimension(0:2*norb)::z1d
  
    !     integer::j,ierr
    
    !     if (errorflag .ne. 0) return 
    !     ierr=0

    !     if(loops.le.2)then
    !         !$omp parallel do schedule(dynamic) &
    !         !$omp & private(j,z1d) &
    !         !$omp & shared(elecs,zstore,an_cr,an2_cr2,haml_diff,cmplt,state)
    !         do j=1, norb
    !             z1d(0:2*norb)=zstore(state)%val(0:2*norb)
    !             z1d(j)=zstore(state)%cos(j)
    !             z1d(j+norb)=zstore(state)%sin(j)*(-1)
    !             call haml_grad_rc(haml_diff(j,:),z1d,zstore,an_cr,an2_cr2,elecs,state,j,cmplt,loops)
    !         end do
    !         !$omp end parallel do
    !     else 
    !         !$omp parallel do schedule(dynamic) num_threads(6)&
    !         !$omp & private(j,z1d) &
    !         !$omp & shared(elecs,zstore,an_cr,an2_cr2,haml_diff,cmplt,state)
    !         do j=1, norb
    !             z1d(0:2*norb)=zstore(state)%val(0:2*norb)
    !             z1d(j)=zstore(state)%cos(j)
    !             z1d(j+norb)=zstore(state)%sin(j)*(-1)
    !             call haml_grad_rc_p(haml_diff(j,:),z1d,zstore,an_cr,an2_cr2,elecs,state,j,cmplt,loops)
    !         end do
    !         !$omp end parallel do
    !     end if 

    !     return

    ! end subroutine haml_grad

    ! subroutine haml_grad_one_elec(haml_diff,zstore,elecs,an_cr,an2_cr2,state,orb,cmplt,loops) 

    !     implicit none
        
    !     real(kind=8),dimension(:,:),intent(inout)::haml_diff 
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     type(elecintrgl),intent(in)::elecs
    !     type(oprts),intent(in)::an_cr,an2_cr2
    !     integer,intent(in)::state,orb,loops
    !     integer,dimension(:),intent(in)::cmplt
    !     real(kind=8),dimension(0:2*norb)::z1d
        
   
    
    !     if (errorflag .ne. 0) return 
        

    !     z1d(0:2*norb)=zstore(state)%val(0:2*norb)
    !     z1d(orb)=zstore(state)%cos(orb)
    !     z1d(orb+norb)=zstore(state)%sin(orb)*(-1)
    !     if(loops.le.2)then
    !         call haml_grad_rc(haml_diff(orb,:),z1d,zstore,an_cr,an2_cr2,elecs,state,orb,cmplt,loops)
    !     else
    !         call haml_grad_rc_p(haml_diff(orb,:),z1d,zstore,an_cr,an2_cr2,elecs,state,orb,cmplt,loops)
    !     end if 

      
      
    !     return

    ! end subroutine haml_grad_one_elec

        !Level 3 routines to calcualte columns of gradient and hessian matrices

!     ! subroutine that finds the gradient w.r.t to a specfic zombie state element 
!     ! for an entire column 
!     function ovrlp_column_grad(z1d,zstore,state,cmplt,loops)

!         implicit none
        
!         type(zombiest),dimension(:),intent(in)::zstore
!         real(kind=8),dimension(0:),intent(in)::z1d
!         integer,intent(in)::state,loops
!         real(kind=8),dimension(ndet)::ovrlp_column_grad
!         integer,dimension(:),intent(in)::cmplt
!         integer::k,j

!         if (errorflag .ne. 0) return
        
!         ovrlp_column_grad=0.0
!         !!$omp parallel do schedule(dynamic) &
!         !!$omp & shared(z1d,zstore,ovrlp_column_grad,state) &
!         !!$omp & private(j)
!         do k=1,loops
!             j=cmplt(k)
!             if(j.ne.state)then 
!                 ovrlp_column_grad(j)=overlap_1(z1d,zstore(j)%val)
!             else
!                 ovrlp_column_grad(j)=0.0
!             end if
!         end do 
!         !!$omp end parallel do
!         return

!     end function ovrlp_column_grad


!     ! Calcualates a column of a hamliltonian Start specifies the row the column
!     ! is started to be calcualted 
!     subroutine haml_grad_rc(hcol,z1d,zstore,an_cr,an2_cr2,elecs,state,orb,cmplt,loops)

!         implicit none
!         real(kind=8),dimension(:),intent(inout)::hcol 
!         type(zombiest),dimension(:),intent(in)::zstore
!         real(kind=8),dimension(0:),intent(in)::z1d
!         type(elecintrgl),intent(in)::elecs
!         type(oprts),intent(in)::an_cr,an2_cr2
!         integer,intent(in)::state,orb,loops
!         integer,dimension(:),intent(in)::cmplt
!         real(kind=8)::h1etot,h2etot
!         integer::j,k
        
        
!         if (errorflag .ne. 0) return
!         !!$omp parallel do private(h1etot,h2etot,j) shared(hcol,zstore,an_cr,an2_cr2,elecs,z1d,state,orb,cmplt)
!         do k=1,loops
!             j=cmplt(k)
!             if(j.ne.state)then
!                 !! Differentiating the bra 1 el
!                 h1etot = haml_vals(z1d,zstore(j)%val,an_cr%ham,elecs%h1ei,elecs%h1_num)
!                 hcol(j)=hcol(j)+h1etot
!                 !Differentiating the bra 2 el
!                 h2etot = haml_vals(z1d,zstore(j)%val,an2_cr2%ham,elecs%h2ei,elecs%h2_num)
!                 hcol(j)=hcol(j)+(0.5*h2etot)
!             else
!                 !Differentiaitn hamiltonian element (a,a) only placed in hamiltonian columm
!                 h1etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an_cr%diff(orb),elecs%h1ei,an_cr%dcnt(0:,orb))
!                 hcol(j)=hcol(j)+h1etot
!                 h2etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an2_cr2%diff(orb),elecs%h2ei,an2_cr2%dcnt(0:,orb))
!                 hcol(j)=hcol(j)+(0.5*h2etot)
!             end if
           
!         end do 
!         !!$omp end parallel do  
!         return
        

!     end subroutine haml_grad_rc
!  ! Calcualates a column of a hamliltonian Start specifies the row the column
!     ! is started to be calcualted 
!     subroutine haml_grad_rc_p(hcol,z1d,zstore,an_cr,an2_cr2,elecs,state,orb,cmplt,loops)

!         implicit none
!         real(kind=8),dimension(:),intent(inout)::hcol 
!         type(zombiest),dimension(:),intent(in)::zstore
!         real(kind=8),dimension(0:),intent(in)::z1d
!         type(elecintrgl),intent(in)::elecs
!         type(oprts),intent(in)::an_cr,an2_cr2
!         integer,intent(in)::state,orb,loops
!         integer,dimension(:),intent(in)::cmplt
!         real(kind=8)::h1etot,h2etot
!         integer::j,k
        
!         if (errorflag .ne. 0) return

!         !$omp parallel do private(h1etot,h2etot,j) shared(hcol,zstore,an_cr,an2_cr2,elecs,z1d,state,orb,cmplt)
!         do k=1,loops
!             j=cmplt(k)
!             if(j.ne.state)then
!                 !! Differentiating the bra 1 el
                
!                 h1etot = haml_vals(z1d,zstore(j)%val,an_cr%ham,elecs%h1ei,elecs%h1_num)
!                 hcol(j)=hcol(j)+h1etot
        
!                 !Differentiating the bra 2 el
!                 h2etot = haml_vals(z1d,zstore(j)%val,an2_cr2%ham,elecs%h2ei,elecs%h2_num)
!                 hcol(j)=hcol(j)+(0.5*h2etot)
!             else
!                 !Differentiaitn hamiltonian element (a,a) only placed in hamiltonian columm
!                 h1etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an_cr%diff(orb),elecs%h1ei,an_cr%dcnt(0:,orb))
!                 hcol(j)=hcol(j)+h1etot
!                 h2etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an2_cr2%diff(orb),elecs%h2ei,an2_cr2%dcnt(0:,orb))
!                 hcol(j)=hcol(j)+(0.5*h2etot)
!             end if
!         end do 
!         !$omp end parallel do  
!         return
        

!     end subroutine haml_grad_rc_p

    ! ! subroutine that finds the gradient of the overlap w.r.t a specified state
    ! subroutine ovrlp_make_hessian(zstore,state,ovrlp_hess,cmplt)

    !     implicit none 
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     real(kind=8),dimension(:,:,:),intent(inout)::ovrlp_hess
    !     integer,intent(in)::state
    !     real(kind=8),dimension(0:2*norb)::z1d,z1dd
    !     integer,dimension(:),intent(in)::cmplt
    !     integer::j,k 

    !     !$omp parallel do shared(zstore,state,ovrlp_hess) private(z1d)
    !     do j=1, norb
    !         z1d(0:2*norb)=zstore(state)%val(0:2*norb)
    !         z1d(j)=zstore(state)%cos(j)
    !         z1d(j+norb)=(-1)*zstore(state)%sin(j)
    !         do k=1,norb
    !             z1dd(0:2*norb)=z1d(0:2*norb)
    !             if(j.ne.k)then
    !                 z1dd(k)=zstore(state)%cos(k)
    !                 z1dd(k+norb)=(-1)*zstore(state)%sin(k)
    !             else 
    !                 z1dd(k)=(-1)*zstore(state)%sin(k)
    !                 z1dd(k+norb)=(-1)*zstore(state)%cos(k)
    !             end if
    !             ovrlp_hess(j,k,:)=ovrlp_column_grad(z1dd,zstore,state,cmplt)
    !         end do
    !     end do
    !     !$omp end parallel do 

    !     return 

    ! end subroutine ovrlp_make_hessian

    ! subroutine haml_hessian(haml_hess,zstore,elecs,an_cr,an2_cr2,state,cmplt) 

    !     implicit none
        
    !     real(kind=8),dimension(:,:,:),intent(inout)::haml_hess 
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     type(elecintrgl),intent(in)::elecs
    !     type(oprts),intent(in)::an_cr,an2_cr2
    !     integer,intent(in)::state
    !     integer,dimension(:),intent(in)::cmplt
    !     real(kind=8),dimension(0:2*norb)::z1d,z1dd
  
    !     integer::j,ierr,k
    
    !     if (errorflag .ne. 0) return 
    !     ierr=0

       
    !     !$omp parallel do &
    !     !$omp & private(j) &
    !     !$omp & shared(elecs,zstore,an_cr,an2_cr2,haml_hess)
    !     do j=1,norb
    !         z1d(0:2*norb)=zstore(state)%val(0:2*norb)
    !         z1d(j)=zstore(state)%cos(j)
    !         z1d(j+norb)=(-1)*zstore(state)%sin(j)
    !         do k=j,norb
    !             z1dd(0:2*norb)=z1d(0:2*norb)
    !             if(j.ne.k)then
    !                 z1dd(k)=zstore(state)%cos(k)
    !                 z1dd(k+norb)=(-1)*zstore(state)%sin(k)
    !             else 
    !                 z1dd(k)=(-1)*zstore(state)%sin(k)
    !                 z1dd(k+norb)=(-1)*zstore(state)%cos(k)
    !             end if
    !             call haml_hess_rc(haml_hess(j,k,:),z1dd,zstore,an_cr,an2_cr2,elecs,state,j,k,cmplt)
    !         end do
    !     end do
    !     !$omp end parallel do

    !     return
        
    ! end subroutine haml_hessian

    ! real(kind=8) function haml_vals_2(z1d,z2d,ops,el,len,neg)

    !     implicit none 
    !     real(kind=8),dimension(0:),intent(in)::z1d,z2d
    !     real(kind=8),dimension(:),intent(in)::el
    !     integer(kind=1),dimension(:),intent(in)::neg
    !     integer,intent(in)::len
    !     type(oprts_2),intent(in)::ops
    !     real(kind=8)::ov
    !     integer::j,k

        
    !     haml_vals_2=0.0
    !     !$omp parallel private(j,k,ov) shared(ops,z1d,z2d)
    !     !$omp do simd reduction(+:haml_vals) 
    !     do j=1,len
    !         ov=1.0
    !         !!$omp do simd reduction(*:ov)
    !         do k=1, norb
    !             ov=ov*((z1d(k)*z2d(ops%alive(k,j))*ops%neg_alive(k,j))+(z1d(k+norb)*z2d(ops%dead(k,j))*ops%neg_dead(k,j))) 
    !         end do
    !         !!$omp end do simd
    !         haml_vals_2=haml_vals_2+(ov*el(j)*neg(j))
    !     end do
    !     !$omp end do simd
    !     !$omp end parallel 
        
    !     return 
      
    ! end function haml_vals_2

     ! ! calculates indvidual hamliltonian elements taking in two Zombie states and a set of 
    ! ! creation and annihilation operations
    ! real(kind=8) function haml_val_grad(z1d,z2d,ops,el,orb)

    !     implicit none 
    !     real(kind=8),dimension(0:),intent(in)::z1d,z2d
    !     real(kind=8),dimension(:),intent(in)::el
    !     integer,intent(in)::orb
    !     type(oprts),intent(in)::ops
    !     real(kind=8)::ov
    !     integer::j,k,len

    !     len=ops%dcnt(0,orb)
    !     haml_val_grad=0.0
    !     !$omp parallel private(j,k,ov) shared(ops,z1d,z2d)
    !     !$omp do simd reduction(+:haml_val_grad) 
    !     do j=1,len
    !         ov=1.0
    !         !!$omp do simd reduction(*:ov)
    !         do k=1, norb
    !             ov=ov*((z1d(k)*z2d(ops%alive_diff(orb,k,(ops%dcnt(j,orb))))*ops%neg_alive_diff(orb,k,(ops%dcnt(j,orb))))&
    !             +(z1d(k+norb)*z2d(ops%dead_diff(orb,k,(ops%dcnt(j,orb))))*ops%neg_dead_diff(orb,k,(ops%dcnt(j,orb))))) 
    !         end do
    !         !!$omp end do simd
    !         haml_val_grad=haml_val_grad+(ov*el(ops%dcnt(j,orb)))
    !     end do
    !     !$omp end do simd
    !     !$omp end parallel 
        
    !     return 
      
    ! end function haml_val_grad

    ! ! calculates indvidual hamliltonian elements taking in two Zombie states and a set of 
    ! ! creation and annihilation operations
    ! real(kind=8) function haml_val_hess(z1d,z2d,ops,el,orb1,orb2)

    !     implicit none 
    !     real(kind=8),dimension(0:),intent(in)::z1d,z2d
    !     real(kind=8),dimension(:),intent(in)::el
    !     integer,intent(in)::orb1,orb2
    !     type(oprts),intent(in)::ops
    !     real(kind=8)::ov
    !     integer::j,k,len

    !     len=ops%hcnt(0,orb1,orb2)
    !     haml_val_hess=0.0
    !     !$omp parallel private(j,k,ov) shared(ops,z1d,z2d)
    !     !$omp do simd reduction(+:haml_val_hess) 
    !     do j=1,len
    !         ov=1.0
    !         !!$omp do simd reduction(*:ov)
    !         do k=1, norb
    !             ov=ov*((z1d(k)*z2d(ops%alive_hess(orb1,orb2,k,(ops%hcnt(j,orb1,orb2))))*&
    !             ops%neg_alive_hess(orb1,orb2,k,(ops%hcnt(j,orb1,orb2))))+&
    !             (z1d(k+norb)*z2d(ops%dead_hess(orb1,orb2,k,(ops%hcnt(j,orb1,orb2))))*&
    !             ops%neg_dead_hess(orb1,orb2,k,(ops%hcnt(j,orb1,orb2))))) 
    !         end do
    !         !!$omp end do simd
    !         haml_val_hess=haml_val_hess+(ov*el(ops%hcnt(j,orb1,orb2)))
    !     end do
    !     !$omp end do simd
    !     !$omp end parallel 
        
    !     return 
      
    ! end function haml_val_hess

    
    ! real(kind=8) function haml_gvals(z1d,z2d,ops,ops2,el,len,orb)

    !     implicit none 
    !     real(kind=8),dimension(0:),intent(in)::z1d,z2d
    !     real(kind=8),dimension(:),intent(in)::el
    !     integer,intent(in)::len,orb
    !     type(oprts),intent(in)::ops,ops2
    !     real(kind=8)::ov
    !     integer::j,k

        
    !     haml_gvals=0.0
    !     !$omp parallel private(j,k,ov) shared(ops,z1d,z2d)
    !     !$omp do simd reduction(+:haml_gvals) 
    !     do j=1,len
    !         ov=1.0
    !         !!$omp do simd reduction(*:ov)
    !         do k=1, norb
    !             ov=ov*((z1d(ops2%alive(k,orb))*ops2%neg_alive(k,orb)*z2d(ops%alive(k,j))*&
    !             ops%neg_alive(k,j))+(z1d(ops2%dead(k,orb))*ops2%neg_dead(k,orb)*z2d(ops%dead(k,j))*ops%neg_dead(k,j)))
                
    !             ! ov=ov*((z1d(ops%alive(k,j+len))*ops%neg_alive(k,j+len)*z2d(ops%alive(k,j))*&
    !             ! ops%neg_alive(k,j))+(z1d(ops%dead(k,j+len))*ops%neg_dead(k,j+len)*z2d(ops%dead(k,j))*ops%neg_dead(k,j))) 
    !         end do
            
    !         haml_gvals=haml_gvals+(ov*el(j))
    !     end do
    !     !$omp end do simd
    !     !$omp end parallel 
        
    !     return 
      
    ! end function haml_gvals

    ! subroutine haml_hess_rc(hcol,z1d,zstore,an_cr,an2_cr2,elecs,state,orb1,orb2,cmplt)

    !     implicit none
    !     real(kind=8),dimension(:),intent(inout)::hcol 
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     real(kind=8),dimension(0:),intent(in)::z1d
    !     type(elecintrgl),intent(in)::elecs
    !     type(oprts),intent(in)::an_cr,an2_cr2
    !     integer,intent(in)::state,orb1,orb2
    !     integer,dimension(:),intent(in)::cmplt
    !     real(kind=8)::h1etot,h2etot
    !     integer::j
        
        
    !     !$omp parallel 
    !     !$omp single
    !     do j=1,ndet
    !         if(cmplt(j).eq.0)then
    !             if(j.ne.state)then
    !                 ! Differentiating the bra 1 el
    !                 !$omp task firstprivate(h1etot,j) shared(hcol,zstore,an_cr,elecs,z1d)
    !                 h1etot = haml_vals(z1d,zstore(j)%val,an_cr%ham,elecs%h1ei,elecs%h1_num)
    !                 !$omp atomic
    !                 hcol(j)=hcol(j)+h1etot
    !                 !$omp end atomic
    !                 !$omp end task
                
    !                 !Differentiating the bra 2 el
    !                 !$omp task firstprivate(h2etot,j) shared(hcol,zstore,an2_cr2,elecs,z1d)
    !                 h2etot = haml_vals(z1d,zstore(j)%val,an2_cr2%ham,elecs%h2ei,elecs%h2_num)
    !                 !$omp atomic
    !                 hcol(j)=hcol(j)+(0.5*h2etot)
    !                 !$omp end atomic
    !                 !$omp end task

    !             else
    !                 !Differentiaitn hamiltonian element (a,a) only placed in hamiltonian column
    !                 !$omp task firstprivate(h1etot,j) shared(hcol,zstore,an_cr,elecs)
    !                 h1etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an_cr%hess(orb1,orb2),&
    !                 elecs%h1ei,an_cr%hcnt(:,orb1,orb2))
    !                 !$omp atomic
    !                 hcol(j)=hcol(j)+h1etot
    !                 !$omp end atomic
    !                 !$omp end task
    !                 !$omp task firstprivate(h2etot,j) shared(hcol,zstore,an2_cr2,elecs)
    !                 h2etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an2_cr2%hess(orb1,orb2),&
    !                 elecs%h2ei,an_cr%hcnt(:,orb1,orb2))
    !                 !$omp atomic
    !                 hcol(j)=hcol(j)+(0.5*h2etot)
    !                 !$omp end atomic
    !                 !$omp end task
    !             end if
    !         end if 
    !     end do 
    !     !$omp end single
    !     !$omp end parallel  

    !     return

    ! end subroutine haml_hess_rc


END MODULE ham
