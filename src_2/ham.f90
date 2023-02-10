MODULE ham 

    use globvars
    use alarrays
    contains

    ! calculates indvidual hamiltonian elements taking in two Zombie states and a set of 
    ! creation and annihilation operations
    real(kind=8) function haml_vals(z1d,z2d,ops,el,len)

        implicit none 
        real(kind=8),dimension(0:),intent(in)::z1d,z2d
        real(kind=8),dimension(:),intent(in)::el
        integer,intent(in)::len
        type(oprts),intent(in)::ops
        real(kind=8)::ov
        integer::j,k

        
        haml_vals=0.0
        !$omp parallel private(j,k,ov) shared(ops,z1d,z2d)
        !$omp do simd reduction(+:haml_vals) 
        do j=1,len
            ov=1.0
            !!$omp do simd reduction(*:ov)
            do k=1, norb
                ov=ov*((z1d(k)*z2d(ops%alive(k,j))*ops%neg_alive(k,j))+(z1d(k+norb)*z2d(ops%dead(k,j))*ops%neg_dead(k,j))) 
            end do
            !!$omp end do simd
            haml_vals=haml_vals+(ov*el(j))
        end do
        !$omp end do simd
        !$omp end parallel 
        
        return 
      
    end function haml_vals

    ! Calcualates a column of a hamiltonian Start specifies the row the column
    ! is started to be calcualted 
    subroutine haml_column(hcol,z1d,zstore,an_cr,an2_cr2,elecs,start)

        implicit none
        real(kind=8),dimension(:),intent(inout)::hcol 
        type(zombiest),dimension(:),intent(in)::zstore
        real(kind=8),dimension(0:),intent(in)::z1d
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::start
        real(kind=8)::h1etot,h2etot
        integer::j
        
        !$omp parallel 
        !$omp single
        do j=1,(ndet-(start-1))
            !$omp task firstprivate(h1etot,j) shared(zstore,hcol,an_cr,an2_cr2,elecs,z1d,start)
            h1etot = haml_vals(z1d,zstore(start+j-1)%val,an_cr,elecs%h1ei,elecs%h1_num)
            !$omp atomic
            hcol(j)=hcol(j)+h1etot
            !$omp end atomic
            !$omp end task
            !$omp task firstprivate(h2etot,j) shared(zstore,hcol,an_cr,an2_cr2,elecs,z1d)
            h2etot = haml_vals(z1d,zstore(start+j-1)%val,an2_cr2,elecs%h2ei,elecs%h2_num)
            !$omp atomic
            hcol(j)=hcol(j)+(0.5*h2etot)
            !$omp end atomic
            !$omp end task
        end do 
        !$omp end single
        !$omp end parallel  
        return

    end subroutine haml_column

    ! hamliltonian calcualtion - calcualtes the whole hamliltonian 
    subroutine haml_make(haml,zstore,elecs,an_cr,an2_cr2,verb) 

        implicit none
        
        real(kind=8),dimension(:,:),intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::verb
        integer::j,ierr
    
        if (errorflag .ne. 0) return 
        ierr=0
    
        !$omp parallel do &
        !$omp & private(j) &
        !$omp & shared(elecs,zstore,an_cr,an2_cr2,haml) 
        do j=1,ndet
            call haml_column(haml(j:,j),zstore(j)%val,zstore,an_cr,an2_cr2,elecs,j) 
            if(verb.eq.1)then
                write(6,"(a,i0,a)") "hamliltonian column ",j, " completed"
            end if 
        end do
        !$omp end parallel do

        do j=1,ndet
            haml(j,:)=haml(:,j)
        end do
        
    end subroutine haml_make

    ! calculates individual overlaps where no creation and annihilation operations are needed
    real(kind=8) function overlap_1(z1d,z2d)

        implicit none
        real(kind=8),dimension(0:)::z1d,z2d
        integer::j
    
    
        overlap_1=1.0
        !$omp parallel do simd reduction(*:overlap_1)
        do j=1,norb
            overlap_1=overlap_1*((z1d(j)*z2d(j))+(z1d(j+norb)*z2d(norb+j)))
        end do
        !$omp end parallel do simd
        

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
        !$omp parallel do &
        !$omp & shared(z1d,zstore,ovrlp_column) &
        !$omp & private(j)
        do j=1,ndet
            if(j.ne.row)then 
                ovrlp_column(j)=overlap_1(z1d,zstore(j)%val)
            else
                ovrlp_column(j)=1.0
            end if 
        end do 
        !$omp end parallel do
        return

    end function ovrlp_column

    !subroutine calcualates whole overlap matrix
    subroutine ovrlp_make(ovrlp,zstore)

        implicit none
        
        real(kind=8),dimension(:,:),intent(inout)::ovrlp
        type(zombiest),dimension(:),intent(in)::zstore
        integer::j,k
    
        if (errorflag .ne. 0) return

        !$omp parallel do &
        !$omp & private(j,k) &
        !$omp & shared(zstore,ovrlp)
        do j=1,ndet
            do k=j,ndet
                if(k.ne.j)then 
                    ovrlp(j,k)=overlap_1(zstore(j)%val,zstore(k)%val); ovrlp(k,j)=ovrlp(j,k)
                else
                    ovrlp(j,k)=1.0
                end if 
            end do
        end do 
        !$omp end parallel do 

    end subroutine ovrlp_make

    ! Subroutine that controls and calcualtes all of the hamiltonian variables 
    subroutine hamgen(haml,zstore,elecs,size,an_cr,an2_cr2,an2_cr2_diff,verb)

        implicit none 

        type(hamiltonian), intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2,an2_cr2_diff
        integer,intent(in)::size,verb
        integer, allocatable,dimension(:)::IPIV1
        real(kind=8),allocatable,dimension(:)::WORK1
        integer::ierr,j
       
        if (errorflag .ne. 0) return
        ierr=0

        call ovrlp_make(haml%ovrlp,zstore)
       
        haml%hjk=haml%ovrlp*elecs%hnuc 
        call haml_make(haml%hjk,zstore,elecs,an_cr,an2_cr2,verb)
        
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

        if(GDflg.eq.'y')then
            ! do j=2,ndet
                ! call gradient_zs(haml,zstore,elecs,an_cr,an2_cr2,j)
                call gradient_zs(haml,zstore,elecs,an_cr,an2_cr2,an2_cr2_diff,2)
            ! end do 
        end if

    end subroutine hamgen

    ! subroutine that finds the gradient w.r.t to a specfic zombie state element 
    ! for an entire column 
    function ovrlp_column_grad(z1d,zstore,state)

        implicit none
        
        type(zombiest),dimension(:),intent(in)::zstore
        real(kind=8),dimension(0:),intent(in)::z1d
        integer,intent(in)::state
        real(kind=8),dimension(ndet)::ovrlp_column_grad
        integer::j
        
        ovrlp_column_grad=0.0
        !$omp parallel do &
        !$omp & shared(z1d,zstore,ovrlp_column_grad,state) &
        !$omp & private(j)
        do j=1,ndet
            if(j.ne.state)then 
                ovrlp_column_grad(j)=overlap_1(z1d,zstore(j)%val)
            else
                ovrlp_column_grad(j)=0.0
            end if 
        end do 
        !$omp end parallel do
        return

    end function ovrlp_column_grad

    ! subroutine that finds the gradient of the overlap w.r.t a specified state
    subroutine ovrlp_make_grad(zstore,state,ovrlp_grad)

        implicit none 
        type(zombiest),dimension(:),intent(in)::zstore
        real(kind=8),dimension(:,:),intent(inout)::ovrlp_grad
        integer,intent(in)::state
        real(kind=8),dimension(0:2*norb)::z1d
        integer::j 

        !$omp parallel do shared(zstore,state,ovrlp_grad) private(z1d)
        do j=1, norb
            z1d=zstore(state)%val
            z1d(j)=zstore(state)%cos(j)
            z1d(j+norb)=(-1)*zstore(state)%sin(j)
            ovrlp_grad(j,:)=ovrlp_column_grad(z1d,zstore,state)
        end do
        !$omp end parallel do 

        return 

    end subroutine ovrlp_make_grad

    

    ! calculates indvidual hamliltonian elements taking in two Zombie states and a set of 
    ! creation and annihilation operations
    real(kind=8) function haml_val_grad(z1d,z2d,ops,el,orb)

        implicit none 
        real(kind=8),dimension(0:),intent(in)::z1d,z2d
        real(kind=8),dimension(:),intent(in)::el
        integer,intent(in)::orb
        type(oprts),intent(in)::ops
        real(kind=8)::ov
        integer::j,k,len

        len=ops%dcnt(0,orb)
        haml_val_grad=0.0
        !$omp parallel private(j,k,ov) shared(ops,z1d,z2d)
        !$omp do simd reduction(+:haml_val_grad) 
        do j=1,len
            ov=1.0
            !!$omp do simd reduction(*:ov)
            do k=1, norb
                ov=ov*((z1d(k)*z2d(ops%alive_diff(orb,k,(ops%dcnt(j,orb))))*ops%neg_alive_diff(orb,k,(ops%dcnt(j,orb))))&
                +(z1d(k+norb)*z2d(ops%dead_diff(orb,k,(ops%dcnt(j,orb))))*ops%neg_dead_diff(orb,k,(ops%dcnt(j,orb))))) 
            end do
            !!$omp end do simd
            haml_val_grad=haml_val_grad+(ov*el(ops%dcnt(j,orb)))
        end do
        !$omp end do simd
        !$omp end parallel 
        
        return 
      
    end function haml_val_grad

    real(kind=8) function haml_gvals(z1d,z2d,ops,ops2,el,len,orb)

        implicit none 
        real(kind=8),dimension(0:),intent(in)::z1d,z2d
        real(kind=8),dimension(:),intent(in)::el
        integer,intent(in)::len,orb
        type(oprts),intent(in)::ops,ops2
        real(kind=8)::ov
        integer::j,k

        
        haml_gvals=0.0
        !$omp parallel private(j,k,ov) shared(ops,z1d,z2d)
        !$omp do simd reduction(+:haml_gvals) 
        do j=1,len
            ov=1.0
            !!$omp do simd reduction(*:ov)
            do k=1, norb
                ov=ov*((z1d(ops2%alive(k,orb))*ops2%neg_alive(k,orb)*z2d(ops%alive(k,j))*&
                ops%neg_alive(k,j))+(z1d(ops2%dead(k,orb))*ops2%neg_dead(k,orb)*z2d(ops%dead(k,j))*ops%neg_dead(k,j)))
                
                ! ov=ov*((z1d(ops%alive(k,j+len))*ops%neg_alive(k,j+len)*z2d(ops%alive(k,j))*&
                ! ops%neg_alive(k,j))+(z1d(ops%dead(k,j+len))*ops%neg_dead(k,j+len)*z2d(ops%dead(k,j))*ops%neg_dead(k,j))) 
            end do
            
            haml_gvals=haml_gvals+(ov*el(j))
        end do
        !$omp end do simd
        !$omp end parallel 
        
        return 
      
    end function haml_gvals

    ! Calcualates a column of a hamliltonian Start specifies the row the column
    ! is started to be calcualted 
    subroutine haml_grad_rc(hcol,z1d,zstore,an_cr,an2_cr2,an2_cr2_diff,elecs,state,orb)

        implicit none
        real(kind=8),dimension(:),intent(inout)::hcol 
        type(zombiest),dimension(:),intent(in)::zstore
        real(kind=8),dimension(0:),intent(in)::z1d
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2,an2_cr2_diff
        integer,intent(in)::state,orb
        real(kind=8)::h1etot,h2etot
        integer::j
        

        !$omp parallel 
        !$omp single
        do j=1,ndet
            if(j.ne.state)then
                ! Differentiating the bra 1 el
                !$omp task firstprivate(h1etot,j) shared(hcol,zstore,an_cr,elecs,z1d)
                h1etot = haml_vals(z1d,zstore(j)%val,an_cr,elecs%h1ei,elecs%h1_num)
                !$omp atomic
                hcol(j)=hcol(j)+h1etot
                !$omp end atomic
                !$omp end task
               
                !Differentiating the bra 2 el
                !$omp task firstprivate(h2etot,j) shared(hcol,zstore,an2_cr2,elecs,z1d)
                h2etot = haml_vals(z1d,zstore(j)%val,an2_cr2,elecs%h2ei,elecs%h2_num)
                !$omp atomic
                hcol(j)=hcol(j)+(0.5*h2etot)
                !$omp end atomic
                !$omp end task

            else
                !Differentiaitn hamiltonian element (a,a) only placed in hamiltonian column
                !$omp task firstprivate(h1etot,j) shared(hcol,zstore,an_cr,elecs)
                h1etot = haml_val_grad(zstore(state)%val,zstore(state)%val,an_cr,elecs%h1ei,orb)
                !$omp atomic
                hcol(j)=hcol(j)+h1etot
                !$omp end atomic
                !$omp end task
                !$omp task firstprivate(h2etot,j) shared(hcol,zstore,an2_cr2,elecs)
                h2etot = haml_val_grad(zstore(state)%val,zstore(state)%val,an2_cr2,elecs%h2ei,orb)
                !$omp atomic
                hcol(j)=hcol(j)+(0.5*h2etot)
                !$omp end atomic
                !$omp end task
            end if
        end do 
        !$omp end single
        !$omp end parallel  
        return

    end subroutine haml_grad_rc

    ! Hamiltonian calcualtion - calcualtes the gradient of ther hamliltonian w.r.t one zombie state 
    subroutine haml_grad(haml_diff,zstore,elecs,an_cr,an2_cr2,an2_cr2_diff,state) 

        implicit none
        
        real(kind=8),dimension(:,:),intent(inout)::haml_diff 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2,an2_cr2_diff
        integer,intent(in)::state
        real(kind=8),dimension(0:2*norb)::z1d
  
        integer::j,ierr
    
        if (errorflag .ne. 0) return 
        ierr=0

       
        !$omp parallel do &
        !$omp & private(j) &
        !$omp & shared(elecs,zstore,an_cr,an2_cr2,haml_diff)
        do j=1,norb
            z1d(0:2*norb)=zstore(state)%val(0:2*norb)
            z1d(j)=zstore(state)%cos(j)
            z1d(j+norb)=zstore(state)%sin(j)*(-1)
            call haml_grad_rc(haml_diff(:,j),z1d,zstore,an_cr,an2_cr2,an2_cr2_diff,elecs,state,j)
        end do
        !$omp end parallel do

        
        
    end subroutine haml_grad


    !subroutine to calcualte gradient of w.r.t to a specific zombie state
    subroutine gradient_zs(haml,zstore,elecs,an_cr,an2_cr2,an2_cr2_diff,state)

        implicit none

        type(hamiltonian), intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2,an2_cr2_diff
        integer,intent(in)::state
        real(kind=8),allocatable,dimension(:,:,:)::temp2 
        integer::ierr,k,l,j,p
        ierr=0

        call ovrlp_make_grad(zstore,state,haml%diff_ovrlp(state,:,:))
        
        ! haml%diff_hjk(state,:,:)=0
        haml%diff_hjk(state,:,:)=haml%diff_ovrlp(state,:,:)*elecs%hnuc  
        call haml_grad(haml%diff_hjk(state,:,:),zstore,elecs,an_cr,an2_cr2,an2_cr2_diff,state)
      
        ! do j=1,norb
        !     print*,haml%diff_hjk(state,j,:)
        ! end do 
        ! print*,"***************************"
       
        
        allocate(temp2(ndet,norb,ndet),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in occupancy diff_inverse temp array allocation . ierr had value ", ierr
            errorflag=1
            return
        end if
        temp2=0.0
    
        do j=1,norb
            do k=1,ndet
                do l=1, ndet
                    if(l.ne.state)then
                        haml%diff_ov_dov(state,k,j,l)=(haml%ovrlp(k,state)*&
                        haml%diff_ovrlp(state,j,l))/abs(haml%ovrlp(k,state))

                        haml%diff_in_dhjk(state,k,j,l)=(haml%inv(k,state)*haml%diff_hjk(state,j,l))

                        temp2(k,j,l)=(haml%inv(k,l))*haml%diff_ovrlp(state,j,l)
                    else
                        do p=1,ndet
                            haml%diff_ov_dov(state,k,j,l)=haml%diff_ov_dov(state,k,j,l)+&
                            (haml%ovrlp(k,p)*haml%diff_ovrlp(state,j,p))/abs(haml%ovrlp(k,p))

                            haml%diff_in_dhjk(state,k,j,l)=haml%diff_in_dhjk(state,k,j,l)+&
                            (haml%inv(k,p)*haml%diff_hjk(state,j,p))

                            temp2(k,j,l)= temp2(k,p,l)+(haml%inv(k,p)*haml%diff_ovrlp(state,j,p))
                        end do 
                    end if 
                end do
            end do 
    
            do k=1, ndet
                do l=1, ndet
                    do p=1,ndet
                        haml%diff_invh(state,k,j,l)=haml%diff_invh(state,k,j,l)+(temp2(k,j,p)*haml%kinvh(p,l))
                    end do 
                    haml%diff_invh(state,k,j,l)=haml%diff_invh(state,k,j,l)*(-1)
                end do  
            end do
        end do


        deallocate(temp2)
        
        return

    end subroutine gradient_zs


END MODULE ham
