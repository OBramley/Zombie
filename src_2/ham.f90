MODULE ham 

    use globvars
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
            ! print*,(ov*ov*el(j))
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
    subroutine hamgen(haml,zstore,elecs,size,an_cr,an2_cr2,verb)

        implicit none 

        type(hamiltonian), intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::size,verb
        ! real(kind=8),allocatable,dimension(:,:)::alive,dead
        integer, allocatable,dimension(:)::IPIV1
        real(kind=8),allocatable,dimension(:)::WORK1
        real(kind=8),dimension(0:2*norb)::z1d,z2d
        integer::ierr,j


        if (errorflag .ne. 0) return
        ierr=0

       
        z1d=zstore(7)%val
        z2d=zstore(5)%val
        print*, haml_vals(z1d,z2d,an_cr,elecs%h1ei,elecs%h1_num)
        print*, haml_vals(z2d,z1d,an_cr,elecs%h1ei,elecs%h1_num)
        print*, haml_vals(z1d,z2d,an2_cr2,elecs%h2ei,elecs%h2_num)
        print*, haml_vals(z2d,z1d,an2_cr2,elecs%h2ei,elecs%h2_num)
        print*,"********************"
        do j=1,norb
            z1d=zstore(7)%val
            z1d(j)=zstore(7)%cos(j)
            z1d(j+norb)=(-1)*zstore(7)%sin(j)
            print*, haml_vals(z1d,z2d,an_cr,elecs%h1ei,elecs%h1_num)
            print*, haml_vals(z2d,z1d,an_cr,elecs%h1ei,elecs%h1_num)
            print*, haml_vals(z1d,z2d,an2_cr2,elecs%h2ei,elecs%h2_num)
            print*, haml_vals(z2d,z1d,an2_cr2,elecs%h2ei,elecs%h2_num)
        end do
        print*,"********************"
        z1d=zstore(7)%val
        z2d=zstore(5)%val
        print*, haml_vals(z1d,z2d,an_cr,elecs%h1ei,elecs%h1_num)
        print*, haml_vals(z2d,z1d,an_cr,elecs%h1ei,elecs%h1_num)
        print*, haml_vals(z1d,z2d,an2_cr2,elecs%h2ei,elecs%h2_num)
        print*, haml_vals(z2d,z1d,an2_cr2,elecs%h2ei,elecs%h2_num)
        print*,"********************"
        do j=1,norb
            z2d=zstore(5)%val
            z2d(j)=zstore(5)%cos(j)
            z2d(j+norb)=(-1)*zstore(5)%sin(j)
            print*, haml_vals(z1d,z2d,an_cr,elecs%h1ei,elecs%h1_num)
            print*, haml_vals(z2d,z1d,an_cr,elecs%h1ei,elecs%h1_num)
            print*, haml_vals(z1d,z2d,an2_cr2,elecs%h2ei,elecs%h2_num)
            print*, haml_vals(z2d,z1d,an2_cr2,elecs%h2ei,elecs%h2_num)
        end do    


        stop
       
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
            call gradient_zs(haml,zstore,elecs,an_cr,an2_cr2,2)
        end if

        ! deallocate(alive,dead,stat=ierr)

    end subroutine hamgen


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

    ! Calcualates a column of a hamliltonian Start specifies the row the column
    ! is started to be calcualted 
    subroutine haml_column_grad(hcol,z1d,zstore,an_cr,an2_cr2,elecs,state,orb)

        implicit none
        real(kind=8),dimension(:),intent(inout)::hcol 
        type(zombiest),dimension(:),intent(in)::zstore
        real(kind=8),dimension(0:),intent(in)::z1d
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::state,orb
        real(kind=8)::h1etot,h2etot
        integer::j
       
        !$omp parallel 
        !$omp single
        do j=1,ndet
            if(j.ne.state)then
                ! !$omp task firstprivate(h1etot,j) shared(hcol,zstore,an_cr,elecs,z1d)
                ! h1etot = haml_vals(z1d,zstore(j)%val,an_cr,elecs%h1ei,elecs%h1_num)
                ! !$omp atomic
                ! hcol(j)=hcol(j)+h1etot
                ! !$omp end atomic
                ! !$omp end task
                !$omp task firstprivate(h2etot,j) shared(hcol,zstore,an2_cr2,elecs,z1d)
                h2etot = haml_vals(z1d,zstore(j)%val,an2_cr2,elecs%h2ei,elecs%h2_num)
                !$omp atomic
                hcol(j)=hcol(j)+(0.5*h2etot)
                !$omp end atomic
                !$omp end task
            else 
                ! !$omp task firstprivate(h1etot,j) shared(hcol,zstore,an_cr,elecs)
                ! h1etot = haml_val_grad(zstore(state)%val,zstore(state)%val,an_cr,elecs%h1ei,orb)
                ! !$omp atomic
                ! hcol(j)=hcol(j)+h1etot
                ! !$omp end atomic
                ! !$omp end task
                ! !$omp task firstprivate(h2etot,j) shared(hcol,zstore,an2_cr2,elecs)
                ! h2etot = haml_val_grad(zstore(state)%val,zstore(state)%val,an2_cr2,elecs%h2ei,orb)
                ! !$omp atomic
                ! hcol(j)=hcol(j)+(0.5*h2etot)
                ! !$omp end atomic
                ! !$omp end task
            end if
        end do 
        !$omp end single
        !$omp end parallel  
        return

    end subroutine haml_column_grad

    ! Hamiltonian calcualtion - calcualtes the gradient of ther hamliltonian w.r.t one zombie state 
    subroutine haml_grad(haml_diff,zstore,elecs,an_cr,an2_cr2,state) 

        implicit none
        
        real(kind=8),dimension(:,:),intent(inout)::haml_diff 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::state
        real(kind=8),dimension(0:2*norb)::z1d
  
        integer::j,ierr
    
        if (errorflag .ne. 0) return 
        ierr=0
    
        !$omp parallel do &
        !$omp & private(j) &
        !$omp & shared(elecs,zstore,an_cr,an2_cr2,haml_diff) 
        do j=1,norb
            z1d=zstore(state)%val
            z1d(j)=zstore(state)%cos(j)
            z1d(j+norb)=zstore(state)%sin(j)*(-1)
            call haml_column_grad(haml_diff(j,:),z1d,zstore,an_cr,an2_cr2,elecs,state,j)
        end do
        !$omp end parallel do

        
        
    end subroutine haml_grad


    !subroutine to calcualte gradient of w.r.t to a specific zombie state
    subroutine gradient_zs(haml,zstore,elecs,an_cr,an2_cr2,state)

        implicit none

        type(hamiltonian), intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::state
        real(kind=8),allocatable,dimension(:,:,:)::temp2 
        integer::ierr,k,l,j,p
        ierr=0

        call ovrlp_make_grad(zstore,state,haml%diff_ovrlp(state,:,:))
        
        haml%diff_hjk(state,:,:)=0
        ! haml%diff_hjk(state,:,:)=haml%diff_ovrlp(state,:,:)*elecs%hnuc 
        call haml_grad(haml%diff_hjk(state,:,:),zstore,elecs,an_cr,an2_cr2,state)
        do j=1, ndet
        print*,haml%diff_hjk(state,:,j)
        end do
        
        allocate(temp2(ndet,norb,ndet),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in occupancy diff_inverse temp array allocation . ierr had value ", ierr
            errorflag=1
            return
        end if
        temp2=0
        
        do k=1, ndet
            do l=1, ndet
                if(l.eq.state)then
                    do j=1,norb
                        do p=1,ndet
                            temp2(k,j,state)= temp2(k,p,state)+(haml%inv(k,p)*haml%diff_ovrlp(state,j,p))
                        end do 
                    end do
                else
                    temp2(k,:,l)=(haml%inv(k,l))*haml%diff_ovrlp(state,:,l)
                end if
            end do
        end do

        do k=1, ndet
            do l=1, ndet
                do j=1,norb
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
