MODULE alarrays

    use globvars


    contains


    ! Routine to allcoate 1&2 electron electron integral matrices
    subroutine allocintgrl(elecs,e1,e2)

        implicit none

        type(elecintrgl), intent(inout)::elecs
        integer,intent(in)::e1,e2
        integer::ierr

        if (errorflag .ne. 0) return
        
        ierr=0
        elecs%h1_num=e1
        elecs%h2_num=e2

        allocate (elecs%h1ei(e1), stat=ierr)
        if(ierr==0) allocate (elecs%h2ei(e2),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in electron integral  allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        ! allocate (elecs%h1ei(norb,norb), stat=ierr)
        ! if(ierr==0) allocate (elecs%h2ei(norb,norb,norb,norb),stat=ierr)
        ! if (ierr/=0) then
        !     write(0,"(a,i0)") "Error in electron integral  allocation. ierr had value ", ierr
        !     errorflag=1
        !     return
        ! end if
        
        elecs%h1ei=0.0d0
        elecs%h2ei=0.0d0
        ! elecs%h1ei(1:norb,1:norb)=0.0d0
        ! elecs%h2ei(1:norb,1:norb,1:norb,1:norb)=0.0d0
        elecs%hnuc= 0.0d0
        
        return
    end subroutine allocintgrl

    ! Routine to deallcoate 1&2 electron electron integral matrices
    subroutine deallocintgrl(elecs)

        implicit none

        type(elecintrgl),intent(inout)::elecs

        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0

        deallocate (elecs%h1ei, stat=ierr)
        if(ierr==0) deallocate (elecs%h2ei,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in electron integral  deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return

    end subroutine deallocintgrl

    ! Routine to allocate set of zombie states

    subroutine alloczs(zstore,nbf)
        implicit none

        type(zombiest),dimension(:),allocatable,intent(inout)::zstore
        integer(kind=16), intent (in) :: nbf

        integer::j,ierr

        if (errorflag .ne. 0) return

        ierr=0

        if(allocated(zstore).eqv..false.)then
            allocate(zstore(nbf),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in zombie set allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if

        do j=1,nbf
            call alloczf(zstore(j))
        end do

        return 

    end subroutine alloczs

    !  Routine to allocate individual zombie states

    subroutine alloczf(zs)
        
        implicit none
        
        type(zombiest),intent(inout)::zs

        integer::ierr

        if (errorflag .ne. 0) return

        ierr = 0

        allocate(zs%alive(norb),stat=ierr)
        if(ierr==0) allocate(zs%dead(norb), stat=ierr)
        if(ierr==0) allocate(zs%sin(norb), stat=ierr)
        if(ierr==0) allocate(zs%cos(norb), stat=ierr)
        if(ierr==0) allocate(zs%phi(norb), stat=ierr)
        if(ierr==0) allocate(zs%val(0:(2*norb)), stat=ierr)
        if(imagflg=='y')then 
            if(ierr==0) allocate(zs%img(norb), stat=ierr)
            zs%img(1:norb)=0.0
        end if
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in Zombie state allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        zs%val=0.0
        zs%alive(1:norb)=1
        zs%dead(1:norb)=1
        zs%phi(1:norb)=0.0d0
        ! zs%cos(1:norb)=(0.0d0,0.0d0)
        ! zs%sin(1:norb)=(0.0d0,0.0d0)
        zs%cos(1:norb)=0.0d0
        zs%sin(1:norb)=0.0d0
        zs%update_num=0
        return
    end subroutine alloczf

    subroutine dealloczs(zstore)
        implicit none

        type(zombiest),dimension(:),allocatable,intent(inout)::zstore

        integer::j,ierr

        if (errorflag .ne. 0) return

        ierr=0

        do j=1,size(zstore)
            call dealloczf(zstore(j))
        end do

        if(ierr==0) deallocate(zstore,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in zombie set deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if


        return 

    end subroutine dealloczs

    subroutine dealloczf(zs)
        type(zombiest),intent(inout)::zs

        integer::ierr

        if (errorflag .ne. 0) return

        ierr = 0

        deallocate(zs%alive,stat=ierr)
        if(ierr==0) deallocate(zs%dead, stat=ierr)
        if(ierr==0) deallocate(zs%sin, stat=ierr)
        if(ierr==0) deallocate(zs%cos, stat=ierr)
        if(ierr==0) deallocate(zs%phi, stat=ierr)
        if(ierr==0) deallocate(zs%val, stat=ierr)
        if(imagflg=='y')then 
            if(ierr==0) deallocate(zs%img, stat=ierr)
            zs%img(1:norb)=0.0
        end if

        if (ierr/=0) then
            write(0,"(a,i0)") "Error in Zombie state deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        
        return
    end subroutine dealloczf

    subroutine alloczs2d(zstore,nbf)
        implicit none

        type(zombiest),dimension(:,:),allocatable,intent(inout)::zstore
        integer, intent (in) :: nbf

        integer::j,k,ierr

        if (errorflag .ne. 0) return

        ierr=0

        if(allocated(zstore).eqv..false.)then
            allocate(zstore(nbf,nbf),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in zombie set allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if
        
        do j=1,nbf
            do k=1, nbf
                call alloczf(zstore(j,k))
            end do
        end do

        return 

    end subroutine alloczs2d

    subroutine dealloczs2d(zstore)
        implicit none

        type(zombiest),dimension(:,:), allocatable, intent(inout)::zstore

        integer::j,k,ierr

        if (errorflag .ne. 0) return

        ierr=0

        do j=1, size(zstore, dim=1)
            do k=1, size(zstore, dim=1)
                call dealloczf(zstore(j,k))
            end do
        end do

       
        deallocate(zstore,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in zombie set allocation. ierr had value ", ierr
            errorflag=1
            return
        end if


        return 

    end subroutine dealloczs2d

    subroutine allocham(ham,size,diff_size)

        implicit none

        type(hamiltonian),intent(inout)::ham 
        integer,intent(in)::size,diff_size
        integer::ierr

        if (errorflag .ne. 0) return

        ierr = 0

        allocate(ham%hjk(size,size), stat=ierr)
        if(ierr==0) allocate(ham%ovrlp(size,size), stat=ierr)
        if(ierr==0) allocate(ham%inv(size,size), stat=ierr)
        if(ierr==0) allocate(ham%kinvh(size,size), stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in Hamiltonian allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        ! ham%hjk(1:size,1:size)=(0.0d0,0.0d0)
        ! ham%ovrlp(1:size,1:size)=(0.0d0,0.0d0)
        ! ham%inv(1:size,1:size)=(0.0d0,0.0d0)
        ham%hjk(1:size,1:size)=0.0d0
        ham%ovrlp(1:size,1:size)=0.0d0
        ham%inv(1:size,1:size)=0.0d0
        if(GDflg.eq.'y')then  
            if(ierr==0) allocate(ham%diff_hjk(size,diff_size,size), stat=ierr)
            if(ierr==0) allocate(ham%diff_ovrlp(size,diff_size,size), stat=ierr)
            ! if(ierr==0) allocate(ham%hess_hjk(size,diff_size,diff_size,size), stat=ierr)
            ! if(ierr==0) allocate(ham%hess_ovrlp(size,diff_size,diff_size,size), stat=ierr)
            if(ierr==0) allocate(ham%diff_invh(size,size,diff_size,size), stat=ierr)
            if(ierr==0) allocate(ham%diff_ov_dov(size,size,diff_size,size), stat=ierr)
            if(ierr==0) allocate(ham%diff_in_dhjk(size,size,diff_size,size), stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in GD Hamiltonian allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
            ham%diff_hjk=0.0
            ham%diff_ovrlp=0.0
            ! ham%hess_hjk=0.0
            ! ham%hess_ovrlp=0.0
            ham%diff_invh=0.0
            ham%diff_ov_dov=0.0
            ham%diff_in_dhjk=0.0
        end if
        ham%gram_num=0
        return
    end subroutine allocham

    subroutine deallocham(ham)

        implicit none

        type(hamiltonian), intent(inout)::ham 

        integer::ierr

        if (errorflag .ne. 0) return

        ierr = 0

        deallocate(ham%Hjk, stat=ierr)
        if(ierr==0) deallocate(ham%ovrlp, stat=ierr)
        if(ierr==0) deallocate(ham%inv, stat=ierr)
        if(ierr==0) deallocate(ham%kinvh, stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in Hamiltonian deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        if(GDflg.eq.'y')then  
            if(ierr==0) deallocate(ham%diff_hjk, stat=ierr)
            if(ierr==0) deallocate(ham%diff_ovrlp, stat=ierr)
            ! if(ierr==0) deallocate(ham%hess_hjk, stat=ierr)
            ! if(ierr==0) deallocate(ham%hess_ovrlp, stat=ierr)
            if(ierr==0) deallocate(ham%diff_invh, stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in GD Hamiltonian deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if

        return

    end subroutine deallocham

    subroutine allocdv(dvecs,length,diff_length)

        implicit none

        type(dvector),intent(inout)::dvecs
        integer, intent(in)::x,length,diff_length
        integer::j,ierr

        if (errorflag .ne. 0) return

        ierr=0
    
        
        allocate(dvecs%d(length),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in d allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        dvecs%d(1:length)=0.0
        dvecs%norm = 0.0d0
        

        if(GDflg.eq.'y')then
            allocate(dvecs%d_diff(length,length,diff_length))
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in d_diff allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
            dvecs%d_diff=0.0
        end if

        return

    end subroutine allocdv

    ! subroutine allocdv(dvecs,x,length,diff_length)

    !     implicit none

    !     type(dvector),intent(inout),allocatable,dimension(:)::dvecs
    !     integer, intent(in)::x,length,diff_length
    !     integer::j,ierr

    !     if (errorflag .ne. 0) return

    !     ierr=0
    !     if(allocated(dvecs).eqv..false.)then
    !         allocate(dvecs(x),stat=ierr)
    !         if (ierr/=0) then
    !             write(0,"(a,i0)") "Error in dvecs allocation. ierr had value ", ierr
    !             errorflag=1
    !             return
    !         end if
    !     end if

    !     do j=1, x
    !         allocate(dvecs(j)%d(length),stat=ierr)
    !         if (ierr/=0) then
    !             write(0,"(a,i0)") "Error in d allocation. ierr had value ", ierr
    !             errorflag=1
    !             return
    !         end if
    !         ! dvecs(j)%d(1:length)=(0.0,0.0)
    !         dvecs(j)%d(1:length)=0.0
    !         dvecs(j)%norm = 0.0d0
    !     end do

    !     if(GDflg.eq.'y')then
    !         allocate(dvecs(1)%d_diff(length,length,diff_length))
    !         if (ierr/=0) then
    !             write(0,"(a,i0)") "Error in d_diff allocation. ierr had value ", ierr
    !             errorflag=1
    !             return
    !         end if
    !         dvecs(1)%d_diff=0.0
    !     end if

    !     return

    ! end subroutine allocdv

    subroutine deallocdv(dvecs)

        implicit none 
        type(dvector),intent(inout),dimension(:)::dvecs
        integer::j,ierr

        if (errorflag .ne. 0) return

        ierr=0

       
        deallocate(dvecs%d,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in d deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        

        if(GDflg.eq.'y')then
            deallocate(dvecs%d_diff)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in d_diff deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if
    
        return
       
    end subroutine deallocdv

    ! subroutine deallocdv(dvecs)

    !     implicit none 
    !     type(dvector),intent(inout),allocatable,dimension(:)::dvecs
    !     integer::j,ierr

    !     if (errorflag .ne. 0) return

    !     ierr=0

    !     do j=1, size(dvecs)
    !         deallocate(dvecs(j)%d,stat=ierr)
    !         if (ierr/=0) then
    !             write(0,"(a,i0)") "Error in d deallocation. ierr had value ", ierr
    !             errorflag=1
    !             return
    !         end if
    !     end do

    !     if(GDflg.eq.'y')then
    !         deallocate(dvecs(1)%d_diff)
    !         if (ierr/=0) then
    !             write(0,"(a,i0)") "Error in d_diff deallocation. ierr had value ", ierr
    !             errorflag=1
    !             return
    !         end if
    !     end if
    
    !     deallocate(dvecs,stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in dvecs deallocation. ierr had value ", ierr
    !         errorflag=1
    !         return
    !     end if

        
        
    !     return
       
    ! end subroutine deallocdv

    subroutine allocgrad(gradients,num,length)

        implicit none

        type(grad),intent(inout)::gradients
        integer, intent(in)::num,length
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0

        
        allocate(gradients%vars(num,length),stat=ierr)
        ! if (ierr==0) allocate(gradients%vars_hess(num,length),stat=ierr)
        if (ierr==0)allocate(gradients%grad_avlb(0:num,num),stat=ierr)
        ! if (ierr==0)allocate(gradients%hessian(num,length,length),stat=ierr)
        ! if (ierr==0)allocate(gradients%prev_mmntm(num,length),stat=ierr)
        ! if (ierr==0)allocate(gradients%hess_sum(num),stat=ierr)
        
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in gradient matrix allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        ! gradients%prev_mmntm=0
        gradients%vars=0
        gradients%grad_avlb=0
        ! gradients%hessian=0
        return
     
    end subroutine allocgrad
   
    subroutine deallocgrad(gradients)

        implicit none

        type(grad),intent(inout)::gradients
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0

        
        deallocate(gradients%vars,stat=ierr)
        if (ierr==0)deallocate(gradients%grad_avlb,stat=ierr)
        ! if (ierr==0) deallocate(gradients%vars_hess,stat=ierr)
        ! if (ierr==0) deallocate(gradients%hessian,stat=ierr)
        ! if (ierr==0)deallocate(gradients%prev_mmntm,stat=ierr)
        ! if (ierr==0)deallocate(gradients%hess_sum,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in gradient matrix deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return
     
    end subroutine deallocgrad


    subroutine value_reset(ham,dvecs,en,size,gradients)

        implicit none
        type(hamiltonian),intent(inout)::ham
        type(dvector), dimension(:),intent(inout):: dvecs
        type(energy),intent(inout):: en
        type(grad),intent(inout)::gradients
        integer,intent(in)::size
        integer::j

        ham%hjk(1:size,1:size)=(0.0d0,0.0d0)
        ham%ovrlp(1:size,1:size)=(0.0d0,0.0d0)
        ham%inv(1:size,1:size)=(0.0d0,0.0d0)
        en%erg=(0.0,0.0)
        if(gramflg.eq."n")then
            dvecs(1)%d(1:size)=(0.0,0.0)
            dvecs(1)%norm = 0.0d0
        else
            do j=1, 1+gramnum
                dvecs(j)%d(1:size)=(0.0,0.0)
                dvecs(j)%norm = 0.0d0
            end do
        end if
        if(GDflg.eq.'y')then 
            ham%diff_hjk=0.0
            ham%diff_ovrlp=0.0
            ham%diff_invh=0.0
            dvecs(1)%d_diff=0.0
            gradients%vars=0
        end if 

    end subroutine value_reset
    
    subroutine alloc_oprts(oper,n)

        implicit none 
        type(oprts),intent(inout)::oper 
        integer,intent(in)::n 
        integer::ierr

        if (errorflag .ne. 0) return
      
        ierr=0
       
        call alloc_oprts_2(oper%ham,n)
        
        if(GDflg.eq.'y')then 
            allocate(oper%diff(norb),stat=ierr)
            ! if(ierr==0) allocate(oper%hess(norb,norb),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in gradient operators allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if
       
        return 

    end subroutine alloc_oprts

    subroutine alloc_oprts_2(oper,n) 

        implicit none 
        type(oprts_2),intent(inout)::oper 
        integer,intent(in)::n 
        integer::ierr,k
        integer(kind=2)::j,norbs

        if (errorflag .ne. 0) return

        ierr=0
        norbs=int(norb,kind=2)
        allocate(oper%alive(norb,n),oper%dead(norb,n),oper%neg_alive(norb,n),oper%neg_dead(norb,n),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in operators deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        oper%neg_alive=1
        oper%neg_dead=1
        oper%alive=0
        oper%dead=0
      
        do k=1,n
            oper%alive(:,k)=[integer(kind=2)::(j,j=1,norbs)]
            oper%dead(:,k)=[integer(kind=2)::((j+norbs),j=1,norbs)]
        end do

        return

    end subroutine alloc_oprts_2

    subroutine dealloc_oprts_2(oper)

        implicit none 
        type(oprts_2),intent(inout)::oper 
        integer::ierr
        
        if (errorflag .ne. 0) return

        ierr=0

        deallocate(oper%alive,oper%dead,oper%neg_alive,oper%neg_dead,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in operators deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

    end subroutine

    subroutine dealloc_oprts(oper)

        implicit none 
        type(oprts),intent(inout)::oper 
        integer::ierr,j,k

        if (errorflag .ne. 0) return

        ierr=0

        call dealloc_oprts_2(oper%ham)
      
        if(GDflg.eq.'y')then 
            do j=1,norb 
                call dealloc_oprts_2(oper%diff(j))
                ! do k=1,norb 
                !     call dealloc_oprts_2(oper%hess(j,k))
                ! end do 
            end do 
            deallocate(oper%diff,oper%dcnt,stat=ierr)
            ! deallocate(oper%diff,oper%hess,oper%dcnt,oper%hcnt,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in gradient operators deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if
       


    end subroutine dealloc_oprts

    subroutine grad_new_alloc(grad_fin,num,e1,e2)

        implicit none
        type(grad),intent(inout)::grad_fin
        integer,intent(in)::e1,e2,num
        integer::ierr
        
        ierr=0
        allocate(grad_fin%one_elec(num,2,e1),stat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in grad_new_setup one_elec array allocation . ierr had value ", ierr
            errorflag=1
            return
        end if
        allocate(grad_fin%two_elec(num,2,e2),stat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in grad_new_setup two_elec array allocation . ierr had value ", ierr
            errorflag=1
            return
        end if
        allocate(grad_fin%ovrlp_div(num),stat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in grad_new_setup ovrlp_div array allocation . ierr had value ", ierr
            errorflag=1
            return
        end if

        return 

    end subroutine grad_new_alloc

    subroutine grad_new_dealloc(grad_fin)

        implicit none
        type(grad),intent(inout)::grad_fin
        integer::ierr
        
        ierr=0
        deallocate(grad_fin%one_elec,stat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in grad_new_setup one_elec array deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if
        deallocate(grad_fin%two_elec,stat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in grad_new_setup two_elec array deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if

        return 

    end subroutine grad_new_dealloc

    subroutine allocgram(ham,num,size)

        implicit none 
        type(hamiltonian),intent(inout)::ham 
        integer,intent(in)::num,size
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0
        ham%gram_num=num
        allocate(ham%gs_ovrlp(num,size,size),stat=ierr)
        if(ierr==0) allocate(ham%gs_kinvh(num,size,size),stat=ierr)
        if(ierr==0) allocate(ham%gs_ovrlp_self(num,size,size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in overlap gram allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return

    end subroutine allocgram

    subroutine deallocgram(ham)

        implicit none 
        type(hamiltonian),intent(inout)::ham 
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0
        deallocate(ham%gs_ovrlp,stat=ierr)
        if(ierr==0) deallocate(ham%gs_kinvh,stat=ierr)
        if(ierr==0) deallocate(ham%gs_ovrlp_self,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in overlap gram allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return

    end subroutine deallocgram

    subroutine allocdvgram(dvec,num,size)

        implicit none 
        type(dvector),intent(inout)::dvec
        integer,intent(in)::num,size
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0

        allocate(dvec%d_gs(num,size),stat=ierr)
        if(ierr/=0)then 
            write(0,"(a,i0)") "Error in dvector gram allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return

    end subroutine allocdvgram

    subroutine deallocdvgram(dvec)

        implicit none 
        type(dvector),intent(inout)::dvec
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0

        deallocate(dvec%d_gs,stat=ierr)
        if(ierr/=0)then 
            write(0,"(a,i0)") "Error in dvector gram deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return

    end subroutine deallocdvgram


END MODULE alarrays

        
        
