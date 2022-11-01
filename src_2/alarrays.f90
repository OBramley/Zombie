MODULE alarrays

    use globvars


    contains


    ! Routine to allcoate 1&2 electron electron integral matrices
    subroutine allocintgrl(elecs)

        implicit none

        type(elecintrgl), intent(inout)::elecs
        
        integer::ierr

        if (errorflag .ne. 0) return
        
        ierr=0
        allocate (elecs%h1ei(norb,norb), stat=ierr)
        if(ierr==0) allocate (elecs%h2ei(norb,norb,norb,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in electron integral  allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        
        elecs%h1ei(1:norb,1:norb)=0.0d0
        elecs%h2ei(1:norb,1:norb,1:norb,1:norb)=0.0d0
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
        integer, intent (in) :: nbf

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
        type(zombiest),intent(inout)::zs

        integer::ierr

        if (errorflag .ne. 0) return

        ierr = 0

        allocate(zs%alive(norb),stat=ierr)
        if(ierr==0) allocate(zs%dead(norb), stat=ierr)
        if(ierr==0) allocate(zs%sin(norb), stat=ierr)
        if(ierr==0) allocate(zs%cos(norb), stat=ierr)
        if(ierr==0) allocate(zs%phi(norb), stat=ierr)
        if(imagflg=='y')then 
            if(ierr==0) allocate(zs%img(norb), stat=ierr)
            zs%img(1:norb)=0.0
        end if
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in Zombie state allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        zs%alive(1:norb)=1
        zs%dead(1:norb)=1
        zs%phi(1:norb)=0.0d0
        zs%cos(1:norb)=(0.0d0,0.0d0)
        zs%sin(1:norb)=(0.0d0,0.0d0)
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
        ham%hjk(1:size,1:size)=(0.0d0,0.0d0)
        ham%ovrlp(1:size,1:size)=(0.0d0,0.0d0)
        ham%inv(1:size,1:size)=(0.0d0,0.0d0)
        if(GDflg.eq.'y')then  
            if(ierr==0) allocate(ham%diff_hjk(size,size,diff_size), stat=ierr)
            if(ierr==0) allocate(ham%diff_ovrlp(size,size,diff_size), stat=ierr)
            if(ierr==0) allocate(ham%diff_invh(size,size,size,diff_size), stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in GD Hamiltonian allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
            ham%diff_hjk(1:size,1:size,1:diff_size)=0.0
            ham%diff_ovrlp(1:size,1:size,1:diff_size)=0.0
            ham%diff_invh(1:size,1:size,1:size,1:diff_size)=0.0
        end if

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
            if(ierr==0) deallocate(ham%diff_invh, stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in GD Hamiltonian deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if

        return

    end subroutine deallocham

    subroutine allocdv(dvecs,x,length,diff_length)

        implicit none

        type(dvector),intent(inout),allocatable,dimension(:)::dvecs
        integer, intent(in)::x,length,diff_length
        integer::j,ierr

        if (errorflag .ne. 0) return

        ierr=0
        if(allocated(dvecs).eqv..false.)then
            allocate(dvecs(x),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in dvecs allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if

        do j=1, x
            allocate(dvecs(j)%d(length),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in d allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
            dvecs(j)%d(1:length)=(0.0,0.0)
            dvecs(j)%norm = 0.0d0
        end do

        if(GDflg.eq.'y')then
            allocate(dvecs(1)%d_diff(length,length,diff_length))
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in d_diff allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
            dvecs(1)%d_diff(1:length,1:length,1:diff_length)=0.0
        end if

        
     
        return

    end subroutine allocdv

    subroutine deallocdv(dvecs)

        implicit none 
        type(dvector),intent(inout),allocatable,dimension(:)::dvecs
        integer::j,ierr

        if (errorflag .ne. 0) return

        ierr=0

        do j=1, size(dvecs)
            deallocate(dvecs(j)%d,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in d deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end do

        if(GDflg.eq.'y')then
            deallocate(dvecs(1)%d_diff)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in d_diff deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if
    
        deallocate(dvecs,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in dvecs deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        
        
        return
       
    end subroutine deallocdv


    subroutine allocerg(en,x)

        implicit none

        type(energy),intent(inout)::en
        integer, intent(in)::x
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0

        
        allocate(en%t(timesteps+1),stat=ierr)
        if(ierr==0) allocate(en%erg(x,timesteps+1))
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in energy allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        en%erg=(0.0,0.0)
        en%t=0.0
        return
     
    end subroutine allocerg

    subroutine deallocerg(en)

        implicit none

        type(energy),intent(inout)::en
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0
        
        deallocate(en%t,stat=ierr)
        if(ierr==0) deallocate(en%erg)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in energy deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return
     
    end subroutine deallocerg

    subroutine allocgrad(gradients,num,length)

        implicit none

        type(grad),intent(inout)::gradients
        integer, intent(in)::num,length
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0

        
        allocate(gradients%vars(num,length),stat=ierr)
        if (ierr==0)allocate(gradients%grad_avlb(num),stat=ierr)
        if (ierr==0)allocate(gradients%prev_mmntm(num,length),stat=ierr)
        
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in gradient matrix allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        gradients%prev_mmntm=0
        gradients%vars=0
        gradients%grad_avlb=0
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
        if (ierr==0)deallocate(gradients%prev_mmntm,stat=ierr)
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
            ham%diff_hjk(1:size,1:size,1:norb)=0.0
            ham%diff_ovrlp(1:size,1:size,1:norb)=0.0
            ham%diff_invh(1:size,1:size,1:size,1:norb)=0.0
            dvecs(1)%d_diff(1:size,1:size,1:norb)=0.0
            gradients%vars=0
        end if 

    end subroutine value_reset
   



END MODULE alarrays

        
        
