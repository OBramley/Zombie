MODULE alarrays

    use globvars


    contains

    ! Routine to allcoate 1&2 electron electron integral matrices
    subroutine allocintgrl(elecs)

        implicit none

        type(elecintrgl),allocatable, intent(inout)::elecs
        
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
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in Zombie state allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        zs%alive(1:norb)=(0.0d0,0.0d0)
        zs%dead(1:norb)=(0.0d0,0.0d0)

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

        type(zombiest),dimension(:,:),allocatable,intent(inout)::zstore

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

    subroutine allocham(ham,size)

        implicit none

        type(hamiltonian),intent(inout)::ham 

        integer::ierr,size

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

    end subroutine deallocham

    subroutine allocdv(dvecs,x,length)

        implicit none

        type(dvector),intent(inout),allocatable,dimension(:)::dvecs
        integer, intent(in)::x,length
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
        end do
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

    
        deallocate(dvecs,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in dvecs deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if
  

       
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



END MODULE alarrays

        
        
