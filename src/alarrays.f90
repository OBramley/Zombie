MODULE alarrays

    use mod_types
    use globvars
    use dnad

    contains


    ! Routine to deallcoate 1&2 electron electron integral matrices
    subroutine deallocintgrl(elecs)

        implicit none

        type(elecintrgl),intent(inout)::elecs

        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0
        deallocate (elecs%integrals, stat=ierr)
        if(ierr==0) deallocate (elecs%alive,stat=ierr)
        if(ierr==0) deallocate (elecs%dead,stat=ierr)
        if(ierr==0) deallocate (elecs%neg_a,stat=ierr)
        if(ierr==0) deallocate (elecs%neg_d,stat=ierr)

        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in electron integral  deallocation. ierr had value ", ierr
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
                write(stderr,"(a,i0)") "Error in zombie set allocation. ierr had value ", ierr
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

       
        allocate(zs%phi(norb), stat=ierr)
        if(ierr==0) allocate(zs%val(0:(2*norb)), stat=ierr)
        ! if(imagflg=='y')then 
        !     if(ierr==0) allocate(zs%img(norb), stat=ierr)
        !     zs%img(1:norb)=0.0
        ! end if
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in Zombie state allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
       
        zs%val=0.0d0
        zs%phi(1:norb)=0.0d0
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
            write(stderr,"(a,i0)") "Error in zombie set deallocation. ierr had value ", ierr
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

       deallocate(zs%val, stat=ierr)
       deallocate(zs%phi, stat=ierr)
       
        ! if(imagflg=='y')then 
        !     if(ierr==0) deallocate(zs%img, stat=ierr)
        ! end if

        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in Zombie state deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        
        return
    end subroutine dealloczf

   
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
            write(stderr,"(a,i0)") "Error in Hamiltonian allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        
        ham%hjk(1:size,1:size)=0.0d0
        ham%ovrlp(1:size,1:size)=0.0d0
        ham%inv(1:size,1:size)=0.0d0
        if(GDflg.eq.'y')then  
            if(ierr==0) allocate(ham%diff_hjk(diff_size,size,size), stat=ierr)
            if(ierr==0) allocate(ham%diff_ovrlp(diff_size,size,size), stat=ierr)
            if (ierr/=0) then
                write(stderr,"(a,i0)") "Error in GD Hamiltonian allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
            ham%diff_hjk=0.0
            ham%diff_ovrlp=0.0
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

        

        deallocate(ham%hjk, stat=ierr)
        if(ierr==0) deallocate(ham%ovrlp, stat=ierr)
        if(ierr==0) deallocate(ham%inv, stat=ierr)
        if(ierr==0) deallocate(ham%kinvh, stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in Hamiltonian deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        if(GDflg.eq.'y')then  
            if(ierr==0) deallocate(ham%diff_hjk, stat=ierr)
            if(ierr==0) deallocate(ham%diff_ovrlp, stat=ierr)
            if (ierr/=0) then
                write(stderr,"(a,i0)") "Error in GD Hamiltonian deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if

        return

    end subroutine deallocham

    subroutine allocdv(dvecs,length)

        implicit none

        type(dvector),intent(inout)::dvecs
        integer, intent(in)::length
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0
    
        allocate(dvecs%d(length),stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in d allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        dvecs%d(1:length)=0.0d0
        dvecs%norm = 0.0d0
        

        return

    end subroutine allocdv

    subroutine deallocdv(dvecs)

        implicit none 
        type(dvector),intent(inout)::dvecs
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0
       
        deallocate(dvecs%d,stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in d deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        
        return
       
    end subroutine deallocdv



    subroutine allocgrad(gradients,num,length)

        implicit none

        type(grad),intent(inout)::gradients
        integer, intent(in)::num,length
        integer::ierr

        if (errorflag .ne. 0) return

        ierr=0

        
        allocate(gradients%vars(num,length),stat=ierr)
        if (ierr==0)allocate(gradients%grad_avlb(num),stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in gradient matrix allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
      
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
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in gradient matrix deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return
     
    end subroutine deallocgrad
    
    subroutine alloc_oprts(oper,n) 

        implicit none 
        type(oprts),intent(inout)::oper 
        integer,intent(in)::n 
        integer::ierr,k
        integer(int16)::j,norbs

        if (errorflag .ne. 0) return

        ierr=0
        norbs=int(norb,kind=int16)
        allocate(oper%alive(norb,n),oper%dead(norb,n),oper%neg_alive(norb,n),oper%neg_dead(norb,n),stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in operators deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        oper%neg_alive=1
        oper%neg_dead=1
        oper%alive=0
        oper%dead=0
      
        do k=1,n
            oper%alive(:,k)=[integer(kind=int16)::(j,j=1,norbs)]
            oper%dead(:,k)=[integer(kind=int16)::((j+norbs),j=1,norbs)]
        end do

        return

    end subroutine alloc_oprts

    subroutine dealloc_oprts(oper)

        implicit none 
        type(oprts),intent(inout)::oper 
        integer::ierr
        
        if (errorflag .ne. 0) return

        ierr=0

        deallocate(oper%alive,oper%dead,oper%neg_alive,oper%neg_dead,stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in operators deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

    end subroutine dealloc_oprts

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
            write(stderr,"(a,i0)") "Error in overlap gram allocation. ierr had value ", ierr
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
            write(stderr,"(a,i0)") "Error in overlap gram allocation. ierr had value ", ierr
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
            write(stderr,"(a,i0)") "Error in dvector gram allocation. ierr had value ", ierr
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
            write(stderr,"(a,i0)") "Error in dvector gram deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return

    end subroutine deallocdvgram

    subroutine alloc_grad_do(grads,size,diff_size)

        implicit none
        type(grad_do),intent(inout)::grads
        integer,intent(in)::size,diff_size
        integer::ierr=0

        if (errorflag .ne. 0) return

        allocate(grads%hjk(size,size), stat=ierr)
        if(ierr==0) allocate(grads%ovrlp(size,size), stat=ierr)
        if(ierr==0) allocate(grads%inv(size,size), stat=ierr)
        if(ierr==0) allocate(grads%kinvh(size,size), stat=ierr)
        if(ierr==0) allocate(grads%diff_hjk_1(diff_size,size), stat=ierr)
        if(ierr==0) allocate(grads%diff_ovrlp_1(diff_size,size), stat=ierr)
        if(ierr==0) allocate(grads%diff_hjk_2(diff_size,size), stat=ierr)
        if(ierr==0) allocate(grads%diff_ovrlp_2(diff_size,size), stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in gradient descent type allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        call alloczf(grads%zom)
        call allocdv(grads%dvec,size)

    end subroutine alloc_grad_do

    subroutine dealloc_grad_do(grads)

        implicit none
        type(grad_do),intent(inout)::grads
        integer::ierr=0

        if (errorflag .ne. 0) return

        call dealloczf(grads%zom)
        call deallocdv(grads%dvec)

        deallocate(grads%hjk, stat=ierr)
        if(ierr==0) deallocate(grads%ovrlp, stat=ierr)
        if(ierr==0) deallocate(grads%inv, stat=ierr)
        if(ierr==0) deallocate(grads%kinvh, stat=ierr)
        if(ierr==0) deallocate(grads%diff_hjk_1, stat=ierr)
        if(ierr==0) deallocate(grads%diff_ovrlp_1, stat=ierr)
        if(ierr==0) deallocate(grads%diff_hjk_2, stat=ierr)
        if(ierr==0) deallocate(grads%diff_ovrlp_2, stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in gradient descent type deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

    end subroutine dealloc_grad_do



END MODULE alarrays

        
        
