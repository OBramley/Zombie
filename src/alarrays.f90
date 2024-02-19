MODULE alarrays

    use mod_types
    use globvars

    contains

    subroutine allocintgrl(elecs,n)

        implicit none

        type(elecintrgl),intent(inout)::elecs
        integer,intent(in)::n
        integer::ierr=0

        if (errorflag .ne. 0) return

        allocate (elecs%integrals(n), stat=ierr)
        allocate (elecs%orbital_choice(norb,n), stat=ierr)
        allocate (elecs%orbital_choice2(0:norb,n), stat=ierr)
        allocate (elecs%orbital_choice3(norb), stat=ierr)

        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in electron integral  deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return

    end subroutine allocintgrl
    ! Routine to deallcoate 1&2 electron electron integral matrices
    subroutine deallocintgrl(elecs)

        implicit none

        type(elecintrgl),intent(inout)::elecs
        integer::ierr=0

        if (errorflag .ne. 0) return

       
        deallocate (elecs%integrals, stat=ierr)
        deallocate (elecs%orbital_choice, stat=ierr)
        deallocate (elecs%orbital_choice2, stat=ierr)
        deallocate (elecs%orbital_choice3, stat=ierr)

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
        integer::ierr=0
        integer::j

        if (errorflag .ne. 0) return

       
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
        integer::ierr=0

        if (errorflag .ne. 0) return

       
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
        zs%phi=0.0d0
        return
    end subroutine alloczf

    subroutine dealloczs(zstore)
        implicit none

        type(zombiest),dimension(:),allocatable,intent(inout)::zstore
        integer::ierr=0
        integer::j

        if (errorflag .ne. 0) return

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
        integer::ierr=0

        if (errorflag .ne. 0) return

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

   
    subroutine allocham(ham,size)

        implicit none

        type(hamiltonian),intent(inout)::ham 
        integer,intent(in)::size
        integer::ierr=0

        if (errorflag .ne. 0) return

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
        return
    end subroutine allocham

    subroutine deallocham(ham)

        implicit none

        type(hamiltonian), intent(inout)::ham 
        integer::ierr=0

        if (errorflag .ne. 0) return

        deallocate(ham%hjk, stat=ierr)
        if(ierr==0) deallocate(ham%ovrlp, stat=ierr)
        if(ierr==0) deallocate(ham%inv, stat=ierr)
        if(ierr==0) deallocate(ham%kinvh, stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in Hamiltonian deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

    
        return

    end subroutine deallocham

    subroutine allocdv(dvecs,length)

        implicit none

        type(dvector),intent(inout)::dvecs
        integer, intent(in)::length
        integer::ierr=0

        if (errorflag .ne. 0) return

        allocate(dvecs%d(length),stat=ierr)
        allocate(dvecs%d_1(length),stat=ierr)
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
        integer::ierr=0

        if (errorflag .ne. 0) return

        deallocate(dvecs%d,stat=ierr)
        deallocate(dvecs%d_1,stat=ierr)
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
        integer::ierr=0

        if (errorflag .ne. 0) return

       
        allocate(gradients%vars(num,length),stat=ierr)
        if (ierr==0)allocate(gradients%grad_avlb(length,num),stat=ierr)
        if (ierr==0)allocate(gradients%ovrlp_grad(length,num,num),stat=ierr)
        if (ierr==0)allocate(gradients%ovrlp_grad_avlb(length,num,num),stat=ierr)
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
        integer::ierr=0

        if (errorflag .ne. 0) return
        
        deallocate(gradients%vars,stat=ierr)
        if (ierr==0)deallocate(gradients%grad_avlb,stat=ierr)
        if (ierr==0)deallocate(gradients%ovrlp_grad,stat=ierr)
        if (ierr==0)deallocate(gradients%ovrlp_grad_avlb,stat=ierr)
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
        integer::k
        integer(int16)::j,norbs
        integer::ierr=0

        if (errorflag .ne. 0) return
     
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
        integer::ierr=0
        
        if (errorflag .ne. 0) return

        deallocate(oper%alive,oper%dead,oper%neg_alive,oper%neg_dead,stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in operators deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

    end subroutine dealloc_oprts

    subroutine allocgram(gram_unit,num,size,length)

        implicit none 
        type(gram),intent(inout)::gram_unit 
        integer,intent(in)::num,size,length
        integer::ierr=0

        if (errorflag .ne. 0) return
        gram_unit%state_num=num
        call alloczs(gram_unit%zstore,size)
        call allocdv(gram_unit%dvec,size)
        call allocham(gram_unit%haml,size)
        if(GDflg=='y')then
            call allocgrad(gram_unit%grads,size,length)
        end if
        if(gram_unit%state_num>1)then
            allocate(gram_unit%wf_ovrlp(num-1,size,size),stat=ierr)
            if(ierr/=0)then 
                write(stderr,"(a,i0)") "Error in Wave function overlap array allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if
       
        return

    end subroutine allocgram

    subroutine deallocgram(gram_unit)

        implicit none 
        type(gram),intent(inout)::gram_unit 
        integer::ierr=0

        if (errorflag .ne. 0) return

        call dealloczs(gram_unit%zstore)
        call deallocdv(gram_unit%dvec)
        call deallocham(gram_unit%haml)
        if(GDflg=='y')then
            call deallocgrad(gram_unit%grads)
        end if
        if(gram_unit%state_num>1)then
            deallocate(gram_unit%wf_ovrlp,stat=ierr)
            if(ierr/=0)then 
                write(stderr,"(a,i0)") "Error in Wave function overlap array deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end if
        

        return

    end subroutine deallocgram

    subroutine allocdvgram(dvec,num,size)

        implicit none 
        type(dvector),intent(inout)::dvec
        integer,intent(in)::num,size
        integer::ierr=0

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

        deallocate(dvec%d_gs,stat=ierr)
        if(ierr/=0)then 
            write(stderr,"(a,i0)") "Error in dvector gram deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return

    end subroutine deallocdvgram

    subroutine alloc_grad_do(grads,size)

        implicit none
        type(grad_do),intent(inout)::grads
        integer,intent(in)::size
        integer::ierr=0

        if (errorflag .ne. 0) return

        allocate(grads%hjk(size,size), stat=ierr)
        if(ierr==0) allocate(grads%ovrlp(size,size), stat=ierr)
        if(ierr==0) allocate(grads%inv(size,size), stat=ierr)
        if(ierr==0) allocate(grads%kinvh(size,size), stat=ierr)
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
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in gradient descent type deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

    end subroutine dealloc_grad_do
    
END MODULE alarrays

        
        
