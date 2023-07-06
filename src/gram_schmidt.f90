MODULE gram_schmidt

    use globvars
    use alarrays

    contains





    subroutine gram_schmidt_control(elect,ndet,an_cr,an2_cr2,nbf)

        implicit none

        type(elecintrgl),intent(in)::elect
        type(oprts),intent(in)::an_cr,an2_cr2
        type(zombiest), dimension(:,:), allocatable:: zstore
        type(dvector),dimension(:),allocatable:: dvecs

        integer::ierr,j,k

        if(errorflag.ne.0) return
        ierr=0


        ! Allocate values 
        

        gramnum
        ! Check run conditions
        if(GDflg.eq."y")then


        else
            allocate(zstore(1,ndet),stat=ierr)
            do j=1,ndet
                call alloczf(zstore(1,j))
            end do           
            if(zomgflg=='y')then
                call genzf(zstore(1,:),ndet)
                do j=1,ndet
                    call zombiewriter(zstore(1,j),j,0)
                end do
                write(6,"(a)") "Zombie states generated"
            else if (zomgflg=='n') then

                
                call read_zombie(zstore)
                write(6,"(a)") "Zombie read in"
            end if
            call flush(6)
            call flush(0)
        end if 
        

        if(propflg=="y")then


            if(GDflg.eq."y")then


            end if 





            if(cleanflg=="n")then
                call deallocdv(dvecs)
                write(6,"(a)") "d-vector deallocated"
                !if(hamgflg=='y')then
                call dealloczs(zstore)
                write(6,"(a)") "Zombie states deallocated"
                call deallocintgrl(elect)
                write(6,"(a)") "Electron integrals deallocated"
                call dealloc_oprts(an_cr)
                call dealloc_oprts(an2_cr2)
                write(6,"(a)") "creation and annihilation operators deallocated"
                !end if
            end if

            call flush(6)
            call flush(0)


        else if((propflg=="n"))then
            if((cleanflg=="y").or.(cleanflg=="f"))then

            else 
                write(6,"(a)") "The program if here has done nothing except read in some values and then deallocate them"
                call deallocintgrl(elect)
                write(6,"(a)") "Electron integrals deallocated"
                call dealloc_oprts(an_cr)
                call dealloc_oprts(an2_cr2)
                write(6,"(a)") "creation and annihilation operators deallocated"
            end if 

        end if 
       


        write(stateno,"(i4.4)")j
                call dvec_writer(dvecs%d,ndet,j)
                call energywriter(erg,"energy_state_"//trim(stateno)//".csv",j)

    end subroutine gram_schmidt_control


    !Gram-Schmidt orthogonalisation 
    subroutine gs(dvecs,haml,diff_state,orb)

        implicit none
        type(dvector), intent(inout),dimension(:)::dvecs
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state,orb
        type(dvector), allocatable,dimension(:)::dvecs_copy
        ! complex(kind=8)::numer,den
        real(kind=8)::numer,den
        ! complex(kind=8),dimension(ndet)::temp
        integer::states,j,k

        if (errorflag .ne. 0) return
    
        states = gramnum+1
        call allocdv(dvecs_copy,ndet,norb)
        
        ! do j=1, states
        !     dvecs_copy(j)%d(:)=dvecs(j)%d(:)
        ! end do
        dvecs_copy=dvecs
        do j=2,states
            do k=1, j-1
                numer=dot_product(dvecs(j)%d,matmul(haml%ovrlp,dvecs_copy(k)%d))
                den=dot_product(dvecs_copy(k)%d,matmul(haml%ovrlp,dvecs_copy(k)%d))
                dvecs_copy(j)%d = dvecs_copy(j)%d - (dvecs_copy(k)%d*(numer/den))
            end do
        end do

        do j=1,states
            call d_norm(dvecs_copy(j),haml,1,diff_state,orb)  
        end do
        dvecs=dvecs_copy
        call deallocdv(dvecs_copy)

        return

    end subroutine gs




End Module gram_schmidt