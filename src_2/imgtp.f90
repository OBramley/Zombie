MODULE imgtp
    use globvars


    contains


    subroutine imgtime_prop(dvecs,en,ham)

        implicit none

        type(dvector),dimension(:),intent(inout)::dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(in)::ham
        integer::j,k,states
        real::db

        do j=1,size(dvecs)
            dvecs(j)%d(j)=(1.0,1.0)
        end do

        states=1
        if(gramflg.eq."y")then
            states=gramnum+1
            call gs(dvecs,ham%ovrlp,states)
        end if

        db=beta/timesteps
        do j=1,timesteps+1
            do k=1,states
                en%t(j)=db*(j-1)
                en%erg(j,k)=ergcalc(ham%hjk,dvecs(k)%d)
            end do
            call timestep(ham%kinvh,dvecs,db,states)
            
            if(gramflg.eq."y")then
                call gs(dvecs,ham%ovrlp,states)
            end if
        
        end do

        return

    end subroutine imgtime_prop

    complex(kind=8) function ergcalc(bham,dvec)

        implicit none
        complex(kind=8),intent(in),dimension(:)::dvec
        complex(kind=8),intent(in),dimension(:,:)::bham
        complex(kind=8),dimension(ndet)::temp
        

        temp=matmul(bham,dvec)
        ergcalc=dot_product(dvec,temp)

        return

    end function ergcalc

    subroutine timestep(kinvh,dvecs,db,states)

        implicit none
        type(dvector),intent(inout),dimension(:)::dvecs
        complex(kind=8),intent(in),dimension(:,:)::kinvh
        integer,intent(in)::states
        real,intent(in)::db
        complex(kind=8),dimension(ndet)::ddot
        integer::j

        do j=1,states
            ddot= -matmul(kinvh,dvecs(j)%d)
            dvecs(j)%d=dvecs(j)%d+(db*ddot)
        end do

        return

    end subroutine timestep

    subroutine gs(dvecs,kinv,states)

        implicit none
        type(dvector), intent(inout),dimension(:)::dvecs
        complex(kind=8),intent(in),dimension(:,:)::kinv
        integer, intent(in)::states

        ! NEED TO WRITE GS ROUTINE

    end subroutine gs



END MODULE imgtp