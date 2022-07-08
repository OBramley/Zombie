MODULE imgtp
    use globvars


    contains





    subroutine imgtime_prop(dvecs,en,ham)

        implicit none

        type(dvector),dimension(:),intent(inout)::dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(in)::ham
        integer::j,ierr
        real::db

        do j=1,size(dvecs)
            dvecs(j)%d(j)=(1.0,1.0)
        end do

        if(gramflg.eq."y")then
            call gs(dvecs,ham%ovrlp)
        end if

        db=beta/timesteps
        do j=1,timesteps+1
        
        end do




    end subroutine imgtime_prop



END MODULE imgtp