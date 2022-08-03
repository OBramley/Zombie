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

        if (errorflag .ne. 0) return

        do j=1,size(dvecs)
            if(imagflg=='n') then
                dvecs(j)%d(j)=(1.0,0.0)
        else if(imagflg=='y') then
            dvecs(j)%d(j)=(1.0,1.0)
        end if
        end do

        states=1
        if(gramflg.eq."y")then
            states=gramnum+1
            call gs(dvecs,ham%ovrlp,states)
        end if

        db=beta/timesteps
    
        do j=1,timesteps+1
            en%t(j)=db*(j-1)
            do k=1,states
                en%erg(k,j)=ergcalc(ham%hjk,dvecs(k)%d)
            end do
            call timestep(ham%kinvh,ham%ovrlp,dvecs,db,states)
            
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
        real(kind=8),dimension(ndet)::temp
        ! complex(kind=8)::result
        real(kind=8)::result
        ! complex(kind=8),dimension(ndet)::temp
        
        if (errorflag .ne. 0) return
        
        ! temp=matmul(REAL(dvec),REAL(bham))
        ! result=dot_product(temp,REAL(dvec))
        temp=matmul(REAL(bham),REAL(dvec))
        result=dot_product(REAL(dvec),temp)
        ergcalc=cmplx(result,0.0,kind=8)
    
        return

    end function ergcalc

    subroutine timestep(kinvh,kover,dvecs,db,states)

        implicit none
        type(dvector),intent(inout),dimension(:)::dvecs
        complex(kind=8),intent(in),dimension(:,:)::kinvh,kover
        integer,intent(in)::states
        real,intent(in)::db
        ! complex(kind=8),dimension(ndet)::ddot
        real(kind=8),dimension(ndet)::ddot
        real(kind=8),dimension(ndet)::temp
        real(kind=8)::norm
        integer::j

        do j=1,states
            ddot= -matmul(REAL(kinvh),REAL(dvecs(j)%d))
            dvecs(j)%d=dvecs(j)%d+cmplx(db*ddot,0.0,kind=8)
            temp=matmul(REAL(kover),REAL(dvecs(j)%d))
            norm=dabs(dot_product(REAL(dvecs(j)%d),temp))
            norm=1/sqrt(norm)
            dvecs(j)%d=norm*dvecs(j)%d
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