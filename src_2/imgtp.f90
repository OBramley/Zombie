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
            !$omp parallel shared(en,j,ham,dvecs) private(k)
            !$omp do
            do k=1,states
                en%erg(k,j)=ergcalc(ham%hjk,dvecs(k)%d)
            end do
            !$omp end do
            !$omp do
            do k=1,states
                call timestep(ham%kinvh,ham%ovrlp,dvecs(k),db)
            end do
            !$omp end do
            !$omp end parallel
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
        complex(kind=8)::result
      
        
        if (errorflag .ne. 0) return
        
        !$omp parallel
        !$omp workshare
        temp=matmul(bham,dvec)
        result=dot_product(dvec,temp)
        ergcalc=result
        !$omp end workshare
        !$omp end parallel
        return

    end function ergcalc

    subroutine timestep(kinvh,kover,dvecs,db)

        implicit none
        type(dvector),intent(inout)::dvecs
        complex(kind=8),intent(in),dimension(:,:)::kinvh,kover
        real,intent(in)::db
        complex(kind=8),dimension(ndet)::ddot,temp
        real(kind=8)::norm

        !$omp parallel 
        !$omp workshare
        ddot= -matmul((kinvh),(dvecs%d))
        dvecs%d=dvecs%d+db*ddot
        temp=matmul(kover,(dvecs%d))
        norm=zabs(dot_product((dvecs%d),temp))
        norm=1/sqrt(norm)
        dvecs%d=norm*dvecs%d
        !$omp end workshare
        !$omp end parallel

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