MODULE imgtp
    use globvars
    use alarrays

    contains

    ! Routine for imaginary time propagation
    subroutine imgtime_prop(dvecs,en,ham)

        implicit none

        type(dvector),dimension(:),intent(inout)::dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(in)::ham
        integer::j,k,states
        real(kind=8)::p
        real::db
        DOUBLE PRECISION, external::ZBQLU01,ZBQLUAB

        if (errorflag .ne. 0) return

        do j=1,size(dvecs)
            if(imagflg=='n') then
                dvecs(j)%d(j)=(1.0,0.0)
                if(zst=='HF') then
                    do k=1, ndet
                    ! k=int(ZBQLUAB(1,ndet))
                        p=ZBQLU01(1)
                        dvecs(j)%d(k)=cmplx(p,0.0,kind=8)
                    end do
                end if
            else if(imagflg=='y') then
                dvecs(j)%d(j)=(1.0,1.0)
             end if
        end do

        states=1
        if(gramflg.eq."y")then
            states=gramnum+1
            call gs(dvecs,ham%ovrlp)
        end if

        db=beta/timesteps

        dvecs(1)%d=(dvecs(1)%d)/sqrt(zabs(dot_product((dvecs(1)%d),matmul(ham%ovrlp,(dvecs(1)%d)))))
    
        do j=1,timesteps+1
            en%t(j)=db*(j-1)
            !$omp parallel shared(en,j,ham,dvecs) private(k)
            !$omp do
            do k=1,states
                en%erg(k,j)=ergcalc(ham%hjk,dvecs(k)%d)
                call timestep(ham%kinvh,ham%ovrlp,dvecs(k),db)
            end do
            !$omp end do
            !$omp end parallel
            if(gramflg.eq."y")then
                call gs(dvecs,ham%ovrlp)
            end if
        end do

        return

    end subroutine imgtime_prop

    ! Calculates the energy
    complex(kind=8) function ergcalc(bham,dvec)

        implicit none

        complex(kind=8),intent(in),dimension(:)::dvec
        complex(kind=8),intent(in),dimension(:,:)::bham
        complex(kind=8)::result
      
        
        if (errorflag .ne. 0) return
        
        !$omp parallel
        !$omp workshare
        result=dot_product(dvec,matmul(bham,dvec))
        ergcalc=result
        !$omp end workshare
        !$omp end parallel
        return

    end function ergcalc

    ! Takes one timestep
    subroutine timestep(kinvh,kover,dvecs,db)

        implicit none

        type(dvector),intent(inout)::dvecs
        complex(kind=8),intent(in),dimension(:,:)::kinvh,kover
        real,intent(in)::db
        complex(kind=8),dimension(ndet)::ddot,temp
        real(kind=8)::norm

        if (errorflag .ne. 0) return

        !$omp parallel 
        !$omp workshare
        ddot= -matmul((kinvh),(dvecs%d))
        dvecs%d=dvecs%d+(db*ddot)
        temp=matmul(kover,(dvecs%d))
        norm=zabs(dot_product((dvecs%d),temp))
        norm=1/sqrt(norm)
        dvecs%d=norm*dvecs%d
        !$omp end workshare
        !$omp end parallel

        return

    end subroutine timestep

    !Gram-Schmidt orthogonalisation 
    subroutine gs(dvecs,kover)

        implicit none
        type(dvector), intent(inout),dimension(:)::dvecs
        complex(kind=8),intent(in),dimension(:,:)::kover
        type(dvector), allocatable,dimension(:)::dvecs_copy
        complex(kind=8)::numer,den,norm
        complex(kind=8),dimension(ndet)::temp
        integer::states,j,k

        if (errorflag .ne. 0) return
    
        states = gramnum+1
        call allocdv(dvecs_copy,states,ndet)
        
        do j=1, states
            dvecs_copy(j)%d(:)=dvecs(j)%d(:)
        end do

        do j=2,states
            do k=1, j-1
                temp=matmul(kover,dvecs_copy(k)%d)
                numer = dot_product(dvecs_copy(j)%d,temp)
                den  = dot_product(dvecs_copy(k)%d,temp)
                dvecs_copy(j)%d = dvecs_copy(j)%d - (dvecs_copy(k)%d*(numer/den))
            end do
        end do

        do j=1,states
            temp=matmul(kover,dvecs_copy(j)%d)
            norm = dot_product(dvecs_copy(j)%d,temp)
            norm = 1/sqrt(norm)
            dvecs_copy(j)%d = dvecs_copy(j)%d*norm
            dvecs(j)%d(:)=dvecs_copy(j)%d(:)
        end do

        call deallocdv(dvecs_copy)

        return

    end subroutine gs



END MODULE imgtp