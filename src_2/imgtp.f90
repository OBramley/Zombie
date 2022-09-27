MODULE imgtp
    use globvars
    use alarrays
    use grad_d
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
            call gs(dvecs,ham)
        else
            call d_norm(dvecs(1),ham,0)
        end if 

        db=beta/timesteps
       
        ! dvecs(1)%d=(dvecs(1)%d)/sqrt(zabs(dot_product((dvecs(1)%d),matmul(ham%ovrlp,(dvecs(1)%d)))))
    
        do j=1,timesteps+1
            en%t(j)=db*(j-1)
            !$omp parallel shared(en,j,ham,dvecs) private(k)
            !$omp do
            do k=1,states
                en%erg(k,j)=ergcalc(ham%hjk,dvecs(k)%d)
                call timestep(ham,dvecs(k),db)
            end do
            !$omp end do
            !$omp end parallel
            if(gramflg.eq."y")then
                call gs(dvecs,ham)
            else
                call d_norm(dvecs(1),ham,1)
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

    subroutine d_norm(dvec,ham,step)

        implicit none
        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        integer,intent(in)::step
        real(kind=8)::norm

        !$omp parallel 
        !$omp workshare
        norm=zabs(dot_product((dvec%d),matmul(ham%ovrlp,(dvec%d))))
        norm=sqrt(norm)
        !$omp end workshare
        !$omp end parallel

        dvec%norm=norm
        if(GDflg.eq.'y')then
            call d_normalise_diff(dvec,ham,step)
        end if
        
        dvec%d=dvec%d/norm
        return
    
    end subroutine d_norm


    ! Takes one timestep
    subroutine timestep(ham,dvec,db) 

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        real,intent(in)::db
        complex(kind=8),dimension(ndet)::ddot
   

        if (errorflag .ne. 0) return
        if(GDflg.eq.'y')then
            call timestep_diff(dvec,ham,db)
        end if
        !$omp parallel 
        !$omp workshare
        ddot= -matmul((ham%kinvh),(dvec%d))
        dvec%d=dvec%d+(db*ddot)
        !$omp end workshare
        !$omp end parallel

        return

    end subroutine timestep

    !Gram-Schmidt orthogonalisation 
    subroutine gs(dvecs,ham)

        implicit none
        type(dvector), intent(inout),dimension(:)::dvecs
        type(hamiltonian),intent(in)::ham
        type(dvector), allocatable,dimension(:)::dvecs_copy
        complex(kind=8)::numer,den
        complex(kind=8),dimension(ndet)::temp
        integer::states,j,k

        if (errorflag .ne. 0) return
    
        states = gramnum+1
        call allocdv(dvecs_copy,states,ndet,norb)
        
        do j=1, states
            dvecs_copy(j)%d(:)=dvecs(j)%d(:)
        end do

        do j=2,states
            do k=1, j-1
                temp=matmul(ham%ovrlp,dvecs_copy(k)%d)
                numer = dot_product(dvecs_copy(j)%d,temp)
                den  = dot_product(dvecs_copy(k)%d,temp)
                dvecs_copy(j)%d = dvecs_copy(j)%d - (dvecs_copy(k)%d*(numer/den))
            end do
        end do

        do j=1,states
            call d_norm(dvecs_copy(j),ham,1)
            ! temp=matmul(kover,dvecs_copy(j)%d)
            ! norm = dot_product(dvecs_copy(j)%d,temp)
            ! norm = 1/sqrt(norm)
            ! dvecs_copy(j)%d = dvecs_copy(j)%d*norm
            dvecs(j)%d=dvecs_copy(j)%d
        end do

        call deallocdv(dvecs_copy)

        return

    end subroutine gs



END MODULE imgtp