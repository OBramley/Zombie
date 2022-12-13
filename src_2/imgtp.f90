MODULE imgtp
    use globvars
    use alarrays
    use grad_d
    contains

    ! Routine for imaginary time propagation
    subroutine imgtime_prop(dvecs,en,haml,diff_state)

        implicit none

        type(dvector),dimension(:),intent(inout)::dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state
        integer::j,k,states
        real(kind=8)::p
        real::db,r
        !DOUBLE PRECISION, external::ZBQLU01,ZBQLUAB

        if (errorflag .ne. 0) return

        do j=1,size(dvecs)
            if(imagflg=='n') then
                dvecs(j)%d(j)=(1.0,0.0)
                if(zst=='HF') then
                    do k=1, ndet
                    ! k=int(ZBQLUAB(1,ndet))
                        call random_number(r)
                        p=r !ZBQLU01(1)
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
            call gs(dvecs,haml,diff_state)
        else
            call d_norm(dvecs(1),haml,0,diff_state)
        end if 
        db=beta/timesteps
       
    
        do j=1,timesteps+1
            en%t(j)=db*(j-1)
            !$omp parallel shared(en,j,haml,dvecs) private(k)
            !$omp do
            do k=1,states
                en%erg(k,j)=ergcalc(haml%hjk,dvecs(k)%d)
                call timestep(haml,dvecs(k),db,diff_state)
            end do
            !$omp end do
            !$omp end parallel
            if(gramflg.eq."y")then
                call gs(dvecs,haml,diff_state)
            else
                call d_norm(dvecs(1),haml,1,diff_state)
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

    subroutine d_norm(dvec,haml,step,diff_state)

        implicit none
        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::step,diff_state
        real(kind=8)::norm

       
        !$omp parallel 
        !$omp workshare
        norm=abs(dot_product((dvec%d),matmul(haml%ovrlp,(dvec%d))))
        norm=sqrt(norm)
        !$omp end workshare
        !$omp end parallel
       
        dvec%norm=norm
        if(GDflg.eq.'y')then
            call d_normalise_diff(dvec,haml,step,diff_state)
        end if
        
        dvec%d=dvec%d/norm
        return
    
    end subroutine d_norm


    ! Takes one timestep
    subroutine timestep(haml,dvec,db,diff_state) 

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        real,intent(in)::db
        integer,intent(in)::diff_state
        complex(kind=8),dimension(ndet)::ddot
   

        if (errorflag .ne. 0) return
        if(GDflg.eq.'y')then
            call timestep_diff(dvec,haml,db,diff_state)
        end if
   
        !$omp parallel 
        !$omp workshare
        ddot= -matmul((haml%kinvh),(dvec%d))
        dvec%d=dvec%d+(db*ddot)
        !$omp end workshare
        !$omp end parallel
     
        return

    end subroutine timestep

    !Gram-Schmidt orthogonalisation 
    subroutine gs(dvecs,haml,diff_state)

        implicit none
        type(dvector), intent(inout),dimension(:)::dvecs
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state
        type(dvector), allocatable,dimension(:)::dvecs_copy
        complex(kind=8)::numer,den
        ! complex(kind=8),dimension(ndet)::temp
        integer::states,j,k

        if (errorflag .ne. 0) return
    
        states = gramnum+1
        call allocdv(dvecs_copy,states,ndet,norb)
        
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
            call d_norm(dvecs_copy(j),haml,1,diff_state)
            ! temp=matmul(kover,dvecs_copy(j)%d)
            ! norm = dot_product(dvecs_copy(j)%d,temp)
            ! norm = 1/sqrt(norm)
            ! dvecs_copy(j)%d = dvecs_copy(j)%d*norm
            
        end do
        dvecs=dvecs_copy
        call deallocdv(dvecs_copy)

        return

    end subroutine gs



END MODULE imgtp