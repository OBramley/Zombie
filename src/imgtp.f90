MODULE imgtp
    use globvars
    use alarrays
    use grad_d
    contains

    ! Routine for imaginary time propagation
    subroutine imgtime_prop(dvecs,en,haml,diff_state,orb)

        implicit none

        type(dvector),dimension(:),intent(inout)::dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state,orb
        integer::j,k,states
        ! real(kind=8)::p
        real::db,r
        !DOUBLE PRECISION, external::ZBQLU01,ZBQLUAB

        if (errorflag .ne. 0) return
    
        do j=1,size(dvecs)
            dvecs(j)%d=0.0
            if(imagflg=='n') then
                dvecs(j)%d(j)=(1.0,0.0)
                if(zst=='HF') then
                    do k=1, ndet
                    ! k=int(ZBQLUAB(1,ndet))
                        call random_number(r)
                        ! dvecs(j)%d(k)=cmplx(p,0.0,kind=8)
                        dvecs(j)%d(k)=r
                    end do
                end if
                
            else if(imagflg=='y') then
                dvecs(j)%d(j)=1.0
                ! dvecs(j)%d(j)=(1.0,1.0)
             end if
        end do
       
       
        states=1
        if(gramflg.eq."y")then
            states=gramnum+1
            call gs(dvecs,haml,diff_state,orb)
        else
            call d_norm(dvecs(1),haml,0,diff_state,orb)
        end if 
        db=beta/timesteps
       
    
        do j=1,timesteps+1
            en%t(j)=db*(j-1)
            !!$omp parallel shared(en,j,haml,dvecs) private(k)
            !!$omp do
            do k=1,states
                en%erg(k,j)=ergcalc(haml%hjk,dvecs(k)%d)
                call timestep(haml,dvecs(k),db,diff_state,orb)
            end do
            !!$omp end do
            !!$omp end parallel
            if(gramflg.eq."y")then
                call gs(dvecs,haml,diff_state,orb)
            else
                call d_norm(dvecs(1),haml,j,diff_state,orb)
            end if
   
        end do

        return

    end subroutine imgtime_prop

    ! Calculates the energy
    real(kind=8) function ergcalc(bham,dvec)

        implicit none

        ! complex(kind=8),intent(in),dimension(:)::dvec
        ! complex(kind=8),intent(in),dimension(:,:)::bham
        real(kind=8),intent(in),dimension(:)::dvec
        real(kind=8),intent(in),dimension(:,:)::bham
        ! complex(kind=8)::result
        real(kind=8)::result
        
        if (errorflag .ne. 0) return
        !$omp parallel
        !$omp workshare
        result=dot_product(dvec,matmul(bham,dvec))
        ergcalc=result
        !$omp end workshare
        !$omp end parallel
   
        return
       
    end function ergcalc

    subroutine d_norm(dvec,haml,step,diff_state,orb)

        implicit none
        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::step,diff_state,orb
        real(kind=8)::norm

       
        !$omp parallel 
        !$omp workshare
        norm=abs(dot_product((dvec%d),matmul(haml%ovrlp,(dvec%d))))
        norm=sqrt(norm)
        !$omp end workshare
        !$omp end parallel
       
        dvec%norm=norm
   
        if(GDflg.eq.'y')then
            call d_normalise_diff(dvec,haml,step,diff_state,orb)
        end if
        
        dvec%d=dvec%d/norm
        return
    
    end subroutine d_norm


    ! Takes one timestep
    subroutine timestep(haml,dvec,db,diff_state,orb) 

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        real,intent(in)::db
        integer,intent(in)::diff_state,orb
        ! complex(kind=8),dimension(ndet)::ddot
        real(kind=8),dimension(ndet)::ddot

   

        if (errorflag .ne. 0) return
        if(GDflg.eq.'y')then
            call timestep_diff(dvec,haml,db,diff_state,orb)
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
            call d_norm(dvecs_copy(j),haml,1,diff_state,orb)  
        end do
        dvecs=dvecs_copy
        call deallocdv(dvecs_copy)

        return

    end subroutine gs


    subroutine imaginary_time_prop2(dvecs,en,haml,diff_state,orb)

        implicit none

        type(dvector),dimension(:),intent(inout)::dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state,orb
        integer::j,k,l
        real::db
        real(kind=8)::norm,result,temp
        real(kind=8),dimension(ndet)::ddot


        if (errorflag .ne. 0) return

        dvecs(1)%d=0.0
        dvecs(1)%d(1)=1.0
        !$acc enter data create(norm,db,result,temp)
        norm=0
       
        !$acc parallel loop gang vector reduction(+:norm) present(haml,dvecs) private (temp,l)   
        do j=1,ndet
            temp=0
            !$acc loop reduction(+:temp)
            do l=1,ndet 
                temp=temp+haml%ovrlp(j,l)*dvecs(1)%d(l)
            end do 
            norm=norm+(temp*dvecs(1)%d(j))
        end do 
       !$acc end parallel loop
        norm = sqrt(abs(norm))
       
        dvecs(1)%norm=norm
       
        call d_normalise_diff(dvecs(1),haml,0,diff_state,orb)
      
        
        dvecs(1)%d=dvecs(1)%d/norm
       
        db=beta/timesteps
       
        do k=1,timesteps+1

            en%t(k)=db*(k-1)
            result=0
            !$acc parallel loop gang vector reduction(+:result) present(haml,dvecs) private(temp,l)   
            do j=1,ndet
                temp=0
                !$acc loop reduction(+:temp)
                do l=1,ndet 
                    temp=temp+haml%hjk(j,l)*dvecs(1)%d(l)
                end do 
                result = result + (dvecs(1)%d(j)*temp)
            end do
            !$acc end parallel loop 
            en%erg(1,k)=result
    
            ddot=0
            call timestep_diff(dvecs(1),haml,db,diff_state,orb)
          
            !$acc parallel loop gang vector present(haml,dvecs) private(temp,l)
            do j=1,ndet 
                temp=0
                 !$acc loop reduction(+:temp)
                do l=1,ndet 
                    temp= temp + haml%kinvh(j,l)*dvecs(1)%d(l)
                end do 
                ddot(j)=temp
            end do
            !$acc end parallel loop 
            dvecs(1)%d=dvecs(1)%d-(db*ddot)
        
            norm=0   
            !$acc parallel loop gang vector reduction(+:norm) present(haml,dvecs) private(temp,l)     
            do j=1,ndet
                temp=0
                !$acc loop reduction(+:temp)
                do l=1,ndet 
                    temp=temp+haml%ovrlp(j,l)*dvecs(1)%d(l)
                end do 
                norm=norm+(temp*dvecs(1)%d(j))
            end do 
            !$acc end parallel loop 
            norm = sqrt(abs(norm))
            dvecs(1)%norm=norm
        
            call d_normalise_diff(dvecs(1),haml,k,diff_state,orb)
               
            dvecs(1)%d=dvecs(1)%d/norm

        end do
        !$acc exit data delete(norm,db,result)

        return

    end subroutine imaginary_time_prop2



END MODULE imgtp