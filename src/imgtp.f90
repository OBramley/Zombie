MODULE imgtp
    use globvars
    use alarrays
    use grad_d
    contains

    ! Routine for imaginary time propagation
    subroutine imgtime_prop(dvecs,en,haml,diff_state,orb)

        implicit none

        type(dvector),intent(inout)::dvecs
        real(kind=8), dimension(:),intent(inout)::erg
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state,orb
        integer::j,k,states
        real::db,r
     
        if (errorflag .ne. 0) return
    
       
        dvecs%d=0.0
        if(imagflg=='n') then
            dvecs(j)%d(j)=(1.0,0.0)
            if(zst=='HF') then
                do k=1, ndet
                    call random_number(r)
                    dvecs(j)%d(k)=r
                end do
            end if 
        else if(imagflg=='y') then
            dvecs%d(1)=1.0
        end if
        
       
       
        states=1
        if(gramflg.eq."y")then
            states=gramnum+1
            call gs(dvecs,haml,diff_state,orb)
        else
            call d_norm(dvecs,haml,0,diff_state,orb)
        end if 
        db=beta/timesteps
       
    
        do j=1,timesteps+1     
            erg(j)=ergcalc(haml%hjk,dvecs%d)
            call timestep(haml,dvecs,db,diff_state,orb)
           
            if(gramflg.eq."y")then
                call gs(dvecs,haml,diff_state,orb)
            else
                call d_norm(dvecs,haml,j,diff_state,orb)
            end if
   
        end do

        return

    end subroutine imgtime_prop

    ! Calculates the energy
    real(kind=8) function ergcalc(bham,dvec)

        implicit none

        real(kind=8),intent(in),dimension(:)::dvec
        real(kind=8),intent(in),dimension(:,:)::bham
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

    subroutine imaginary_time_prop2(dvecs,erg,haml,size,diff_state,orb)

        implicit none

        type(dvector),intent(inout)::dvecs
        real(kind=8), dimension(:),intent(inout)::erg
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state,orb,size
        integer::j,k,l,g
        real(kind=8)::norm,result,temp,db
        real(kind=8),dimension(size)::ddot


        if (errorflag .ne. 0) return
        ! Imaginary time Propagation if no gram-schmidt orthogonalization is needed
        if(haml%gram_num.eq.0)then
            dvecs%d=0.0
            dvecs%d(1)=1.0
        
            norm=0
        
        
            do j=1,size
                temp=0
                do l=1,size 
                    temp=temp+haml%ovrlp(j,l)*dvecs%d(l)
                end do 
                norm=norm+(temp*dvecs%d(j))
            end do 
        
            norm = sqrt(abs(norm))
        
            dvecs%norm=norm
        
            call d_normalise_diff(dvecs,haml,0,diff_state,orb)
        
            
            dvecs%d=dvecs%d/norm
        
            db=beta/timesteps
        
            do k=1,timesteps+1
                result=0
            
                do j=1,size
                    temp=0
                
                    do l=1,size 
                        temp=temp+haml%hjk(j,l)*dvecs%d(l)
                    end do 
                    result = result + (dvecs%d(j)*temp)
                end do
            
                erg(k)=result
        
                ddot=0
                call timestep_diff(dvecs,haml,db,diff_state,orb)
            
            
                do j=1,size 
                    temp=0
            
                    do l=1,size 
                        temp= temp + haml%kinvh(j,l)*dvecs%d(l)
                    end do 
                    ddot(j)=temp
                end do
            
                dvecs%d=dvecs%d-(db*ddot)
            
                norm=0   
        
                do j=1,size
                    temp=0
            
                    do l=1,size 
                        temp=temp+haml%ovrlp(j,l)*dvecs%d(l)
                    end do 
                    norm=norm+(temp*dvecs%d(j))
                end do 
        
                norm = sqrt(abs(norm))
                dvecs%norm=norm
            
                call d_normalise_diff(dvecs,haml,k,diff_state,orb)
                
                dvecs%d=dvecs%d/norm

            end do
        else !Imaginary time propagation with GSO
            
            ! Set up gs dvectors
            do k=1, haml%gram_num
                dvecs%d_gs=0.0
                dvecs%d_gs(k,k)=1.0
            end do
            dvecs%d(haml%gram_num+1)=1.0
            !Orthogonalise
            call gs(dvecs,haml)
            
            ! Normalise vectors
            do k=1, haml%gram_num
                norm=0
                do j=1,size
                    temp=0
                    do l=1,size 
                        temp=temp+haml%gs_ovrlp(k,j,l)*dvecs%d_gs(k,l)
                    end do 
                    norm=norm+(temp*dvecs%d_gs(k,j))
                end do 
            
                norm = sqrt(abs(norm))    
                dvecs%d_gs(k,:)=dvecs%d_gs(k,:)/norm
    
            end do
            !Normalise dvector of interest
            do j=1,size
                temp=0
                do l=1,size 
                    temp=temp+haml%ovrlp(j,l)*dvecs%d(l)
                end do 
                norm=norm+(temp*dvecs%d(j))
            end do 
        
            norm = sqrt(abs(norm))
        
            dvecs%norm=norm
            
            !Propagate derivative of dvector
            call d_normalise_diff(dvecs,haml,0,diff_state,orb)
        
            db=beta/timesteps
           
            !Begin time steps
            do k=1,timesteps+1
                
                ! Calculate energy of system being looked at
                result=0
                do j=1,size
                    temp=0
                    do l=1,size 
                        temp=temp+haml%hjk(j,l)*dvecs%d(l)
                    end do 
                    result = result + (dvecs%d(j)*temp)
                end do
            
                erg(k)=result
                ! Take a time step for dvector and propagate derivative
                ddot=0
                call timestep_diff(dvecs,haml,db,diff_state,orb)
            
                do j=1,size 
                    temp=0
                    do l=1,size 
                        temp= temp + haml%kinvh(j,l)*dvecs%d(l)
                    end do 
                    ddot(j)=temp
                end do
            
                dvecs%d=dvecs%d-(db*ddot)

                ! Make gs vectors take time step
                do g=1, haml%gram_num
                    ddot=0
                    do j=1,size 
                        temp=0
                        do l=1,size 
                            temp= temp + haml%gs_kinvh(g,j,l)*dvecs%d_gs(g,l)
                        end do 
                        ddot(j)=temp
                    end do
                    dvecs%d_gs(g,:)=dvecs%d_gs(g,:)-(db*ddot)
                end do 

                ! Orthogonalise 
                call gs(dvecs,haml)
                
                !Normnalise gs dvectors
                do g=1, haml%gram_num
                    norm=0
                    do j=1,size
                        temp=0
                        do l=1,size 
                            temp=temp+haml%gs_ovrlp(g,j,l)*dvecs%d_gs(g,l)
                        end do 
                        norm=norm+(temp*dvecs%d_gs(g,j))
                    end do 
                    norm = sqrt(abs(norm))    
                    dvecs%d_gs(g,:)=dvecs%d_gs(g,:)/norm
                end do
                
                !Normalise dvector
                do j=1,size
                    temp=0
                    do l=1,size 
                        temp=temp+haml%ovrlp(j,l)*dvecs%d(l)
                    end do 
                    norm=norm+(temp*dvecs%d(j))
                end do 
            
                norm = sqrt(abs(norm))
            
                dvecs%norm=norm
                
                ! Propagate derivative 
                call d_normalise_diff(dvecs,haml,k,diff_state,orb)
    
               
            end do
        end if 
      
        return

    end subroutine imaginary_time_prop2



END MODULE imgtp