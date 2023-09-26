MODULE imgtp
    
    use mod_types
    use globvars
    use alarrays
   
    contains


    subroutine imaginary_time_erg(values,size)

        implicit none
        type(grad_do),intent(inout)::values
        integer,intent(in)::size
        integer::j,k,l
        real(wp)::norm,result,temp
        real(wp),dimension(size)::ddot
        real(wp)::db
        
        if (errorflag .ne. 0) return
    
        values%dvec%d_1=0.0d0
        values%dvec%d_1(1)=1.0d0
    
        db=beta/timesteps

        do k=1,timesteps+1
            norm=0.0d0

            do j=1,size
                temp=0.0d0
                do l=1,size 
                    temp=temp+values%ovrlp(j,l)*values%dvec%d_1(l)
                end do 
                norm=norm+(temp*values%dvec%d_1(j))
            end do 
        
            values%dvec%norm = sqrt(abs(norm))
        
            values%dvec%d=values%dvec%d_1/values%dvec%norm

            ddot=0.0d0
                
            do j=1,size 
                temp=0.0d0
                do l=1,size 
                    temp= temp + values%kinvh(j,l)*values%dvec%d(l)
                end do 
                ddot(j)=temp
            end do
        
            values%dvec%d_1=values%dvec%d-(db*ddot)
        end do

        norm=0.0d0

        do j=1,size
            temp=0.0d0
            do l=1,size 
                temp=temp+values%ovrlp(j,l)*values%dvec%d_1(l)
            end do 
            norm=norm+(temp*values%dvec%d_1(j))
        end do 
        values%dvec%d_o_d=sign_d_o_d(norm)
        values%dvec%norm = sqrt(abs(norm))
    
        values%dvec%d=values%dvec%d_1/values%dvec%norm

        result=0.0d0
        do j=1,size
            temp=0.0d0
            do l=1,size 
                temp=temp+values%hjk(j,l)*values%dvec%d(l)
            end do 
            result = result + (values%dvec%d(j)*temp)
        end do
    
        values%erg=result

    end subroutine  imaginary_time_erg

    ! Routine for imaginary time propagation
    ! subroutine imgtime_prop(dvecs,erg,haml)

    !     implicit none

    !     type(dvector),intent(inout)::dvecs
    !     type(dual), dimension(:),intent(inout)::erg
    !     type(hamiltonian),intent(in)::haml
    !     integer::j,k,states
    !     real(wp)::db,r
     
    !     if (errorflag .ne. 0) return
    
       
    !     dvecs%d=0.0d0
    !     if(imagflg=='n') then
    !         dvecs%d(j)=1.0d0
    !         if(zst=='HF') then
    !             do k=1, ndet
    !                 call random_number(r)
    !                 dvecs%d(k)=r
    !             end do
    !         end if 
    !     ! else if(imagflg=='y') then
    !     !     dvecs%d(1)=1.0
    !     end if
        
       
    !     states=1
    !     if(gramflg.eq."y")then
    !         states=gramnum+1
    !         ! call gs(dvecs,haml)
    !     else
    !         call d_norm(dvecs,haml)
    !     end if 
    !     db=beta/timesteps
       
    
    !     do j=1,timesteps+1     
    !         erg(j)=ergcalc(haml%hjk,dvecs%d)
    !         call timestep(haml,dvecs,db)
           
    !         if(gramflg.eq."y")then
    !             ! call gs(dvecs,haml)
    !         else
    !             call d_norm(dvecs,haml)
    !         end if
   
    !     end do

    !     return

    ! end subroutine imgtime_prop

    ! Calculates the energy
    function ergcalc(bham,dvec) result(result)

        implicit none

        real(wp),intent(in),dimension(:)::dvec
        real(wp),intent(in),dimension(:,:)::bham
        real(wp)::result
        real(wp)::temp
        integer::j,l
      
        result=0.0d0
            
        do j=1,ndet
            temp=0.0d0
            do l=1,ndet 
                temp=temp+bham(j,l)*dvec(l)
            end do 
            result = result + (dvec(j)*temp)
        end do

        ! !$omp parallel
        ! !$omp workshare
        ! result=dot_product(dvec,matmul(bham,dvec))
        ! !$omp end workshare
        ! !$omp end parallel
   
        return
       
    end function ergcalc

    ! subroutine d_norm(dvec,haml)

    !     implicit none
    !     type(dvector),intent(inout)::dvec
    !     type(hamiltonian),intent(in)::haml
    !     type(dual)::norm

       
    !     !$omp parallel 
    !     !$omp workshare
    !     norm=abs(dot_product((dvec%d),matmul(haml%ovrlp,(dvec%d))))
    !     dvec%norm=sqrt(norm)
    !     !$omp end workshare
    !     !$omp end parallel
       
    !     dvec%d=dvec%d/norm

    !     return
    
    ! end subroutine d_norm


    ! Takes one timestep
    ! subroutine timestep(haml,dvec,db) 

    !     implicit none

    !     type(dvector),intent(inout)::dvec
    !     type(hamiltonian),intent(in)::haml
    !     real(wp),intent(in)::db
    !     type(dual),dimension(ndet)::ddot

   

    !     if (errorflag .ne. 0) return
  
    !     !$omp parallel 
    !     !$omp workshare
    !     ddot= -matmul((haml%kinvh),(dvec%d))
    !     dvec%d=dvec%d+(db*ddot)
    !     !$omp end workshare
    !     !$omp end parallel
     
    !     return

    ! end subroutine timestep

    subroutine imaginary_time_prop2(dvecs,erg,haml,size)

        implicit none

        type(dvector),intent(inout)::dvecs
        real(wp), dimension(:),intent(inout)::erg
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::size
        integer::j,k,l!,g
        real(wp)::norm,result,temp
        real(wp)::db
        real(wp),dimension(size)::ddot


        if (errorflag .ne. 0) return
        ! Imaginary time Propagation if no gram-schmidt orthogonalization is needed
        if(haml%gram_num.eq.0)then
            dvecs%d_1=0.0d0
            dvecs%d_1(1)=1.0d0
        
            db=beta/timesteps
        
            do k=1,timesteps+1
                norm=0.0d0
        
        
                do j=1,size
                    temp=0.0d0
                    do l=1,size 
                        temp=temp+haml%ovrlp(j,l)*dvecs%d_1(l)
                    end do 
                    norm=norm+(temp*dvecs%d_1(j))
                end do 
                dvecs%norm = sqrt(abs(norm))
                
                dvecs%d=dvecs%d_1/dvecs%norm
        
                result=0.0d0
            
                do j=1,size
                    temp=0.0d0
                    do l=1,size 
                        temp=temp+haml%hjk(j,l)*dvecs%d(l)
                    end do 
                    result = result + (dvecs%d(j)*temp)
                end do
            
                erg(k)=result
               
        
                ddot=0.0d0
                
                do j=1,size 
                    temp=0.0d0
                    do l=1,size 
                        temp= temp + haml%kinvh(j,l)*dvecs%d(l)
                    end do 
                    ddot(j)=temp
                end do
            
                dvecs%d_1=dvecs%d-(db*ddot)
            
            end do

            norm=0.0d0
            do j=1,size
                temp=0.0d0
                do l=1,size 
                    temp=temp+haml%ovrlp(j,l)*dvecs%d_1(l)
                end do 
                norm=norm+(temp*dvecs%d_1(j))
            end do 
            dvecs%d_o_d=sign_d_o_d(norm)
            dvecs%norm = sqrt(abs(norm))
            
            dvecs%d=dvecs%d_1/dvecs%norm
            
        else !Imaginary time propagation with GSO
            
            ! ! Set up gs dvectors
            ! do k=1, haml%gram_num
            !     dvecs%d_gs=0.0
            !     dvecs%d_gs(k,k)=1.0
            ! end do
            ! dvecs%d(haml%gram_num+1)=1.0d0
            ! !Orthogonalise
            ! call gs(dvecs,haml)
            
            ! ! Normalise vectors
            ! do k=1, haml%gram_num
            !     norm=0.0d0
            !     do j=1,size
            !         temp=0.0d0
            !         do l=1,size 
            !             temp=temp+haml%gs_ovrlp(k,j,l)*dvecs%d_gs(k,l)
            !         end do 
            !         norm=norm+(temp*dvecs%d_gs(k,j))
            !     end do 
            
            !     norm = sqrt(abs(norm))    
            !     dvecs%d_gs(k,:)=dvecs%d_gs(k,:)/norm
    
            ! end do
            ! !Normalise dvector of interest
            ! do j=1,size
            !     temp=0.0d0
            !     do l=1,size 
            !         temp=temp+haml%ovrlp(j,l)*dvecs%d(l)
            !     end do 
            !     norm=norm+(temp*dvecs%d(j))
            ! end do 
        
            ! norm = sqrt(abs(norm))
        
            ! dvecs%norm=norm
        
            ! db=beta/timesteps
           
            ! !Begin time steps
            ! do k=1,timesteps+1
                
            !     ! Calculate energy of system being looked at
            !     result=0.0d0
            !     do j=1,size
            !         temp=0.0d0
            !         do l=1,size 
            !             temp=temp+haml%hjk(j,l)*dvecs%d(l)
            !         end do 
            !         result = result + (dvecs%d(j)*temp)
            !     end do
            
            !     erg(k)=result
            !     ! Take a time step for dvector and propagate derivative
            !     ddot=0.0d0

            !     do j=1,size 
            !         temp=0.0d0
            !         do l=1,size 
            !             temp= temp + haml%kinvh(j,l)*dvecs%d(l)
            !         end do 
            !         ddot(j)=temp
            !     end do
            
            !     dvecs%d=dvecs%d-(db*ddot)

            !     ! Make gs vectors take time step
            !     do g=1, haml%gram_num
            !         ddot=0.0d0
            !         do j=1,size 
            !             temp=0.0d0
            !             do l=1,size 
            !                 temp= temp + haml%gs_kinvh(g,j,l)*dvecs%d_gs(g,l)
            !             end do 
            !             ddot(j)=temp
            !         end do
            !         dvecs%d_gs(g,:)=dvecs%d_gs(g,:)-(db*ddot)
            !     end do 

            !     ! Orthogonalise 
            !     call gs(dvecs,haml)
                
            !     !Normnalise gs dvectors
            !     do g=1, haml%gram_num
            !         norm=0.0d0
            !         do j=1,size
            !             temp=0.0d0
            !             do l=1,size 
            !                 temp=temp+haml%gs_ovrlp(g,j,l)*dvecs%d_gs(g,l)
            !             end do 
            !             norm=norm+(temp*dvecs%d_gs(g,j))
            !         end do 
            !         norm = sqrt(abs(norm))    
            !         dvecs%d_gs(g,:)=dvecs%d_gs(g,:)/norm
            !     end do
                
            !     !Normalise dvector
            !     do j=1,size
            !         temp=0.0d0
            !         do l=1,size 
            !             temp=temp+haml%ovrlp(j,l)*dvecs%d(l)
            !         end do 
            !         norm=norm+(temp*dvecs%d(j))
            !     end do 
            
            !     norm = sqrt(abs(norm))
            
            !     dvecs%norm=norm
                
            ! end do
        end if 
      
        return

    end subroutine imaginary_time_prop2



END MODULE imgtp