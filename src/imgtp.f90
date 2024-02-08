MODULE imgtp
    
    use mod_types
    use globvars
    use alarrays
   
    contains


    subroutine imaginary_time_erg(values,size)

        implicit none
        type(grad_do),intent(inout)::values
        integer,intent(in)::size
        integer::k
        real(wp)::norm,db
        real(wp),dimension(size)::ddot,temp
        
        if (errorflag .ne. 0) return
    
        values%dvec%d_1=0.0d0
        values%dvec%d_1(1)=1.0d0
    
        db=beta/timesteps

        do k=1,timesteps+1
            call DGEMV("N",ndet,ndet,1.d0,values%ovrlp,ndet,values%dvec%d_1,1,0.d0,temp,1)
            norm=dot_product(temp,values%dvec%d_1)
            values%dvec%norm = sqrt(abs(norm))
            values%dvec%d=values%dvec%d_1/values%dvec%norm

            call DGEMV("N",ndet,ndet,db,values%kinvh,ndet,values%dvec%d,1,0.d0,ddot,1)
            values%dvec%d_1=values%dvec%d-ddot
        end do

        call DGEMV("N",ndet,ndet,1.d0,values%ovrlp,ndet,values%dvec%d_1,1,0.d0,temp,1)
        norm=dot_product(temp,values%dvec%d_1)

        values%dvec%d_o_d=sign_d_o_d(norm)
        values%dvec%norm = sqrt(abs(norm))
    
        values%dvec%d=values%dvec%d_1/values%dvec%norm
        call DGEMV("N",ndet,ndet,1.d0,values%hjk,ndet,values%dvec%d,1,0.d0,temp,1)
        values%erg=dot_product(temp,values%dvec%d)

    end subroutine  imaginary_time_erg

    ! Calculates the energy
    function ergcalc(bham,dvec) result(result)

        implicit none

        real(wp),intent(in),dimension(:)::dvec
        real(wp),intent(in),dimension(:,:)::bham
        real(wp)::result
        real(wp),dimension(ndet)::temp
       
        call DGEMV("N",ndet,ndet,1.d0,bham,ndet,dvec,1,0.d0,temp,1)
        result=dot_product(temp,dvec)

        return
       
    end function ergcalc

    subroutine imaginary_time_prop2(dvecs,erg,haml,size)

        implicit none

        type(dvector),intent(inout)::dvecs
        real(wp), dimension(:),intent(inout)::erg
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::size
        integer::k!,g
        real(wp)::norm
        real(wp)::db
        real(wp),dimension(size)::ddot,temp


        if (errorflag .ne. 0) return
        ! Imaginary time Propagation if no gram-schmidt orthogonalization is needed
        if(haml%gram_num.eq.0)then
            dvecs%d_1=0.0d0
            dvecs%d_1(1)=1.0d0
        
            db=beta/timesteps
        
            do k=1,timesteps+1
            
                call DGEMV("N",ndet,ndet,1.d0,haml%ovrlp,ndet,dvecs%d_1,1,0.d0,temp,1)
                norm=dot_product(temp,dvecs%d_1)
               
                dvecs%norm = sqrt(abs(norm))
                dvecs%d=dvecs%d_1/dvecs%norm
    
                call DGEMV("N",ndet,ndet,1.d0,haml%hjk,ndet,dvecs%d,1,0.d0,temp,1)
            
                erg(k)=dot_product(temp,dvecs%d)
                
                call DGEMV("N",ndet,ndet,db,haml%kinvh,ndet,dvecs%d,1,0.d0,ddot,1)
                dvecs%d_1=dvecs%d-ddot
            
            end do
            call DGEMV("N",ndet,ndet,1.d0,haml%ovrlp,ndet,dvecs%d_1,1,0.d0,temp,1)
            norm=dot_product(temp,dvecs%d_1)

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