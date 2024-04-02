MODULE imgtp
    
    use mod_types
    use globvars
    use alarrays
    
    interface imaginary_time
        module procedure imaginary_time_erg, imaginary_time_prop2, &
            imaginary_time_prop2_gs,imaginary_time_wavefunction_gs,imaginary_time_wavefunction_erg_gs
    end interface imaginary_time
   
    contains

    subroutine imaginary_time_erg(values,size)

        implicit none
        type(grad_do),intent(inout)::values
        integer,intent(in)::size
        integer::k
        real(wp)::norm,db
        real(wp),dimension(size)::ddot,temp
        
        if (errorflag .ne. 0) return
    
        values%dvecs%d_1=0.0d0
        values%dvecs%d_1(1)=1.0d0
    
        db=beta/timesteps

        do k=1,timesteps+1
            call DGEMV("N",size,size,1.d0,values%ovrlp,size,values%dvecs%d_1,1,0.d0,temp,1)
            norm=dot_product(temp,values%dvecs%d_1)
            values%dvecs%norm = sqrt(abs(norm))
            values%dvecs%d=values%dvecs%d_1/values%dvecs%norm

            call DGEMV("N",size,size,db,values%kinvh,ndet,values%dvecs%d,1,0.d0,ddot,1)
            values%dvecs%d_1=values%dvecs%d-ddot
        end do

        call DGEMV("N",size,size,1.d0,values%ovrlp,size,values%dvecs%d_1,1,0.d0,temp,1)
        norm=dot_product(temp,values%dvecs%d_1)

        values%dvecs%d_o_d=sign_d_o_d(norm)
        values%dvecs%norm = sqrt(abs(norm))
    
        values%dvecs%d=values%dvecs%d_1/values%dvecs%norm
        call DGEMV("N",size,size,1.d0,values%hjk,size,values%dvecs%d,1,0.d0,temp,1)
        values%erg=dot_product(temp,values%dvecs%d)

    end subroutine  imaginary_time_erg

    ! Calculates the energy
    function ergcalc(bham,dvecs) result(result)

        implicit none

        real(wp),intent(in),dimension(:)::dvecs
        real(wp),intent(in),dimension(:,:)::bham
        real(wp)::result
        real(wp),dimension(ndet)::temp
       
        call DGEMV("N",ndet,ndet,1.d0,bham,ndet,dvecs,1,0.d0,temp,1)
        result=dot_product(temp,dvecs)

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
   
        dvecs%d_1=0.0d0
        dvecs%d_1(1)=1.0d0
    
        db=beta/timesteps
    
        do k=1,timesteps+1
        
            call DGEMV("N",size,size,1.d0,haml%ovrlp,size,dvecs%d_1,1,0.d0,temp,1)
            norm=dot_product(temp,dvecs%d_1)
            
            dvecs%norm = sqrt(abs(norm))
            dvecs%d=dvecs%d_1/dvecs%norm

            call DGEMV("N",size,size,1.d0,haml%hjk,size,dvecs%d,1,0.d0,temp,1)
        
            erg(k)=dot_product(temp,dvecs%d)
            
            call DGEMV("N",size,size,db,haml%kinvh,size,dvecs%d,1,0.d0,ddot,1)
            dvecs%d_1=dvecs%d-ddot
        
        end do
        call DGEMV("N",size,size,1.d0,haml%ovrlp,size,dvecs%d_1,1,0.d0,temp,1)
        norm=dot_product(temp,dvecs%d_1)

        dvecs%d_o_d=sign_d_o_d(norm)
        dvecs%norm = sqrt(abs(norm))
        dvecs%d=dvecs%d_1/dvecs%norm
            
        return

    end subroutine imaginary_time_prop2


    subroutine imaginary_time_prop2_gs(dvecs,erg,haml,size)

        implicit none

        type(dvector),intent(inout)::dvecs
        real(wp), dimension(:,:),intent(inout)::erg
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::size
        integer::k,g
        real(wp)::norm
        real(wp)::db
        real(wp),dimension(size)::ddot,temp


        dvecs%d_1=0.0d0
        dvecs%d_1(1)=1.0d0
        dvecs%d_gs=0.0d0
        do k=1, gramnum
            dvecs%d_gs(k,k+1)=1.0d0
        end do
        call gs_dvector(dvecs,haml%ovrlp)
        db=beta/timesteps
    
        do k=1,timesteps+1
        
            call DGEMV("N",size,size,1.d0,haml%ovrlp,size,dvecs%d_1,1,0.d0,temp,1)
            norm=dot_product(temp,dvecs%d_1)
            
            dvecs%norm = sqrt(abs(norm))
            dvecs%d=dvecs%d_1/dvecs%norm
            do g=1, gramnum
                call DGEMV("N",size,size,1.d0,haml%ovrlp,size,dvecs%d_gs(g,:),1,0.d0,temp,1)
                norm=dot_product(temp,dvecs%d_gs(g,:))
                dvecs%d_gs(g,:)=dvecs%d_gs(g,:)/sqrt(abs(norm))
            end do
            

            call DGEMV("N",size,size,1.d0,haml%hjk,size,dvecs%d,1,0.d0,temp,1)
            erg(1,k)=dot_product(temp,dvecs%d)
            do g=1, gramnum
                call DGEMV("N",size,size,1.d0,haml%hjk,size,dvecs%d_gs(g,:),1,0.d0,temp,1)
                erg(g+1,k)=dot_product(temp,dvecs%d_gs(g,:))
            end do
            
            call DGEMV("N",size,size,db,haml%kinvh,size,dvecs%d,1,0.d0,ddot,1)
            dvecs%d_1=dvecs%d-ddot
            do g=1, gramnum
                call DGEMV("N",size,size,db,haml%kinvh,size,dvecs%d_gs(g,:),1,0.d0,ddot,1)
                dvecs%d_gs(g,:)=dvecs%d_gs(g,:)-ddot
            end do
            call gs_dvector(dvecs,haml%ovrlp)
        
        end do
        call DGEMV("N",size,size,1.d0,haml%ovrlp,size,dvecs%d_1,1,0.d0,temp,1)
        norm=dot_product(temp,dvecs%d_1)

        dvecs%d_o_d=sign_d_o_d(norm)
        dvecs%norm = sqrt(abs(norm))
        dvecs%d=dvecs%d_1/dvecs%norm

        return

    end subroutine imaginary_time_prop2_gs

    subroutine gs_dvector(dvecs,ovrlp)
        implicit none
        type(dvector),intent(inout)::dvecs
        real(wp),dimension(:,:),intent(in)::ovrlp
        real(wp)::numer,denom
        real(wp),dimension(ndet)::temp
        integer::j,k

        
        do j=1,gramnum
            numer=dot_product(dvecs%d_gs(j,:),matmul(ovrlp,dvecs%d_1))
            denom=dot_product(dvecs%d_1,matmul(ovrlp,dvecs%d_1))
            temp=dvecs%d_1*numer/denom
            if(j>1)then
                do k=1,j-1
                    numer=dot_product(dvecs%d_gs(k,:),matmul(ovrlp,dvecs%d_gs(j,:)))
                    denom=dot_product(dvecs%d_gs(k,:),matmul(ovrlp,dvecs%d_gs(k,:)))
                    temp = temp + dvecs%d_gs(k,:)*numer/denom
                end do
            end if
            dvecs%d_gs(j,:) = dvecs%d_gs(j,:)-temp
        end do

        return 

    end subroutine gs_dvector

    subroutine gs_wavefunction_vector(dvecs,gs_numer1,gs_numer2,size,state)
        implicit none
        type(dvector),intent(inout)::dvecs
        real(wp),dimension(:,:),intent(in)::gs_numer1,gs_numer2
        integer,intent(in)::size,state
        real(wp),dimension(size)::temp
        integer::j
        temp=0.0d0
        do j=1,state-1 
            temp=temp+dot_product(dvecs%d_1,gs_numer1(j,:))*gs_numer2(j,:)
        end do
        dvecs%d_1=dvecs%d_1-temp
        return 

    end subroutine gs_wavefunction_vector 

    subroutine imaginary_time_wavefunction_gs(gramstore,erg,size,state)

        implicit none

        type(gram),dimension(:)::gramstore
        real(wp), dimension(:),intent(inout)::erg
        integer,intent(in)::size,state
        integer::k,j
        real(wp)::norm
        real(wp)::db
        real(wp),dimension(size)::ddot,temp
        real(wp),dimension(state-1,size)::gs_numer1,gs_numer2
        
        if(errorflag.ne.0) return

        do j=1,state-1 
            gs_numer1(j,:)=matmul(gramstore(j)%dvecs%d,gramstore(state)%wf_ovrlp(j,:,:))
            gs_numer2(j,:)=matmul(gramstore(state)%haml%inv,gs_numer1(j,:))/gramstore(j)%d_ovrlp_d
        end do 

        gramstore(state)%dvecs%d_1=0.0d0
        gramstore(state)%dvecs%d_1(1)=1.0d0
    
        db=beta/timesteps
    
        do k=1,timesteps+1
            call  gs_wavefunction_vector(gramstore(state)%dvecs,gs_numer1,gs_numer2,size,state)

            call DGEMV("N",size,size,1.d0,gramstore(state)%haml%ovrlp,size,gramstore(state)%dvecs%d_1,1,0.d0,temp,1)
            norm=dot_product(temp,gramstore(state)%dvecs%d_1)
            
            gramstore(state)%dvecs%norm = sqrt(abs(norm))
            gramstore(state)%dvecs%d=gramstore(state)%dvecs%d_1/gramstore(state)%dvecs%norm
        
            call DGEMV("N",size,size,1.d0,gramstore(state)%haml%hjk,size,gramstore(state)%dvecs%d,1,0.d0,temp,1)
            erg(k)=dot_product(temp,gramstore(state)%dvecs%d)
          
            call DGEMV("N",size,size,db,gramstore(state)%haml%kinvh,size,gramstore(state)%dvecs%d,1,0.d0,ddot,1)
            gramstore(state)%dvecs%d_1=gramstore(state)%dvecs%d-ddot
           
        end do
      
            call DGEMV("N",size,size,1.d0,gramstore(state)%haml%ovrlp,size,gramstore(state)%dvecs%d_1,1,0.d0,temp,1)
            norm=dot_product(temp,gramstore(state)%dvecs%d_1)

            gramstore(state)%dvecs%d_o_d=sign_d_o_d(norm)
            gramstore(state)%dvecs%norm = sqrt(abs(norm))
            gramstore(state)%dvecs%d=gramstore(state)%dvecs%d_1/gramstore(state)%dvecs%norm
       

        return

    end subroutine imaginary_time_wavefunction_gs

    subroutine imaginary_time_wavefunction_erg_gs(values,gramstore,size,state)

        implicit none
        type(grad_do),intent(inout)::values
        type(gram),dimension(:)::gramstore
        integer,intent(in)::size,state
        integer::k,j
        real(wp)::norm
        real(wp)::db
        real(wp),dimension(size)::ddot,temp
        real(wp),dimension(state-1,size)::gs_numer1,gs_numer2
        
        if(errorflag.ne.0) return

        do j=1,state-1 
            gs_numer1(j,:)=matmul(gramstore(j)%dvecs%d,values%wf_ovrlp(j,:,:))
            gs_numer2(j,:)=matmul(values%inv,gs_numer1(j,:))/gramstore(j)%d_ovrlp_d
        end do 

        values%dvecs%d_1=0.0d0
        values%dvecs%d_1(1)=1.0d0
    
        db=beta/timesteps
    
        do k=1,timesteps+1
            call  gs_wavefunction_vector(values%dvecs,gs_numer1,gs_numer2,size,state)

            call DGEMV("N",size,size,1.d0,values%ovrlp,size,values%dvecs%d_1,1,0.d0,temp,1)
            norm=dot_product(temp,values%dvecs%d_1)
            values%dvecs%norm = sqrt(abs(norm))
            values%dvecs%d=values%dvecs%d_1/values%dvecs%norm

            call DGEMV("N",size,size,db,values%kinvh,ndet,values%dvecs%d,1,0.d0,ddot,1)
            values%dvecs%d_1=values%dvecs%d-ddot
           
        end do
      
        call DGEMV("N",size,size,1.d0,values%ovrlp,size,values%dvecs%d_1,1,0.d0,temp,1)
        norm=dot_product(temp,values%dvecs%d_1)

        values%dvecs%d_o_d=sign_d_o_d(norm)
        values%dvecs%norm = sqrt(abs(norm))
    
        values%dvecs%d=values%dvecs%d_1/values%dvecs%norm
        call DGEMV("N",size,size,1.d0,values%hjk,size,values%dvecs%d,1,0.d0,temp,1)
        values%erg=dot_product(temp,values%dvecs%d)
       

        return

    end subroutine imaginary_time_wavefunction_erg_gs



    

END MODULE imgtp