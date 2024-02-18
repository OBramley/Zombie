MODULE imgtp
    
    use mod_types
    use globvars
    use alarrays
    
    interface imaginary_time
        module procedure imaginary_time_erg, imaginary_time_prop2, imaginary_time_prop2_gs
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
    
        values%dvec%d_1=0.0d0
        values%dvec%d_1(1)=1.0d0
    
        db=beta/timesteps

        do k=1,timesteps+1
            call DGEMV("N",size,size,1.d0,values%ovrlp,size,values%dvec%d_1,1,0.d0,temp,1)
            norm=dot_product(temp,values%dvec%d_1)
            values%dvec%norm = sqrt(abs(norm))
            values%dvec%d=values%dvec%d_1/values%dvec%norm

            call DGEMV("N",size,size,db,values%kinvh,ndet,values%dvec%d,1,0.d0,ddot,1)
            values%dvec%d_1=values%dvec%d-ddot
        end do

        call DGEMV("N",size,size,1.d0,values%ovrlp,size,values%dvec%d_1,1,0.d0,temp,1)
        norm=dot_product(temp,values%dvec%d_1)

        values%dvec%d_o_d=sign_d_o_d(norm)
        values%dvec%norm = sqrt(abs(norm))
    
        values%dvec%d=values%dvec%d_1/values%dvec%norm
        call DGEMV("N",size,size,1.d0,values%hjk,size,values%dvec%d,1,0.d0,temp,1)
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
        integer::j,k

        
        do j=1,gramnum
            numer=dot_product(dvecs%d_gs(j,:),matmul(ovrlp,dvecs%d_1))
            denom=dot_product(dvecs%d_1,matmul(ovrlp,dvecs%d_1))
            dvecs%d_gs(j,:) = dvecs%d_gs(j,:)-dvecs%d*numer/denom
            if(j>1)then
                do k=1,j-1
                    numer=dot_product(dvecs%d_gs(k,:),matmul(ovrlp,dvecs%d_gs(j,:)))
                    denom=dot_product(dvecs%d_gs(k,:),matmul(ovrlp,dvecs%d_gs(k,:)))
                    dvecs%d_gs(j,:) = dvecs%d_gs(j,:)-dvecs%d_gs(k,:)*numer/denom
                end do
            end if
        end do

        return 

    end subroutine gs_dvector

END MODULE imgtp