module comparison

    use globvars
    use alarrays
    use imgtp
    use ham
    use operators
    contains

    subroutine compare(zstore,erg,elect,haml,ergchange,achange,orbital,iters)
        
        implicit none
    
        type(zombiest), dimension(:),intent(inout):: zstore
        type(elecintrgl),intent(in)::elect
        type(hamiltonian),intent(inout)::haml
        real(kind=8),intent(inout)::erg
        integer,intent(in)::orbital,iters
        type(zombiest), dimension(:),allocatable::zstoretemp
        type(hamiltonian)::hamltemp
        type(dvector), dimension(:), allocatable:: dvecstemp
        type(energy)::entemp
        real(kind=8),dimension(:,:),intent(inout)::achange
        real(kind=8),dimension(:),intent(inout)::ergchange
        real(kind=8),dimension(iters)::result1,result2
        real(kind=8):: dummy,temperg,ac,upper,lower, ul,ll
        DOUBLE PRECISION, external::ZBQLUAB, ZBQLU01
        integer:: p,j,k,updwn,start,max,check
        


        call alloczs(zstoretemp,ndet)
        call allocham(hamltemp,ndet)
        call allocdv(dvecstemp,1,ndet)
        call allocerg(entemp,1)
        ac=achange(1,orbital)
        print*,"Starting coeffcients"
        print*,REAL(zstore(2)%alive(orbital))
        print*,REAL(zstore(2)%dead(orbital))
        print*,ac,erg
        print*,"Improvement loop begun"
        
        max=100
        

        k=0
        check=0
        do while(k.lt.iters) 
            zstoretemp=zstore
                ! if(updwn.eq.1) then
                !     dummy=-100
                !     p=0
                !     do while(dummy.lt.ac)
                !         !$omp critical
                !         dummy= lower+((upper-lower)*ZBQLU01(1))
                !         !$omp end critical
                !         if(dummy.ge.ul)then
                !             dummy=-100
                !         end if
                !         p=p+1
                !         if(p.eq.max)then
                !             exit 
                !         end if
                !     end do
                !     if(p.eq.max)then
                !         print*,"Biasing at max"
                !         result1(k:)=erg
                !         result2(k:)=ac
                !         exit 
                !     end if
                !     ! print*,p
                ! else if(updwn.eq.0) then
                !     dummy=100
                !     p=0
                !     do while(dummy.gt.ac)
                !         !$omp critical
                !         dummy= lower+((upper-lower)*ZBQLU01(1))
                !         !$omp end critical
                !         if(dummy.le.ll)then
                !             dummy=100
                !         end if
                !         p=p+1
                !         if(p.eq.max)then
                !             exit 
                !         end if
                !     end do
                !     if(p.eq.max)then
                !         print*,"Biasing at max"
                !         result1(k:)=erg
                !         result2(k:)=ac
                !         exit 
                !     end if
                !     ! print*,p
                ! end if

            !$omp critical
            dummy= 0.5*pirl*ZBQLU01(1)
            ! dummy= lower+((upper-lower)*ZBQLU01(1))  !0.5*pirl*ZBQLU01(1)
            !$omp end critical
           
            zstoretemp(2)%alive(orbital)=cmplx(sin(dummy),0.0d0,kind=8)
            zstoretemp(2)%dead(orbital)=cmplx(cos(dummy),0.0d0,kind=8)
         
            hamltemp=haml
            hamltemp%ovrlp(1,2)=overlap(zstoretemp(1),zstoretemp(2))
            hamltemp%ovrlp(2,1)=hamltemp%ovrlp(1,2)
            hamltemp%inv(1,1)=hamltemp%ovrlp(1,1)/(((hamltemp%ovrlp(1,1))**2)-((hamltemp%ovrlp(1,2))**2))
            hamltemp%inv(1,1)=hamltemp%inv(2,2)
            hamltemp%inv(1,2)=-hamltemp%ovrlp(1,2)/(((hamltemp%ovrlp(1,1))**2)-((hamltemp%ovrlp(1,2))**2))
            hamltemp%inv(2,1)=hamltemp%inv(1,2)
            hamltemp%hjk(1,2)=hamval(zstoretemp(1),zstoretemp(2),elect,hamltemp%ovrlp(1,2))
            hamltemp%hjk(2,1)=hamltemp%hjk(1,2)
            hamltemp%hjk(2,2)=hamval(zstoretemp(2),zstoretemp(2),elect,hamltemp%ovrlp(2,2))
            hamltemp%kinvh=matmul(hamltemp%inv,hamltemp%hjk)
            

            ! call hamgen(hamltemp,zstoretemp,elect,ndet)
            call imgtime_prop(dvecstemp,entemp,hamltemp)
            temperg=REAL(entemp%erg(1,timesteps+1))
       
            if(temperg.lt.erg)then
                erg=temperg
                print*,"Accept new energ is",erg
                k=k+1
                check=0
                zstore=zstoretemp
                ac=dummy
                haml=hamltemp
                achange(1+k,orbital)=ac
                ergchange((orbital*iters)-iters+1+k)=erg
    

            else 
                check=check+1
                if(check.eq.max)then
                    achange(2+k:,orbital)=achange(1+k,orbital)
                    ergchange((orbital*iters)-iters+2+k:(orbital*iters)+1)=ergchange((orbital*iters)-iters+1+k)
                    exit 
                end if
            end if   
            

            dvecstemp(1)%d(1:ndet)=(0.0,0.0)
            entemp%t(1:timesteps+1)=(0.0)
            entemp%erg(1,1:timesteps+1)=(0.0,0.0)
        end do
        print*,"Final coeffcients"
        print*,REAL(zstore(2)%alive(orbital))
        print*,REAL(zstore(2)%dead(orbital))
        print*,ac,erg
        
        call deallocdv(dvecstemp)
        call deallocerg(entemp)
        call deallocham(hamltemp)
        call dealloczs(zstoretemp)
        return
    
    end subroutine compare

end module comparison