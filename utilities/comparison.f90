module comparison

    use globvars
    use imgtp
    use ham

    contains

    subroutine compare(zstore,zstoretemp,dvecstemp,en,entemp,elect,hamltemp,ergchange,achange,orbital,iters)
        use globvars
        use ham
        use imgtp
        
        implicit none
    
        type(zombiest), dimension(:),intent(inout):: zstore,zstoretemp
        type(dvector), dimension(:), intent(inout):: dvecstemp
        type(energy),intent(inout):: en, entemp
        type(elecintrgl),intent(in)::elect
        type(hamiltonian),intent(inout)::hamltemp
        integer,intent(in)::orbital,iters
        real(kind=8),dimension(:,:),intent(inout)::achange
        real(kind=8),dimension(:),intent(inout)::ergchange
        real(kind=8),dimension(iters)::result1,result2
        real(kind=8):: dummy,dummy1,erg,temperg,ac,at
        DOUBLE PRECISION, external::ZBQLUAB, ZBQLU01
        integer:: j,k,l,p, istat,ierr
        
        erg=REAL(en%erg(1,timesteps))
        ac=achange(1,orbital)
        print*,erg,ac 
        do k=1, iters
            zstoretemp=zstore
            !$omp critical
            dummy=0.5*pirl*ZBQLU01(1)
            !$omp end critical    
            zstoretemp(2)%alive(orbital)=cmplx(sin(dummy),0.0d0,kind=8)
            zstoretemp(2)%dead(orbital)=cmplx(cos(dummy),0.0d0,kind=8)
            call hamgen(hamltemp,zstoretemp,elect,ndet)
            call imgtime_prop(dvecstemp,entemp,hamltemp)
            temperg=REAL(entemp%erg(1,timesteps))
            print*,temperg
            print*,dummy
            if(temperg<erg)then
                erg=temperg
                ! print*,j,"improve"
                ! do p=1, ndet
                zstore=zstoretemp
                ! end do
                en=entemp
                ac=dummy
                ! en%erg(1,1:timesteps+1)=entemp%erg(1,1:timesteps+1)
                !ergchange((orbital*iters)-iters+k+1)=REAL(en%erg(1,timesteps))
                !achange(orbital,k+1)=dummy
            !else 
            !    achange(orbital,k+1)=achange(orbital,k)
            !    ergchange((orbital*iters)-iters+k+1)=ergchange((orbital*iters)-iters+k)
            end if
            result1(k)=erg
            result2(k)=ac
            ! ergchange((orbital*iters)-iters+k+1)=erg
            ! achange(k+1,orbital)=ac

            dvecstemp(1)%d(1:ndet)=(0.0,0.0)
            dvecstemp(1)%d(1)=(1.0,0.0)
            hamltemp%hjk(1:ndet,1:ndet)=(0.0d0,0.0d0)
            hamltemp%ovrlp(1:ndet,1:ndet)=(0.0d0,0.0d0)
            hamltemp%inv(1:ndet,1:ndet)=(0.0d0,0.0d0)
            entemp%t(1:timesteps+1)=(0.0)
            entemp%erg(1,1:timesteps+1)=(0.0,0.0)
        end do

        ergchange((orbital*iters)-iters+2:(orbital*iters)+1)=result1
        achange(2:,orbital)=result2


        return
    
    end subroutine compare

end module comparison