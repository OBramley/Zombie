module comparison

    use globvars
    use alarrays
    use imgtp
    use ham
    use operators
    contains

    subroutine compare(zstore,erg,elect,haml,ergchange,achange,orbital,iters,pass)
        
        implicit none
    
        type(zombiest), dimension(:),intent(inout):: zstore
        type(elecintrgl),intent(in)::elect
        type(hamiltonian),intent(inout)::haml
        real(kind=8),intent(inout)::erg
        integer,intent(in)::orbital,iters,pass
        type(zombiest), dimension(:),allocatable::zstoretemp
        type(hamiltonian)::hamltemp
        type(dvector), dimension(:), allocatable:: dvecstemp
        type(energy)::entemp
        real(kind=8),dimension(:,:),intent(inout)::achange
        real(kind=8),dimension(:),intent(inout)::ergchange
        real(kind=8):: dummy,temperg,ac,ul,ll
        DOUBLE PRECISION, external::ZBQLU01
        integer:: k,max,check,thresh,upper,lower,dummyforce,j,p,abort
        
 

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
        thresh=max/10
        ul=0.5*pirl
        ll=0
        k=0+(iters*pass)
        check=0
        upper=0
        lower=0
        dummyforce=0
        p=0
        abort=0
        dummy=3*pirl

        do while(k.lt.iters) 
            zstoretemp=zstore
            ! j=0
            ! do while((dummy.lt.ll).and.(dummy.gt.ul))
                ! j=j+1
            !$omp critical
            dummy= 0.5*pirl*ZBQLU01(1)
            !$omp end critical
                ! if(j.eq.200000000)then
                    ! abort=abort+2
                    ! print*,"maxed out"
                    ! exit 
                ! end if
            ! end do

            if(dummyforce.eq.-1)then
                ! print*,"dummyforce -1"
                ! check=0
                !dummy less than ac
                j=0
                do while((dummy.gt.ac))!.and.(dummy.lt.ll))
                    j=j+1
                    !$omp critical
                    dummy= 0.5*pirl*ZBQLU01(1)
                    !$omp end critical
                    ! if(modulo(j,1000000).eq.0)then
                    !     print*,"j is",j
                    ! end if
                    if(j.eq.100000000)then
                        print*,"maxed out"
                        abort=abort+1
                        p=max-1
                        exit 
                    end if
                end do
                p=p+1
                ! print*,"p has value",p
            end if
            if(dummyforce.eq.2)then
                ! print*,"dummyforce 2"
                ! check=0
                ! p=p+upper
                !dummy greater than ac
                j=0
                do while((dummy.lt.ac))!.and.(dummy.gt.ul))
                    j=j+1
                    !$omp critical
                    dummy= 0.5*pirl*ZBQLU01(1)
                    !$omp end critical
                    ! if(modulo(j,1000000).eq.0)then
                    !     print*,"j is",j
                    ! end if
                    if(j.eq.100000000)then
                        print*,"maxed out"
                        abort=abort+1
                        p=max-1
                        exit 
                    end if
                end do
                p=p+1
                ! print*,"p has value",p
            end if

            if(p.eq.max)then
                ! print*,"Still no luck"
                if(abort.ge.2)then
                    print*,"final abort"
                    achange(2+k:,orbital)=achange(1+k,orbital)
                    ergchange((orbital*iters)-iters+2+k:(orbital*iters)+1)=ergchange((orbital*iters)-iters+1+k)
                    exit
                end if 

                if(dummyforce.eq.-1)then 
                    dummyforce=2
                    check=upper
                    p=0
                    abort=abort+1
                else if (dummyforce.eq.2)then
                    dummyforce=-1
                    abort=abort+1
                    check=lower
                    p=0
                end if  
            end if
           
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
                lower=0
                upper=0
                p=0
                abort=0
                dummyforce=0
            else 
                ! print*,"reject"
                check=check+1
                if(check.eq.max)then
                    ! print*,"Clocking out"
                    achange(2+k:,orbital)=achange(1+k,orbital)
                    ergchange((orbital*iters)-iters+2+k:(orbital*iters)+1)=ergchange((orbital*iters)-iters+1+k)
                    exit 
                end if
                if(dummy.gt.ac)then
                    upper=upper+1
                else if(dummy.lt.ac)then
                    lower=lower+1
                end if
                if(dummyforce.eq.0)then
                    if((upper.ge.thresh).or.(lower.ge.thresh))then
                        ! print*,"upper lower = ",upper,lower
                        if((upper.eq.lower).and.(upper.eq.(5*thresh)))then
                            ! print*,"Equal upper and lower = 5xthresh"
                            achange(2+k:,orbital)=achange(1+k,orbital)
                            ergchange((orbital*iters)-iters+2+k:(orbital*iters)+1)=ergchange((orbital*iters)-iters+1+k)
                            exit
                        end if
                        if(upper.gt.lower)then
                            ! print*,"forcing dummy lower"
                            dummyforce=-1
                            check=0
                        else if(upper.lt.lower)then
                            ! print*,"forcing dummy higher"
                            dummyforce=2
                            check=0
                        end if
                    end if
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

    subroutine random_ord(n,array) !result(array)
        
        implicit none

        integer, intent(in)::n
        integer, dimension(:),allocatable, intent(out)::array
        integer::j,k,l,temp,ierr
        real::r
        allocate(array(n),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in mixed allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        do j=1, n
            array(j)=j
        end do
        ! print*,array(1:n)
        do j=1,2
            do k=1,n
                call random_number(r)
                l=1+floor((n+1)*r)
                temp=array(l)
                array(l)=array(k)
                array(k)=temp
            end do
        end do
        return
    end subroutine random_ord




end module comparison