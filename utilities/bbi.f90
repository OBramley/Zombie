program bbi

    use globvars
    use readpars
    use alarrays
    use electrons
    use zom
    use ham
    use imgtp
    use outputs
    use comparison

    implicit none

    ! interface
    !     subroutine compare(zstore,zstoretemp,dvecstemp,en,entemp,elect,hamltemp,ergchange,achange,orbital)
    !         use globvars
    !         use imgtp
    !         use ham
    !         type(zombiest), dimension(:),intent(inout):: zstore,zstoretemp
    !         type(dvector), dimension(:), intent(inout):: dvecstemp
    !         type(energy),intent(inout):: en, entemp
    !         type(elecintrgl),intent(in)::elect
    !         type(hamiltonian),intent(inout)::hamltemp
    !         integer,intent(in)::orbital
    !         real(kind=8),dimension(:,:),intent(inout)::ergchange,achange
    !     end subroutine compare
    ! end interface

    type(zombiest), dimension(:), allocatable:: zstore,zstoretemp
    type(dvector), dimension(:), allocatable:: dvecs, dvecstemp
    type(energy):: en, entemp
    type(elecintrgl)::elect
    type(hamiltonian)::haml, hamltemp
    real(kind=8):: starttime, stoptime, runtime, erg, dummy,dummy1
    DOUBLE PRECISION, external::ZBQLUAB, ZBQLU01, ZBQLNOR
    character(LEN=100) :: CWD
    real(kind=8),dimension(:,:),allocatable::achange
    real(kind=8),dimension(:),allocatable::ergchange
    integer:: j,k,l,p, istat,ierr,iters
    integer(kind=8):: randseed

    real(kind=8)::mu(19),sig(19)
    real(kind=8),dimension(38)::val


    call CPU_TIME(starttime) !used to calculate the runtime, which is output at the end
    write(6,"(a)") " ________________________________________________________________ "
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|             Zombie State bias improve Program v1.10            |"
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|________________________________________________________________|"
    write(6,"(a)") ""
    write(6,"(a)") ""
    write(6,"(a)") ""

    ierr=0
    istat=0
    call initialise
    call readrunconds

    call allocintgrl(elect)
    call electronintegrals(elect)
    
    iters=10
    
    allocate(ergchange((iters*norb)+1),stat=ierr)
    if(ierr==0)allocate(achange(iters+1,norb),stat=ierr)
    if (ierr/=0) then
        write(0,"(a,i0)") "Error in results allocation. ierr had value ", ierr
        errorflag=1
        return
    end if
    ergchange(:)=0.0
    achange(:,:)=0.0

    open(unit=570, file="/dev/urandom", access="stream", &
    form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
        read(570) randseed    ! This takes the random seed from the true-random bin. If
        close(570)            ! the urandom bin does not exist the random seed is set
    else                       ! to zero which forces the date to be used
        randseed=0
    end if

    randseed = abs(randseed)    ! Negative seed values seem to cause instability

    call ZBQLINI(randseed,0)   ! Generates the seed value using the UCL random library
    write(6,"(a)") "Random seed set"

    call alloczs(zstore,ndet)
    if(zomgflg=='y') then
        call musig(mu,sig)
        do j=1, ndet
            !$omp critical
            do k=1,norb/2
                val(2*k-1)=2*pirl*ZBQLNOR(mu(k),sig(k))
                val(2*k)=2*pirl*ZBQLNOR(mu(k),sig(k))
            end do
            !$omp end critical
            do k=1, norb
                zstore(j)%alive(k)=cmplx(sin(val(k)),0.0d0,kind=8)
                zstore(j)%dead(K)=cmplx(cos(val(k)),0.0d0,kind=8)
            end do
        end do
        if(rhf_1=='y') then
            zstore(1)%alive(1:nel)=(1.0d0,0.0d0)
            zstore(1)%dead(1:nel)=(0.0d0,0.0d0)
            zstore(1)%alive((nel+1):norb)=(0.0d0,0.0d0)
            zstore(1)%dead((nel+1):norb)=(1.0d0,0.0d0)
        end if 
        ! call genzf(zstore,ndet)
        do j=1,ndet
            call zombiewriter(zstore(j),j)
        end do
    else 
        call read_zombie(zstore)
    end if

    write(6,"(a)") "Initial Zombie states generated"
    do j=1,norb
        achange(1,j)=(dasin(REAL(zstore(2)%alive(j))))!/(2*pirl))
        ! print*,achange(1,j)
    end do

    ! do j=1,norb
    !     print*,achange(1,j)
    ! end do
   
    
    call allocham(haml,ndet)
    call hamgen(haml,zstore,elect,ndet)
    print*,REAL(haml%hjk(2,1))
    call matrixwriter(haml%hjk,ndet,"ham.csv")
    call matrixwriter(haml%ovrlp,ndet,"ovrl.csv")
    call matrixwriter(haml%inv,ndet,"inv.csv")
    write(6,"(a)") "Initial Hamiltonian generated"
    call allocdv(dvecs,1,ndet)
    call allocerg(en,1)
    call imgtime_prop(dvecs,en,haml)

    erg=REAL(en%erg(1,timesteps+1))
    ergchange(1)=erg
    print*, "Starting energy, ", erg


    
    do j=1, norb
        call compare(zstore,erg,elect,haml,ergchange,achange,j,iters)
        write(6,"(a)") "Orbital completed"
        ! stop
    end do
    write(6,"(a)") "run finished"
    do j=1,ndet
        call zombiewriter_c(zstore(j),j)
    end do
    call matrixwriter(haml%hjk,ndet,"finalham.csv")
    call matrixwriter(haml%ovrlp,ndet,"finalovrl.csv")
    open(unit=200,file="ergchange.csv",status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening matrix file. ierr had value ", ierr
            errorflag=1
            return
        end if
        
        do j=1, (iters*norb)+1
            write(200,'(*(e25.17e3 :", "))') ergchange(j)
        end do
    close(200)

    open(unit=201,file="achange.csv",status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening matrix file. ierr had value ", ierr
            errorflag=1
            return
        end if
        
        do j=1, iters+1
            write(201,'(*(e25.17e3 :", "))') ((achange(j,k)),k=1,norb)
        end do
    close(201)

    call deallocerg(en)
    call deallocham(haml)
    call dealloczs(zstore)
    call deallocdv(dvecs)
    
    deallocate(ergchange,stat=ierr)
    if(ierr==0)deallocate(achange,stat=ierr)
    if (ierr/=0) then
        write(0,"(a,i0)") "Error in results deallocation. ierr had value ", ierr
        errorflag=1
        return
    end if

    call CPU_TIME(stoptime)
    runtime = stoptime-starttime
    call getcwd(CWD)

    if (errorflag.eq.0) write(6,"(a,a)") 'Successfully Executed Zombie states Program in ', trim(CWD)
    if (errorflag.ne.0) write(6,"(a,a)") 'Unsuccessfully Executed Zombie states  Program in ', trim(CWD)
  
    if (runtime/3600.0d0 .gt. 1.0d0)then
        runtime = runtime/3600.0d0
        write(6,"(a,es12.5,a)") 'Time taken : ', runtime, ' hours'
    else if (runtime/60.0d0 .gt. 1.0d0)then
        runtime = runtime/60.0d0
        write(6,"(a,es12.5,a)") 'Time taken : ', runtime , ' mins'
    else
        write(6,"(a,es12.5,a)") 'Time taken : ', runtime, ' seconds'
    end if


    call flush(6)
    call flush(0)

    stop

end program 
