program d_check

    use globvars
    use readpars
    use alarrays
    use electrons
    use zom
    use ham
    use imgtp

    implicit none


    type(zombiest), dimension(:), allocatable:: zstore
    type(dvector), dimension(:), allocatable:: dvecs
    type(energy):: en
    type(elecintrgl)::elect
    type(hamiltonian)::haml
    real(kind=8):: starttime, stoptime, runtime
    character(LEN=100) :: CWD
    real(kind=8),dimension(:,:,:),allocatable::results
    integer:: j,k,l, istat,ierr,iters
    integer(kind=8):: randseed
    character(LEN=4)::orbital


    call CPU_TIME(starttime) !used to calculate the runtime, which is output at the end
    write(6,"(a)") " ________________________________________________________________ "
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|               Zombie State Analysis Program v1.00              |"
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
    iters=1000
    allocate(results(norb,ndet*iters,3),stat=ierr)
    if (ierr/=0) then
        write(0,"(a,i0)") "Error in results allocation. ierr had value ", ierr
        errorflag=1
        return
    end if
    results(:,:,:)=0.0

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

    !$omp parallel shared(elect,results,randseed) private(j,k,l,zstore,haml,dvecs,en,ierr,errorflag) 
    !$omp do
    do l=1, iters
        print*,"Iteration ",l, " started"
        call alloczs(zstore,ndet)
        call genzf(zstore,ndet)
        ! print*,"here",l
        call allocham(haml,ndet)
        call hamgen(haml,zstore,elect,ndet)
        ! print*,"here",l
        call allocdv(dvecs,1,ndet)
        call allocerg(en,1)
        ! print*,"here",l
        call imgtime_prop(dvecs,en,haml)
        ! print*,"here",l
        do j=1,norb
            results(j,((ndet*l+1)-ndet):ndet*l,2)=REAL(en%erg(1,timesteps))
            do k=1, ndet
                results(j,(ndet*l-ndet)+k,1)=REAL(zstore(k)%alive(j))
                results(j,(ndet*l-ndet)+k,3)=REAL(dvecs(1)%d(k))
            end do
        end do
        ! print*,"here",l
        call deallocerg(en)
        call deallocham(haml)
        call dealloczs(zstore)
        call deallocdv(dvecs)
        ! print*,"here",l

        print*,"Iteration ",l, " finished"
    end do
    !$omp end do
    !$omp end parallel

    call deallocintgrl(elect)
    
    do j=1,norb
        write(orbital,"(i4.4)")j
        open(unit=100+j,file="electron"//trim(orbital)//".csv",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening result file. ierr had value ", ierr
            errorflag=1
        end if
        do k=1,(ndet*iters)
            write(100+j,'(*(e25.17e3 :", "))') (results(j,k,l),l=1,3)
        end do
        close(100+j)
    end do

    deallocate(results)

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

end program d_check