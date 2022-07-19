program MainZombie
    
    
    use globvars
    use readpars
    use alarrays
    use electrons
    use ham
    use outputs
    use imgtp
    use clean

    
    
    implicit none
    
    ! Private variables
    type(zombiest), dimension(:), allocatable:: zstore, cstore
    type(dvector), dimension(:), allocatable:: dvecs, dvec_clean
    type(energy):: en, en_clean
    type(elecintrgl),allocatable::elect
    type(hamiltonian)::haml, clean_haml
    integer:: j, k, n, m, istat, clean_ndet
    complex(kind=8)::clean_norm, clean_erg
    character(LEN=2)::stateno
    character(LEN=100) :: CWD

    ! Public variables
    real(kind=8):: starttime, stoptime, runtime
    integer(kind=8):: randseed

    call CPU_TIME(starttime) !used to calculate the runtime, which is output at the end
    write(6,"(a)") " ________________________________________________________________ "
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|               Zombie State Simulation Program v2.00            |"
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|________________________________________________________________|"
    write(6,"(a)") ""
    write(6,"(a)") ""
    write(6,"(a)") ""

    call initialise
    call readrunconds
   


    open(unit=570, file="/dev/urandom", access="stream", &
    form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
        read(570) ranseed    ! This takes the random seed from the true-random bin. If
        close(570)           ! the urandom bin does not exist the random seed is set
    else                   ! to zero which forces the date to be used
        ranseed=0
    end if

    ranseed = abs(ranseed)    ! Negative seed values seem to cause instability

    call ZBQLINI(ranseed,0)   ! Generates the seed value using the UCL random library

    ! generate 1 and 2 electron integrals
    call allocintgrl(elect)
    call electronintegrals(elecs)
    write(6,"(a)") "1 & 2 electron integrals successfully generated"

    ! generate zombie states
    call alloczs(zstore,ndet)
    if(zomgflg=='y')then
        call genzf(zstore,ndet)
        write(6,"(a)") " Zombie states generated"
    else if (zomgflg=='n') then
        write(6,"(a)") " Need to write read in routine"
    end if
    ! generate Hamiltonian and overlap

    call allocham(haml,ndet)
    if(hamgflg=='y')then
        call hamgen(haml,zstore,elecs,ndet)
        call matrixwriter(ham%hjk,ndet,"ham.csv")
        call matrixwriter(ham%ovrlp,ndet,"ovlp.csv")
        write(6,"(a)") "Hamiltonian successfully generated"
    else if (hamgflg=='n')then
        write(6,"(a)") " Need to write read in routine"
    end if

    if(propflg=='y') then
        ! Imaginary time propagation
        if(gramflg.eq."n")then
            call allocdv(dvecs,1,ndet)
            call allocerg(en,1)
        else if(gramflg.eq."y")then
            call allocdv(dvecs,1+gramnum,ndet)
            call allocerg(en,1+gramnum)
        else
            write(0,) "Error in gramflg setting. This should have been caught ", ierr
                errorflag=1
        end if

        write(6,"(a)") "Imaginary time propagation started"
        call imgtime_prop(dvecs,en,haml)
        write(6,"(a)") "Imaginary time propagation finished"

        if(gram.eq."n")then
            call energywriter(en%t,en%erg,"energy.csv",0)
        else if(gramflg.eq."y")then
            do j=1, 1+gramnum
                write(stateno,"(i4.4)")j
                call energywriter(en(j)%t,en(j)%erg,"energy_state_"//trim(stateno)//".csv",j)
            end do
        end if

        call deallocham(haml)
        write(6,"(a)") "Hamiltonian deallocated"

        if(cleanflg=="y")then
            call clean_setup(cstore,nel,clean_haml,en,clean_ndet)
            write(6,"(a)") "Cleaning hamiltonian generated"
            call allocdv(dvec_clean,1,clean_ndet)
            call cleaner(zstore,cstore,dvecs,dvec_clean,clean_ndet,clean_norm)
            clean_erg=ergcalc(clean_haml,dvec_clean)
            write(6,"(a)") "Cleaning process complete"
            call allocerg(en_clean,1)
            en_clean%t(1:timestep+1)=en%t(1:timestep+1)
            en_clean%erg(1:timestep+1)=clean_erg/clean_norm
            call energywriter(en_clean%t,en_clean%erg,"clean_energy.csv",99)
            call deallocerg(en_clean)
            call deallocdv(dvec_clean)
            call deallocham(clean_haml)
            call dealloczs(cstore)
            write(6,"(a)") "Cleaning dealocated"
        end if

        call deallocerg(en)
        call deallocdv(dvecs)
        call dealloczs(zstore)
        call deallocintgrl(elect)
    else 
        call deallocham(haml)
        call dealloczs(zstore)
        call deallocintgrl(elect)
    end if

    write(6,"(a)") "All values deallocated"

    if (errorflag .ne. 0) then
        write(6,"(a)") "Program terminated early."
        write(6,"(a,i0)") "errorflag value is ", errorflag
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


end program MainZombie



