program MainZombie
    
    
    use globvars
    use alarrays
    use electrons
    use ham
    use outputs
    use imgtp

    
    
    implicit none
    
    ! Private variables
    type(zombiest), dimension(:), allocatable:: zstore
    type(dvector), dimension(:), allocatable:: dvecs
    type(energy):: en
    type(elecintrgl),allocatable::elect
    type(hamiltonian)::haml
    integer:: j, k, n, m, istat 

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

    ! Maybe point to initialise Random flag

    ! Read in run conditions


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

    ! generate zombie states
    call alloczs(zstore,ndet)

    ! generate Hamiltonian and overlap

    call allocham(haml,ndet)
    call hamgen(haml,zstore,elecs,ndet)
    call matrixwriter(ham%hjk,ndet,"ham.csv")
    call matrixwriter(ham%ovrlp,ndet,"ovlp.csv")

    
    ! Imaginary time propagation
    if(gram.eq."n")then
        call allocdv(dvecs,1)
        call allocerg(en,1)
    else if(gramflg.eq."y")then
        call allocdv(dvecs,1+gramnum)
        call allocerg(en,1+gramnum)
    else
        write(0,) "Error in gramflg setting. This should have been caught ", ierr
            errorflag=1
    end if

    call imgtime_prop(dvecs,en,haml)






end program MainZombie



