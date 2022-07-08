program MainZombie
    
    
    use globvars
    use alarrays
    use electrons
    use ham
    
    
    implicit none
    
    ! Private variables
    type(zombiest), dimension(:), allocatable:: zstore
    type(dvector), dimension(:), allocatable:: dvecs
    type(energy):: en
    type(elecintrgl),allocatable::elect
    integer:: j, k, n, m 

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

    ! generate 1 and 2 electron integrals
    call allocintgrl(elect)
    call electronintegrals(elecs)


    !written

    ! generate zombie states
    call alloczs(zstore,ndet)

    ! generate Hamiltonian and overlap

    call allocham(ham,ndet)
    call hamgen(ham,zstore,elecs,ndet)
    
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






end program MainZombie



