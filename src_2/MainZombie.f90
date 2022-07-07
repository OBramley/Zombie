program MainZombie
    
    
    use globvars
    use alarrays
    use electrons
    
    
    implicit none
    
    ! Private variables
    type(zombiest), dimension(:), allocatable:: zstore
    type(dvector), dimension(:), allocatable:: dvecs
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
    
    ! Imaginary time propagation






end program MainZombie



