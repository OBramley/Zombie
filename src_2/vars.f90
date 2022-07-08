MODULE globvars



    ! Type defining the zombie state
    type zombiest
        complex(kind=8), dimension(:), allocatable::alive
        complex(kind=8), dimension(:), allocatable::dead
    end type zombiest

    ! Type defining the Hamiltonian matrix and the overlap matrix
    type hamiltonian
        complex(kind=8), dimension(:,:), allocatable::Hjk
        complex(kind=8), dimension(:,:), allocatable::ovrlp
    end type hamiltonian

    type dvector
        complex(kind=8), dimension(:), allocatable::d
    end type dvector

    ! Type defining the 1&2 electron integrals
    type elecintrgl
        real(kind=8), dimension(:,:), allocatable::h1ei
        real(kind=8), dimension(:,:,:,:), allocatable::h2ei
        real(kind=8) :: hnuc
    end type elecintrgl


    integer::ndet       ! Number of Zombie states
    integer::norb       ! Number of spin orbitals
    integer::nel        ! Number of electrons in molecule
    real(kind=8)::spin  ! Spin of the molecule
    integer::beta       ! Distance proagated in imaginary time
    integer::timesteps  ! Number of time steps

    character(LEN=1)::zomgflg    ! Flag to generate zombie states or not
    character(LEN=1)::hamgflg    ! Flag to generate Hamiltonian or not
    character(LEN=1)::cleanflg   ! Flag to determine if cleaning is
    character(LEN=1)::gram       ! Flag to determine if gram schmidt orthogolnalisation should be carried out
    integer::gramnum        ! Number of additional states to be generated for GS orthogonalisation
    character(LEN=1)::rhf_1      ! Flag to decide if the first zombie state should be st as the RHF determinant

    real(kind=8)::pirl      ! pi
    real(kind=8) :: sqrtpi  ! square root of pi
    complex(kind=8)::i      ! The imagianry unit

    integer:: errorflag      ! Error flag

    contains

    subroutine initialise
        implicit none
        
        ! Initialises all the globabl variables to zero and sets the imaginary unit 

        ndet=0
        norb=0
        nel=0
        spin=0
        beta=0
        timesteps=0
        gramnum=0

        errorflag=0

        sqrtpi = 1.7724538509055160272981674833411451827975494561223871d0
        pirl = sqrtpi**2.0d0

        i = (0.0d0,1.0d0)

    end subroutine initialise

END MODULE globvars





