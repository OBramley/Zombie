MODULE globvars



    ! Type defining the zombie state
    type zombiest
        complex(kind=8), dimension(:), allocatable::alive
        complex(kind=8), dimension(:), allocatable::dead
        real(kind=8),dimension(:),allocatable::diffalive
        real(kind=8),dimension(:),allocatable::diffdead
    end type zombiest

    ! Type defining the Hamiltonian matrix and the overlap matrix
    type hamiltonian
        complex(kind=8), dimension(:,:), allocatable::hjk
        complex(kind=8), dimension(:,:), allocatable::ovrlp
        complex(kind=8), dimension(:,:), allocatable::inv
        complex(kind=8), dimension(:,:), allocatable::kinvh

        real(kind=8), dimension(:,:,:), allocatable::diff_hjk
        real(kind=8), dimension(:,:,:), allocatable::diff_ovrlp
        real(kind=8), dimension(:,:,:), allocatable::diff_inv 
    end type hamiltonian

    type dvector
        complex(kind=8), dimension(:), allocatable::d
        real(kind=8), dimension(:,:),allocatable::d_diff
    end type dvector

    type energy
        real(kind=8), dimension(:),allocatable::t
        complex(kind=8), dimension(:,:),allocatable::erg 
    end type energy

    ! Type defining the 1&2 electron integrals
    type elecintrgl
        real(kind=8), dimension(:,:), allocatable::h1ei
        real(kind=8), dimension(:,:,:,:), allocatable::h2ei
        real(kind=8) :: hnuc
    end type elecintrgl

    type grad 
        real(kind=8),dimension(:,:), allocatable::vars
    end type grad

    integer::ndet       ! Number of Zombie states
    integer::norb       ! Number of spin orbitals
    integer::nel        ! Number of electrons in molecule
    real(kind=8)::spin  ! Spin of the molecule
    real::beta       ! Distance proagated in imaginary time
    integer::timesteps  ! Number of time steps
    character(LEN=2)::zst !Type of zombie state to be generated

    character(LEN=1)::GDflg      ! Flag to decide if to use Gradient descent 
    character(LEN=1)::zomgflg    ! Flag to generate zombie states or not
    character(LEN=1)::hamgflg    ! Flag to generate Hamiltonian or not
    character(LEN=1)::propflg    ! Flag to propagate in imaginary time or not
    character(LEN=1)::cleanflg   ! Flag to determine if cleaning is to occur
    character(LEN=1)::gramflg    ! Flag to determine if gram schmidt orthogolnalisation should be carried out
    integer::gramnum             ! Number of additional states to be generated for GS orthogonalisation
    character(LEN=1)::rhf_1      ! Flag to decide if the first zombie state should be st as the RHF determinant
    character(LEN=1)::imagflg    ! Flag to decide if zombie states are imaginary or real
    integer::bb_improv           ! Flag to improve the biasing (not yet implemented)
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





