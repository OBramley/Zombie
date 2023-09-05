MODULE globvars

    use mod_types
    use dnad
   
    implicit none 

    integer, parameter::wp = dp ! Sets the precision of the program to double precision


    ! Type defining the zombie state
    type zombiest
        type(dual),dimension(:),allocatable::phi
        ! real(wp),dimension(:),allocatable::img
        type(dual),dimension(:),allocatable::val
    end type zombiest

   
    ! Type defining the Hamiltonian matrix and the overlap matrix
    type hamiltonian
        type(dual), dimension(:,:), allocatable::hjk
        type(dual), dimension(:,:), allocatable::ovrlp
        type(dual), dimension(:,:), allocatable::inv
        type(dual), dimension(:,:), allocatable::kinvh
        real(wp), dimension(:,:,:), allocatable::diff_hjk !ZS to be differentiated,zs is bra or ket, orbital, bra/ket pairing
        real(wp), dimension(:,:,:), allocatable::diff_ovrlp!ZS to be differntiated, orbital, bra/ket pairing
        ! real(wp), dimension(:,:,:,:),allocatable::diff_ov_dov
        ! real(wp), dimension(:,:,:,:),allocatable::diff_in_dhjk
        !The inverse of the overlap matrix multiplied by the hamiltonian
        !diff_invh(j,l,k,m) j specifies the dependnce on zs_j, l,k give the

        real(wp), dimension(:,:,:), allocatable::gs_ovrlp
        real(wp), dimension(:,:,:), allocatable::gs_kinvh
        real(wp), dimension(:,:,:), allocatable::gs_ovrlp_self
        integer::gram_num
    end type hamiltonian

    type dvector
        integer::n
        type(dual), dimension(:), allocatable::d
        real(wp), dimension(:,:),allocatable::d_gs
        type(dual):: norm
    end type dvector

    ! Type defining the 1&2 electron integrals
    type elecintrgl
        integer::h1_num
        integer::h2_num
        real(wp), dimension(:), allocatable::h1ei
        real(wp), dimension(:), allocatable::h2ei
        real(wp) :: hnuc
        
    end type elecintrgl

    type oprts_2
        integer(int16),dimension(:,:), allocatable::alive,dead
        integer(int8),dimension(:,:), allocatable::neg_alive,neg_dead
    end type oprts_2

    type oprts
        type(oprts_2)::ham 
        type(oprts_2),allocatable,dimension(:)::diff
        integer,dimension(:,:), allocatable::dcnt
    end type oprts

    
    type grad 
        real(wp),dimension(:,:), allocatable::vars
        real(wp):: prev_erg
        real(wp):: current_erg
        integer,dimension(:),allocatable::grad_avlb
    end type grad

    
    integer::ndet       ! Number of Zombie states
    integer::norb       ! Number of spin orbitals
    integer::nel        ! Number of electrons in molecule
    real(wp)::spin  ! Spin of the molecule
    real::beta       ! Distance proagated in imaginary time
    integer::timesteps  ! Number of time steps
    character(LEN=2)::zst !Type of zombie state to be generated

    character(LEN=1)::rstrtflg   ! Flag to restart program
    character(LEN=1)::GDflg      ! Flag to decide if to use Gradient descent 
    character(LEN=1)::zomgflg    ! Flag to generate zombie states or not
    character(LEN=1)::hamgflg    ! Flag to generate Hamiltonian or not
    character(LEN=1)::propflg    ! Flag to propagate in imaginary time or not
    character(LEN=1)::cleanflg   ! Flag to determine if cleaning is to occur
    character(LEN=1)::gramflg    ! Flag to determine if gram schmidt orthogolnalisation should be carried out
    character(LEN=1)::GPUflg     ! GPU flag
    integer::gramnum             ! Number of additional states to be generated for GS orthogonalisation
    character(LEN=1)::rhf_1      ! Flag to decide if the first zombie state should be st as the RHF determinant
    character(LEN=1)::imagflg    ! Flag to decide if zombie states are imaginary or real
    real(wp), parameter::pirl  = ATAN(1.0d0)*4.0d0      ! pi
    real(wp), parameter :: sqrtpi = sqrt(pirl)  ! square root of pi
    complex(wp)::i      ! The imagianry unit

    integer:: errorflag      ! Error flag

    integer::num_devices     !number of devices 
    integer::max_threads
    integer::max_teams
    integer::threadpteam

    

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

        i = (0.0d0,1.0d0)

    end subroutine initialise

END MODULE globvars





