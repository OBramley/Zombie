MODULE globvars

    use mod_types

   
    implicit none 

    integer, parameter::wp = dp ! Sets the precision of the program to double precision

   
    ! Type defining the zombie state
    type zombiest
        real(wp),dimension(:),allocatable::phi
        ! real(wp),dimension(:),allocatable::img
        real(wp),dimension(:),allocatable::val
        real(wp)::num
        real(wp)::num_strt
    end type zombiest

    interface val_set
        module procedure val_set_single
        module procedure val_set_whole
    end interface val_set

   
    ! Type defining the Hamiltonian matrix and the overlap matrix
    type hamiltonian
        real(wp), dimension(:,:), allocatable::hjk
        real(wp), dimension(:,:), allocatable::ovrlp
        real(wp), dimension(:,:), allocatable::inv
        real(wp), dimension(:,:), allocatable::kinvh
    end type hamiltonian

    type dvector
        integer::n
        real(wp), dimension(:), allocatable::d
        real(wp), dimension(:), allocatable::d_1
        real(wp), dimension(:,:),allocatable::d_gs
        real(wp):: norm
        integer(int8):: d_o_d
    end type dvector

    ! Type defining the 1&2 electron integrals
    type elecintrgl
        integer::num
        real(wp), dimension(:), allocatable::integrals
        integer,dimension(:,:),allocatable::orbital_choice
        integer,dimension(:,:),allocatable::orbital_choice2
        integer,dimension(:),allocatable::orbital_choice3
        real(wp) :: hnuc
    end type elecintrgl

    type oprts
        integer(int16),dimension(:,:), allocatable::alive,dead
        integer(int8),dimension(:,:), allocatable::neg_alive,neg_dead
    end type oprts
    
    type grad_do
        real(dp), dimension(:,:), allocatable::hjk
        real(dp), dimension(:,:), allocatable::ovrlp
        real(dp), dimension(:,:), allocatable::inv
        real(dp), dimension(:,:), allocatable::kinvh
        type(zombiest)::zom
        type(dvector)::dvec
        real(dp)::erg
    end type grad_do

    type grad 
        real(wp),dimension(:,:), allocatable::vars
        real(wp):: prev_erg
        real(wp):: current_erg
        integer,dimension(:,:),allocatable::grad_avlb
        real(wp),dimension(:,:,:),allocatable::ovrlp_grad
        integer,dimension(:,:,:),allocatable::ovrlp_grad_avlb
    end type grad

    type gram
        integer::state_num
        type(zombiest),dimension(:),allocatable::zstore
        type(dvector)::dvec
        type(hamiltonian)::haml
        type(grad)::grads
        real(wp),dimension(:,:,:),allocatable::wf_ovrlp
    end type gram 

    type trial_data
        real(wp)::dead
        real(wp)::alive
        real(wp)::new_ham
        real(wp)::old_ham
        integer::orbital
        integer::diag
        integer::rhf
        integer::z1
        integer::z2
        real(wp),dimension(:),allocatable::input_features
    end type trial_data

    type neural_network_layer
        real(wp),dimension(:,:),allocatable::weights
        real(wp),dimension(:),allocatable::biases
        real(wp),dimension(:),allocatable::hidden_layer
        real(wp)::output
    end type neural_network_layer

    
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
    integer::gramwave            ! Number of wave functions to be generated for GS orthogonalisation
    character(LEN=1)::rhf_1      ! Flag to decide if the first zombie state should be st as the RHF determinant
    character(LEN=1)::imagflg    ! Flag to decide if zombie states are imaginary or real
    integer(wp):: randseed       ! Random seed
    real(wp), parameter::pirl  = ATAN(1.0d0)*4.0d0      ! pi
    real(wp), parameter :: sqrtpi = sqrt(pirl)  ! square root of pi
    complex(wp)::i      ! The imagianry unit

    integer:: errorflag      ! Error flag
    

    !!!!!! Gradient descent variables!!!!!!
    real(wp)::lr  ! learning rate
    real(wp)::lr_alpha ! learning rate reduction
    integer::epoc_max ! maximum number of epochs
    integer::lr_loop_max    ! maximum number of learning rates
    integer::ndet_increase ! number of determinants to increase by
    integer::blind_clone_num !number of steps before blind cloning
    integer::ndet_max ! maximum number of determinants
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer::num_devices     !number of devices 
    integer::max_threads
    integer::max_teams
    integer::threadpteam

   
    
    contains

    subroutine val_set_whole(this)
        implicit none
        class(zombiest),intent(inout)::this
        this%val(0)=0.0d0
        
      
        this%val(1:norb)=sin(this%phi)
        this%val(1+norb:2*norb)=cos(this%phi)
    
        return

    end subroutine val_set_whole

    subroutine val_set_single(this,n)
        implicit none
        class(zombiest),intent(inout)::this
        integer,intent(in)::n

        this%val(n)=sin(this%phi(n))
        this%val(n+norb)=cos(this%phi(n))
       
    
        return
    end subroutine val_set_single

    function sign_d_o_d(x) result(s)

        implicit none 
        real(wp),intent(in)::x
        integer(int8)::s

        if(x>0) then
            s=-1
        else if(x<0) then
            s=1
        else
            s=0
        end if

    end function sign_d_o_d

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
        gramwave=0
        lr=0  
        lr_alpha=0 
        epoc_max=0 
        lr_loop_max=0    
        ndet_increase=0 
        blind_clone_num=0
        ndet_max=0 
        rstrtflg='n'
        errorflag=0

        i = (0.0d0,1.0d0)

    end subroutine initialise

END MODULE globvars





