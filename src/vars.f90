MODULE globvars

    use mod_types

   
    implicit none 

    integer, parameter::wp = dp ! Sets the precision of the program to double precision

    type :: dual ! make this private will create difficulty to use the
                        ! original write/read commands, hence x and dx are
                        ! variables which can be accessed using D%x and D%dx in
                        ! other units using this module in which D is defined
                        ! as type(dual).
        sequence
        real(wp) :: x        ! functional value
        real(wp), allocatable, dimension(:):: dx  ! derivative

    end type dual

    ! Type defining the zombie state
    type zombiest
        type(dual),dimension(:),allocatable::phi
        ! real(wp),dimension(:),allocatable::img
        type(dual),dimension(:),allocatable::val
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
    
        real(wp), dimension(:,:,:), allocatable::gs_ovrlp
        real(wp), dimension(:,:,:), allocatable::gs_kinvh
        real(wp), dimension(:,:,:), allocatable::gs_ovrlp_self
        integer::gram_num
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
        integer(int16),dimension(:,:), allocatable::alive
        integer(int16),dimension(:,:), allocatable::dead
        integer(int8),dimension(:,:), allocatable::neg_a
        integer(int8),dimension(:,:), allocatable::neg_d
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

    elemental subroutine dx_zero(d)
        type(dual), intent(inout) :: d
        d%dx = 0.0d0

    end subroutine dx_zero

 !-----------------------------------------
    ! COS of dual numbers
    ! <res, dres> = cos(<u, du>) = <cos(u), -sin(u) * du>
    !----------------------------------------
    elemental function cos_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res
        
        allocate(res%dx(size(u%dx)))
        res%x = cos(u%x)
        res%dx = -sin(u%x) * u%dx

    end function cos_d

    !-----------------------------------------
    ! SIN of dual numbers
    ! <res, dres> = sin(<u, du>) = <sin(u), cos(u) * du>
    !----------------------------------------
    elemental function sin_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        allocate(res%dx(size(u%dx)))
        res%x = sin(u%x)
        res%dx = cos(u%x) * u%dx

    end function sin_d

    subroutine val_set_whole(this)
        implicit none
        class(zombiest),intent(inout)::this
      
        this%val(0)%x=0.0d0
  
        this%val(1:norb)=sin_d(this%phi)
        this%val(1+norb:2*norb)=cos_d(this%phi)
    
        return

    end subroutine val_set_whole

    subroutine val_set_single(this,n)
        implicit none
        class(zombiest),intent(inout)::this
        integer,intent(in)::n

     
        this%val(n)=sin_d(this%phi(n))
        this%val(n+norb)=cos_d(this%phi(n))
       
    
        return
    end subroutine val_set_single

    function sign_d_o_d(x) result(s)

        implicit none 
        real(wp),intent(in)::x
        integer(int8)::s

        if(x>0) then
            s=1
        else if(x<0) then
            s=-1
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

        errorflag=0

        i = (0.0d0,1.0d0)

    end subroutine initialise

END MODULE globvars





