MODULE devvars

    use globvars

    ! Type defining the zombie state
    type zombiest_dev
        ! complex(kind=8), dimension(:), allocatable,device::sin
        ! complex(kind=8), dimension(:), allocatable,device::cos
        real(kind=8), dimension(:), allocatable,device::sin
        real(kind=8), dimension(:), allocatable,device::cos
        real(kind=8),dimension(:),allocatable,device::phi
        real(kind=8),dimension(:),allocatable,device::img
        real(kind=8),dimension(:),allocatable,device::val
        integer(kind=1),dimension(:),allocatable,device::dead
        integer(kind=1),dimension(:),allocatable,device::alive
        integer::update_num
    end type zombiest_dev

    ! Type defining the Hamiltonian matrix and the overlap matrix
    type hamiltonian_dev
        real(kind=8), dimension(:,:), allocatable,device::hjk
        real(kind=8), dimension(:,:), allocatable,device::ovrlp
        real(kind=8), dimension(:,:), allocatable,device::inv
        real(kind=8), dimension(:,:), allocatable,device::kinvh
        ! complex(kind=8), dimension(:,:), allocatable,device::hjk
        ! complex(kind=8), dimension(:,:), allocatable,device::ovrlp
        ! complex(kind=8), dimension(:,:), allocatable,device::inv
        ! complex(kind=8), dimension(:,:), allocatable,device::kinvh

        real(kind=8), dimension(:,:,:), allocatable,device::diff_hjk !ZS to be differentiated,zs is bra or ket, orbital, bra/ket pairing
        real(kind=8), dimension(:,:,:), allocatable,device::diff_ovrlp!ZS to be differntiated, orbital, bra/ket pairing
        real(kind=8), dimension(:,:,:,:), allocatable,device::hess_hjk !ZS to be differentiated,zs is bra or ket, orbital, bra/ket pairing
        real(kind=8), dimension(:,:,:,:), allocatable,device::hess_ovrlp!ZS to be differntiated, orbital, bra/ket pairing
        real(kind=8), dimension(:,:,:,:), allocatable,device::diff_invh
        real(kind=8), dimension(:,:,:,:),allocatable,device::diff_ov_dov
        real(kind=8), dimension(:,:,:,:),allocatable,device::diff_in_dhjk
        !The inverse of the overlap matrix multiplied by the hamiltonian
        !diff_invh(j,l,k,m) j specifies the dependnce on zs_j, l,k give the 
    end type hamiltonian_dev

    type dvector_dev
        !complex(kind=8), dimension(:), allocatable,device::d
        real(kind=8), dimension(:), allocatable,device::d
        real(kind=8), dimension(:,:,:),allocatable,device::d_diff
        ! d_diff(k,j,m) strucutred k specifies the position in the vector d
        ! j specifies the dependence on ZS j. m corresponds to the coeffcient m within 
        ! zombie stae j
        real(kind=8):: norm
    end type dvector_dev

    type energy_dev
        real(kind=8), dimension(:),allocatable,device::t
        real(kind=8), dimension(:,:),allocatable,device::erg 
        ! complex(kind=8), dimension(:,:),allocatable,device::erg 
    end type energy_dev

    ! Type defining the 1&2 electron integrals
    type elecintrgl_dev
        integer,device::h1_num
        integer,device::h2_num
        ! real(kind=8), dimension(:,:), allocatable,device::h1ei
        ! real(kind=8), dimension(:,:,:,:), allocatable,device::h2ei
        real(kind=8), dimension(:), allocatable,device::h1ei
        real(kind=8), dimension(:), allocatable,device::h2ei
        real(kind=8) :: hnuc
        
    end type elecintrgl_dev

    type oprts_2_dev
        integer(kind=2),dimension(:,:), allocatable,device::alive,dead
        integer(kind=1),dimension(:,:), allocatable,device::neg_alive,neg_dead
    end type oprts_2_dev

    type oprts_dev
        type(oprts_2)::ham 
        type(oprts_2),allocatable,dimension(:)::diff
        type(oprts_2),allocatable,dimension(:,:)::hess
        integer,dimension(:,:), allocatable::dcnt
        integer,dimension(:,:,:), allocatable::hcnt
    end type oprts_dev

    type grad_dev 
        real(kind=8),dimension(:,:), allocatable,device::vars
        real(kind=8),dimension(:,:), allocatable,device::vars_hess
        real(kind=8),dimension(:,:,:), allocatable,device::hessian
        real(kind=8):: prev_erg
        real(kind=8):: current_erg
        integer,dimension(:,:),allocatable,device::grad_avlb
        ! real(kind=8),dimension(:,:),allocatable,device::prev_mmntm
        ! real(kind=8),dimension(:),allocatable,device::hess_sum
    end type grad_dev

    
    integer,device,constant::ndet_dev       ! Number of Zombie states
    integer,device,constant::norb_dev       ! Number of spin orbitals
    integer,device,constant::nel_dev        ! Number of electrons in molecule
    real(kind=8),device,constant::spin_dev  ! Spin of the molecule
    real,device,constant::beta_dev       ! Distance proagated in imaginary time
    integer,device,constant::timesteps_dev  ! Number of time steps
    real(kind=8),device,constant::pirl      ! pi
    integer,device:: errorflag_dev 

    contains

    subroutine initialise_dev
        
        implicit none
        
        ! Initialises all the globabl variables to zero and sets the imaginary unit 

        ndet_dev=ndet
        norb_dev=norb
        nel_dev=nel
        spin_dev=spin
        beta_dev=beta
        timesteps_dev=timesteps
        errorflag_dev=errorflag
        pirl_dev = pirl
       

    end subroutine initialise_dev

END MODULE devvars





