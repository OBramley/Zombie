!******************************************************************************
!* dual Number Automatic Differentiation (DNAD) of Fortran Codes
!*-----------------------------------------------------------------------------
!* COPYRIGHT (c) Joshua Hodson, All rights reserved, you are free to copy,
!* modify, or translate this code to other languages such as c/c++. This is a
!* fork of the original Fortran DNAD module developed by Dr. Wenbin Yu. See
!* original copyright information below. You can download the original version
!* at https://cdmhub.org/resources/374
!*
!* COPYRIGHT (c) Wenbin Yu, All rights reserved, you are free to copy,
!* modify or translate this code to other languages such as c/c++. If
!* you find a bug please let me know through wenbinyu.heaven@gmail.com. If
!* you added new functions and want to share with others, please let me know
!* too. You are welcome to share your successful stories with us through
!* http://groups.google.com/group/hifi-comp.
!******************************************************************************
!* Acknowledgements
!*-----------------------------------------------------------------------------
!* The development of DNAD is supported, in part, by the Chief Scientist
!* Innovative Research Fund at AFRL/RB WPAFB, and by Department of Army
!* SBIR (Topic A08-022) through Advanced Dynamics Inc. The views and
!* conclusions contained herein are those of the authors and should not be
!* interpreted as necessarily representing the official policies or
!* endorsement, either expressed or implied, of the funding agency.
!*
!* Additional development of DNAD has been supported under a Department of
!* Energy (DOE) Nuclear Energy University Program (NEUP) Graduate Fellowship.
!* Any opinions, findings, conclusions or recommendations expressed in this
!* publication are those of the authors and do not necessarily reflect the
!* views of the Department of Energy Office of Nuclear Energy.
!******************************************************************************
!* Citation
!*-----------------------------------------------------------------------------
!* Your citation of the following two papers is appreciated:
!* Yu, W. and Blair, M.: "DNAD, a Simple Tool for Automatic Differentiation of
!* Fortran Codes Using dual Numbers," Computer Physics Communications, vol.
!* 184, 2013, pp. 1446-1452.
!*
!* Spall, R. and Yu, W.: "Imbedded dual-Number Automatic Differentiation for
!* CFD Sensitivity Analysis," Journal of Fluids Engineering, vol. 135, 2013,
!* 014501.
!******************************************************************************
!* Quick Start Guide
!*-----------------------------------------------------------------------------
!* To integrate DNAD into an existing Fortran program, do the following:
!*
!*   1. Include the DNAD module in the source files by adding "use dnadmod" to
!*      the beginning of all modules, global functions, and global subroutines
!*      that include definitions of floating-point variables.
!*   2. Redefine all floating-point variables as type(dual). This can be done
!*      using precompiler directives so that the integration can be turned on
!*      or off at compile-time, eliminating the need for maintaining two
!*      separate code bases for the same project.
!*   3. All I/O involving floating-point variables will need to be examined.
!*      A method will need to be determined for inputting and outputting
!*      derivative values. This customization is typically unique for each
!*      piece of software and needs to be determined on a case-by-case basis.
!*   4. When compiling DNAD, use the compiler option "-Dndv=#", where # is the
!*      number of design variables desired. This sizes the derivative array
!*      that is stored with each floating point number.
!*   5. When compiling DNAD, use compiler options to specify precision. If no
!*      compiler options are specified, DNAD will default to single-precision
!*      floating-point arithmetic. Most popular Fortran compilers provide
!*      options for specifying precision at compile-time so that it does not
!*      have to be hard-coded into the source code. For example, use the
!*      "-fdefault-real-8" compiler in gfortran or the "-r8" compiler option
!*      with Intel Fortran to compile DNAD as double-precision.
!*   6. Modify the compilation process for the target software to include the
!*      DNAD module in the resulting executable or library.
!******************************************************************************
!* Change Log
!*-----------------------------------------------------------------------------
!*  2022-05-05  Nicholas Brady
!*  - added overloading of integers and reals to dual powers: pow_id, pow_rd
!*  - as well as for complex dual numbers: pow_i_dc, pow_r_dc, pow_c_dc
!*
!*  2022-03-24  Nicholas Brady
!*  - added type dual_complex (these have various physical applications)
!*  - extended the overloading of many intrinsic functions (but not all) to
!*      include type dual_complex
!*  - defined erf_c to calculate the error function of a complex number
!*
!*  2021-05-29  Nicholas Brady
!*  - overloaded hyperbolic, inverse hyperbolic, and error functions (intrinsic)
!*      sinh, cosh, tanh, asinh, acosh, atanh, erf
!*
!*  2016-04-29  Joshua Hodson
!*  - Updated copyright, acknowledgments, and quick start guide.
!*  - Removed overloads for single-precision reals.
!*  - Added tan, dtan, atan, and atan2 intrinsic function overloads.
!*  - Removed macro for precision and defined all floating-point variables as
!*    default real. Compiler options can now be used to set precision.
!*  - Added checks for undefined derivatives when only constants are used in
!*    the calculation (i.e. all partial derivatives are zero). This limits the
!*    perpetuation of NaN values in the code.
!*  - Combined the header and source files into a single file.
!*
!*  2015-07-29  Joshua Hodson
!*  - Added maxloc intrinsic function overload.
!*  - Converted UPPERCASE to lowercase for readability.
!*  - Added macros for defining precision and number of design variables.
!*  - Renamed module from dual_Num_Auto_Diff to dnadmod
!*  - Renamed dual number type from dual_NUM to dual
!*  - Renamed components of dual number type from (xp_ad_, xp_ad_) to (x, dx)
!*
!*  2014-06-05  Wenbin Yu
!*  - Forked from original DNAD repository, see https://cdmhub.org/resources/374
!*
!******************************************************************************

! Number of design variables (default = 1)
! #ifndef ndv
! #define ndv 1
! #endif

module dnad
    
    use mod_types, only: wp=>dp, int16, int8 
    use dual_set

    implicit none
    integer, PARAMETER :: dual_size2=dual_size*2

    private

    real(wp) :: negative_one = -1.0d0
    real(wp), parameter :: PI = ATAN(1.0d0)*4.0      ! pi - Geometric constant
   
   

    ! type,public :: dual_complex ! make this private will create difficulty to use the
    !                     ! original write/read commands, hence z and dz are
    !                     ! variables which can be accessed using D%z and D%dz in
    !                     ! other units using this module in which D is defined
    !                     ! as type(dual_complex).
    !     sequence
    !     complex(wp) :: z        ! functional value
    !     complex(wp):: dz(dual_size)  ! derivative
    ! end type dual_complex

    ! type,public :: dual2 ! make this private will create difficulty to use the
    !                     ! original write/read commands, hence x and dx are
    !                     ! variables which can be accessed using D%x and D%dx in
    !                     ! other units using this module in which D is defined
    !                     ! as type(dual).
    !     sequence
    !     real(wp) :: x        ! functional value
    !     real(wp), dimension(dual_size2):: dx  ! derivative

    ! end type dual2

    ! type,public :: dual_complex2 ! make this private will create difficulty to use the
    !                     ! original write/read commands, hence z and dz are
    !                     ! variables which can be accessed using D%z and D%dz in
    !                     ! other units using this module in which D is defined
    !                     ! as type(dual_complex).
    !     sequence
    !     complex(wp) :: z        ! functional value
    !     complex(wp):: dz(dual_size2)  ! derivative
    ! end type dual_complex2


!******** Interfaces for operator overloading
    public assignment (=)
    interface assignment (=)
        module procedure assign_di ! dual=integer, elemental
        module procedure assign_dr ! dual=real(wp) elemental
        module procedure assign_id ! integer=dual, elemental
        module procedure assign_rd ! real(wp)=dual, elemental
        
        ! module procedure assign_dc_i    ! dual_complex = integer
        ! module procedure assign_dc_r    ! dual_complex = real
        ! module procedure assign_dc_c    ! dual_complex = complex
        ! ! module procedure assign_c_dc    ! complex = dual_complex
        
        module procedure assign_di2 ! dual=integer, elemental
        module procedure assign_dr2 ! dual=real(wp) elemental
        module procedure assign_id2 ! integer=dual, elemental
        module procedure assign_rd2 ! real(wp)=dual, elemental
        
        ! module procedure assign_dc_i2    ! dual_complex = integer
        ! module procedure assign_dc_r2    ! dual_complex = real
        ! module procedure assign_dc_c2    ! dual_complex = complex
        ! module procedure assign_c_dc2    ! complex = dual_complex
    end interface

    public operator (+)
    interface operator (+)
        module procedure add_d  ! +dual number, elemental
        module procedure add_dd ! dual + dual, elemental
        module procedure add_di ! dual + integer, elemental
        module procedure add_dr ! dual + real, elemental
        module procedure add_id ! integer + dual, elemental
        module procedure add_rd ! real + dual, elemental

        ! module procedure add_dc     ! +dual_complex number, elemental
        ! module procedure add_dc_dc  ! dual_complex + dual_complex
        ! module procedure add_c_dc   ! complex + dual_complex
        ! module procedure add_dc_c   ! dual_complex + complex
        ! module procedure add_i_dc   ! integer + dual_complex
        ! module procedure add_dc_i   ! dual_complex + integer
        ! module procedure add_r_dc   ! real + dual_complex
        ! module procedure add_dc_r   ! dual_complex + real

        module procedure add_d2  ! +dual number, elemental
        module procedure add_dd2 ! dual + dual, elemental
        module procedure add_di2 ! dual + integer, elemental
        module procedure add_dr2 ! dual + real, elemental
        module procedure add_id2 ! integer + dual, elemental
        module procedure add_rd2 ! real + dual, elemental
        
        ! module procedure add_dc2     ! +dual_complex number, elemental
        ! module procedure add_dc_dc2  ! dual_complex + dual_complex
        ! module procedure add_c_dc2   ! complex + dual_complex
        ! module procedure add_dc_c2   ! dual_complex + complex
        ! module procedure add_i_dc2   ! integer + dual_complex
        ! module procedure add_dc_i2   ! dual_complex + integer
        ! module procedure add_r_dc2   ! real + dual_complex
        ! module procedure add_dc_r2   ! dual_complex + real
    end interface

    public operator (-)
    interface operator (-)
        module procedure minus_d  ! negate a dual number,elemental
        module procedure minus_dd ! dual -dual,elemental
        module procedure minus_di ! dual-integer,elemental
        module procedure minus_dr ! dual-real,elemental
        module procedure minus_id ! integer-dual,elemental
        module procedure minus_rd ! real-dual,elemental

        ! module procedure minus_dc     ! -dual_complex number, elemental
        ! module procedure minus_dc_dc  ! dual_complex - dual_complex
        ! module procedure minus_c_dc   ! complex - dual_complex
        ! module procedure minus_dc_c   ! dual_complex - complex
        ! module procedure minus_i_dc   ! integer - dual_complex
        ! module procedure minus_dc_i   ! dual_complex - integer
        ! module procedure minus_r_dc   ! real - dual_complex
        ! module procedure minus_dc_r   ! dual_complex - real

        module procedure minus_d2  ! negate a dual number,elemental
        module procedure minus_dd2 ! dual -dual,elemental
        module procedure minus_di2 ! dual-integer,elemental
        module procedure minus_dr2 ! dual-real,elemental
        module procedure minus_id2 ! integer-dual,elemental
        module procedure minus_rd2 ! real-dual,elemental

        ! module procedure minus_dc2     ! -dual_complex number, elemental
        ! module procedure minus_dc_dc2  ! dual_complex - dual_complex
        ! module procedure minus_c_dc2   ! complex - dual_complex
        ! module procedure minus_dc_c2   ! dual_complex - complex
        ! module procedure minus_i_dc2   ! integer - dual_complex
        ! module procedure minus_dc_i2   ! dual_complex - integer
        ! module procedure minus_r_dc2   ! real - dual_complex
        ! module procedure minus_dc_r2   ! dual_complex - real
    end interface

    public operator (*)
    interface operator (*)
        module procedure mult_dd ! dual*dual, elemental
        module procedure mult_di ! dual*integer,elemental
        module procedure mult_di16 ! dual*int16,elemental
        module procedure mult_di8 ! dual*int8,elemental
        module procedure mult_dr ! dual*real,elemental
        module procedure mult_id ! integer*dual,elemental
        module procedure mult_rd ! real*dual,elemental

        ! module procedure mult_dc_dc ! dual_complex*dual_complex, elemental
        ! module procedure mult_dc_i ! dual_complex*integer,elemental
        ! module procedure mult_dc_r ! dual_complex*real,elemental
        ! module procedure mult_i_dc ! integer*dual,elemental
        ! module procedure mult_r_dc ! real*dual_complex,elemental
        ! module procedure mult_dc_c ! dual_complex*complex,elemental
        ! module procedure mult_c_dc ! complex*dual_complex,elemental

        module procedure mult_dd2 ! dual*dual, elemental
        module procedure mult_di2 ! dual*integer,elemental
        module procedure mult_di162 ! dual*int16,elemental
        module procedure mult_di82 ! dual*int8,elemental
        module procedure mult_dr2 ! dual*real,elemental
        module procedure mult_id2 ! integer*dual,elemental
        module procedure mult_rd2 ! real*dual,elemental

        ! module procedure mult_dc_dc2 ! dual_complex*dual_complex, elemental
        ! module procedure mult_dc_i2 ! dual_complex*integer,elemental
        ! module procedure mult_dc_r2 ! dual_complex*real,elemental
        ! module procedure mult_i_dc2 ! integer*dual,elemental
        ! module procedure mult_r_dc2 ! real*dual_complex,elemental
        ! module procedure mult_dc_c2 ! dual_complex*complex,elemental
        ! module procedure mult_c_dc2 ! complex*dual_complex,elemental
    end interface

    public operator (/)
    interface operator (/)
        module procedure div_dd ! dual/dual,elemental
        module procedure div_di ! dual/integer, elemental
        module procedure div_dr ! dual/real,emental
        module procedure div_id ! integer/dual, elemental
        module procedure div_rd ! real/dual, elemental

        ! module procedure div_dc_dc ! dual_complex/dual_complex,elemental
        ! module procedure div_dc_i ! dual_complex/integer, elemental
        ! module procedure div_dc_r ! dual_complex/real,emental
        ! module procedure div_i_dc ! integer/dual, elemental
        ! module procedure div_r_dc ! real/dual_complex, elemental
        ! module procedure div_dc_c ! dual_complex/complex,emental
        ! module procedure div_c_dc ! complex/dual_complex, elemental

        module procedure div_dd2 ! dual/dual,elemental
        module procedure div_di2 ! dual/integer, elemental
        module procedure div_dr2 ! dual/real,emental
        module procedure div_id2 ! integer/dual, elemental
        module procedure div_rd2 ! real/dual, elemental

        ! module procedure div_dc_dc2 ! dual_complex/dual_complex,elemental
        ! module procedure div_dc_i2 ! dual_complex/integer, elemental
        ! module procedure div_dc_r2 ! dual_complex/real,emental
        ! module procedure div_i_dc2 ! integer/dual, elemental
        ! module procedure div_r_dc2 ! real/dual_complex, elemental
        ! module procedure div_dc_c2 ! dual_complex/complex,emental
        ! module procedure div_c_dc2 ! complex/dual_complex, elemental
    end interface

    public operator (**)
    interface operator (**)
        module procedure pow_di ! dual number to an integer power,elemental
        module procedure pow_dr ! dual number to a real power, elemental
        module procedure pow_dd ! dual number to a dual power, elemental
        module procedure pow_id ! integer to a dual power, elemental
        module procedure pow_rd ! real to a dual power, elemental

        ! module procedure pow_dc_i ! dual_complex number to an integer power,elemental
        ! module procedure pow_dc_r ! dual_complex number to a real power, elemental
        ! module procedure pow_dc_c ! dual_complex number to a complex power, elemental
        ! module procedure pow_dc_dc ! dual_complex number to a dual_complex power, elemental
        ! module procedure pow_i_dc ! integer to a dual_complex power, elemental
        ! module procedure pow_r_dc ! real to a dual_complex power, elemental
        ! module procedure pow_c_dc ! complex number to a dual_complex power, elemental

        module procedure pow_di2 ! dual number to an integer power,elemental
        module procedure pow_dr2 ! dual number to a real power, elemental
        module procedure pow_dd2 ! dual number to a dual power, elemental
        module procedure pow_id2 ! integer to a dual power, elemental
        module procedure pow_rd2 ! real to a dual power, elemental

        ! module procedure pow_dc_i2 ! dual_complex number to an integer power,elemental
        ! module procedure pow_dc_r2 ! dual_complex number to a real power, elemental
        ! module procedure pow_dc_c2 ! dual_complex number to a complex power, elemental
        ! module procedure pow_dc_dc2 ! dual_complex number to a dual_complex power, elemental
        ! module procedure pow_i_dc2 ! integer to a dual_complex power, elemental
        ! module procedure pow_r_dc2 ! real to a dual_complex power, elemental
        ! module procedure pow_c_dc2 ! complex number to a dual_complex power, elemental
    end interface

    public operator (==)
    interface operator (==)
        module procedure eq_dd ! compare two dual numbers, elemental
        module procedure eq_di ! compare a dual and an integer, elemental
        module procedure eq_dr ! compare a dual and a real, elemental
        module procedure eq_id ! compare integer with a dual number, elemental
        module procedure eq_rd ! compare a real with a dual number, elemental

        ! module procedure eq_dc_dc ! compare two dual_complex numbers, elemental
        ! module procedure eq_dc_i ! compare a dual_complex and an integer, elemental
        ! module procedure eq_dc_r ! compare a dual_complex and a real, elemental
        ! module procedure eq_i_dc ! compare integer with a dual_complex number, elemental
        ! module procedure eq_r_dc ! compare a real with a dual_complex number, elemental

        module procedure eq_dd2 ! compare two dual numbers, elemental
        module procedure eq_di2 ! compare a dual and an integer, elemental
        module procedure eq_dr2 ! compare a dual and a real, elemental
        module procedure eq_id2 ! compare integer with a dual number, elemental
        module procedure eq_rd2 ! compare a real with a dual number, elemental

        ! module procedure eq_dc_dc2 ! compare two dual_complex numbers, elemental
        ! module procedure eq_dc_i2 ! compare a dual_complex and an integer, elemental
        ! module procedure eq_dc_r2 ! compare a dual_complex and a real, elemental
        ! module procedure eq_i_dc2 ! compare integer with a dual_complex number, elemental
        ! module procedure eq_r_dc2 ! compare a real with a dual_complex number, elemental
    end interface

    public operator (<=)
    interface operator (<=)
        module procedure le_dd ! compare two dual numbers, elemental
        module procedure le_di ! compare a dual and an integer, elemental
        module procedure le_dr ! compare a dual and a real,elemental
        module procedure le_id ! compare integer with a dual number, elemental
        module procedure le_rd ! compare a real with a dual number, elemental

        module procedure le_dd2 ! compare two dual numbers, elemental
        module procedure le_di2 ! compare a dual and an integer, elemental
        module procedure le_dr2 ! compare a dual and a real,elemental
        module procedure le_id2 ! compare integer with a dual number, elemental
        module procedure le_rd2 ! compare a real with a dual number, elemental

        ! <= cannot be done for complex numbers
    end interface

    public operator (<)
    interface operator (<)
        module procedure lt_dd ! compare two dual numbers, elemental
        module procedure lt_di ! compare a dual and an integer, elemental
        module procedure lt_dr ! compare dual with a real, elemental
        module procedure lt_id ! compare integer with a dual number, elemental
        module procedure lt_rd ! compare a real with a dual number, elemental
        module procedure lt_dd2 ! compare two dual numbers, elemental
        module procedure lt_di2 ! compare a dual and an integer, elemental
        module procedure lt_dr2 ! compare dual with a real, elemental
        module procedure lt_id2 ! compare integer with a dual number, elemental
        module procedure lt_rd2 ! compare a real with a dual number, elemental
        ! < cannot be done for complex numbers
    end interface

    public operator (>=)
    interface operator (>=)
        module procedure ge_dd ! compare two dual numbers, elemental
        module procedure ge_di ! compare dual with integer, elemental
        module procedure ge_dr ! compare dual with a real number, elemental
        module procedure ge_id ! compare integer with a dual number, elemental
        module procedure ge_rd ! compare a real with a dual number, elemental
        module procedure ge_dd2 ! compare two dual numbers, elemental
        module procedure ge_di2 ! compare dual with integer, elemental
        module procedure ge_dr2 ! compare dual with a real number, elemental
        module procedure ge_id2 ! compare integer with a dual number, elemental
        module procedure ge_rd2 ! compare a real with a dual number, elemental

        ! >= cannot be done for complex numbers
    end interface

    public operator (>)
    interface operator (>)
        module procedure gt_dd ! compare two dual numbers, elemental
        module procedure gt_di ! compare a dual and an integer, elemental
        module procedure gt_dr ! compare dual with a real, elemental
        module procedure gt_id ! compare integer with a dual number, elemental
        module procedure gt_rd ! compare a real with a dual number, elemental
        module procedure gt_dd2 ! compare two dual numbers, elemental
        module procedure gt_di2 ! compare a dual and an integer, elemental
        module procedure gt_dr2 ! compare dual with a real, elemental
        module procedure gt_id2 ! compare integer with a dual number, elemental
        module procedure gt_rd2 ! compare a real with a dual number, elemental
        

        ! > cannot be done for complex numbers
    end interface

    public operator (/=)
    interface operator (/=)
        module procedure ne_dd ! compare two dual numbers, elemental
        module procedure ne_di ! compare a dual and an integer, elemental
        module procedure ne_dr ! compare dual with a real, elemental
        module procedure ne_id ! compare integer with a dual number, elemental
        module procedure ne_rd ! compare a real with a dual number, elemental

        ! module procedure ne_dc_dc ! compare two dual_complex numbers, elemental
        ! module procedure ne_dc_i ! compare a dual_complex and an integer, elemental
        ! module procedure ne_dc_r ! compare dual_complex with a real, elemental
        ! module procedure ne_i_dc ! compare integer with a dual_complex number, elemental
        ! module procedure ne_r_dc ! compare a real with a dual_complex number, elemental
        module procedure ne_dd2 ! compare two dual numbers, elemental
        module procedure ne_di2 ! compare a dual and an integer, elemental
        module procedure ne_dr2 ! compare dual with a real, elemental
        module procedure ne_id2 ! compare integer with a dual number, elemental
        module procedure ne_rd2 ! compare a real with a dual number, elemental

        ! module procedure ne_dc_dc2 ! compare two dual_complex numbers, elemental
        ! module procedure ne_dc_i2 ! compare a dual_complex and an integer, elemental
        ! module procedure ne_dc_r2 ! compare dual_complex with a real, elemental
        ! module procedure ne_i_dc2 ! compare integer with a dual_complex number, elemental
        ! module procedure ne_r_dc2 ! compare a real with a dual_complex number, elemental
    end interface


!------------------------------------------------
! Interfaces for intrinsic functions overloading
!------------------------------------------------
    public abs
    interface abs
        module procedure abs_d  ! absolute value of a dual number, elemental
        ! module procedure abs_dc
        module procedure abs_d2  ! absolute value of a dual number, elemental
        ! module procedure abs_dc2
    end interface

    public dabs
    interface dabs
        module procedure abs_d  ! same as abs, used for some old fortran commands
        ! module procedure abs_dc
        module procedure abs_d2  ! same as abs, used for some old fortran commands
        ! module procedure abs_dc2
    end interface

    public acos
    interface acos
        module procedure acos_d ! arccosine of a dual number, elemental
        ! module procedure acos_dc ! arccosine of a dual number, elemental
        module procedure acos_d2 ! arccosine of a dual number, elemental
        ! module procedure acos_dc2 ! arccosine of a dual number, elemental
    end interface

    public asin
    interface asin
        module procedure asin_d ! arcsine of a dual number, elemental
        ! module procedure asin_dc ! arcsine of a dual number, elemental

        module procedure asin_d2 ! arcsine of a dual number, elemental
        ! module procedure asin_dc2 ! arcsine of a dual number, elemental
    end interface

    public atan
    interface atan
        module procedure atan_d ! arctan of a dual number, elemental
        ! module procedure atan_dc ! arctan of a dual number, elemental
        module procedure atan_d2 ! arctan of a dual number, elemental
        ! module procedure atan_dc2 ! arctan of a dual number, elemental
         
    end interface

    public atan2
    interface atan2
        module procedure atan2_d ! arctan of a dual number, elemental
        module procedure atan2_d2 ! arctan of a dual number, elemental

    end interface

    public cos
    interface cos
        module procedure cos_d ! cosine of a dual number, elemental
        ! module procedure cos_dc ! cosine of a compulex_dual number, elemental
        module procedure cos_d2 ! cosine of a dual number, elemental
        ! module procedure cos_dc2 ! cosine of a compulex_dual number, elemental
    end interface

    public dcos
    interface dcos
        module procedure cos_d ! cosine of a dual number, elemental
        ! module procedure cos_dc ! cosine of a compulex_dual number, elemental
        module procedure cos_d2 ! cosine of a dual number, elemental
        ! module procedure cos_dc2 ! cosine of a compulex_dual number, elemental
    end interface

    public dot_product
    interface dot_product
        module procedure dot_product_dd ! dot product two dual number vectors
        module procedure exp_d2 ! exponential of a dual number, elemental
        ! module procedure exp_dc2 ! exponential of a dual number, elemental
    end interface

    public exp
    interface exp
        module procedure exp_d ! exponential of a dual number, elemental
        ! module procedure exp_dc ! exponential of a dual number, elemental
        module procedure exp_d2 ! exponential of a dual number, elemental
        ! module procedure exp_dc2 ! exponential of a dual number, elemental
    end interface

    public int
    interface int
        module procedure int_d ! integer part of a dual number, elemental
        module procedure int_d2 ! integer part of a dual number, elemental

    end interface

    public log
    interface log
        module procedure log_d ! log of a dual number, elemental
        ! module procedure log_dc ! log of a dual_complex number, elemental
        module procedure log_d2 ! log of a dual number, elemental
        ! module procedure log_dc2 ! log of a dual_complex number, elemental
    end interface

    public log10
    interface log10
        module procedure log10_d ! log10 of a dual number, elemental
        ! module procedure log10_dc ! log10 of a dual_complex number, elemental
        module procedure log10_d2 ! log10 of a dual number, elemental
        ! module procedure log10_dc2 ! log10 of a dual_complex number, elemental
    end interface

    public matmul
    interface matmul
        module procedure matmul_dd ! multiply two dual matrices
        module procedure matmul_dv ! multiply a dual matrix with a dual vector
        module procedure matmul_vd ! multiply a dual vector with a dual matrix
        module procedure matmul_dd2 ! multiply two dual matrices
        module procedure matmul_dv2 ! multiply a dual matrix with a dual vector
        module procedure matmul_vd2 ! multiply a dual vector with a dual matrix
    end interface


    public max
    interface max
        module procedure max_dd ! max of from two to four dual numbers, elemental
        module procedure max_di ! max of a dual number and an integer, elemental
        module procedure max_dr ! max of a dual number and a real, elemental
        module procedure max_rd ! max of a real,and a dual number,  elemental
        module procedure max_dd2 ! max of from two to four dual numbers, elemental
        module procedure max_di2 ! max of a dual number and an integer, elemental
        module procedure max_dr2 ! max of a dual number and a real, elemental
        module procedure max_rd2 ! max of a real,and a dual number,  elemental
    end interface

    public dmax1
    interface dmax1
        module procedure max_dd ! max of from two to four dual numbers, elemental
        module procedure max_dd2 ! max of from two to four dual numbers, elemental
    end interface

    public maxval
    interface maxval
        module procedure maxval_d ! maxval of a dual number vector
        module procedure maxval_d2 ! maxval of a dual number vector
    end interface

    public min
    interface min
        module procedure min_dd ! min of from two to four dual numbers, elemental
        module procedure min_dr ! min of a dual and a real, elemental
        module procedure min_dd2 ! min of from two to four dual numbers, elemental
        module procedure min_dr2 ! min of a dual and a real, elemental
    end interface

    public dmin1
    interface dmin1
        module procedure min_dd ! min of from two to four dual numbers, elemental
        module procedure min_dd2 ! min of from two to four dual numbers, elemental
    end interface

    public minval
    interface minval
        module procedure minval_d ! obtain the maxval  of a dual number vectgor
        module procedure minval_d2 ! obtain the maxval  of a dual number vectgor
    end interface

    public nint
    interface nint
        module procedure nint_d ! nearest integer to the argument, elemental
        module procedure nint_d2 ! nearest integer to the argument, elemental
    end interface

    public sign
    interface  sign
        module procedure  sign_dd ! sign(a,b) with two dual numbers, elemental
        module procedure  sign_rd ! sign(a,b) with a real and a dual, elemental
        module procedure  sign_dd2 ! sign(a,b) with two dual numbers, elemental
        module procedure  sign_rd2 ! sign(a,b) with a real and a dual, elemental
    end interface

    public sin
    interface sin
        module procedure sin_d ! obtain sine of a dual number, elemental
        ! module procedure sin_dc ! obtain sine of a dual_complex number, elemental
        module procedure sin_d2 ! obtain sine of a dual number, elemental
        ! module procedure sin_dc2 ! obtain sine of a dual_complex number, elemental
    end interface

    public dsin
    interface dsin
        module procedure sin_d ! obtain sine of a dual number, elemental
        ! module procedure sin_dc ! obtain sine of a dual_complex number, elemental
        module procedure sin_d2 ! obtain sine of a dual number, elemental
        ! module procedure sin_dc2 ! obtain sine of a dual_complex number, elemental
    end interface

    public tan
    interface tan
        module procedure tan_d ! obtain sine of a dual number, elemental
        ! module procedure tan_dc ! obtain sine of a dual_complex number, elemental
        module procedure tan_d2 ! obtain sine of a dual number, elemental
        ! module procedure tan_dc2 ! obtain sine of a dual_complex number, elemental
    end interface


    public dtan
    interface dtan
        module procedure tan_d ! obtain sine of a dual number, elemental
        ! module procedure tan_dc ! obtain sine of a dual_complex number, elemental
        module procedure tan_d2 ! obtain sine of a dual number, elemental
        ! module procedure tan_dc2 ! obtain sine of a dual_complex number, elemental
    end interface

    public sqrt
    interface sqrt
        module procedure sqrt_d ! obtain the sqrt of a dual number, elemental
        ! module procedure sqrt_dc ! obtain the sqrt of a dual_complex number, elemental
            module procedure sqrt_d2 ! obtain the sqrt of a dual number, elemental
                ! module procedure sqrt_dc2 ! obtain the sqrt of a dual_complex number, elemental
    end interface

    public sum
    interface sum
        module procedure sum_d ! sum a dual array
        ! module procedure sum_dc ! sum a dual_complex array
        module procedure sum_d2 ! sum a dual array
        ! module procedure sum_dc2 ! sum a dual_complex array
    end interface

    public maxloc
    interface maxloc
        module procedure maxloc_d ! location of max in a dual array
        module procedure maxloc_d2 ! location of max in a dual array
    end interface

    public sinh
    interface sinh
        module procedure sinh_d ! obtain sinh of a dual number, elemental
        ! module procedure sinh_dc ! obtain sinh of a dual number, elemental
        module procedure sinh_d2 ! obtain sinh of a dual number, elemental
        ! module procedure sinh_dc2 ! obtain sinh of a dual number, elemental
    end interface

    public cosh
    interface cosh
        module procedure cosh_d ! obtain cosh of a dual number, elemental
        ! module procedure cosh_dc ! obtain cosh of a dual number, elemental
            module procedure cosh_d2 ! obtain cosh of a dual number, elemental
                ! module procedure cosh_dc2 ! obtain cosh of a dual number, elemental
    end interface

    public tanh
    interface tanh
        module procedure tanh_d ! obtain tanh of a dual number, elemental
        ! module procedure tanh_dc ! obtain tanh of a dual number, elemental
            module procedure tanh_d2 ! obtain tanh of a dual number, elemental
                ! module procedure tanh_dc2 ! obtain tanh of a dual number, elemental
    end interface

    public asinh
    interface asinh
        module procedure asinh_d ! obtain asinh of a dual number, elemental
        ! module procedure asinh_dc ! obtain asinh of a dual number, elemental
        module procedure asinh_d2 ! obtain asinh of a dual number, elemental
            ! module procedure asinh_dc2 ! obtain asinh of a dual number, elemental
    end interface

    public acosh
    interface acosh
        module procedure acosh_d ! obtain acosh of a dual number, elemental
        ! module procedure acosh_dc ! obtain acosh of a dual number, elemental
        module procedure acosh_d2 ! obtain acosh of a dual number, elemental
        ! module procedure acosh_dc2 ! obtain acosh of a dual number, elemental
    end interface

    public atanh
    interface atanh
        module procedure atanh_d ! obtain atanh of a dual number, elemental
        ! module procedure atanh_dc ! obtain atanh of a dual number, elemental
        module procedure atanh_d2 ! obtain atanh of a dual number, elemental
        ! module procedure atanh_dc2 ! obtain atanh of a dual number, elemental
    end interface

    public erf
    interface erf
        module procedure erf_d      ! erf of a dual number 
        ! module procedure erf_c      ! erf of a complex number
        ! module procedure erf_dc     ! erf of a dual_complex number
        module procedure erf_d2      ! erf of a dual number 
        ! module procedure erf_dc2     ! erf of a dual_complex number
    end interface

    public dx_zero
    interface dx_zero
        module procedure zero_out_dual2 ! set the dx part of a dual2 number to zero
        module procedure zero_out_dual ! set the dx part of a dual number to zero
    end interface

    public dual_2_dual2
    public dual2_2_dual        
    public typ2_2_typ1

contains

elemental function typ2_2_typ1(typ2) result(typ1)
    !$acc routine seq
    type(dual2),intent(in) :: typ2
    type(dual2):: typ1
  
    typ1%x = typ2%x
    typ1%dx(1:dual_size) = typ2%dx(1+dual_size:2*dual_size)
    typ1%dx(dual_size+1:) = 0.0d0
   
end function typ2_2_typ1

!********* Set of functions to map dual to dual2
elemental function dual_2_dual2(d,typ) result(d2)
    type(dual), intent(in) :: d
    type(dual2):: d2
    integer,intent(in) :: typ

    select case(typ)
        case(1)
            d2%x = d%x
            d2%dx(1:dual_size) = d%dx
            d2%dx(dual_size+1:) = 0.0d0
        case(2)
            d2%x = d%x
            d2%dx(1:dual_size) = 0.0d0
            d2%dx(dual_size+1:) = d%dx
    end select
    
end function dual_2_dual2

!******** Set of functions to map dual2 to dual

function dual2_2_dual(d2,typ) result(d)
    type(dual2),intent(in) :: d2(:)
    type(dual), dimension(size(d2)) :: d
    integer,intent(in) :: typ

    select case(typ)
        case(1)
            d = dual2_2_dual_1(d2)
        case(2)
            d = dual2_2_dual_2(d2)
    end select
    
   
end function dual2_2_dual

elemental function dual2_2_dual_1(d2) result(d)
    type(dual) :: d
    type(dual2), intent(in) :: d2

    d%x = d2%x
    d%dx(1:dual_size) = d2%dx(1:dual_size)
   
end function dual2_2_dual_1

elemental function dual2_2_dual_2(d2) result(d)
    type(dual) :: d
    type(dual2), intent(in) :: d2

    d%x = d2%x
    d%dx(1:dual_size) = d2%dx(1+dual_size:dual_size2)

end function dual2_2_dual_2

!********* Function to zero out the second part of dual and dual2
elemental subroutine zero_out_dual2(d2)
    type(dual2), intent(inout) :: d2
    d2%dx = 0.0d0

end subroutine zero_out_dual2

elemental subroutine zero_out_dual(d)
    type(dual), intent(inout) :: d
    d%dx = 0.0d0

end subroutine zero_out_dual

!*********Begin: functions/subroutines for overloading operators

!******* Begin: (=)
!---------------------

    !-----------------------------------------
    ! dual = integer
    ! <u, du> = <j 0>
    !-----------------------------------------
    elemental subroutine assign_di(u, j)
         type(dual), intent(out) :: u
         integer, intent(in) :: j

         u%x =j  ! This is faster than direct assignment
         u%dx = 0.0d0

    end subroutine assign_di

    !-----------------------------------------
    ! dual = real(double)
    ! <u, du> = <r, 0>
    !-----------------------------------------
    elemental subroutine assign_dr(u, r)
        !$acc routine seq
        type(dual), intent(out) :: u
        real(wp), intent(in) :: r
       
        u%x = real(r,kind=wp)
        u%dx = 0.0d0

    end subroutine assign_dr


    !-----------------------------------------
    ! integer = dual
    ! j = <u, du>
    !-----------------------------------------
    elemental subroutine assign_id(j, v)
         type(dual), intent(in) :: v
         integer, intent(out) :: j

         j = int(v%x)

    end subroutine assign_id

    !-----------------------------------------
    ! real = dual
    ! j = <u, du>
    !-----------------------------------------
    elemental subroutine assign_rd(j, v)
         type(dual), intent(in) :: v
         real(wp), intent(out) :: j

         j = v%x

    end subroutine assign_rd

!******* end: (=)
!---------------------


!******* Begin: (+)
!---------------------

    !-----------------------------------------
    ! Unary positive
    ! <res, dres> = +<u, du>
    !-----------------------------------------
    elemental function add_d(u) result(res)
         type(dual), intent(in) :: u
         type(dual) :: res

         res = u  ! Faster than assigning component wise

    end function add_d


    !-----------------------------------------
    ! dual + dual
    ! <res, dres> = <u, du> + <v, dv> = <u + v, du + dv>
    !-----------------------------------------
    elemental function add_dd(u, v) result(res)
         type(dual), intent(in) :: u, v
         type(dual) :: res

         res%x = u%x + v%x
         res%dx = u%dx + v%dx

    end function add_dd


    !-----------------------------------------
    ! dual + integer
    ! <res, dres> = <u, du> + j = <u + j du>
    !-----------------------------------------
    elemental function add_di(u, j) result(res)
         type(dual), intent(in) :: u
         integer, intent(in) :: j
         type(dual) :: res

         res%x = real(j,kind=wp) + u%x
         res%dx = u%dx

    end function add_di


    !-----------------------------------------
    ! dual + double
    ! <res, dres> = <u, du> + <r, 0> = <u + r, du>
    !-----------------------------------------
    elemental function add_dr(u, r) result(res)
      !$acc routine seq
        type(dual), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual) :: res

        res%x = r + u%x
        res%dx = u%dx

    end function add_dr


    !-----------------------------------------
    ! integer + dual
    ! <res, dres> = <j 0> + <v, dv> = <i + v, dv>
    !-----------------------------------------
    elemental function add_id(j, v) result(res)
        integer, intent(in) :: j
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = real(j,kind=wp) + v%x
        res%dx = v%dx

    end function add_id


    !-----------------------------------------
    ! double + dual
    ! <res, dres> = <r, 0> + <v, dv> = <r + v, dv>
    !-----------------------------------------
    elemental function add_rd(r, v) result(res)
        real(wp), intent(in) :: r
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = r + v%x
        res%dx = v%dx

    end function add_rd

!******* end: (+)
!---------------------


!******* Begin: (-)
!---------------------

    !-------------------------------------------------
    ! negate a dual
    ! <res, dres> = -<u, du>
    !-------------------------------------------------
    elemental function minus_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = -u%x
        res%dx = -u%dx

    end function minus_d


    !-------------------------------------------------
    ! dual - dual
    ! <res, dres> = <u, du> - <v, dv> = <u - v, du - dv>
    !-------------------------------------------------
    elemental function minus_dd(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        res%x = u%x - v%x
        res%dx = u%dx - v%dx

    end function minus_dd

    !-------------------------------------------------
    ! dual - integer
    ! <res, dres> = <u, du> - j = <u - j du>
    !-------------------------------------------------
    elemental function minus_di(u, j) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: j
        type(dual) :: res

        res%x = u%x - real(j,kind=wp)
        res%dx = u%dx

    end function minus_di


    !-------------------------------------------------
    ! dual - double
    ! <res, dres> = <u, du> - r = <u - r, du>
    !-------------------------------------------------
    elemental function minus_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real(wp),intent(in) :: r
        type(dual) :: res

        res%x = u%x - r
        res%dx = u%dx

    end function minus_dr


    !-------------------------------------------------
    ! integer - dual
    ! <res, dres> = j - <v, dv> = <i - v, -dv>
    !-------------------------------------------------
    elemental function minus_id(j,v) result(res)
        integer, intent(in) :: j
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = real(j,kind=wp) - v%x
        res%dx = -v%dx

    end function minus_id


    !-------------------------------------------------
    ! double - dual
    ! <res, dres> = r - <v, dv> = <r - v, -dv>
    !-------------------------------------------------
    elemental function minus_rd(r, v) result(res)
         real(wp), intent(in) :: r
         type(dual), intent(in) :: v
         type(dual) :: res

        res%x = r - v%x
        res%dx = -v%dx

    end function minus_rd

!******* end: (-)
!---------------------


!******* BEGIN: (*)
!---------------------

    !----------------------------------------
    ! dual * dual
    ! <res, dres> = <u, du> * <v, dv> = <u * v, u * dv + v * du>
    !----------------------------------------
    elemental function mult_dd(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        res%x = u%x * v%x
        res%dx = u%x * v%dx + v%x * u%dx

    end function mult_dd


    !-----------------------------------------
    ! dual * integer
    ! <res, dres> = <u, du> * j = <u * j du * i>
    !-----------------------------------------
    elemental function mult_di(u, j) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: j
        type(dual) :: res

        real(wp):: r

        r = real(j,kind=wp)
        res%x = r * u%x
        res%dx = r * u%dx

    end function mult_di

    !-----------------------------------------
    ! dual * integer
    ! <res, dres> = <u, du> * j = <u * j du * i>
    !-----------------------------------------
    elemental function mult_di16(u, j) result(res)
        type(dual), intent(in) :: u
        integer(int16), intent(in) :: j
        type(dual) :: res

        real(wp):: r

        r = real(j,kind=wp)
        res%x = r * u%x
        res%dx = r * u%dx

    end function mult_di16

            !-----------------------------------------
    ! dual * integer
    ! <res, dres> = <u, du> * j = <u * j du * i>
    !-----------------------------------------
    elemental function mult_di8(u, j) result(res)
        type(dual), intent(in) :: u
        integer(int8), intent(in) :: j
        type(dual) :: res

        real(wp):: r

        r = real(j,kind=wp)
        res%x = r * u%x
        res%dx = r * u%dx

    end function mult_di8

    !-----------------------------------------
    ! dual * double
    ! <res, dres> = <u, du> * r = <u * r, du * r>
    !----------------------------------------
    elemental function mult_dr(u, r) result(res)
        !$acc routine seq
        type(dual), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual) :: res

        res%x = u%x * r
        res%dx = u%dx * r

    end function mult_dr


    !-----------------------------------------
    ! integer * dual
    ! <res, dres> = j * <v, dv> = <i * v, j * dv>
    !-----------------------------------------
    elemental function mult_id(j,v) result(res)
        integer, intent(in) :: j
        type(dual), intent(in) :: v
        type(dual) :: res

        real(wp):: r

        r = real(j,kind=wp)
        res%x = r * v%x
        res%dx = r * v%dx

    end function mult_id


    !-----------------------------------------
    ! double * dual
    ! <res, dres> = r * <v, dv> = <r * v, r * dv>
    !-----------------------------------------
    elemental function mult_rd(r, v) result(res)
        real(wp), intent(in) :: r
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = r * v%x
        res%dx = r * v%dx

    end function mult_rd

!******* end: (*)
!---------------------


!******* BEGIN: (/)
!---------------------

    !-----------------------------------------
    ! dual / dual
    ! <res, dres> = <u, du> / <v, dv> = <u / v, du / v - u * dv / v^2>
    !-----------------------------------------
    elemental function div_dd(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        real(wp):: inv

        inv = 1.0d0 / v%x
        res%x = u%x * inv
        res%dx = (u%dx - res%x * v%dx) * inv

    end function div_dd


    !-----------------------------------------
    ! dual / integer
    ! <res, dres> = <u, du> / j = <u / j du / i>
    !-----------------------------------------
    elemental function div_di(u, j) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: j
        type(dual) :: res

        real(wp):: inv

        inv = 1.0d0 / real(j,kind=wp)
        res%x = u%x * inv
        res%dx = u%dx * inv

    end function div_di


    !-----------------------------------------
    ! dual / double
    ! <res, dres> = <u, du> / r = <u / r, du / r>
    !----------------------------------------
    elemental function div_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual):: res

        real(wp):: inv

        inv = 1.0d0 / r
        res%x = u%x * inv
        res%dx = u%dx * inv

    end function div_dr


    !-----------------------------------------
    ! integer / dual
    ! <res, dres> = j / <v, dv> = <i / v, -i / v^2 * du>
    !-----------------------------------------
    elemental function div_id(j,v) result(res)
        integer, intent(in) :: j
        type(dual), intent(in) :: v
        type(dual) :: res

        real(wp):: inv

        inv = 1.0d0 / v%x
        res%x = real(j,kind=wp) * inv
        res%dx = -res%x * inv * v%dx

    end function div_id


    !-----------------------------------------
    ! double / dual
    ! <res, dres> = r / <u, du> = <r / u, -r / u^2 * du>
    !-----------------------------------------
    elemental function div_rd(r, v) result(res)
        real(wp), intent(in) :: r
        type(dual), intent(in) :: v
        type(dual) :: res

        real(wp):: inv

        inv = 1.0d0 / v%x
        res%x = r * inv
        res%dx = -res%x * inv * v%dx

    end function div_rd

!******* end: (/)
!---------------------

!******* BEGIN: (**)
!---------------------

    !-----------------------------------------
    ! power(dual, integer)
    ! <res, dres> = <u, du> ^ j = <u ^ j j * u ^ (j- 1) * du>
    !-----------------------------------------
    elemental function pow_di(u, j) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: j
        type(dual) :: res

        real(wp):: pow_x

        pow_x = u%x ** (j- 1)
        res%x = u%x * pow_x
        res%dx = real(j,kind=wp) * pow_x * u%dx

    end function pow_di

    !-----------------------------------------
    ! power(dual, double)
    ! <res, dres> = <u, du> ^ r = <u ^ r, r * u ^ (r - 1) * du>
    !-----------------------------------------
    elemental function pow_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual) :: res

        real(wp):: pow_x

        pow_x = u%x ** (r - 1.0d0)
        res%x = u%x * pow_x
        res%dx = r * pow_x * u%dx

    end function pow_dr

    !-----------------------------------------
    ! POWER dual number to a dual power
    ! <res, dres> = <u, du> ^ <v, dv>
    !     = <u ^ v, u ^ v * (v / u * du + Log(u) * dv)>
    !-----------------------------------------
    elemental function pow_dd(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        res%x = u%x ** v%x
        res%dx = res%x * (v%x / u%x * u%dx + log(u%x) * v%dx)

    end function pow_dd

    !-----------------------------------------
    ! POWER integer to a dual power
    ! <res, dres> = j ^ <u, du>
    !     = <i ^ u, j ^ u * Log(j) * du)>
    !-----------------------------------------
    elemental function pow_id(j,u) result(res)
        integer, intent(in) :: j
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = j ** u%x
        res%dx = res%x * log(real(j,kind=wp)) * u%dx

    end function pow_id

    !-----------------------------------------
    ! POWER real to a dual power
    ! <res, dres> = r ^ <u, du>
    !     = <r ^ u, r ^ u * Log(r) * du)>
    !-----------------------------------------
    elemental function pow_rd(r, u) result(res)
        real(wp), intent(in) :: r
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = r ** u%x
        res%dx = res%x * log(r) * u%dx

    end function pow_rd    

!******* end: (**)
!---------------------


!******* BEGIN: (==)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_dd(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs, rhs
         logical :: res

         res = (lhs%x == rhs%x)

    end function eq_dd


    !-----------------------------------------
    ! compare a dual with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_di(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs
         integer, intent(in) :: rhs
         logical :: res

         res = (lhs%x == real(rhs))

    end function eq_di


    !-----------------------------------------
    ! compare a dual number with a real number,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real(wp), intent(in) :: rhs
        logical::res

        res = (lhs%x == rhs)

    end function eq_dr


    !-----------------------------------------
    ! compare an integer with a dual,
    ! simply compare the functional value.
    !----------------------------------------
    elemental function eq_id(lhs, rhs) result(res)
         integer, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs == rhs%x)

    end function eq_id


    !-----------------------------------------
    ! compare a real with a dual,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_rd(lhs, rhs) result(res)
         real(wp), intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs == rhs%x)

    end function eq_rd

!******* end: (==)
!---------------------


!******* BEGIN: (<=)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function le_dd(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs, rhs
         logical :: res

         res = (lhs%x <= rhs%x)

    end function le_dd


    !-----------------------------------------
    ! compare a dual with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_di(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs
         integer, intent(in) :: rhs
         logical :: res

         res = (lhs%x <= rhs)

    end function le_di


    !-----------------------------------------
    ! compare a dual number with a real number,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_dr(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs
         real(wp), intent(in) :: rhs
         logical :: res

         res = (lhs%x <= rhs)

    end function le_dr


    !-----------------------------------------
    ! compare a dual number with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_id(j,rhs) result(res)
         integer, intent(in) :: j
         type(dual), intent(in) :: rhs
         logical :: res

         res = (j<= rhs%x)

    end function le_id


    !-----------------------------------------
    ! compare a real with a dual,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_rd(lhs, rhs) result(res)
         real(wp), intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs <= rhs%x)

    end function le_rd

!******* end: (<=)
!---------------------

!******* BEGIN: (<)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function lt_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x < rhs%x)

    end function lt_dd

    !-----------------------------------------
    ! compare a dual with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function lt_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x < rhs)

    end function lt_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !----------------------------------------
    elemental function lt_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real(wp), intent(in) :: rhs
        logical :: res

        res = (lhs%x < rhs)

    end function lt_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function lt_id(j,rhs) result(res)
         integer, intent(in) :: j
         type(dual), intent(in) :: rhs
         logical :: res

         res = (j< rhs%x)

    end function lt_id


    !-----------------------------------------
    ! compare a real with a dual
    !----------------------------------------
    elemental function lt_rd(lhs, rhs) result(res)
         real(wp), intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs < rhs%x)

    end function lt_rd

!******* end: (<)
!---------------------

!******* BEGIN: (>=)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !----------------------------------------
    elemental function ge_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x >= rhs%x)

    end function ge_dd


    !-----------------------------------------
    ! compare a dual with an integer
    !-----------------------------------------
    elemental function ge_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x >= rhs)

    end function ge_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ge_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real(wp), intent(in) :: rhs
        logical :: res

        res = (lhs%x >= rhs)

    end function ge_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function ge_id(j,rhs) result(res)
        integer, intent(in) :: j
        type(dual), intent(in) :: rhs
        logical :: res

        res = (j>= rhs%x)

    end function ge_id


    !-----------------------------------------
    ! compare a real with a dual
    !-----------------------------------------
    elemental function ge_rd(lhs, rhs) result(res)
         real(wp), intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs >= rhs%x)

    end function ge_rd

!******* end: (>=)
!---------------------

!******* BEGIN: (>)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function gt_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x > rhs%x)

    end function gt_dd


    !-----------------------------------------
    ! compare a dual with an integer
    !-----------------------------------------
    elemental function gt_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x > rhs)

    end function gt_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function gt_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real(wp), intent(in) :: rhs
        logical :: res

        res = (lhs%x > rhs)

    end function gt_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function gt_id(j,rhs) result(res)
        integer, intent(in) :: j
        type(dual), intent(in) :: rhs
        logical :: res

        res = (j> rhs%x)

    end function gt_id


    !-----------------------------------------
    ! compare a real with a dual
    !-----------------------------------------
    elemental function gt_rd(lhs, rhs) result(res)
         real(wp), intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs > rhs%x)

    end function gt_rd

!******* end: (>)
!---------------------

!******* BEGIN: (/=)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ne_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x /= rhs%x)

    end function ne_dd


    !-----------------------------------------
    ! compare a dual with an integer
    !-----------------------------------------
    elemental function ne_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x /= rhs)

    end function ne_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ne_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real(wp), intent(in) :: rhs
        logical :: res

        res = (lhs%x /= rhs)

    end function ne_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function ne_id(j,rhs) result(res)
        integer, intent(in) :: j
        type(dual), intent(in) :: rhs
        logical :: res

        res = (j/= rhs%x)

    end function ne_id


    !-----------------------------------------
    ! compare a real with a dual
    !-----------------------------------------
    elemental function ne_rd(lhs, rhs) result(res)
        real(wp), intent(in) :: lhs
        type(dual), intent(in) :: rhs
        logical :: res

        res = (lhs /= rhs%x)

    end function ne_rd

!******* end: (/=)
!---------------------

    !---------------------------------------------------
    ! Absolute value of dual numbers
    ! <res, dres> = abs(<u, du>) = <abs(u), du * sign(u)>
    !---------------------------------------------------
    elemental function abs_d(u) result(res)
         type(dual), intent(in) :: u
         type(dual) :: res
         integer :: j

         if(u%x > 0) then
            res%x = u%x
            res%dx = u%dx
         else if (u%x < 0) then
            res%x = -u%x
            res%dx = -u%dx
         else
            res%x = 0.0d0
            do j = 1, dual_size
                if (u%dx(j) .eq. 0.0d0) then
                    res%dx(j) = 0.0d0
                else
                    res%dx(j) = set_NaN()
                end if
            end do
         endif

    end function abs_d


    !-----------------------------------------
    ! ACOS of dual numbers
    ! <res, dres> = acos(<u, du>) = <acos(u), -du / sqrt(1 - u^2)>
    !----------------------------------------
    elemental function acos_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = acos(u%x)
        if (u%x == 1.0d0 .or. u%x == -1.0d0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = -u%dx / sqrt(1.0d0 - u%x**2)
        end if

    end function acos_d


    !-----------------------------------------
    ! ASIN of dual numbers
    ! <res, dres> = asin(<u, du>) = <asin(u), du / sqrt(1 - u^2)>
    !----------------------------------------
    elemental function asin_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = asin(u%x)
        if (u%x == 1.0d0 .or. u%x == -1.0d0) then
            res%dx = set_NaN()  ! Undefined derivative
        else
            res%dx = u%dx / sqrt(1.0d0 - u%x**2)
        end if

    end function asin_d


    !-----------------------------------------
    ! ATAN of dual numbers
    ! <res, dres> = atan(<u, du>) = <atan(u), du / (1 + u^2)>
    !----------------------------------------
    elemental function atan_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = atan(u%x)
        res%dx = u%dx / (1.0d0 + u%x**2)

    end function atan_d


    !-----------------------------------------
    ! ATAN2 of dual numbers
    ! <res, dres> = atan2(<u, du>, <v, dv>)
    !             = <atan2(u, v), v / (u^2 + v^2) * du - u / (u^2 + v^2) * dv>
    !----------------------------------------
    elemental function atan2_d(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        real(wp):: usq_plus_vsq

        res%x = atan2(u%x, v%x)

        usq_plus_vsq = u%x**2 + v%x**2
        res%dx = v%x / usq_plus_vsq * u%dx - u%x / usq_plus_vsq * v%dx

    end function atan2_d


    !-----------------------------------------
    ! COS of dual numbers
    ! <res, dres> = cos(<u, du>) = <cos(u), -sin(u) * du>
    !----------------------------------------
    elemental function cos_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = cos(u%x)
        res%dx = -sin(u%x) * u%dx

    end function cos_d


    !-----------------------------------------
    ! DOT PRODUCT two dual number vectors
    ! <res, dres> = <u, du> . <v, dv> = <u . v, u . dv + v . du>
    !-----------------------------------------
    function dot_product_dd(u, v) result(res)
        type(dual), intent(in) :: u(:), v(:)
        type(dual) :: res

        integer :: j

        res%x = dot_product(u%x, v%x)
        do j = 1, dual_size
            res%dx(j) = dot_product(u%x, v%dx(j)) + dot_product(v%x, u%dx(j))
        end do

    end function dot_product_dd


    !-----------------------------------------
    ! EXPONENTIAL OF dual numbers
    ! <res, dres> = exp(<u, du>) = <exp(u), exp(u) * du>
    !-----------------------------------------
    elemental function exp_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        real(wp):: exp_x

        exp_x = exp(u%x)
        res%x = exp_x
        res%dx = u%dx * exp_x

    end function exp_d


    !-----------------------------------------
    ! Convert dual to integer
    ! j = int(<u, du>) = int(u)
    !----------------------------------------
    elemental function int_d(u) result(res)
         type(dual), intent(in) :: u
         integer :: res

         res = int(u%x)

    end function int_d


    !-----------------------------------------
    ! LOG OF dual numbers,defined for u%x>0 only
    ! the error control should be done in the original code
    ! in other words, if u%x<=0, it is not possible to obtain LOG.
    ! <res, dres> = log(<u, du>) = <log(u), du / u>
    !----------------------------------------
    elemental function log_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        real(wp):: inv

        inv = 1.0d0 / u%x
        res%x = log(u%x)
        res%dx = u%dx * inv

    end function log_d


    !-----------------------------------------
    ! LOG10 OF dual numbers,defined for u%x>0 only
    ! the error control should be done in the original code
    ! in other words, if u%x<=0, it is not possible to obtain LOG.
    ! <res, dres> = log10(<u, du>) = <log10(u), du / (u * log(10))>
    ! LOG<u,up>=<LOG(u),up/u>
    !----------------------------------------
    elemental function log10_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        real(wp):: inv

        inv = 1.0d0 / (u%x * log(10.0d0))
        res%x = log10(u%x)
        res%dx = u%dx * inv

    end function log10_d


    !-----------------------------------------
    ! MULTIPLY two dual number matrices
    ! <res, dres> = <u, du> . <v, dv> = <u . v, du . v + u . dv>
    !----------------------------------------
    function matmul_dd(u,v) result(res)
        type(dual), intent(in) :: u(:,:), v(:,:)
        type(dual) :: res(size(u,1), size(v,2))

        integer :: j

        res%x = matmul(u%x, v%x)
        do j = 1, dual_size
            res%dx(j) = matmul(u%dx(j), v%x) + matmul(u%x, v%dx(j))
        end do

    end function matmul_dd


    !-----------------------------------------
    ! MULTIPLY a dual number matrix with a dual number
    ! vector
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    function matmul_dv(u, v) result(res)
        type(dual), intent(in) :: u(:,:), v(:)
        type(dual) :: res(size(u,1))
        integer :: j

        res%x = matmul(u%x, v%x)
        do j = 1, dual_size
            res%dx(j) = matmul(u%dx(j), v%x) + matmul(u%x, v%dx(j))
        end do

    end function matmul_dv


    !-----------------------------------------
    ! MULTIPLY a dual vector with a  dual matrix
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    function matmul_vd(u, v) result(res)
        type(dual), intent(in) :: u(:), v(:,:)
        type(dual) :: res(size(v, 2))
        integer::j

        res%x = matmul(u%x, v%x)
        do j = 1, dual_size
            res%dx(j) = matmul(u%dx(j), v%x) + matmul(u%x, v%dx(j))
        end do

    end function matmul_vd

    !-----------------------------------------
    ! Obtain the max of 2 to 5 dual numbers
    !----------------------------------------
    elemental function max_dd(val1, val2, val3, val4,val5) result(res)
        type(dual), intent(in) :: val1, val2
        type(dual), intent(in), optional :: val3, val4,val5
        type(dual) :: res

        if (val1%x > val2%x) then
            res = val1
        else
            res = val2
        endif
        if(present(val3))then
           if(res%x < val3%x) res = val3
        endif
        if(present(val4))then
           if(res%x < val4%x) res = val4
        endif
        if(present(val5))then
           if(res%x < val5%x) res = val5
        endif

    end function max_dd


    !-----------------------------------------
    ! Obtain the max of a dual number and an integer
    !----------------------------------------
    elemental function max_di(u, j) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: j
        type(dual) :: res

        if (u%x > j) then
            res = u
        else
            res = j
        endif

    end function max_di

    !-----------------------------------------
    ! Obtain the max of a dual number and a real number
    !----------------------------------------
    elemental function max_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual) :: res

        if (u%x > r) then
            res = u
        else
            res%x = r
        endif

    end function max_dr


    !---------------------------------------------------
    ! Obtain the max of a real and a dual
    !---------------------------------------------------
     elemental function max_rd(n, u) result(res)
        real(wp), intent(in) :: n
        type(dual), intent(in) :: u
        type(dual) :: res

        if (u%x > n) then
            res = u
        else
            res%x = n
        endif

    end function max_rd


    !-----------------------------------------
    ! Obtain the max value of vector u
    !----------------------------------------
    function maxval_d(u) result(res)
        type(dual), intent(in) :: u(:)
        integer :: iloc(1)
        type(dual) :: res

        iloc=maxloc(u%x)
        res=u(iloc(1))

    end function maxval_d


    !-----------------------------------------
    ! Obtain the min of 2 to 4 dual numbers
    !----------------------------------------
    elemental function min_dd(val1, val2, val3, val4) result(res)
        type(dual), intent(in) :: val1, val2
        type(dual), intent(in), optional :: val3, val4
        type(dual) :: res

        if (val1%x < val2%x) then
            res = val1
        else
            res = val2
        endif
        if(present(val3))then
           if(res%x > val3%x) res = val3
        endif
        if(present(val4))then
           if(res%x > val4%x) res = val4
        endif

    end function min_dd


    !-----------------------------------------
    ! Obtain the min of a dual and a double
    !----------------------------------------
    elemental function min_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual) :: res

        if (u%x < r) then
            res = u
        else
            res%x = r
        endif

    end function min_dr


  !-----------------------------------------
    ! Obtain the min value of vector u
    !----------------------------------------
    function minval_d(u) result(res)
        type(dual), intent(in) :: u(:)
        integer :: iloc(1)
        type(dual) :: res

        iloc=minloc(u%x)
        res=u(iloc(1))

    end function minval_d


    !------------------------------------------------------
    !Returns the nearest integer to u%x, ELEMENTAL
    !------------------------------------------------------
    elemental function nint_d(u) result(res)
        type(dual), intent(in) :: u
        integer :: res

        res=nint(u%x)

    end function nint_d


    !----------------------------------------------------------------
    ! SIGN(a,b) with two dual numbers as inputs,
    ! the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
    !----------------------------------------------------------------
    elemental function sign_dd(val1, val2) result(res)
        type(dual), intent(in) :: val1, val2
        type(dual) :: res

        if (val2%x < 0.0d0) then
            res = -abs(val1)
        else
            res =  abs(val1)
        endif

     end function sign_dd


    !----------------------------------------------------------------
    ! SIGN(a,b) with one real and one dual number as inputs,
    ! the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
    !----------------------------------------------------------------
    elemental function sign_rd(val1, val2) result(res)
        real(wp), intent(in) :: val1
        type(dual), intent(in) :: val2
        type(dual) :: res

        if (val2%x < 0.0d0) then
            res%x = -abs(val1)
        else
            res%x = abs(val1)
        endif

     end function sign_rd


    !-----------------------------------------
    ! SIN of dual numbers
    ! <res, dres> = sin(<u, du>) = <sin(u), cos(u) * du>
    !----------------------------------------
    elemental function sin_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = sin(u%x)
        res%dx = cos(u%x) * u%dx

    end function sin_d


    !-----------------------------------------
    ! TAN of dual numbers
    ! <res, dres> = tan(<u, du>) = <tan(u), du / cos(u)^2>
    !----------------------------------------
    elemental function tan_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = tan(u%x)
        res%dx = u%dx / cos(u%x)**2

    end function tan_d


    !-----------------------------------------
    ! SQRT of dual numbers
    ! <res, dres> = sqrt(<u, du>) = <sqrt(u), du / (2 * sqrt(u))>
    !----------------------------------------
    elemental function sqrt_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res
        integer :: j

        res%x = sqrt(u%x)

        if (res%x .ne. 0.0d0) then
            res%dx = 0.5 * u%dx / res%x
        else
            do j = 1, dual_size
                if (u%dx(j) .eq. 0.0d0) then
                    res%dx(j) = 0.0d0
                else
                    res%dx(j) = set_NaN()
                end if
            end do
        end if

    end function sqrt_d


    !-----------------------------------------
    ! Sum of a dual array
    !-----------------------------------------
    function sum_d(u) result(res)
        type(dual), intent(in) :: u(:)
        type(dual) :: res
        integer :: j

        res%x = sum(u%x)
        do j = 1, dual_size
            res%dx(j) = sum(u%dx(j))
        end do

    end function sum_d


    !-----------------------------------------
    ! Find the location of the max value in an
    ! array of dual numbers
    !-----------------------------------------
    function maxloc_d(array) result(ind)
        type(dual), intent(in) :: array(:)
        integer :: ind(1)

        ind = maxloc(array%x)

    end function maxloc_d


    elemental function set_NaN() result(res)
        real(wp):: res

        res = sqrt(negative_one)

    end function set_NaN


    !-----------------------------------------
    ! Hyperbolic functions: sinh, cosh, tanh
    ! and their inverses: asinh, acosh, atanh
    !-----------------------------------------
    !-----------------------------------------
    ! SINH OF dual numbers
    ! <res, dres> = sinh(<u, du>) = <sinh(u), cosh(u) * du>
    !-----------------------------------------
    elemental function sinh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = sinh(u%x)
        res%dx = u%dx * cosh(u%x)

    end function sinh_d

    !-----------------------------------------
    ! COSH OF dual numbers
    ! <res, dres> = cosh(<u, du>) = <cosh(u), sinh(u) * du>
    !-----------------------------------------
    elemental function cosh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = cosh(u%x)
        res%dx = u%dx * sinh(u%x)

    end function cosh_d

    !-----------------------------------------
    ! TANH OF dual numbers
    ! <res, dres> = tanh(<u, du>) = <tanh(u), 1.0d0/cosh(u)**2 * du>
    !-----------------------------------------
    elemental function tanh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = tanh(u%x)
        res%dx = u%dx * 1.0d0/cosh(u%x)**2

    end function tanh_d

    !-----------------------------------------
    ! ASINH OF dual numbers
    ! <res, dres> = asinh(<u, du>) = <asinh(u), 1/sqrt(u**2 + 1) * du>
    !-----------------------------------------
    elemental function asinh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = asinh(u%x)
        res%dx = u%dx * 1.0d0/sqrt(u%x**2 + 1.0d0)

    end function asinh_d

    !-----------------------------------------
    ! ACOSH OF dual numbers
    ! <res, dres> = acosh(<u, du>) = <acosh(u), 1/sqrt(u**2 - 1) * du>
    !-----------------------------------------
    elemental function acosh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = acosh(u%x)
        if (u%x <= 1.0d0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = u%dx * 1.0d0/sqrt(u%x**2 - 1.0d0)
        end if

    end function acosh_d

    !-----------------------------------------
    ! ATAHN OF dual numbers
    ! <res, dres> = atanh(<u, du>) = <atanh(u), 1/(1 - u**2) * du>
    !-----------------------------------------
    elemental function atanh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = atanh(u%x)
        if (abs(u%x) >= 1.0d0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = u%dx * 1.0d0/(1.0d0 - u%x**2)
        end if

    end function atanh_d

    !-----------------------------------------
    ! ERF OF dual numbers
    ! <res, dres> = ERF(<u, du>) = <ERF(u), 2/sqrt(Pj)*exp(-u**2) * du>
    !-----------------------------------------
    elemental function erf_d(u) result(res)
      type(dual), intent(in) :: u
      type(dual) :: res

      res%x  = erf(u%x)
      res%dx = u%dx * 2.0/sqrt(PI) * exp(-u%x**2)

    end function erf_d







! ! ******************************************************************************
! ! COMPLEX DUAL NUMBERS
! ! ******************************************************************************

! !*********Begin: functions/subroutines for overloading operators

! !******* Begin: (=)
! !---------------------

!     !-----------------------------------------
!     ! dual_complex = integer
!     ! <u, du> = <j 0>
!     !-----------------------------------------
!     elemental subroutine assign_dc_i(u, j)
!          type(dual_complex), intent(out) :: u
!          integer, intent(in) :: j

!          u%z = complex(real(j,kind=wp), 0.0d0)  ! This is faster than direct assignment
!          u%dz = complex(0.0d0, 0.0d0)

!     end subroutine assign_dc_i


!     !-----------------------------------------
!     ! dual_complex = real(double)
!     ! <u, du> = <r, 0>
!     !-----------------------------------------
!     elemental subroutine assign_dc_r(u, r)
!         type(dual_complex), intent(out) :: u
!         real(wp), intent(in) :: r

!         u%z = complex(r, 0.0d0)
!         u%dz = complex(0.0d0, 0.0d0)

!     end subroutine assign_dc_r

!     !-----------------------------------------
!     ! dual_complex = complex
!     ! <u, du> = <r, 0>
!     !-----------------------------------------
!     elemental subroutine assign_dc_c(u, r)
!         type(dual_complex), intent(out) :: u
!         complex(wp), intent(in) :: r

!         u%z = r
!         u%dz = complex(0.0d0, 0.0d0)

!     end subroutine assign_dc_c


!     !-----------------------------------------
!     ! integer = dual        Is there a situation where complex --> integer makes sense?
!     ! j = <u, du>
!     !-----------------------------------------
!     ! elemental subroutine assign_id(j,v)
!     !      type(dual_complex), intent(in) :: v
!     !      integer, intent(out) :: j
!     !
!     !      j = int(v%z)
!     !
!     ! end subroutine assign_id

!     !-----------------------------------------
!     ! complex = dual_complex        Is there a situation where complex --> integer makes sense?
!     ! j = <u, du>
!     !-----------------------------------------
!     elemental subroutine assign_c_dc(u, v)
!          type(dual_complex), intent(in) :: v
!          complex(wp), intent(out) :: u

!          u = v%z

!     end subroutine assign_c_dc


! !******* end: (=)
! !---------------------

! !******* Begin: (+)
! !---------------------

!     !-----------------------------------------
!     ! Unary positive
!     ! <res, dres> = +<u, du>
!     !-----------------------------------------
!     elemental function add_dc(u) result(res)
!          type(dual_complex), intent(in) :: u
!          type(dual_complex) :: res

!          res = u  ! Faster than assigning component wise

!     end function add_dc

!     !-----------------------------------------
!     ! dual_complex + dual_complex
!     ! <res, dres> = <u, du> + <v, dv> = <u + v, du + dv>
!     !-----------------------------------------
!     elemental function add_dc_dc(u, v) result(res)
!        type(dual_complex), intent(in) :: u, v
!        type(dual_complex) :: res

!        res%z = u%z + v%z
!        res%dz = u%dz + v%dz

!     end function add_dc_dc

!     !-----------------------------------------
!     ! complex + dual_complex
!     ! <res, dres> = <r, 0> + <v, dv> = <r + v, dv>
!     !-----------------------------------------
!     elemental function add_c_dc(r, v) result(res)
!         complex(wp), intent(in) :: r
!         type(dual_complex), intent(in) :: v
!         type(dual_complex) :: res

!         res%z = r + v%z
!         res%dz = v%dz

!     end function add_c_dc

!     !-----------------------------------------
!     ! dual_complex + complex
!     ! <res, dres> = <u, du> + <r, 0> = <u + r, du>
!     !-----------------------------------------
!     elemental function add_dc_c(u, r) result(res)
!       type(dual_complex), intent(in) :: u
!       complex(wp), intent(in) :: r
!       type(dual_complex) :: res

!       res%z = r + u%z
!       res%dz = u%dz

!     end function add_dc_c


!   !-----------------------------------------
!   ! dual_complex + integer
!   ! <res, dres> = <u, du> + j = <u + j du>
!   !-----------------------------------------
!   elemental function add_dc_i(u, j) result(res)
!        type(dual_complex), intent(in) :: u
!        integer, intent(in) :: j
!        type(dual_complex) :: res

!        res%z = real(j,kind=wp) + u%z
!        res%dz = u%dz

!   end function add_dc_i

!   !-----------------------------------------
!   ! integer + dual_complex
!   ! <res, dres> = <u, du> + j = <u + j du>
!   !-----------------------------------------
!   elemental function add_i_dc(j,u) result(res)
!        type(dual_complex), intent(in) :: u
!        integer, intent(in) :: j
!        type(dual_complex) :: res

!        res%z = real(j,kind=wp) + u%z
!        res%dz = u%dz

!   end function add_i_dc

!   !-----------------------------------------
!   ! dual_complex + integer
!   ! <res, dres> = <u, du> + j = <u + j du>
!   !-----------------------------------------
!   elemental function add_dc_r(u, r) result(res)
!        type(dual_complex), intent(in) :: u
!        real(wp), intent(in) :: r
!        type(dual_complex) :: res

!        res%z = r + u%z
!        res%dz = u%dz

!   end function add_dc_r

!   !-----------------------------------------
!   ! integer + dual_complex
!   ! <res, dres> = <u, du> + j = <u + j du>
!   !-----------------------------------------
!   elemental function add_r_dc(r, u) result(res)
!        type(dual_complex), intent(in) :: u
!        real(wp), intent(in) :: r
!        type(dual_complex) :: res

!        res%z = r + u%z
!        res%dz = u%dz

!   end function add_r_dc
! !
! ! !******* end: (+)
! ! !---------------------
! !
! !
! !******* Begin: (-)
! !---------------------

!   !-------------------------------------------------
!   ! negate a dual_complex
!   ! <res, dres> = -<u, du>
!   !-------------------------------------------------
!   elemental function minus_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = -u%z
!       res%dz = -u%dz

!   end function minus_dc


!   !-------------------------------------------------
!   ! dual_complex - dual_complex
!   ! <res, dres> = <u, du> - <v, dv> = <u - v, du - dv>
!   !-------------------------------------------------
!   elemental function minus_dc_dc(u, v) result(res)
!       type(dual_complex), intent(in) :: u, v
!       type(dual_complex) :: res

!       res%z = u%z - v%z
!       res%dz = u%dz - v%dz

!   end function minus_dc_dc

!   !-------------------------------------------------
!   ! dual_complex - integer
!   ! <res, dres> = <u, du> - j = <u - j du>
!   !-------------------------------------------------
!   elemental function minus_dc_i(u, j) result(res)
!       type(dual_complex), intent(in) :: u
!       integer, intent(in) :: j
!       type(dual_complex) :: res

!       res%z = u%z - real(j,kind=wp)
!       res%dz = u%dz

!   end function minus_dc_i

!   !-------------------------------------------------
!   ! dual_complex - real
!   ! <res, dres> = <u, du> - r = <u - r, du>
!   !-------------------------------------------------
!   elemental function minus_dc_r(u, r) result(res)
!       type(dual_complex), intent(in) :: u
!       real(wp), intent(in) :: r
!       type(dual_complex) :: res

!       res%z = u%z - r
!       res%dz = u%dz

!   end function minus_dc_r


!   !-------------------------------------------------
!   ! dual_complex - complex
!   ! <res, dres> = <u, du> - r = <u - r, du>
!   !-------------------------------------------------
!   elemental function minus_dc_c(u, c) result(res)
!       type(dual_complex), intent(in) :: u
!       complex(wp), intent(in) :: c
!       type(dual_complex) :: res

!       res%z = u%z - c
!       res%dz = u%dz

!   end function minus_dc_c


!   !-------------------------------------------------
!   ! integer - dual_complex
!   ! <res, dres> = j - <v, dv> = <i - v, -dv>
!   !-------------------------------------------------
!   elemental function minus_i_dc(j,v) result(res)
!       integer, intent(in) :: j
!       type(dual_complex), intent(in) :: v
!       type(dual_complex) :: res

!       res%z = real(j,kind=wp) - v%z
!       res%dz = -v%dz

!   end function minus_i_dc

!   !-------------------------------------------------
!   ! real - dual_complex
!   ! <res, dres> = r - <v, dv> = <r - v, -dv>
!   !-------------------------------------------------
!   elemental function minus_r_dc(r, v) result(res)
!       real(wp), intent(in) :: r
!       type(dual_complex), intent(in) :: v
!       type(dual_complex) :: res

!       res%z = r - v%z
!       res%dz = -v%dz

!   end function minus_r_dc


!   !-------------------------------------------------
!   ! complex - dual_complex
!   ! <res, dres> = c - <v, dv> = <c - v, -dv>
!   !-------------------------------------------------
!   elemental function minus_c_dc(c, v) result(res)
!        complex(wp), intent(in) :: c
!        type(dual_complex), intent(in) :: v
!        type(dual_complex) :: res

!       res%z = c - v%z
!       res%dz = -v%dz

!   end function minus_c_dc
! !
! ! !******* end: (-)
! ! !---------------------
! !
! !
! !******* BEGIN: (*)
! !---------------------

!   !----------------------------------------
!   ! dual_complex * dual_complex
!   ! <res, dres> = <u, du> * <v, dv> = <u * v, u * dv + v * du>
!   !----------------------------------------
!   elemental function mult_dc_dc(u, v) result(res)
!       type(dual_complex), intent(in) :: u, v
!       type(dual_complex) :: res

!       res%z = u%z * v%z
!       res%dz = u%z * v%dz + v%z * u%dz

!   end function mult_dc_dc


!   !-----------------------------------------
!   ! dual_complex * integer
!   ! <res, dres> = <u, du> * j = <u * j du * i>
!   !-----------------------------------------
!   elemental function mult_dc_i(u, j) result(res)
!       type(dual_complex), intent(in) :: u
!       integer, intent(in) :: j
!       type(dual_complex) :: res

!       real(wp):: r

!       r = real(j,kind=wp)
!       res%z = r * u%z
!       res%dz = r * u%dz

!   end function mult_dc_i

!   !-----------------------------------------
!   ! dual_complex * real
!   ! <res, dres> = <u, du> * r = <u * r, du * r>
!   !----------------------------------------
!   elemental function mult_dc_r(u, r) result(res)
!       type(dual_complex), intent(in) :: u
!       real(wp), intent(in) :: r
!       type(dual_complex) :: res

!       res%z = u%z * r
!       res%dz = u%dz * r

!   end function mult_dc_r


!   !-----------------------------------------
!   ! integer * dual_complex
!   ! <res, dres> = j * <v, dv> = <i * v, j * dv>
!   !-----------------------------------------
!   elemental function mult_i_dc(j,v) result(res)
!       integer, intent(in) :: j
!       type(dual_complex), intent(in) :: v
!       type(dual_complex) :: res

!       real(wp):: r

!       r = real(j,kind=wp)
!       res%z = r * v%z
!       res%dz = r * v%dz

!   end function mult_i_dc


!   !-----------------------------------------
!   ! double * dual_complex
!   ! <res, dres> = r * <v, dv> = <r * v, r * dv>
!   !-----------------------------------------
!   elemental function mult_r_dc(r, v) result(res)
!       real(wp), intent(in) :: r
!       type(dual_complex), intent(in) :: v
!       type(dual_complex) :: res

!       res%z = r * v%z
!       res%dz = r * v%dz

!   end function mult_r_dc


!   !-----------------------------------------
!   ! complex * dual_complex
!   ! <res, dres> = c * <v, dv> = <c * v, c * dv>
!   !-----------------------------------------
!   elemental function mult_c_dc(c, v) result(res)
!       complex(wp), intent(in) :: c
!       type(dual_complex), intent(in) :: v
!       type(dual_complex) :: res

!       res%z = c * v%z
!       res%dz = c * v%dz

!   end function mult_c_dc


!   !-----------------------------------------
!   ! dual_complex * complex
!   ! <res, dres> = c * <v, dv> = <c * v, c * dv>
!   !-----------------------------------------
!   elemental function mult_dc_c(v, c) result(res)
!       complex(wp), intent(in) :: c
!       type(dual_complex), intent(in) :: v
!       type(dual_complex) :: res

!       res%z = c * v%z
!       res%dz = c * v%dz

!   end function mult_dc_c
! !
! ! !******* end: (*)
! ! !---------------------
! !
! !
! ! !******* BEGIN: (/)
! ! !---------------------
! !
!   !-----------------------------------------
!   ! dual_complex / dual_complex
!   ! <res, dres> = <u, du> / <v, dv> = <u / v, du / v - u * dv / v^2>
!   !-----------------------------------------
!   elemental function div_dc_dc(u, v) result(res)
!       type(dual_complex), intent(in) :: u, v
!       type(dual_complex) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / v%z
!       res%z = u%z * inv
!       res%dz = (u%dz - res%z * v%dz) * inv

!   end function div_dc_dc
! !
! !
!   !-----------------------------------------
!   ! dual_complex / integer
!   ! <res, dres> = <u, du> / j = <u / j du / i>
!   !-----------------------------------------
!   elemental function div_dc_i(u, j) result(res)
!       type(dual_complex), intent(in) :: u
!       integer, intent(in) :: j
!       type(dual_complex) :: res

!       real(wp):: inv

!       inv = 1.0d0 / real(j,kind=wp)
!       res%z = u%z * inv
!       res%dz = u%dz * inv

!   end function div_dc_i
! !
! !
!   !-----------------------------------------
!   ! dual_complex / double
!   ! <res, dres> = <u, du> / r = <u / r, du / r>
!   !----------------------------------------
!   elemental function div_dc_r(u, r) result(res)
!       type(dual_complex), intent(in) :: u
!       real(wp), intent(in) :: r
!       type(dual_complex):: res

!       real(wp):: inv

!       inv = 1.0d0 / r
!       res%z = u%z * inv
!       res%dz = u%dz * inv

!   end function div_dc_r
! !
! !
!   !-----------------------------------------
!   ! integer / dual_complex
!   ! <res, dres> = j / <v, dv> = <i / v, -i / v^2 * du>
!   !-----------------------------------------
!   elemental function div_i_dc(j,v) result(res)
!       integer, intent(in) :: j
!       type(dual_complex), intent(in) :: v
!       type(dual_complex) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / v%z
!       res%z = real(j,kind=wp) * inv
!       res%dz = -res%z * inv * v%dz

!   end function div_i_dc


!   !-----------------------------------------
!   ! double / dual_complex
!   ! <res, dres> = r / <u, du> = <r / u, -r / u^2 * du>
!   !-----------------------------------------
!   elemental function div_r_dc(r, v) result(res)
!       real(wp), intent(in) :: r
!       type(dual_complex), intent(in) :: v
!       type(dual_complex) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / v%z
!       res%z = r * inv
!       res%dz = -res%z * inv * v%dz

!   end function div_r_dc


!   !-----------------------------------------
!   ! complex / dual_complex
!   ! <res, dres> = c / <u, du> = <c / u, -r / u^2 * du>
!   !-----------------------------------------
!   elemental function div_c_dc(c, v) result(res)
!       complex(wp), intent(in) :: c
!       type(dual_complex), intent(in) :: v
!       type(dual_complex) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / v%z
!       res%z = c * inv
!       res%dz = -res%z * inv * v%dz

!   end function div_c_dc


!   !-----------------------------------------
!   ! dual_complex / complex
!   ! <res, dres> = c / <u, du> = <c / u, -c / u^2 * du>
!   !-----------------------------------------
!   elemental function div_dc_c(v, c) result(res)
!       complex(wp), intent(in) :: c
!       type(dual_complex), intent(in) :: v
!       type(dual_complex) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / c
!       res%z = v%z * inv
!       res%dz = v%dz * inv

!   end function div_dc_c

! !******* end: (/)
! !---------------------

! !******* BEGIN: (**)
! !---------------------

!   !-----------------------------------------
!   ! power(dual_complex, integer)
!   ! <res, dres> = <u, du> ^ j = <u ^ j j * u ^ (j- 1) * du>
!   !-----------------------------------------
!   elemental function pow_dc_i(u, j) result(res)
!       type(dual_complex), intent(in) :: u
!       integer, intent(in) :: j
!       type(dual_complex) :: res

!       complex(wp):: pow_x

!       pow_x = u%z ** (j- 1)
!       res%z = u%z * pow_x
!       res%dz = real(j,kind=wp) * pow_x * u%dz

!   end function pow_dc_i

!   !-----------------------------------------
!   ! power(dual_complex, double)
!   ! <res, dres> = <u, du> ^ r = <u ^ r, r * u ^ (r - 1) * du>
!   !-----------------------------------------
!   elemental function pow_dc_r(u, r) result(res)
!       type(dual_complex), intent(in) :: u
!       real(wp), intent(in) :: r
!       type(dual_complex) :: res

!       complex(wp):: pow_x

!       pow_x = u%z ** (r - 1.0d0)
!       res%z = u%z * pow_x
!       res%dz = r * pow_x * u%dz

!   end function pow_dc_r

!   !-----------------------------------------
!   ! power(dual_complex, complex)
!   ! <res, dres> = <u, du> ^ c = <u ^ c, c * u ^ (c - 1) * du>
!   !-----------------------------------------
!   elemental function pow_dc_c(u, c) result(res)
!       type(dual_complex), intent(in) :: u
!       complex(wp), intent(in) :: c
!       type(dual_complex) :: res

!       complex(wp) :: pow_x

!       pow_x = u%z ** (c - 1.0d0)
!       res%z = u%z * pow_x
!       res%dz = c * pow_x * u%dz

!   end function pow_dc_c

!   !-----------------------------------------
!   ! POWER dual_complex numbers to a dual_complex power
!   ! <res, dres> = <u, du> ^ <v, dv>
!   !     = <u ^ v, u ^ v * (v / u * du + Log(u) * dv)>
!   !-----------------------------------------
!   elemental function pow_dc_dc(u, v) result(res)
!       type(dual_complex), intent(in) :: u, v
!       type(dual_complex) :: res

!       res%z = u%z ** v%z
!       res%dz = res%z * (v%z / u%z * u%dz + log(u%z) * v%dz)

!   end function pow_dc_dc

!   !-----------------------------------------
!   ! POWER integer numbers to a dual_complex power
!   ! <res, dres> = j ^ <u, du>
!   !     = <i ^ u, j ^ u * Log(j) * du>
!   !-----------------------------------------
!   elemental function pow_i_dc(j,u) result(res)
!       integer, intent(in)   :: j
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = j ** u%z
!       res%dz = res%z * log(real(j,kind=wp)) * u%dz

!   end function pow_i_dc

!   !-----------------------------------------
!   ! POWER real numbers to a dual_complex power
!   ! <res, dres> = r ^ <u, du>
!   !     = <r ^ u, r ^ u * Log(r) * du)>
!   !-----------------------------------------
!   elemental function pow_r_dc(r, u) result(res)
!       real(wp), intent(in) :: r
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = r ** u%z
!       res%dz = res%z * log(r) * u%dz

!   end function pow_r_dc

!   !-----------------------------------------
!   ! POWER complex number to a dual_complex power
!   ! <res, dres> = c ^ <u, du>
!   !     = <c ^ u, c ^ u * Log(c) * du)>
!   !-----------------------------------------
!   elemental function pow_c_dc(c, u) result(res)
!       complex(wp), intent(in) :: c
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = c ** u%z
!       res%dz = res%z * log(c) * u%dz

!   end function pow_c_dc

! !******* end: (**)
! !---------------------


! !******* BEGIN: (==)
! !---------------------
!   !-----------------------------------------
!   ! compare two dual_complex numbers,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_dc_dc(lhs, rhs) result(res)
!        type(dual_complex), intent(in) :: lhs, rhs
!        logical :: res

!        res = (lhs%z == rhs%z)

!   end function eq_dc_dc


!   !-----------------------------------------
!   ! compare a dual_complex with an integer,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_dc_i(lhs, rhs) result(res)
!        type(dual_complex), intent(in) :: lhs
!        integer, intent(in) :: rhs
!        logical :: res

!        res = (lhs%z == complex(real(rhs), 0.0d0))

!   end function eq_dc_i


!   !-----------------------------------------
!   ! compare a dual_complex number with a real number,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_dc_r(lhs, rhs) result(res)
!       type(dual_complex), intent(in) :: lhs
!       real(wp), intent(in) :: rhs
!       logical::res

!       res = (lhs%z == complex(rhs, 0.0d0))

!   end function eq_dc_r


!   !-----------------------------------------
!   ! compare an integer with a dual_complex,
!   ! simply compare the functional value.
!   !----------------------------------------
!   elemental function eq_i_dc(lhs, rhs) result(res)
!        integer, intent(in) :: lhs
!        type(dual_complex), intent(in) :: rhs
!        logical :: res

!        res = (complex(real(lhs), 0.0d0) == rhs%z)

!   end function eq_i_dc


!   !-----------------------------------------
!   ! compare a real with a dual_complex,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_r_dc(lhs, rhs) result(res)
!        real(wp), intent(in) :: lhs
!        type(dual_complex), intent(in) :: rhs
!        logical :: res

!        res = (complex(lhs, 0.0d0) == rhs%z)

!   end function eq_r_dc


!   !-----------------------------------------
!   ! compare a complex number with a dual_complex,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_c_dc(lhs, rhs) result(res)
!        real(wp), intent(in) :: lhs
!        type(dual_complex), intent(in) :: rhs
!        logical :: res

!        res = (lhs == rhs%z)

!   end function eq_c_dc


!   !-----------------------------------------
!   ! compare a dual_complex with a comlex,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_dc_c(lhs, rhs) result(res)
!        type(dual_complex), intent(in) :: lhs
!        complex(wp), intent(in) :: rhs
!        logical :: res

!        res = (lhs%z == rhs)

!   end function eq_dc_c

! !******* end: (==)
! !---------------------


! !******* can't compare complex numbers
! ! (<=)
! ! (<)
! ! (>=)
! ! (>)
! !
! !******* BEGIN: (/=)
! !---------------------
!   !-----------------------------------------
!   ! compare two dual_complex numbers, simply compare
!   ! the functional value.
!   !-----------------------------------------
!   elemental function ne_dc_dc(lhs, rhs) result(res)
!       type(dual_complex), intent(in) :: lhs, rhs
!       logical :: res

!       res = (lhs%z /= rhs%z)

!   end function ne_dc_dc

!   !-----------------------------------------
!   ! compare a dual_complex number with a complex number,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function ne_dc_c(lhs, rhs) result(res)
!       type(dual_complex), intent(in) :: lhs
!       complex(wp), intent(in) :: rhs
!       logical :: res

!       res = (lhs%z /= rhs)

!   end function ne_dc_c

!   !-----------------------------------------
!   ! compare a complex with a dual_complex
!   !-----------------------------------------
!   elemental function ne_c_dc(lhs, rhs) result(res)
!       complex(wp), intent(in) :: lhs
!       type(dual_complex), intent(in) :: rhs
!       logical :: res

!       res = (lhs /= rhs%z)

!   end function ne_c_dc


!   !-----------------------------------------
!   ! compare a dual_complex with an integer,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function ne_dc_i(lhs, rhs) result(res)
!        type(dual_complex), intent(in) :: lhs
!        integer, intent(in) :: rhs
!        logical :: res

!        res = (lhs%z /= complex(real(rhs), 0.0d0))

!   end function ne_dc_i


!   !-----------------------------------------
!   ! compare a dual_complex number with a real number,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function ne_dc_r(lhs, rhs) result(res)
!       type(dual_complex), intent(in) :: lhs
!       real(wp), intent(in) :: rhs
!       logical::res

!       res = (lhs%z /= complex(rhs, 0.0d0))

!   end function ne_dc_r


!   !-----------------------------------------
!   ! compare an integer with a dual_complex,
!   ! simply compare the functional value.
!   !----------------------------------------
!   elemental function ne_i_dc(lhs, rhs) result(res)
!        integer, intent(in) :: lhs
!        type(dual_complex), intent(in) :: rhs
!        logical :: res

!        res = (complex(real(lhs), 0.0d0) /= rhs%z)

!   end function ne_i_dc


!   !-----------------------------------------
!   ! compare a real with a dual_complex,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function ne_r_dc(lhs, rhs) result(res)
!        real(wp), intent(in) :: lhs
!        type(dual_complex), intent(in) :: rhs
!        logical :: res

!        res = (complex(lhs, 0.0d0) /= rhs%z)

!   end function ne_r_dc

! !******* end: (/=)
! !---------------------
! !
!   !---------------------------------------------------
!   ! Absolute value of dual_complex numbers
!   ! <res, dres> = abs(<u, du>) = <abs(u), du * u / |u|>
!   !---------------------------------------------------
!   elemental function abs_dc(u) result(res)
!        type(dual_complex), intent(in) :: u
!        type(dual_complex) :: res

!        res%z = complex(abs(u%z), 0.)
!        res%dz = u%dz * u%z / abs(u%z)

!   end function abs_dc
! !
! !
!   !-----------------------------------------
!   ! ACOS of dual_complex numbers
!   ! <res, dres> = acos(<u, du>) = <acos(u), -du / sqrt(1 - u^2)>
!   !----------------------------------------
!   elemental function acos_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = acos(u%z)
!       if (u%z == complex(1.0d0, 0.0d0) .or. u%z == complex(-1.0d0, 0.0d0)) then
!           res%dz = set_Nan()  ! Undefined derivative
!       else
!           res%dz = -u%dz / sqrt(1.0d0 - u%z**2)
!       end if

!   end function acos_dc


!   !-----------------------------------------
!   ! ASIN of dual_complex numbers
!   ! <res, dres> = asin(<u, du>) = <asin(u), du / sqrt(1 - u^2)>
!   !----------------------------------------
!   elemental function asin_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = asin(u%z)
!       if (u%z == complex(1.0d0, 0.0d0) .or. u%z == complex(-1.0d0, 0.0d0)) then
!           res%dz = set_NaN()  ! Undefined derivative
!       else
!           res%dz = u%dz / sqrt(1.0d0 - u%z**2)
!       end if

!   end function asin_dc


!   !-----------------------------------------
!   ! ATAN of dual_complex numbers
!   ! <res, dres> = atan(<u, du>) = <atan(u), du / (1 + u^2)>
!   !----------------------------------------
!   elemental function atan_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = atan(u%z)
!       res%dz = u%dz / (1.0d0 + u%z**2)

!   end function atan_dc


!   !-----------------------------------------
!   ! ATAN2 of dual_complex numbers
!   ! <res, dres> = atan2(<u, du>, <v, dv>)
!   !             = <atan2(u, v), v / (u^2 + v^2) * du - u / (u^2 + v^2) * dv>
!   !----------------------------------------
!   ! elemental function atan2_dc(u, v) result(res)
!   !     type(dual_complex), intent(in) :: u, v
!   !     type(dual_complex) :: res
!   !
!   !     complex(wp) :: usq_plus_vsq
!   !
!   !     res%z = atan2(u%z, v%z)
!   !
!   !     usq_plus_vsq = u%z**2 + v%z**2
!   !     res%dz = v%z / usq_plus_vsq * u%dz - u%z / usq_plus_vsq * v%dz
!   !
!   ! end function atan2_dc
! !
! !
!   !-----------------------------------------
!   ! COS of dual_complex numbers
!   ! <res, dres> = cos(<u, du>) = <cos(u), -sin(u) * du>
!   !----------------------------------------
!   elemental function cos_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = cos(u%z)
!       res%dz = -sin(u%z) * u%dz

!   end function cos_dc
! !
! !
! !   !-----------------------------------------
! !   ! DOT PRODUCT two dual_complex number vectors
! !   ! <res, dres> = <u, du> . <v, dv> = <u . v, u . dv + v . du>
! !   !-----------------------------------------
! !   function dot_product_dc_d(u, v) result(res)
! !       type(dual_complex), intent(in) :: u(:), v(:)
! !       type(dual_complex) :: res
! !
! !       integer :: j
! !
! !       res%z = dot_product(u%z, v%z)
! !       do j = 1, dual_size
! !           res%dz(j) = dot_product(u%z, v%dz(j)) + dot_product(v%z, u%dz(j))
! !       end do
! !
! !   end function dot_product_dc_d
! !
! !
!   !-----------------------------------------
!   ! EXPONENTIAL OF dual_complex numbers
!   ! <res, dres> = exp(<u, du>) = <exp(u), exp(u) * du>
!   !-----------------------------------------
!   elemental function exp_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       complex(wp) :: exp_x

!       exp_x = exp(u%z)
!       res%z = exp_x
!       res%dz = u%dz * exp_x

!   end function exp_dc
! !
! !
! !   !-----------------------------------------
! !   ! Convert dual_complex to integer
! !   ! j = int(<u, du>) = int(u)
! !   !----------------------------------------
! !   elemental function int_dc(u) result(res)
! !        type(dual_complex), intent(in) :: u
! !        integer :: res
! !
! !        res = int(u%z)
! !
! !   end function int_dc
! !
! !
!   !-----------------------------------------
!   ! LOG OF dual_complex numbers,defined for u%z>0 only
!   ! the error control should be done in the original code
!   ! in other words, if u%z<=0, it is not possible to obtain LOG.
!   ! <res, dres> = log(<u, du>) = <log(u), du / u>
!   !----------------------------------------
!   elemental function log_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / u%z
!       res%z = log(u%z)
!       res%dz = u%dz * inv

!   end function log_dc


!   !-----------------------------------------
!   ! LOG10 OF dual_complex numbers,defined for u%z>0 only
!   ! the error control should be done in the original code
!   ! in other words, if u%z<=0, it is not possible to obtain LOG.
!   ! <res, dres> = log10(<u, du>) = <log10(u), du / (u * log(10))>
!   ! LOG<u,up>=<LOG(u),up/u>
!   !----------------------------------------
!   elemental function log10_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / (u%z * log(10.0d0))
!       res%z = log(u%z) / log(10.0d0)
!       res%dz = u%dz * inv

!   end function log10_dc
! !
! !
! !   !-----------------------------------------
! !   ! MULTIPLY two dual_complex number matrices
! !   ! <res, dres> = <u, du> . <v, dv> = <u . v, du . v + u . dv>
! !   !----------------------------------------
! !   function matmul_dc_d(u,v) result(res)
! !       type(dual_complex), intent(in) :: u(:,:), v(:,:)
! !       type(dual_complex) :: res(size(u,1), size(v,2))
! !
! !       integer :: j
! !
! !       res%z = matmul(u%z, v%z)
! !       do j = 1, dual_size
! !           res%dz(j) = matmul(u%dz(j), v%z) + matmul(u%z, v%dz(j))
! !       end do
! !
! !   end function matmul_dc_d
! !
! !
! !   !-----------------------------------------
! !   ! MULTIPLY a dual_complex number matrix with a dual_complex number
! !   ! vector
! !   !
! !   ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
! !   !----------------------------------------
! !   function matmul_dc_v(u, v) result(res)
! !       type(dual_complex), intent(in) :: u(:,:), v(:)
! !       type(dual_complex) :: res(size(u,1))
! !       integer :: j
! !
! !       res%z = matmul(u%z, v%z)
! !       do j = 1, dual_size
! !           res%dz(j) = matmul(u%dz(j), v%z) + matmul(u%z, v%dz(j))
! !       end do
! !
! !   end function matmul_dc_v
! !
! !
! !   !-----------------------------------------
! !   ! MULTIPLY a dual_complex vector with a  dual_complex matrix
! !   !
! !   ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
! !   !----------------------------------------
! !   function matmul_vd(u, v) result(res)
! !       type(dual_complex), intent(in) :: u(:), v(:,:)
! !       type(dual_complex) :: res(size(v, 2))
! !       integer::i
! !
! !       res%z = matmul(u%z, v%z)
! !       do j = 1, dual_size
! !           res%dz(j) = matmul(u%dz(j), v%z) + matmul(u%z, v%dz(j))
! !       end do
! !
! !   end function matmul_vd
! !
! !   !-----------------------------------------
! !   ! Obtain the max of 2 to 5 dual_complex numbers
! !   !----------------------------------------
! !   elemental function max_dc_d(val1, val2, val3, val4,val5) result(res)
! !       type(dual_complex), intent(in) :: val1, val2
! !       type(dual_complex), intent(in), optional :: val3, val4,val5
! !       type(dual_complex) :: res
! !
! !       if (val1%z > val2%z) then
! !           res = val1
! !       else
! !           res = val2
! !       endif
! !       if(present(val3))then
! !          if(res%z < val3%z) res = val3
! !       endif
! !       if(present(val4))then
! !          if(res%z < val4%z) res = val4
! !       endif
! !       if(present(val5))then
! !          if(res%z < val5%z) res = val5
! !       endif
! !
! !   end function max_dc_d
! !
! !
! !   !-----------------------------------------
! !   ! Obtain the max of a dual_complex number and an integer
! !   !----------------------------------------
! !   elemental function max_dc_i(u, j) result(res)
! !       type(dual_complex), intent(in) :: u
! !       integer, intent(in) :: j
! !       type(dual_complex) :: res
! !
! !       if (u%z > j) then
! !           res = u
! !       else
! !           res = i
! !       endif
! !
! !   end function max_dc_i
! !
! !   !-----------------------------------------
! !   ! Obtain the max of a dual_complex number and a real number
! !   !----------------------------------------
! !   elemental function max_dc_r(u, r) result(res)
! !       type(dual_complex), intent(in) :: u
! !       real(wp), intent(in) :: r
! !       type(dual_complex) :: res
! !
! !       if (u%z > r) then
! !           res = u
! !       else
! !           res = r
! !       endif
! !
! !   end function max_dc_r
! !
! !
! !   !---------------------------------------------------
! !   ! Obtain the max of a real and a dual_complex
! !   !---------------------------------------------------
! !    elemental function max_rd(n, u) result(res)
! !       real(wp), intent(in) :: n
! !       type(dual_complex), intent(in) :: u
! !       type(dual_complex) :: res
! !
! !       if (u%z > n) then
! !           res = u
! !       else
! !           res = n
! !       endif
! !
! !   end function max_rd
! !
! !
! !   !-----------------------------------------
! !   ! Obtain the max value of vector u
! !   !----------------------------------------
! !   function maxval_dc(u) result(res)
! !       type(dual_complex), intent(in) :: u(:)
! !       integer :: iloc(1)
! !       type(dual_complex) :: res
! !
! !       iloc=maxloc(u%z)
! !       res=u(iloc(1))
! !
! !   end function maxval_dc
! !
! !
! !   !-----------------------------------------
! !   ! Obtain the min of 2 to 4 dual_complex numbers
! !   !----------------------------------------
! !   elemental function min_dc_d(val1, val2, val3, val4) result(res)
! !       type(dual_complex), intent(in) :: val1, val2
! !       type(dual_complex), intent(in), optional :: val3, val4
! !       type(dual_complex) :: res
! !
! !       if (val1%z < val2%z) then
! !           res = val1
! !       else
! !           res = val2
! !       endif
! !       if(present(val3))then
! !          if(res%z > val3%z) res = val3
! !       endif
! !       if(present(val4))then
! !          if(res%z > val4%z) res = val4
! !       endif
! !
! !   end function min_dc_d
! !
! !
! !   !-----------------------------------------
! !   ! Obtain the min of a dual_complex and a double
! !   !----------------------------------------
! !   elemental function min_dc_r(u, r) result(res)
! !       type(dual_complex), intent(in) :: u
! !       real(wp), intent(in) :: r
! !       type(dual_complex) :: res
! !
! !       if (u%z < r) then
! !           res = u
! !       else
! !           res = r
! !       endif
! !
! !   end function min_dc_r
! !
! !
! ! !-----------------------------------------
! !   ! Obtain the min value of vector u
! !   !----------------------------------------
! !   function minval_dc(u) result(res)
! !       type(dual_complex), intent(in) :: u(:)
! !       integer :: iloc(1)
! !       type(dual_complex) :: res
! !
! !       iloc=minloc(u%z)
! !       res=u(iloc(1))
! !
! !   end function minval_dc
! !
! !
! !   !------------------------------------------------------
! !   !Returns the nearest integer to u%z, ELEMENTAL
! !   !------------------------------------------------------
! !   elemental function nint_dc(u) result(res)
! !       type(dual_complex), intent(in) :: u
! !       integer :: res
! !
! !       res=nint(u%z)
! !
! !   end function nint_dc
! !
! !
! !   !----------------------------------------------------------------
! !   ! SIGN(a,b) with two dual_complex numbers as inputs,
! !   ! the result will be |a| if b%z>=0, -|a| if b%z<0,ELEMENTAL
! !   !----------------------------------------------------------------
! !   elemental function sign_dc_d(val1, val2) result(res)
! !       type(dual_complex), intent(in) :: val1, val2
! !       type(dual_complex) :: res
! !
! !       if (val2%z < 0.0d0) then
! !           res = -abs(val1)
! !       else
! !           res =  abs(val1)
! !       endif
! !
! !    end function sign_dc_d
! !
! !
! !   !----------------------------------------------------------------
! !   ! SIGN(a,b) with one real and one dual_complex number as inputs,
! !   ! the result will be |a| if b%z>=0, -|a| if b%z<0,ELEMENTAL
! !   !----------------------------------------------------------------
! !   elemental function sign_rd(val1, val2) result(res)
! !       real(wp), intent(in) :: val1
! !       type(dual_complex), intent(in) :: val2
! !       type(dual_complex) :: res
! !
! !       if (val2%z < 0.0d0) then
! !           res = -abs(val1)
! !       else
! !           res = abs(val1)
! !       endif
! !
! !    end function sign_rd
! !
! !
! !   !-----------------------------------------
! !   ! SIN of dual_complex numbers
! !   ! <res, dres> = sin(<u, du>) = <sin(u), cos(u) * du>
! !   !----------------------------------------
!   elemental function sin_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = sin(u%z)
!       res%dz = cos(u%z) * u%dz

!   end function sin_dc
! !
! !
!   !-----------------------------------------
!   ! TAN of dual_complex numbers
!   ! <res, dres> = tan(<u, du>) = <tan(u), du / cos(u)^2>
!   !----------------------------------------
!   elemental function tan_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = tan(u%z)
!       res%dz = u%dz / cos(u%z)**2

!   end function tan_dc
! !
! !
!   !-----------------------------------------
!   ! SQRT of dual_complex numbers
!   ! <res, dres> = sqrt(<u, du>) = <sqrt(u), du / (2 * sqrt(u))>
!   !----------------------------------------
!   elemental function sqrt_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res
!       integer :: j

!       res%z = sqrt(u%z)

!       if (res%z /= complex(0.0d0, 0.0d0)) then
!           res%dz = 0.5 * u%dz / res%z
!       else
!           do j = 1, dual_size
!               if (u%dz(j) == complex(0.0d0, 0.0d0)) then
!                   res%dz(j) = complex(0.0d0, 0.0d0)
!               else
!                   res%dz(j) = set_NaN()
!               end if
!           end do
!       end if

!   end function sqrt_dc


!   !-----------------------------------------
!   ! Sum of a dual_complex array
!   !-----------------------------------------
!   function sum_dc(u) result(res)
!       type(dual_complex), intent(in) :: u(:)
!       type(dual_complex) :: res
!       integer :: j

!       res%z = sum(u%z)
!       do j = 1, dual_size
!           res%dz(j) = sum(u%dz(j))
!       end do

!   end function sum_dc
! !
! !
! !   !-----------------------------------------
! !   ! Find the location of the max value in an
! !   ! array of dual_complex numbers
! !   !-----------------------------------------
! !   function maxloc_dc(array) result(ind)
! !       type(dual_complex), intent(in) :: array(:)
! !       integer :: ind(1)
! !
! !       ind = maxloc(array%z)
! !
! !   end function maxloc_dc
! !
! !
! !   elemental function set_NaN() result(res)
! !       real(wp):: res
! !
! !       res = sqrt(negative_one)
! !
! !   end function set_NaN
! !
! !
!   !-----------------------------------------
!   ! Hyperbolic functions: sinh, cosh, tanh
!   ! and their inverses: asinh, acosh, atanh
!   !-----------------------------------------
!   !-----------------------------------------
!   ! SINH OF dual_complex numbers
!   ! <res, dres> = sinh(<u, du>) = <sinh(u), cosh(u) * du>
!   !-----------------------------------------
!   elemental function sinh_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = sinh(u%z)
!       res%dz = u%dz * cosh(u%z)

!   end function sinh_dc

!   !-----------------------------------------
!   ! COSH OF dual_complex numbers
!   ! <res, dres> = cosh(<u, du>) = <cosh(u), sinh(u) * du>
!   !-----------------------------------------
!   elemental function cosh_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = cosh(u%z)
!       res%dz = u%dz * sinh(u%z)

!   end function cosh_dc

!   !-----------------------------------------
!   ! TANH OF dual_complex numbers
!   ! <res, dres> = tanh(<u, du>) = <tanh(u), 1.0d0/cosh(u)**2 * du>
!   !-----------------------------------------
!   elemental function tanh_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = tanh(u%z)
!       res%dz = u%dz * 1.0d0/cosh(u%z)**2

!   end function tanh_dc

!   !-----------------------------------------
!   ! ASINH OF dual_complex numbers
!   ! <res, dres> = asinh(<u, du>) = <asinh(u), 1/sqrt(u**2 + 1) * du>
!   !-----------------------------------------
!   elemental function asinh_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = asinh(u%z)
!       res%dz = u%dz * 1.0d0/sqrt(u%z**2 + 1.0d0)

!   end function asinh_dc

!   !-----------------------------------------
!   ! ACOSH OF dual_complex numbers
!   ! <res, dres> = acosh(<u, du>) = <acosh(u), 1/sqrt(u**2 - 1) * du>
!   !-----------------------------------------
!   elemental function acosh_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = acosh(u%z)
!       if (u%z == complex(1.0d0, 0.0d0)) then
!           res%dz = set_Nan()  ! Undefined derivative (∞)
!       else
!           res%dz = u%dz * 1.0d0/sqrt(u%z**2 - 1.0d0)
!       end if


!   end function acosh_dc

!   !-----------------------------------------
!   ! ATAHN OF dual_complex numbers
!   ! <res, dres> = atanh(<u, du>) = <atanh(u), 1/(1 - u**2) * du>
!   !-----------------------------------------
!   elemental function atanh_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z = atanh(u%z)
!       if (u%z == complex(1.0d0, 0.0d0) .or. u%z == complex(-1.0d0, 0.0d0)) then
!           res%dz = set_Nan()  ! Undefined derivative
!       else
!           res%dz = u%dz * 1.0d0/(1.0d0 - u%z**2)
!       end if

!   end function atanh_dc


!   ! -----------------------------------------
!   ! ERF OF complex number (see Abramowitz & Stegun: Handbook of Mathematical Functions)
!   ! https://math.stackexchange.com/questions/712434/erfaib-error-function-separate-into-real-and-imaginary-part
!   ! https://personal.math.ubc.ca/~cbm/aands/
!   !
!   ! erf(x + iy) = erf(x) + e⁻ˣ² / (2πx) [(1 - cos(2xy)) + i sin(2xy)] +
!   !                 2/π e⁻ˣ² ∑ₐ₌₁^(∞) [e^(-a²/4)/(a² + 4x²) ( fₐ(x,y) + i gₐ(x,y) ) ] + ϵ(x,y)
!   ! fₐ(x,y) = 2x(1 - cos(2xy) cosh(ay)) + a sin(2xy) sinh(ay)
!   ! gₐ(x,y) = 2x sin(2xy) cosh(ay) + a cos(2xy) sinh(ay)

!   ! for improved numerical stability and accuracy, push e⁻ˣ² * e^(-a²/4) inside of fₐ(x,y) and gₐ(x,y)
!   ! i.e. cosh(ay) * e⁻ˣ² * e^(-a²/4) = 0.5*(exp(ay - x² - a²/4) + exp(-ay - x² - a²/4))
!   ! i.e. sinh(ay) * e⁻ˣ² * e^(-a²/4) = 0.5*(exp(ay - x² - a²/4) - exp(-ay - x² - a²/4))
!   ! this is important because sinh(ay) --> ∞ as a --> ∞
!   ! but e^(-a²/4) --> 0 as a --> ∞
!   ! -----------------------------------------
!     elemental function erf_c(c) result(res)
!         complex(wp), intent(in) :: c
!         complex(wp) :: res

!         real(wp):: x, y
!         complex(wp), parameter :: j=(0,1)
        
       
!         integer :: k

!         x = real(c)
!         y = aimag(c)

!         ! res = erf(x) + exp(-x**2) / (2.*PI*x)*((1. - cos(2*x*y)) + i*sin(2*x*y))
!         res = erf(x) + 1. / (2.*PI*x)*((1. - cos(2*x*y)) + j*sin(2*x*y))

!         do k = 1, 100   ! inefficiently defined
!             ! res = res + (2./PI * exp(-x**2)) * (exp(-k**2/4.0)/(k**2 + 4*x**2)) * (f_k(k,x,y) + i*g_k(k,x,y))
!             res = res + 2./PI * (1./(k**2 + 4*x**2)) * (f_k(k,x,y) + j*g_k(k,x,y))
!         end do


!         contains
!             elemental function f_k(k, x, y) result(res)
!                 implicit none

!                 integer, intent(in) :: k
!                 real(wp), intent(in) :: x, y
!                 real(wp):: res
!                 ! 0.5*(exp(k*y - x**2 - k**2/4.0) + exp(-k*y - x**2 - k**2/4.0)) = cosh(k*y) * exp(-x**2 - k**2/4.0)
!                 ! 0.5*(exp(k*y - x**2 - k**2/4.0) - exp(-k*y - x**2 - k**2/4.0)) = sinh(k*y) * exp(-x**2 - k**2/4.0)

!                 ! res = 2.*x*(1. - cos(2*x*y)*cosh(k*y)) + k*sin(2*x*y)*sinh(k*y)
!                 ! res = ( 2.*x*(1. - cos(2*x*y)*cosh(k*y)) + k*sin(2*x*y)*sinh(k*y) ) * exp(-x**2 - k**2/4.0)
!                 ! res = ( 2.*x*(1. - cos(2*x*y)*( 0.5*(exp(k*y)+exp(-k*y)) )   ) + &
!                 !     & k*sin(2*x*y)*( 0.5*(exp(k*y)-exp(-k*y)) ) ) * exp(-x**2 - k**2/4.0)
!                 res = 2.*x*(exp(-x**2 - k**2/4.0) &
!                     & - cos(2*x*y)*( 0.5*(exp(k*y - x**2 - k**2/4.0) + exp(-k*y - x**2 - k**2/4.0)) ) ) + &
!                     & k*sin(2*x*y)*( 0.5*(exp(k*y - x**2 - k**2/4.0) - exp(-k*y - x**2 - k**2/4.0)) )


!             end function f_k

!             elemental function g_k(k, x, y) result(res)
!                 implicit none

!                 integer, intent(in) :: k
!                 real(wp), intent(in) :: x, y
!                 real(wp):: res
!                 ! 0.5*(exp(k*y - x**2 - k**2/4.0) + exp(-k*y - x**2 - k**2/4.0)) = cosh(k*y) * exp(-x**2 - k**2/4.0)
!                 ! 0.5*(exp(k*y - x**2 - k**2/4.0) - exp(-k*y - x**2 - k**2/4.0)) = sinh(k*y) * exp(-x**2 - k**2/4.0)

!                 ! res = 2*x*sin(2*x*y)*cosh(k*y) + k*cos(2*x*y)*sinh(k*y)
!                 ! res = ( 2*x*sin(2*x*y)*cosh(k*y) + k*cos(2*x*y)*sinh(k*y) ) * exp(-x**2 - k**2/4.0)
!                 ! res = ( 2*x*sin(2*x*y)*( 0.5*(exp(k*y)+exp(-k*y)) ) + &
!                 !     & k*cos(2*x*y)*( 0.5*(exp(k*y)-exp(-k*y)) ) ) * exp(-x**2 - k**2/4.0)
!                 res = 2*x*sin(2*x*y)*( 0.5*(exp(k*y - x**2 - k**2/4.0) + exp(-k*y - x**2 - k**2/4.0)) ) + &
!                     &   k*cos(2*x*y)*( 0.5*(exp(k*y - x**2 - k**2/4.0) - exp(-k*y - x**2 - k**2/4.0)) )
!             end function g_k

!     end function erf_c

!     ! -----------------------------------------
!     ! ERF OF dual_complex numbers
!     ! <res, dres> = ERF(<u, du>) = <ERF(u), 2/sqrt(PI)*exp(-u**2) * du>
!     ! -----------------------------------------
!     elemental function erf_dc(u) result(res)
!       type(dual_complex), intent(in) :: u
!       type(dual_complex) :: res

!       res%z  = erf(u%z)
!       res%dz = u%dz * 2.0/sqrt(PI) * exp(-u%z**2)

!     end function erf_dc


    !*********Begin: functions/subroutines for overloading operators

!******* Begin: (=)
!---------------------

    !-----------------------------------------
    ! dual = integer
    ! <u, du> = <j 0>
    !-----------------------------------------
    elemental subroutine assign_di2(u, j)
         type(dual2), intent(out) :: u
         integer, intent(in) :: j

         u%x =j  ! This is faster than direct assignment
         u%dx = 0.0d0

    end subroutine assign_di2



    !-----------------------------------------
    ! dual = real(double)
    ! <u, du> = <r, 0>
    !-----------------------------------------
    elemental subroutine assign_dr2(u, r)
        !$acc routine seq
        type(dual2), intent(out) :: u
        real(wp), intent(in) :: r
       
        u%x = real(r,kind=wp)
        u%dx = 0.0d0

    end subroutine assign_dr2


    !-----------------------------------------
    ! integer = dual2
    ! j = <u, du>
    !-----------------------------------------
    elemental subroutine assign_id2(j, v)
         type(dual2), intent(in) :: v
         integer, intent(out) :: j

         j = int(v%x)

    end subroutine assign_id2

    !-----------------------------------------
    ! real = dual2
    ! j = <u, du>
    !-----------------------------------------
    elemental subroutine assign_rd2(j, v)
         type(dual2), intent(in) :: v
         real(wp), intent(out) :: j

         j = real(v%x)

    end subroutine assign_rd2

!******* end: (=)
!---------------------


!******* Begin: (+)
!---------------------

    !-----------------------------------------
    ! Unary positive
    ! <res, dres> = +<u, du>
    !-----------------------------------------
    elemental function add_d2(u) result(res)
         type(dual2), intent(in) :: u
         type(dual2) :: res

         res = u  ! Faster than assigning component wise

    end function add_d2


    !-----------------------------------------
    ! dual2 + dual2
    ! <res, dres> = <u, du> + <v, dv> = <u + v, du + dv>
    !-----------------------------------------
    elemental function add_dd2(u, v) result(res)
        !$acc routine seq
         type(dual2), intent(in) :: u, v
         type(dual2) :: res

         res%x = u%x + v%x
         res%dx = u%dx + v%dx

    end function add_dd2


    !-----------------------------------------
    ! dual2 + integer
    ! <res, dres> = <u, du> + j = <u + j du>
    !-----------------------------------------
    elemental function add_di2(u, j) result(res)
         type(dual2), intent(in) :: u
         integer, intent(in) :: j
         type(dual2) :: res

         res%x = real(j,kind=wp) + u%x
         res%dx = u%dx

    end function add_di2


    !-----------------------------------------
    ! dual2 + double
    ! <res, dres> = <u, du> + <r, 0> = <u + r, du>
    !-----------------------------------------
    elemental function add_dr2(u, r) result(res)
      !$acc routine seq
        type(dual2), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual2) :: res

        res%x = r + u%x
        res%dx = u%dx

    end function add_dr2


    !-----------------------------------------
    ! integer + dual2
    ! <res, dres> = <j 0> + <v, dv> = <i + v, dv>
    !-----------------------------------------
    elemental function add_id2(j, v) result(res)
        integer, intent(in) :: j
        type(dual2), intent(in) :: v
        type(dual2) :: res

        res%x = real(j,kind=wp) + v%x
        res%dx = v%dx

    end function add_id2


    !-----------------------------------------
    ! double + dual2
    ! <res, dres> = <r, 0> + <v, dv> = <r + v, dv>
    !-----------------------------------------
    elemental function add_rd2(r, v) result(res)
        real(wp), intent(in) :: r
        type(dual2), intent(in) :: v
        type(dual2) :: res

        res%x = r + v%x
        res%dx = v%dx

    end function add_rd2

!******* end: (+)
!---------------------


!******* Begin: (-)
!---------------------

    !-------------------------------------------------
    ! negate a dual2
    ! <res, dres> = -<u, du>
    !-------------------------------------------------
    elemental function minus_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = -u%x
        res%dx = -u%dx

    end function minus_d2


    !-------------------------------------------------
    ! dual2 - dual2
    ! <res, dres> = <u, du> - <v, dv> = <u - v, du - dv>
    !-------------------------------------------------
    elemental function minus_dd2(u, v) result(res)
        type(dual2), intent(in) :: u, v
        type(dual2) :: res

        res%x = u%x - v%x
        res%dx = u%dx - v%dx

    end function minus_dd2

    !-------------------------------------------------
    ! dual2 - integer
    ! <res, dres> = <u, du> - j = <u - j du>
    !-------------------------------------------------
    elemental function minus_di2(u, j) result(res)
        type(dual2), intent(in) :: u
        integer, intent(in) :: j
        type(dual2) :: res

        res%x = u%x - real(j,kind=wp)
        res%dx = u%dx

    end function minus_di2


    !-------------------------------------------------
    ! dual2 - double
    ! <res, dres> = <u, du> - r = <u - r, du>
    !-------------------------------------------------
    elemental function minus_dr2(u, r) result(res)
        type(dual2), intent(in) :: u
        real(wp),intent(in) :: r
        type(dual2) :: res

        res%x = u%x - r
        res%dx = u%dx

    end function minus_dr2


    !-------------------------------------------------
    ! integer - dual2
    ! <res, dres> = j - <v, dv> = <i - v, -dv>
    !-------------------------------------------------
    elemental function minus_id2(j,v) result(res)
        integer, intent(in) :: j
        type(dual2), intent(in) :: v
        type(dual2) :: res

        res%x = real(j,kind=wp) - v%x
        res%dx = -v%dx

    end function minus_id2


    !-------------------------------------------------
    ! double - dual2
    ! <res, dres> = r - <v, dv> = <r - v, -dv>
    !-------------------------------------------------
    elemental function minus_rd2(r, v) result(res)
         real(wp), intent(in) :: r
         type(dual2), intent(in) :: v
         type(dual2) :: res

        res%x = r - v%x
        res%dx = -v%dx

    end function minus_rd2

!******* end: (-)
!---------------------


!******* BEGIN: (*)
!---------------------

    !----------------------------------------
    ! dual2 * dual2
    ! <res, dres> = <u, du> * <v, dv> = <u * v, u * dv + v * du>
    !----------------------------------------
    elemental function mult_dd2(u, v) result(res)
        type(dual2), intent(in) :: u, v
        type(dual2) :: res

        res%x = u%x * v%x
        res%dx = u%x * v%dx + v%x * u%dx

    end function mult_dd2


    !-----------------------------------------
    ! dual2 * integer
    ! <res, dres> = <u, du> * j = <u * j du * i>
    !-----------------------------------------
    elemental function mult_di2(u, j) result(res)
        type(dual2), intent(in) :: u
        integer, intent(in) :: j
        type(dual2) :: res

        real(wp):: r

        r = real(j,kind=wp)
        res%x = r * u%x
        res%dx = r * u%dx

    end function mult_di2

    !-----------------------------------------
    ! dual * integer
    ! <res, dres> = <u, du> * j = <u * j du * i>
    !-----------------------------------------
    elemental function mult_di162(u, j) result(res)
        type(dual2), intent(in) :: u
        integer(int16), intent(in) :: j
        type(dual2) :: res

        real(wp):: r

        r = real(j,kind=wp)
        res%x = r * u%x
        res%dx = r * u%dx

    end function mult_di162

            !-----------------------------------------
    ! dual * integer
    ! <res, dres> = <u, du> * j = <u * j du * i>
    !-----------------------------------------
    elemental function mult_di82(u, j) result(res)
        type(dual2), intent(in) :: u
        integer(int8), intent(in) :: j
        type(dual2) :: res

        real(wp):: r

        r = real(j,kind=wp)
        res%x = r * u%x
        res%dx = r * u%dx

    end function mult_di82

    !-----------------------------------------
    ! dual2 * double
    ! <res, dres> = <u, du> * r = <u * r, du * r>
    !----------------------------------------
    elemental function mult_dr2(u, r) result(res)
        !$acc routine seq
        type(dual2), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual2) :: res

        res%x = u%x * r
        res%dx = u%dx * r

    end function mult_dr2


    !-----------------------------------------
    ! integer * dual2
    ! <res, dres> = j * <v, dv> = <i * v, j * dv>
    !-----------------------------------------
    elemental function mult_id2(j,v) result(res)
        integer, intent(in) :: j
        type(dual2), intent(in) :: v
        type(dual2) :: res

        real(wp):: r

        r = real(j,kind=wp)
        res%x = r * v%x
        res%dx = r * v%dx

    end function mult_id2


    !-----------------------------------------
    ! double * dual2
    ! <res, dres> = r * <v, dv> = <r * v, r * dv>
    !-----------------------------------------
    elemental function mult_rd2(r, v) result(res)
        real(wp), intent(in) :: r
        type(dual2), intent(in) :: v
        type(dual2) :: res

        res%x = r * v%x
        res%dx = r * v%dx

    end function mult_rd2

!******* end: (*)
!---------------------


!******* BEGIN: (/)
!---------------------

    !-----------------------------------------
    ! dual2 / dual2
    ! <res, dres> = <u, du> / <v, dv> = <u / v, du / v - u * dv / v^2>
    !-----------------------------------------
    elemental function div_dd2(u, v) result(res)
        type(dual2), intent(in) :: u, v
        type(dual2) :: res

        real(wp):: inv

        inv = 1.0d0 / v%x
        res%x = u%x * inv
        res%dx = (u%dx - res%x * v%dx) * inv

    end function div_dd2


    !-----------------------------------------
    ! dual2 / integer
    ! <res, dres> = <u, du> / j = <u / j du / i>
    !-----------------------------------------
    elemental function div_di2(u, j) result(res)
        type(dual2), intent(in) :: u
        integer, intent(in) :: j
        type(dual2) :: res

        real(wp):: inv

        inv = 1.0d0 / real(j,kind=wp)
        res%x = u%x * inv
        res%dx = u%dx * inv

    end function div_di2


    !-----------------------------------------
    ! dual2 / double
    ! <res, dres> = <u, du> / r = <u / r, du / r>
    !----------------------------------------
    elemental function div_dr2(u, r) result(res)
        type(dual2), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual2):: res

        real(wp):: inv

        inv = 1.0d0 / r
        res%x = u%x * inv
        res%dx = u%dx * inv

    end function div_dr2


    !-----------------------------------------
    ! integer / dual2
    ! <res, dres> = j / <v, dv> = <i / v, -i / v^2 * du>
    !-----------------------------------------
    elemental function div_id2(j,v) result(res)
        integer, intent(in) :: j
        type(dual2), intent(in) :: v
        type(dual2) :: res

        real(wp):: inv

        inv = 1.0d0 / v%x
        res%x = real(j,kind=wp) * inv
        res%dx = -res%x * inv * v%dx

    end function div_id2


    !-----------------------------------------
    ! double / dual2
    ! <res, dres> = r / <u, du> = <r / u, -r / u^2 * du>
    !-----------------------------------------
    elemental function div_rd2(r, v) result(res)
        real(wp), intent(in) :: r
        type(dual2), intent(in) :: v
        type(dual2) :: res

        real(wp):: inv

        inv = 1.0d0 / v%x
        res%x = r * inv
        res%dx = -res%x * inv * v%dx

    end function div_rd2

!******* end: (/)
!---------------------

!******* BEGIN: (**)
!---------------------

    !-----------------------------------------
    ! power(dual2, integer)
    ! <res, dres> = <u, du> ^ j = <u ^ j j * u ^ (j- 1) * du>
    !-----------------------------------------
    elemental function pow_di2(u, j) result(res)
        type(dual2), intent(in) :: u
        integer, intent(in) :: j
        type(dual2) :: res

        real(wp):: pow_x

        pow_x = u%x ** (j- 1)
        res%x = u%x * pow_x
        res%dx = real(j,kind=wp) * pow_x * u%dx

    end function pow_di2

    !-----------------------------------------
    ! power(dual2, double)
    ! <res, dres> = <u, du> ^ r = <u ^ r, r * u ^ (r - 1) * du>
    !-----------------------------------------
    elemental function pow_dr2(u, r) result(res)
        type(dual2), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual2) :: res

        real(wp):: pow_x

        pow_x = u%x ** (r - 1.0d0)
        res%x = u%x * pow_x
        res%dx = r * pow_x * u%dx

    end function pow_dr2

    !-----------------------------------------
    ! POWER dual2 number to a dual2 power
    ! <res, dres> = <u, du> ^ <v, dv>
    !     = <u ^ v, u ^ v * (v / u * du + Log(u) * dv)>
    !-----------------------------------------
    elemental function pow_dd2(u, v) result(res)
        type(dual2), intent(in) :: u, v
        type(dual2) :: res

        res%x = u%x ** v%x
        res%dx = res%x * (v%x / u%x * u%dx + log(u%x) * v%dx)

    end function pow_dd2

    !-----------------------------------------
    ! POWER integer to a dual2 power
    ! <res, dres> = j ^ <u, du>
    !     = <i ^ u, j ^ u * Log(j) * du)>
    !-----------------------------------------
    elemental function pow_id2(j,u) result(res)
        integer, intent(in) :: j
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = j ** u%x
        res%dx = res%x * log(real(j,kind=wp)) * u%dx

    end function pow_id2

    !-----------------------------------------
    ! POWER real to a dual2 power
    ! <res, dres> = r ^ <u, du>
    !     = <r ^ u, r ^ u * Log(r) * du)>
    !-----------------------------------------
    elemental function pow_rd2(r, u) result(res)
        real(wp), intent(in) :: r
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = r ** u%x
        res%dx = res%x * log(r) * u%dx

    end function pow_rd2    

!******* end: (**)
!---------------------


!******* BEGIN: (==)
!---------------------
    !-----------------------------------------
    ! compare two dual2 numbers,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_dd2(lhs, rhs) result(res)
         type(dual2), intent(in) :: lhs, rhs
         logical :: res

         res = (lhs%x == rhs%x)

    end function eq_dd2


    !-----------------------------------------
    ! compare a dual2 with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_di2(lhs, rhs) result(res)
         type(dual2), intent(in) :: lhs
         integer, intent(in) :: rhs
         logical :: res

         res = (lhs%x == real(rhs))

    end function eq_di2


    !-----------------------------------------
    ! compare a dual2 number with a real number,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_dr2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs
        real(wp), intent(in) :: rhs
        logical::res

        res = (lhs%x == rhs)

    end function eq_dr2


    !-----------------------------------------
    ! compare an integer with a dual2,
    ! simply compare the functional value.
    !----------------------------------------
    elemental function eq_id2(lhs, rhs) result(res)
         integer, intent(in) :: lhs
         type(dual2), intent(in) :: rhs
         logical :: res

         res = (lhs == rhs%x)

    end function eq_id2


    !-----------------------------------------
    ! compare a real with a dual2,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_rd2(lhs, rhs) result(res)
         real(wp), intent(in) :: lhs
         type(dual2), intent(in) :: rhs
         logical :: res

         res = (lhs == rhs%x)

    end function eq_rd2

!******* end: (==)
!---------------------


!******* BEGIN: (<=)
!---------------------
    !-----------------------------------------
    ! compare two dual2 numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function le_dd2(lhs, rhs) result(res)
         type(dual2), intent(in) :: lhs, rhs
         logical :: res

         res = (lhs%x <= rhs%x)

    end function le_dd2


    !-----------------------------------------
    ! compare a dual2 with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_di2(lhs, rhs) result(res)
         type(dual2), intent(in) :: lhs
         integer, intent(in) :: rhs
         logical :: res

         res = (lhs%x <= rhs)

    end function le_di2


    !-----------------------------------------
    ! compare a dual2 number with a real number,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_dr2(lhs, rhs) result(res)
         type(dual2), intent(in) :: lhs
         real(wp), intent(in) :: rhs
         logical :: res

         res = (lhs%x <= rhs)

    end function le_dr2


    !-----------------------------------------
    ! compare a dual2 number with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_id2(j,rhs) result(res)
         integer, intent(in) :: j
         type(dual2), intent(in) :: rhs
         logical :: res

         res = (j<= rhs%x)

    end function le_id2


    !-----------------------------------------
    ! compare a real with a dual2,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_rd2(lhs, rhs) result(res)
         real(wp), intent(in) :: lhs
         type(dual2), intent(in) :: rhs
         logical :: res

         res = (lhs <= rhs%x)

    end function le_rd2

!******* end: (<=)
!---------------------

!******* BEGIN: (<)
!---------------------
    !-----------------------------------------
    ! compare two dual2 numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function lt_dd2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x < rhs%x)

    end function lt_dd2

    !-----------------------------------------
    ! compare a dual2 with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function lt_di2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x < rhs)

    end function lt_di2


    !-----------------------------------------
    ! compare a dual2 number with a real number, simply compare
    ! the functional value.
    !----------------------------------------
    elemental function lt_dr2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs
        real(wp), intent(in) :: rhs
        logical :: res

        res = (lhs%x < rhs)

    end function lt_dr2


    !-----------------------------------------
    ! compare a dual2 number with an integer
    !-----------------------------------------
    elemental function lt_id2(j,rhs) result(res)
         integer, intent(in) :: j
         type(dual2), intent(in) :: rhs
         logical :: res

         res = (j< rhs%x)

    end function lt_id2


    !-----------------------------------------
    ! compare a real with a dual2
    !----------------------------------------
    elemental function lt_rd2(lhs, rhs) result(res)
         real(wp), intent(in) :: lhs
         type(dual2), intent(in) :: rhs
         logical :: res

         res = (lhs < rhs%x)

    end function lt_rd2

!******* end: (<)
!---------------------

!******* BEGIN: (>=)
!---------------------
    !-----------------------------------------
    ! compare two dual2 numbers, simply compare
    ! the functional value.
    !----------------------------------------
    elemental function ge_dd2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x >= rhs%x)

    end function ge_dd2


    !-----------------------------------------
    ! compare a dual2 with an integer
    !-----------------------------------------
    elemental function ge_di2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x >= rhs)

    end function ge_di2


    !-----------------------------------------
    ! compare a dual2 number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ge_dr2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs
        real(wp), intent(in) :: rhs
        logical :: res

        res = (lhs%x >= rhs)

    end function ge_dr2


    !-----------------------------------------
    ! compare a dual2 number with an integer
    !-----------------------------------------
    elemental function ge_id2(j,rhs) result(res)
        integer, intent(in) :: j
        type(dual2), intent(in) :: rhs
        logical :: res

        res = (j>= rhs%x)

    end function ge_id2


    !-----------------------------------------
    ! compare a real with a dual2
    !-----------------------------------------
    elemental function ge_rd2(lhs, rhs) result(res)
         real(wp), intent(in) :: lhs
         type(dual2), intent(in) :: rhs
         logical :: res

         res = (lhs >= rhs%x)

    end function ge_rd2

!******* end: (>=)
!---------------------

!******* BEGIN: (>)
!---------------------
    !-----------------------------------------
    ! compare two dual2 numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function gt_dd2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x > rhs%x)

    end function gt_dd2


    !-----------------------------------------
    ! compare a dual2 with an integer
    !-----------------------------------------
    elemental function gt_di2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x > rhs)

    end function gt_di2


    !-----------------------------------------
    ! compare a dual2 number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function gt_dr2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs
        real(wp), intent(in) :: rhs
        logical :: res

        res = (lhs%x > rhs)

    end function gt_dr2


    !-----------------------------------------
    ! compare a dual2 number with an integer
    !-----------------------------------------
    elemental function gt_id2(j,rhs) result(res)
        integer, intent(in) :: j
        type(dual2), intent(in) :: rhs
        logical :: res

        res = (j> rhs%x)

    end function gt_id2


    !-----------------------------------------
    ! compare a real with a dual2
    !-----------------------------------------
    elemental function gt_rd2(lhs, rhs) result(res)
         real(wp), intent(in) :: lhs
         type(dual2), intent(in) :: rhs
         logical :: res

         res = (lhs > rhs%x)

    end function gt_rd2

!******* end: (>)
!---------------------

!******* BEGIN: (/=)
!---------------------
    !-----------------------------------------
    ! compare two dual2 numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ne_dd2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x /= rhs%x)

    end function ne_dd2


    !-----------------------------------------
    ! compare a dual2 with an integer
    !-----------------------------------------
    elemental function ne_di2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x /= rhs)

    end function ne_di2


    !-----------------------------------------
    ! compare a dual2 number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ne_dr2(lhs, rhs) result(res)
        type(dual2), intent(in) :: lhs
        real(wp), intent(in) :: rhs
        logical :: res

        res = (lhs%x /= rhs)

    end function ne_dr2


    !-----------------------------------------
    ! compare a dual2 number with an integer
    !-----------------------------------------
    elemental function ne_id2(j,rhs) result(res)
        integer, intent(in) :: j
        type(dual2), intent(in) :: rhs
        logical :: res

        res = (j/= rhs%x)

    end function ne_id2


    !-----------------------------------------
    ! compare a real with a dual2
    !-----------------------------------------
    elemental function ne_rd2(lhs, rhs) result(res)
        real(wp), intent(in) :: lhs
        type(dual2), intent(in) :: rhs
        logical :: res

        res = (lhs /= rhs%x)

    end function ne_rd2

!******* end: (/=)
!---------------------

    !---------------------------------------------------
    ! Absolute value of dual2 numbers
    ! <res, dres> = abs(<u, du>) = <abs(u), du * sign(u)>
    !---------------------------------------------------
    elemental function abs_d2(u) result(res)
         type(dual2), intent(in) :: u
         type(dual2) :: res
         integer :: j

         if(u%x > 0) then
            res%x = u%x
            res%dx = u%dx
         else if (u%x < 0) then
            res%x = -u%x
            res%dx = -u%dx
         else
            res%x = 0.0d0
            do j = 1, dual_size2
                if (u%dx(j) .eq. 0.0d0) then
                    res%dx(j) = 0.0d0
                else
                    res%dx(j) = set_NaN()
                end if
            end do
         endif

    end function abs_d2


    !-----------------------------------------
    ! ACOS of dual2 numbers
    ! <res, dres> = acos(<u, du>) = <acos(u), -du / sqrt(1 - u^2)>
    !----------------------------------------
    elemental function acos_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = acos(u%x)
        if (u%x == 1.0d0 .or. u%x == -1.0d0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = -u%dx / sqrt(1.0d0 - u%x**2)
        end if

    end function acos_d2


    !-----------------------------------------
    ! ASIN of dual2 numbers
    ! <res, dres> = asin(<u, du>) = <asin(u), du / sqrt(1 - u^2)>
    !----------------------------------------
    elemental function asin_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = asin(u%x)
        if (u%x == 1.0d0 .or. u%x == -1.0d0) then
            res%dx = set_NaN()  ! Undefined derivative
        else
            res%dx = u%dx / sqrt(1.0d0 - u%x**2)
        end if

    end function asin_d2


    !-----------------------------------------
    ! ATAN of dual2 numbers
    ! <res, dres> = atan(<u, du>) = <atan(u), du / (1 + u^2)>
    !----------------------------------------
    elemental function atan_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = atan(u%x)
        res%dx = u%dx / (1.0d0 + u%x**2)

    end function atan_d2


    !-----------------------------------------
    ! ATAN2 of dual2 numbers
    ! <res, dres> = atan2(<u, du>, <v, dv>)
    !             = <atan2(u, v), v / (u^2 + v^2) * du - u / (u^2 + v^2) * dv>
    !----------------------------------------
    elemental function atan2_d2(u, v) result(res)
        type(dual2), intent(in) :: u, v
        type(dual2) :: res

        real(wp):: usq_plus_vsq

        res%x = atan2(u%x, v%x)

        usq_plus_vsq = u%x**2 + v%x**2
        res%dx = v%x / usq_plus_vsq * u%dx - u%x / usq_plus_vsq * v%dx

    end function atan2_d2


    !-----------------------------------------
    ! COS of dual2 numbers
    ! <res, dres> = cos(<u, du>) = <cos(u), -sin(u) * du>
    !----------------------------------------
    elemental function cos_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = cos(u%x)
        res%dx = -sin(u%x) * u%dx

    end function cos_d2


    !-----------------------------------------
    ! DOT PRODUCT two dual2 number vectors
    ! <res, dres> = <u, du> . <v, dv> = <u . v, u . dv + v . du>
    !-----------------------------------------
    function dot_product_dd2(u, v) result(res)
        type(dual2), intent(in) :: u(:), v(:)
        type(dual2) :: res

        integer :: j

        res%x = dot_product(u%x, v%x)
        do j = 1, dual_size2
            res%dx(j) = dot_product(u%x, v%dx(j)) + dot_product(v%x, u%dx(j))
        end do

    end function dot_product_dd2


    !-----------------------------------------
    ! EXPONENTIAL OF dual2 numbers
    ! <res, dres> = exp(<u, du>) = <exp(u), exp(u) * du>
    !-----------------------------------------
    elemental function exp_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        real(wp):: exp_x

        exp_x = exp(u%x)
        res%x = exp_x
        res%dx = u%dx * exp_x

    end function exp_d2


    !-----------------------------------------
    ! Convert dual2 to integer
    ! j = int(<u, du>) = int(u)
    !----------------------------------------
    elemental function int_d2(u) result(res)
         type(dual2), intent(in) :: u
         integer :: res

         res = int(u%x)

    end function int_d2


    !-----------------------------------------
    ! LOG OF dual2 numbers,defined for u%x>0 only
    ! the error control should be done in the original code
    ! in other words, if u%x<=0, it is not possible to obtain LOG.
    ! <res, dres> = log(<u, du>) = <log(u), du / u>
    !----------------------------------------
    elemental function log_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        real(wp):: inv

        inv = 1.0d0 / u%x
        res%x = log(u%x)
        res%dx = u%dx * inv

    end function log_d2


    !-----------------------------------------
    ! LOG10 OF dual2 numbers,defined for u%x>0 only
    ! the error control should be done in the original code
    ! in other words, if u%x<=0, it is not possible to obtain LOG.
    ! <res, dres> = log10(<u, du>) = <log10(u), du / (u * log(10))>
    ! LOG<u,up>=<LOG(u),up/u>
    !----------------------------------------
    elemental function log10_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        real(wp):: inv

        inv = 1.0d0 / (u%x * log(10.0d0))
        res%x = log10(u%x)
        res%dx = u%dx * inv

    end function log10_d2


    !-----------------------------------------
    ! MULTIPLY two dual2 number matrices
    ! <res, dres> = <u, du> . <v, dv> = <u . v, du . v + u . dv>
    !----------------------------------------
    function matmul_dd2(u,v) result(res)
        type(dual2), intent(in) :: u(:,:), v(:,:)
        type(dual2) :: res(size(u,1), size(v,2))

        integer :: j

        res%x = matmul(u%x, v%x)
        do j = 1, dual_size2
            res%dx(j) = matmul(u%dx(j), v%x) + matmul(u%x, v%dx(j))
        end do

    end function matmul_dd2


    !-----------------------------------------
    ! MULTIPLY a dual2 number matrix with a dual2 number
    ! vector
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    function matmul_dv2(u, v) result(res)
        type(dual2), intent(in) :: u(:,:), v(:)
        type(dual2) :: res(size(u,1))
        integer :: j

        res%x = matmul(u%x, v%x)
        do j = 1, dual_size2
            res%dx(j) = matmul(u%dx(j), v%x) + matmul(u%x, v%dx(j))
        end do

    end function matmul_dv2


    !-----------------------------------------
    ! MULTIPLY a dual2 vector with a  dual2 matrix
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    function matmul_vd2(u, v) result(res)
        type(dual2), intent(in) :: u(:), v(:,:)
        type(dual2) :: res(size(v, 2))
        integer::j

        res%x = matmul(u%x, v%x)
        do j = 1, dual_size2
            res%dx(j) = matmul(u%dx(j), v%x) + matmul(u%x, v%dx(j))
        end do

    end function matmul_vd2

    !-----------------------------------------
    ! Obtain the max of 2 to 5 dual2 numbers
    !----------------------------------------
    elemental function max_dd2(val1, val2, val3, val4,val5) result(res)
        type(dual2), intent(in) :: val1, val2
        type(dual2), intent(in), optional :: val3, val4,val5
        type(dual2) :: res

        if (val1%x > val2%x) then
            res = val1
        else
            res = val2
        endif
        if(present(val3))then
           if(res%x < val3%x) res = val3
        endif
        if(present(val4))then
           if(res%x < val4%x) res = val4
        endif
        if(present(val5))then
           if(res%x < val5%x) res = val5
        endif

    end function max_dd2


    !-----------------------------------------
    ! Obtain the max of a dual2 number and an integer
    !----------------------------------------
    elemental function max_di2(u, j) result(res)
        type(dual2), intent(in) :: u
        integer, intent(in) :: j
        type(dual2) :: res

        if (u%x > j) then
            res = u
        else
            res = j
        endif

    end function max_di2

    !-----------------------------------------
    ! Obtain the max of a dual2 number and a real number
    !----------------------------------------
    elemental function max_dr2(u, r) result(res)
        type(dual2), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual2) :: res

        if (u%x > r) then
            res = u
        else
            res%x = r
        endif

    end function max_dr2


    !---------------------------------------------------
    ! Obtain the max of a real and a dual2
    !---------------------------------------------------
     elemental function max_rd2(n, u) result(res)
        real(wp), intent(in) :: n
        type(dual2), intent(in) :: u
        type(dual2) :: res

        if (u%x > n) then
            res = u
        else
            res%x = n
        endif

    end function max_rd2


    !-----------------------------------------
    ! Obtain the max value of vector u
    !----------------------------------------
    function maxval_d2(u) result(res)
        type(dual2), intent(in) :: u(:)
        integer :: iloc(1)
        type(dual2) :: res

        iloc=maxloc(u%x)
        res=u(iloc(1))

    end function maxval_d2


    !-----------------------------------------
    ! Obtain the min of 2 to 4 dual2 numbers
    !----------------------------------------
    elemental function min_dd2(val1, val2, val3, val4) result(res)
        type(dual2), intent(in) :: val1, val2
        type(dual2), intent(in), optional :: val3, val4
        type(dual2) :: res

        if (val1%x < val2%x) then
            res = val1
        else
            res = val2
        endif
        if(present(val3))then
           if(res%x > val3%x) res = val3
        endif
        if(present(val4))then
           if(res%x > val4%x) res = val4
        endif

    end function min_dd2


    !-----------------------------------------
    ! Obtain the min of a dual2 and a double
    !----------------------------------------
    elemental function min_dr2(u, r) result(res)
        type(dual2), intent(in) :: u
        real(wp), intent(in) :: r
        type(dual2) :: res

        if (u%x < r) then
            res = u
        else
            res%x = r
        endif

    end function min_dr2


  !-----------------------------------------
    ! Obtain the min value of vector u
    !----------------------------------------
    function minval_d2(u) result(res)
        type(dual2), intent(in) :: u(:)
        integer :: iloc(1)
        type(dual2) :: res

        iloc=minloc(u%x)
        res=u(iloc(1))

    end function minval_d2


    !------------------------------------------------------
    !Returns the nearest integer to u%x, ELEMENTAL
    !------------------------------------------------------
    elemental function nint_d2(u) result(res)
        type(dual2), intent(in) :: u
        integer :: res

        res=nint(u%x)

    end function nint_d2


    !----------------------------------------------------------------
    ! SIGN(a,b) with two dual2 numbers as inputs,
    ! the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
    !----------------------------------------------------------------
    elemental function sign_dd2(val1, val2) result(res)
        type(dual2), intent(in) :: val1, val2
        type(dual2) :: res

        if (val2%x < 0.0d0) then
            res = -abs(val1)
        else
            res =  abs(val1)
        endif

     end function sign_dd2


    !----------------------------------------------------------------
    ! SIGN(a,b) with one real and one dual2 number as inputs,
    ! the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
    !----------------------------------------------------------------
    elemental function sign_rd2(val1, val2) result(res)
        real(wp), intent(in) :: val1
        type(dual2), intent(in) :: val2
        type(dual2) :: res

        if (val2%x < 0.0d0) then
            res%x = -abs(val1)
        else
            res%x = abs(val1)
        endif

     end function sign_rd2


    !-----------------------------------------
    ! SIN of dual2 numbers
    ! <res, dres> = sin(<u, du>) = <sin(u), cos(u) * du>
    !----------------------------------------
    elemental function sin_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = sin(u%x)
        res%dx = cos(u%x) * u%dx

    end function sin_d2


    !-----------------------------------------
    ! TAN of dual2 numbers
    ! <res, dres> = tan(<u, du>) = <tan(u), du / cos(u)^2>
    !----------------------------------------
    elemental function tan_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = tan(u%x)
        res%dx = u%dx / cos(u%x)**2

    end function tan_d2


    !-----------------------------------------
    ! SQRT of dual2 numbers
    ! <res, dres> = sqrt(<u, du>) = <sqrt(u), du / (2 * sqrt(u))>
    !----------------------------------------
    elemental function sqrt_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res
        integer :: j

        res%x = sqrt(u%x)

        if (res%x .ne. 0.0d0) then
            res%dx = 0.5 * u%dx / res%x
        else
            do j = 1, dual_size2
                if (u%dx(j) .eq. 0.0d0) then
                    res%dx(j) = 0.0d0
                else
                    res%dx(j) = set_NaN()
                end if
            end do
        end if

    end function sqrt_d2


    !-----------------------------------------
    ! Sum of a dual2 array
    !-----------------------------------------
    function sum_d2(u) result(res)
        type(dual2), intent(in) :: u(:)
        type(dual2) :: res
        integer :: j

        res%x = sum(u%x)
        do j = 1, dual_size2
            res%dx(j) = sum(u%dx(j))
        end do

    end function sum_d2


    !-----------------------------------------
    ! Find the location of the max value in an
    ! array of dual2 numbers
    !-----------------------------------------
    function maxloc_d2(array) result(ind)
        type(dual2), intent(in) :: array(:)
        integer :: ind(1)

        ind = maxloc(array%x)

    end function maxloc_d2



    !-----------------------------------------
    ! Hyperbolic functions: sinh, cosh, tanh
    ! and their inverses: asinh, acosh, atanh
    !-----------------------------------------
    !-----------------------------------------
    ! SINH OF dual2 numbers
    ! <res, dres> = sinh(<u, du>) = <sinh(u), cosh(u) * du>
    !-----------------------------------------
    elemental function sinh_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = sinh(u%x)
        res%dx = u%dx * cosh(u%x)

    end function sinh_d2

    !-----------------------------------------
    ! COSH OF dual2 numbers
    ! <res, dres> = cosh(<u, du>) = <cosh(u), sinh(u) * du>
    !-----------------------------------------
    elemental function cosh_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = cosh(u%x)
        res%dx = u%dx * sinh(u%x)

    end function cosh_d2

    !-----------------------------------------
    ! TANH OF dual2 numbers
    ! <res, dres> = tanh(<u, du>) = <tanh(u), 1.0d0/cosh(u)**2 * du>
    !-----------------------------------------
    elemental function tanh_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = tanh(u%x)
        res%dx = u%dx * 1.0d0/cosh(u%x)**2

    end function tanh_d2

    !-----------------------------------------
    ! ASINH OF dual2 numbers
    ! <res, dres> = asinh(<u, du>) = <asinh(u), 1/sqrt(u**2 + 1) * du>
    !-----------------------------------------
    elemental function asinh_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = asinh(u%x)
        res%dx = u%dx * 1.0d0/sqrt(u%x**2 + 1.0d0)

    end function asinh_d2

    !-----------------------------------------
    ! ACOSH OF dual2 numbers
    ! <res, dres> = acosh(<u, du>) = <acosh(u), 1/sqrt(u**2 - 1) * du>
    !-----------------------------------------
    elemental function acosh_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = acosh(u%x)
        if (u%x <= 1.0d0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = u%dx * 1.0d0/sqrt(u%x**2 - 1.0d0)
        end if

    end function acosh_d2

    !-----------------------------------------
    ! ATAHN OF dual2 numbers
    ! <res, dres> = atanh(<u, du>) = <atanh(u), 1/(1 - u**2) * du>
    !-----------------------------------------
    elemental function atanh_d2(u) result(res)
        type(dual2), intent(in) :: u
        type(dual2) :: res

        res%x = atanh(u%x)
        if (abs(u%x) >= 1.0d0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = u%dx * 1.0d0/(1.0d0 - u%x**2)
        end if

    end function atanh_d2

    !-----------------------------------------
    ! ERF OF dual2 numbers
    ! <res, dres> = ERF(<u, du>) = <ERF(u), 2/sqrt(Pj)*exp(-u**2) * du>
    !-----------------------------------------
    elemental function erf_d2(u) result(res)
      type(dual2), intent(in) :: u
      type(dual2) :: res

      res%x  = erf(u%x)
      res%dx = u%dx * 2.0/sqrt(PI) * exp(-u%x**2)

    end function erf_d2







! ! ******************************************************************************
! ! COMPLEX DUAL NUMBERS
! ! ******************************************************************************

! !*********Begin: functions/subroutines for overloading operators

! !******* Begin: (=)
! !---------------------

!     !-----------------------------------------
!     ! dual_complex2 = integer
!     ! <u, du> = <j 0>
!     !-----------------------------------------
!     elemental subroutine assign_dc_i2(u, j)
!          type(dual_complex2), intent(out) :: u
!          integer, intent(in) :: j

!          u%z = complex(real(j,kind=wp), 0.0d0)  ! This is faster than direct assignment
!          u%dz = complex(0.0d0, 0.0d0)

!     end subroutine assign_dc_i2


!     !-----------------------------------------
!     ! dual_complex2 = real(double)
!     ! <u, du> = <r, 0>
!     !-----------------------------------------
!     elemental subroutine assign_dc_r2(u, r)
!         type(dual_complex2), intent(out) :: u
!         real(wp), intent(in) :: r

!         u%z = complex(r, 0.0d0)
!         u%dz = complex(0.0d0, 0.0d0)

!     end subroutine assign_dc_r2

!     !-----------------------------------------
!     ! dual_complex2 = complex
!     ! <u, du> = <r, 0>
!     !-----------------------------------------
!     elemental subroutine assign_dc_c2(u, r)
!         type(dual_complex2), intent(out) :: u
!         complex(wp), intent(in) :: r

!         u%z = r
!         u%dz = complex(0.0d0, 0.0d0)

!     end subroutine assign_dc_c2


!     !-----------------------------------------
!     ! integer = dual        Is there a situation where complex --> integer makes sense?
!     ! j = <u, du>
!     !-----------------------------------------
!     ! elemental subroutine assign_id2(j,v)
!     !      type(dual_complex2), intent(in) :: v
!     !      integer, intent(out) :: j
!     !
!     !      j = int(v%z)
!     !
!     ! end subroutine assign_id2

!     !-----------------------------------------
!     ! complex = dual_complex2        Is there a situation where complex --> integer makes sense?
!     ! j = <u, du>
!     !-----------------------------------------
!     elemental subroutine assign_c_dc2(u, v)
!          type(dual_complex2), intent(in) :: v
!          complex(wp), intent(out) :: u

!          u = v%z

!     end subroutine assign_c_dc2


! !******* end: (=)
! !---------------------

! !******* Begin: (+)
! !---------------------

!     !-----------------------------------------
!     ! Unary positive
!     ! <res, dres> = +<u, du>
!     !-----------------------------------------
!     elemental function add_dc2(u) result(res)
!          type(dual_complex2), intent(in) :: u
!          type(dual_complex2) :: res

!          res = u  ! Faster than assigning component wise

!     end function add_dc2

!     !-----------------------------------------
!     ! dual_complex2 + dual_complex2
!     ! <res, dres> = <u, du> + <v, dv> = <u + v, du + dv>
!     !-----------------------------------------
!     elemental function add_dc_dc2(u, v) result(res)
!        type(dual_complex2), intent(in) :: u, v
!        type(dual_complex2) :: res

!        res%z = u%z + v%z
!        res%dz = u%dz + v%dz

!     end function add_dc_dc2

!     !-----------------------------------------
!     ! complex + dual_complex2
!     ! <res, dres> = <r, 0> + <v, dv> = <r + v, dv>
!     !-----------------------------------------
!     elemental function add_c_dc2(r, v) result(res)
!         complex(wp), intent(in) :: r
!         type(dual_complex2), intent(in) :: v
!         type(dual_complex2) :: res

!         res%z = r + v%z
!         res%dz = v%dz

!     end function add_c_dc2

!     !-----------------------------------------
!     ! dual_complex2 + complex
!     ! <res, dres> = <u, du> + <r, 0> = <u + r, du>
!     !-----------------------------------------
!     elemental function add_dc_c2(u, r) result(res)
!       type(dual_complex2), intent(in) :: u
!       complex(wp), intent(in) :: r
!       type(dual_complex2) :: res

!       res%z = r + u%z
!       res%dz = u%dz

!     end function add_dc_c2


!   !-----------------------------------------
!   ! dual_complex2 + integer
!   ! <res, dres> = <u, du> + j = <u + j du>
!   !-----------------------------------------
!   elemental function add_dc_i2(u, j) result(res)
!        type(dual_complex2), intent(in) :: u
!        integer, intent(in) :: j
!        type(dual_complex2) :: res

!        res%z = real(j,kind=wp) + u%z
!        res%dz = u%dz

!   end function add_dc_i2

!   !-----------------------------------------
!   ! integer + dual_complex2
!   ! <res, dres> = <u, du> + j = <u + j du>
!   !-----------------------------------------
!   elemental function add_i_dc2(j,u) result(res)
!        type(dual_complex2), intent(in) :: u
!        integer, intent(in) :: j
!        type(dual_complex2) :: res

!        res%z = real(j,kind=wp) + u%z
!        res%dz = u%dz

!   end function add_i_dc2

!   !-----------------------------------------
!   ! dual_complex2 + integer
!   ! <res, dres> = <u, du> + j = <u + j du>
!   !-----------------------------------------
!   elemental function add_dc_r2(u, r) result(res)
!        type(dual_complex2), intent(in) :: u
!        real(wp), intent(in) :: r
!        type(dual_complex2) :: res

!        res%z = r + u%z
!        res%dz = u%dz

!   end function add_dc_r2

!   !-----------------------------------------
!   ! integer + dual_complex2
!   ! <res, dres> = <u, du> + j = <u + j du>
!   !-----------------------------------------
!   elemental function add_r_dc2(r, u) result(res)
!        type(dual_complex2), intent(in) :: u
!        real(wp), intent(in) :: r
!        type(dual_complex2) :: res

!        res%z = r + u%z
!        res%dz = u%dz

!   end function add_r_dc2
! !
! ! !******* end: (+)
! ! !---------------------
! !
! !
! !******* Begin: (-)
! !---------------------

!   !-------------------------------------------------
!   ! negate a dual_complex2
!   ! <res, dres> = -<u, du>
!   !-------------------------------------------------
!   elemental function minus_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = -u%z
!       res%dz = -u%dz

!   end function minus_dc2


!   !-------------------------------------------------
!   ! dual_complex2 - dual_complex2
!   ! <res, dres> = <u, du> - <v, dv> = <u - v, du - dv>
!   !-------------------------------------------------
!   elemental function minus_dc_dc2(u, v) result(res)
!       type(dual_complex2), intent(in) :: u, v
!       type(dual_complex2) :: res

!       res%z = u%z - v%z
!       res%dz = u%dz - v%dz

!   end function minus_dc_dc2

!   !-------------------------------------------------
!   ! dual_complex2 - integer
!   ! <res, dres> = <u, du> - j = <u - j du>
!   !-------------------------------------------------
!   elemental function minus_dc_i2(u, j) result(res)
!       type(dual_complex2), intent(in) :: u
!       integer, intent(in) :: j
!       type(dual_complex2) :: res

!       res%z = u%z - real(j,kind=wp)
!       res%dz = u%dz

!   end function minus_dc_i2

!   !-------------------------------------------------
!   ! dual_complex2 - real
!   ! <res, dres> = <u, du> - r = <u - r, du>
!   !-------------------------------------------------
!   elemental function minus_dc_r2(u, r) result(res)
!       type(dual_complex2), intent(in) :: u
!       real(wp), intent(in) :: r
!       type(dual_complex2) :: res

!       res%z = u%z - r
!       res%dz = u%dz

!   end function minus_dc_r2


!   !-------------------------------------------------
!   ! dual_complex2 - complex
!   ! <res, dres> = <u, du> - r = <u - r, du>
!   !-------------------------------------------------
!   elemental function minus_dc_c2(u, c) result(res)
!       type(dual_complex2), intent(in) :: u
!       complex(wp), intent(in) :: c
!       type(dual_complex2) :: res

!       res%z = u%z - c
!       res%dz = u%dz

!   end function minus_dc_c2


!   !-------------------------------------------------
!   ! integer - dual_complex2
!   ! <res, dres> = j - <v, dv> = <i - v, -dv>
!   !-------------------------------------------------
!   elemental function minus_i_dc2(j,v) result(res)
!       integer, intent(in) :: j
!       type(dual_complex2), intent(in) :: v
!       type(dual_complex2) :: res

!       res%z = real(j,kind=wp) - v%z
!       res%dz = -v%dz

!   end function minus_i_dc2

!   !-------------------------------------------------
!   ! real - dual_complex2
!   ! <res, dres> = r - <v, dv> = <r - v, -dv>
!   !-------------------------------------------------
!   elemental function minus_r_dc2(r, v) result(res)
!       real(wp), intent(in) :: r
!       type(dual_complex2), intent(in) :: v
!       type(dual_complex2) :: res

!       res%z = r - v%z
!       res%dz = -v%dz

!   end function minus_r_dc2


!   !-------------------------------------------------
!   ! complex - dual_complex2
!   ! <res, dres> = c - <v, dv> = <c - v, -dv>
!   !-------------------------------------------------
!   elemental function minus_c_dc2(c, v) result(res)
!        complex(wp), intent(in) :: c
!        type(dual_complex2), intent(in) :: v
!        type(dual_complex2) :: res

!       res%z = c - v%z
!       res%dz = -v%dz

!   end function minus_c_dc2
! !
! ! !******* end: (-)
! ! !---------------------
! !
! !
! !******* BEGIN: (*)
! !---------------------

!   !----------------------------------------
!   ! dual_complex2 * dual_complex2
!   ! <res, dres> = <u, du> * <v, dv> = <u * v, u * dv + v * du>
!   !----------------------------------------
!   elemental function mult_dc_dc2(u, v) result(res)
!       type(dual_complex2), intent(in) :: u, v
!       type(dual_complex2) :: res

!       res%z = u%z * v%z
!       res%dz = u%z * v%dz + v%z * u%dz

!   end function mult_dc_dc2


!   !-----------------------------------------
!   ! dual_complex2 * integer
!   ! <res, dres> = <u, du> * j = <u * j du * i>
!   !-----------------------------------------
!   elemental function mult_dc_i2(u, j) result(res)
!       type(dual_complex2), intent(in) :: u
!       integer, intent(in) :: j
!       type(dual_complex2) :: res

!       real(wp):: r

!       r = real(j,kind=wp)
!       res%z = r * u%z
!       res%dz = r * u%dz

!   end function mult_dc_i2

!   !-----------------------------------------
!   ! dual_complex2 * real
!   ! <res, dres> = <u, du> * r = <u * r, du * r>
!   !----------------------------------------
!   elemental function mult_dc_r2(u, r) result(res)
!       type(dual_complex2), intent(in) :: u
!       real(wp), intent(in) :: r
!       type(dual_complex2) :: res

!       res%z = u%z * r
!       res%dz = u%dz * r

!   end function mult_dc_r2


!   !-----------------------------------------
!   ! integer * dual_complex2
!   ! <res, dres> = j * <v, dv> = <i * v, j * dv>
!   !-----------------------------------------
!   elemental function mult_i_dc2(j,v) result(res)
!       integer, intent(in) :: j
!       type(dual_complex2), intent(in) :: v
!       type(dual_complex2) :: res

!       real(wp):: r

!       r = real(j,kind=wp)
!       res%z = r * v%z
!       res%dz = r * v%dz

!   end function mult_i_dc2


!   !-----------------------------------------
!   ! double * dual_complex2
!   ! <res, dres> = r * <v, dv> = <r * v, r * dv>
!   !-----------------------------------------
!   elemental function mult_r_dc2(r, v) result(res)
!       real(wp), intent(in) :: r
!       type(dual_complex2), intent(in) :: v
!       type(dual_complex2) :: res

!       res%z = r * v%z
!       res%dz = r * v%dz

!   end function mult_r_dc2


!   !-----------------------------------------
!   ! complex * dual_complex2
!   ! <res, dres> = c * <v, dv> = <c * v, c * dv>
!   !-----------------------------------------
!   elemental function mult_c_dc2(c, v) result(res)
!       complex(wp), intent(in) :: c
!       type(dual_complex2), intent(in) :: v
!       type(dual_complex2) :: res

!       res%z = c * v%z
!       res%dz = c * v%dz

!   end function mult_c_dc2


!   !-----------------------------------------
!   ! dual_complex2 * complex
!   ! <res, dres> = c * <v, dv> = <c * v, c * dv>
!   !-----------------------------------------
!   elemental function mult_dc_c2(v, c) result(res)
!       complex(wp), intent(in) :: c
!       type(dual_complex2), intent(in) :: v
!       type(dual_complex2) :: res

!       res%z = c * v%z
!       res%dz = c * v%dz

!   end function mult_dc_c2
! !
! ! !******* end: (*)
! ! !---------------------
! !
! !
! ! !******* BEGIN: (/)
! ! !---------------------
! !
!   !-----------------------------------------
!   ! dual_complex2 / dual_complex2
!   ! <res, dres> = <u, du> / <v, dv> = <u / v, du / v - u * dv / v^2>
!   !-----------------------------------------
!   elemental function div_dc_dc2(u, v) result(res)
!       type(dual_complex2), intent(in) :: u, v
!       type(dual_complex2) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / v%z
!       res%z = u%z * inv
!       res%dz = (u%dz - res%z * v%dz) * inv

!   end function div_dc_dc2
! !
! !
!   !-----------------------------------------
!   ! dual_complex2 / integer
!   ! <res, dres> = <u, du> / j = <u / j du / i>
!   !-----------------------------------------
!   elemental function div_dc_i2(u, j) result(res)
!       type(dual_complex2), intent(in) :: u
!       integer, intent(in) :: j
!       type(dual_complex2) :: res

!       real(wp):: inv

!       inv = 1.0d0 / real(j,kind=wp)
!       res%z = u%z * inv
!       res%dz = u%dz * inv

!   end function div_dc_i2
! !
! !
!   !-----------------------------------------
!   ! dual_complex2 / double
!   ! <res, dres> = <u, du> / r = <u / r, du / r>
!   !----------------------------------------
!   elemental function div_dc_r2(u, r) result(res)
!       type(dual_complex2), intent(in) :: u
!       real(wp), intent(in) :: r
!       type(dual_complex2):: res

!       real(wp):: inv

!       inv = 1.0d0 / r
!       res%z = u%z * inv
!       res%dz = u%dz * inv

!   end function div_dc_r2
! !
! !
!   !-----------------------------------------
!   ! integer / dual_complex2
!   ! <res, dres> = j / <v, dv> = <i / v, -i / v^2 * du>
!   !-----------------------------------------
!   elemental function div_i_dc2(j,v) result(res)
!       integer, intent(in) :: j
!       type(dual_complex2), intent(in) :: v
!       type(dual_complex2) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / v%z
!       res%z = real(j,kind=wp) * inv
!       res%dz = -res%z * inv * v%dz

!   end function div_i_dc2


!   !-----------------------------------------
!   ! double / dual_complex2
!   ! <res, dres> = r / <u, du> = <r / u, -r / u^2 * du>
!   !-----------------------------------------
!   elemental function div_r_dc2(r, v) result(res)
!       real(wp), intent(in) :: r
!       type(dual_complex2), intent(in) :: v
!       type(dual_complex2) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / v%z
!       res%z = r * inv
!       res%dz = -res%z * inv * v%dz

!   end function div_r_dc2


!   !-----------------------------------------
!   ! complex / dual_complex2
!   ! <res, dres> = c / <u, du> = <c / u, -r / u^2 * du>
!   !-----------------------------------------
!   elemental function div_c_dc2(c, v) result(res)
!       complex(wp), intent(in) :: c
!       type(dual_complex2), intent(in) :: v
!       type(dual_complex2) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / v%z
!       res%z = c * inv
!       res%dz = -res%z * inv * v%dz

!   end function div_c_dc2


!   !-----------------------------------------
!   ! dual_complex2 / complex
!   ! <res, dres> = c / <u, du> = <c / u, -c / u^2 * du>
!   !-----------------------------------------
!   elemental function div_dc_c2(v, c) result(res)
!       complex(wp), intent(in) :: c
!       type(dual_complex2), intent(in) :: v
!       type(dual_complex2) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / c
!       res%z = v%z * inv
!       res%dz = v%dz * inv

!   end function div_dc_c2

! !******* end: (/)
! !---------------------

! !******* BEGIN: (**)
! !---------------------

!   !-----------------------------------------
!   ! power(dual_complex2, integer)
!   ! <res, dres> = <u, du> ^ j = <u ^ j j * u ^ (j- 1) * du>
!   !-----------------------------------------
!   elemental function pow_dc_i2(u, j) result(res)
!       type(dual_complex2), intent(in) :: u
!       integer, intent(in) :: j
!       type(dual_complex2) :: res

!       complex(wp):: pow_x

!       pow_x = u%z ** (j- 1)
!       res%z = u%z * pow_x
!       res%dz = real(j,kind=wp) * pow_x * u%dz

!   end function pow_dc_i2

!   !-----------------------------------------
!   ! power(dual_complex2, double)
!   ! <res, dres> = <u, du> ^ r = <u ^ r, r * u ^ (r - 1) * du>
!   !-----------------------------------------
!   elemental function pow_dc_r2(u, r) result(res)
!       type(dual_complex2), intent(in) :: u
!       real(wp), intent(in) :: r
!       type(dual_complex2) :: res

!       complex(wp):: pow_x

!       pow_x = u%z ** (r - 1.0d0)
!       res%z = u%z * pow_x
!       res%dz = r * pow_x * u%dz

!   end function pow_dc_r2

!   !-----------------------------------------
!   ! power(dual_complex2, complex)
!   ! <res, dres> = <u, du> ^ c = <u ^ c, c * u ^ (c - 1) * du>
!   !-----------------------------------------
!   elemental function pow_dc_c2(u, c) result(res)
!       type(dual_complex2), intent(in) :: u
!       complex(wp), intent(in) :: c
!       type(dual_complex2) :: res

!       complex(wp) :: pow_x

!       pow_x = u%z ** (c - 1.0d0)
!       res%z = u%z * pow_x
!       res%dz = c * pow_x * u%dz

!   end function pow_dc_c2

!   !-----------------------------------------
!   ! POWER dual_complex2 numbers to a dual_complex2 power
!   ! <res, dres> = <u, du> ^ <v, dv>
!   !     = <u ^ v, u ^ v * (v / u * du + Log(u) * dv)>
!   !-----------------------------------------
!   elemental function pow_dc_dc2(u, v) result(res)
!       type(dual_complex2), intent(in) :: u, v
!       type(dual_complex2) :: res

!       res%z = u%z ** v%z
!       res%dz = res%z * (v%z / u%z * u%dz + log(u%z) * v%dz)

!   end function pow_dc_dc2

!   !-----------------------------------------
!   ! POWER integer numbers to a dual_complex2 power
!   ! <res, dres> = j ^ <u, du>
!   !     = <i ^ u, j ^ u * Log(j) * du>
!   !-----------------------------------------
!   elemental function pow_i_dc2(j,u) result(res)
!       integer, intent(in)   :: j
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = j ** u%z
!       res%dz = res%z * log(real(j,kind=wp)) * u%dz

!   end function pow_i_dc2

!   !-----------------------------------------
!   ! POWER real numbers to a dual_complex2 power
!   ! <res, dres> = r ^ <u, du>
!   !     = <r ^ u, r ^ u * Log(r) * du)>
!   !-----------------------------------------
!   elemental function pow_r_dc2(r, u) result(res)
!       real(wp), intent(in) :: r
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = r ** u%z
!       res%dz = res%z * log(r) * u%dz

!   end function pow_r_dc2

!   !-----------------------------------------
!   ! POWER complex number to a dual_complex2 power
!   ! <res, dres> = c ^ <u, du>
!   !     = <c ^ u, c ^ u * Log(c) * du)>
!   !-----------------------------------------
!   elemental function pow_c_dc2(c, u) result(res)
!       complex(wp), intent(in) :: c
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = c ** u%z
!       res%dz = res%z * log(c) * u%dz

!   end function pow_c_dc2

! !******* end: (**)
! !---------------------


! !******* BEGIN: (==)
! !---------------------
!   !-----------------------------------------
!   ! compare two dual_complex2 numbers,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_dc_dc2(lhs, rhs) result(res)
!        type(dual_complex2), intent(in) :: lhs, rhs
!        logical :: res

!        res = (lhs%z == rhs%z)

!   end function eq_dc_dc2


!   !-----------------------------------------
!   ! compare a dual_complex2 with an integer,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_dc_i2(lhs, rhs) result(res)
!        type(dual_complex2), intent(in) :: lhs
!        integer, intent(in) :: rhs
!        logical :: res

!        res = (lhs%z == complex(real(rhs), 0.0d0))

!   end function eq_dc_i2


!   !-----------------------------------------
!   ! compare a dual_complex2 number with a real number,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_dc_r2(lhs, rhs) result(res)
!       type(dual_complex2), intent(in) :: lhs
!       real(wp), intent(in) :: rhs
!       logical::res

!       res = (lhs%z == complex(rhs, 0.0d0))

!   end function eq_dc_r2


!   !-----------------------------------------
!   ! compare an integer with a dual_complex2,
!   ! simply compare the functional value.
!   !----------------------------------------
!   elemental function eq_i_dc2(lhs, rhs) result(res)
!        integer, intent(in) :: lhs
!        type(dual_complex2), intent(in) :: rhs
!        logical :: res

!        res = (complex(real(lhs), 0.0d0) == rhs%z)

!   end function eq_i_dc2


!   !-----------------------------------------
!   ! compare a real with a dual_complex2,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_r_dc2(lhs, rhs) result(res)
!        real(wp), intent(in) :: lhs
!        type(dual_complex2), intent(in) :: rhs
!        logical :: res

!        res = (complex(lhs, 0.0d0) == rhs%z)

!   end function eq_r_dc2


!   !-----------------------------------------
!   ! compare a complex number with a dual_complex2,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_c_dc2(lhs, rhs) result(res)
!        real(wp), intent(in) :: lhs
!        type(dual_complex2), intent(in) :: rhs
!        logical :: res

!        res = (lhs == rhs%z)

!   end function eq_c_dc2


!   !-----------------------------------------
!   ! compare a dual_complex2 with a comlex,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function eq_dc_c2(lhs, rhs) result(res)
!        type(dual_complex2), intent(in) :: lhs
!        complex(wp), intent(in) :: rhs
!        logical :: res

!        res = (lhs%z == rhs)

!   end function eq_dc_c2

! !******* end: (==)
! !---------------------


! !******* can't compare complex numbers
! ! (<=)
! ! (<)
! ! (>=)
! ! (>)
! !
! !******* BEGIN: (/=)
! !---------------------
!   !-----------------------------------------
!   ! compare two dual_complex2 numbers, simply compare
!   ! the functional value.
!   !-----------------------------------------
!   elemental function ne_dc_dc2(lhs, rhs) result(res)
!       type(dual_complex2), intent(in) :: lhs, rhs
!       logical :: res

!       res = (lhs%z /= rhs%z)

!   end function ne_dc_dc2

!   !-----------------------------------------
!   ! compare a dual_complex2 number with a complex number,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function ne_dc_c2(lhs, rhs) result(res)
!       type(dual_complex2), intent(in) :: lhs
!       complex(wp), intent(in) :: rhs
!       logical :: res

!       res = (lhs%z /= rhs)

!   end function ne_dc_c2

!   !-----------------------------------------
!   ! compare a complex with a dual_complex2
!   !-----------------------------------------
!   elemental function ne_c_dc2(lhs, rhs) result(res)
!       complex(wp), intent(in) :: lhs
!       type(dual_complex2), intent(in) :: rhs
!       logical :: res

!       res = (lhs /= rhs%z)

!   end function ne_c_dc2


!   !-----------------------------------------
!   ! compare a dual_complex2 with an integer,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function ne_dc_i2(lhs, rhs) result(res)
!        type(dual_complex2), intent(in) :: lhs
!        integer, intent(in) :: rhs
!        logical :: res

!        res = (lhs%z /= complex(real(rhs), 0.0d0))

!   end function ne_dc_i2


!   !-----------------------------------------
!   ! compare a dual_complex2 number with a real number,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function ne_dc_r2(lhs, rhs) result(res)
!       type(dual_complex2), intent(in) :: lhs
!       real(wp), intent(in) :: rhs
!       logical::res

!       res = (lhs%z /= complex(rhs, 0.0d0))

!   end function ne_dc_r2


!   !-----------------------------------------
!   ! compare an integer with a dual_complex2,
!   ! simply compare the functional value.
!   !----------------------------------------
!   elemental function ne_i_dc2(lhs, rhs) result(res)
!        integer, intent(in) :: lhs
!        type(dual_complex2), intent(in) :: rhs
!        logical :: res

!        res = (complex(real(lhs), 0.0d0) /= rhs%z)

!   end function ne_i_dc2


!   !-----------------------------------------
!   ! compare a real with a dual_complex2,
!   ! simply compare the functional value.
!   !-----------------------------------------
!   elemental function ne_r_dc2(lhs, rhs) result(res)
!        real(wp), intent(in) :: lhs
!        type(dual_complex2), intent(in) :: rhs
!        logical :: res

!        res = (complex(lhs, 0.0d0) /= rhs%z)

!   end function ne_r_dc2

! !******* end: (/=)
! !---------------------
! !
!   !---------------------------------------------------
!   ! Absolute value of dual_complex2 numbers
!   ! <res, dres> = abs(<u, du>) = <abs(u), du * u / |u|>
!   !---------------------------------------------------
!   elemental function abs_dc2(u) result(res)
!        type(dual_complex2), intent(in) :: u
!        type(dual_complex2) :: res

!        res%z = complex(abs(u%z), 0.)
!        res%dz = u%dz * u%z / abs(u%z)

!   end function abs_dc2
! !
! !
!   !-----------------------------------------
!   ! ACOS of dual_complex2 numbers
!   ! <res, dres> = acos(<u, du>) = <acos(u), -du / sqrt(1 - u^2)>
!   !----------------------------------------
!   elemental function acos_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = acos(u%z)
!       if (u%z == complex(1.0d0, 0.0d0) .or. u%z == complex(-1.0d0, 0.0d0)) then
!           res%dz = set_Nan()  ! Undefined derivative
!       else
!           res%dz = -u%dz / sqrt(1.0d0 - u%z**2)
!       end if

!   end function acos_dc2


!   !-----------------------------------------
!   ! ASIN of dual_complex2 numbers
!   ! <res, dres> = asin(<u, du>) = <asin(u), du / sqrt(1 - u^2)>
!   !----------------------------------------
!   elemental function asin_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = asin(u%z)
!       if (u%z == complex(1.0d0, 0.0d0) .or. u%z == complex(-1.0d0, 0.0d0)) then
!           res%dz = set_NaN()  ! Undefined derivative
!       else
!           res%dz = u%dz / sqrt(1.0d0 - u%z**2)
!       end if

!   end function asin_dc2


!   !-----------------------------------------
!   ! ATAN of dual_complex2 numbers
!   ! <res, dres> = atan(<u, du>) = <atan(u), du / (1 + u^2)>
!   !----------------------------------------
!   elemental function atan_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = atan(u%z)
!       res%dz = u%dz / (1.0d0 + u%z**2)

!   end function atan_dc2


!   !-----------------------------------------
!   ! ATAN2 of dual_complex2 numbers
!   ! <res, dres> = atan2(<u, du>, <v, dv>)
!   !             = <atan2(u, v), v / (u^2 + v^2) * du - u / (u^2 + v^2) * dv>
!   !----------------------------------------
!   ! elemental function atan2_dc2(u, v) result(res)
!   !     type(dual_complex2), intent(in) :: u, v
!   !     type(dual_complex2) :: res
!   !
!   !     complex(wp) :: usq_plus_vsq
!   !
!   !     res%z = atan2(u%z, v%z)
!   !
!   !     usq_plus_vsq = u%z**2 + v%z**2
!   !     res%dz = v%z / usq_plus_vsq * u%dz - u%z / usq_plus_vsq * v%dz
!   !
!   ! end function atan2_dc2
! !
! !
!   !-----------------------------------------
!   ! COS of dual_complex2 numbers
!   ! <res, dres> = cos(<u, du>) = <cos(u), -sin(u) * du>
!   !----------------------------------------
!   elemental function cos_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = cos(u%z)
!       res%dz = -sin(u%z) * u%dz

!   end function cos_dc2
! !
! !
! !   !-----------------------------------------
! !   ! DOT PRODUCT two dual_complex2 number vectors
! !   ! <res, dres> = <u, du> . <v, dv> = <u . v, u . dv + v . du>
! !   !-----------------------------------------
! !   function dot_product_dc_d2(u, v) result(res)
! !       type(dual_complex2), intent(in) :: u(:), v(:)
! !       type(dual_complex2) :: res
! !
! !       integer :: j
! !
! !       res%z = dot_product(u%z, v%z)
! !       do j = 1, dual_size2
! !           res%dz(j) = dot_product(u%z, v%dz(j)) + dot_product(v%z, u%dz(j))
! !       end do
! !
! !   end function dot_product_dc_d
! !
! !
!   !-----------------------------------------
!   ! EXPONENTIAL OF dual_complex2 numbers
!   ! <res, dres> = exp(<u, du>) = <exp(u), exp(u) * du>
!   !-----------------------------------------
!   elemental function exp_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       complex(wp) :: exp_x

!       exp_x = exp(u%z)
!       res%z = exp_x
!       res%dz = u%dz * exp_x

!   end function exp_dc2
! !
! !
! !   !-----------------------------------------
! !   ! Convert dual_complex2 to integer
! !   ! j = int(<u, du>) = int(u)
! !   !----------------------------------------
! !   elemental function int_dc2(u) result(res)
! !        type(dual_complex2), intent(in) :: u
! !        integer :: res
! !
! !        res = int(u%z)
! !
! !   end function int_dc2
! !
! !
!   !-----------------------------------------
!   ! LOG OF dual_complex2 numbers,defined for u%z>0 only
!   ! the error control should be done in the original code
!   ! in other words, if u%z<=0, it is not possible to obtain LOG.
!   ! <res, dres> = log(<u, du>) = <log(u), du / u>
!   !----------------------------------------
!   elemental function log_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / u%z
!       res%z = log(u%z)
!       res%dz = u%dz * inv

!   end function log_dc2


!   !-----------------------------------------
!   ! LOG10 OF dual_complex2 numbers,defined for u%z>0 only
!   ! the error control should be done in the original code
!   ! in other words, if u%z<=0, it is not possible to obtain LOG.
!   ! <res, dres> = log10(<u, du>) = <log10(u), du / (u * log(10))>
!   ! LOG<u,up>=<LOG(u),up/u>
!   !----------------------------------------
!   elemental function log10_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       complex(wp) :: inv

!       inv = 1.0d0 / (u%z * log(10.0d0))
!       res%z = log(u%z) / log(10.0d0)
!       res%dz = u%dz * inv

!   end function log10_dc2
! !
! !
! !   !-----------------------------------------
! !   ! MULTIPLY two dual_complex2 number matrices
! !   ! <res, dres> = <u, du> . <v, dv> = <u . v, du . v + u . dv>
! !   !----------------------------------------
! !   function matmul_dc_d2(u,v) result(res)
! !       type(dual_complex2), intent(in) :: u(:,:), v(:,:)
! !       type(dual_complex2) :: res(size(u,1), size(v,2))
! !
! !       integer :: j
! !
! !       res%z = matmul(u%z, v%z)
! !       do j = 1, dual_size2
! !           res%dz(j) = matmul(u%dz(j), v%z) + matmul(u%z, v%dz(j))
! !       end do
! !
! !   end function matmul_dc_d
! !
! !
! !   !-----------------------------------------
! !   ! MULTIPLY a dual_complex2 number matrix with a dual_complex2 number
! !   ! vector
! !   !
! !   ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
! !   !----------------------------------------
! !   function matmul_dc_v(u, v) result(res)
! !       type(dual_complex2), intent(in) :: u(:,:), v(:)
! !       type(dual_complex2) :: res(size(u,1))
! !       integer :: j
! !
! !       res%z = matmul(u%z, v%z)
! !       do j = 1, dual_size2
! !           res%dz(j) = matmul(u%dz(j), v%z) + matmul(u%z, v%dz(j))
! !       end do
! !
! !   end function matmul_dc_v
! !
! !
! !   !-----------------------------------------
! !   ! MULTIPLY a dual_complex2 vector with a  dual_complex2 matrix
! !   !
! !   ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
! !   !----------------------------------------
! !   function matmul_vd2(u, v) result(res)
! !       type(dual_complex2), intent(in) :: u(:), v(:,:)
! !       type(dual_complex2) :: res(size(v, 2))
! !       integer::i
! !
! !       res%z = matmul(u%z, v%z)
! !       do j = 1, dual_size2
! !           res%dz(j) = matmul(u%dz(j), v%z) + matmul(u%z, v%dz(j))
! !       end do
! !
! !   end function matmul_vd2
! !
! !   !-----------------------------------------
! !   ! Obtain the max of 2 to 5 dual_complex2 numbers
! !   !----------------------------------------
! !   elemental function max_dc_d2(val1, val2, val3, val4,val5) result(res)
! !       type(dual_complex2), intent(in) :: val1, val2
! !       type(dual_complex2), intent(in), optional :: val3, val4,val5
! !       type(dual_complex2) :: res
! !
! !       if (val1%z > val2%z) then
! !           res = val1
! !       else
! !           res = val2
! !       endif
! !       if(present(val3))then
! !          if(res%z < val3%z) res = val3
! !       endif
! !       if(present(val4))then
! !          if(res%z < val4%z) res = val4
! !       endif
! !       if(present(val5))then
! !          if(res%z < val5%z) res = val5
! !       endif
! !
! !   end function max_dc_d
! !
! !
! !   !-----------------------------------------
! !   ! Obtain the max of a dual_complex2 number and an integer
! !   !----------------------------------------
! !   elemental function max_dc_i2(u, j) result(res)
! !       type(dual_complex2), intent(in) :: u
! !       integer, intent(in) :: j
! !       type(dual_complex2) :: res
! !
! !       if (u%z > j) then
! !           res = u
! !       else
! !           res = i
! !       endif
! !
! !   end function max_dc_i2
! !
! !   !-----------------------------------------
! !   ! Obtain the max of a dual_complex2 number and a real number
! !   !----------------------------------------
! !   elemental function max_dc_r2(u, r) result(res)
! !       type(dual_complex2), intent(in) :: u
! !       real(wp), intent(in) :: r
! !       type(dual_complex2) :: res
! !
! !       if (u%z > r) then
! !           res = u
! !       else
! !           res = r
! !       endif
! !
! !   end function max_dc_r2
! !
! !
! !   !---------------------------------------------------
! !   ! Obtain the max of a real and a dual_complex2
! !   !---------------------------------------------------
! !    elemental function max_rd2(n, u) result(res)
! !       real(wp), intent(in) :: n
! !       type(dual_complex2), intent(in) :: u
! !       type(dual_complex2) :: res
! !
! !       if (u%z > n) then
! !           res = u
! !       else
! !           res = n
! !       endif
! !
! !   end function max_rd2
! !
! !
! !   !-----------------------------------------
! !   ! Obtain the max value of vector u
! !   !----------------------------------------
! !   function maxval_dc(u) result(res)
! !       type(dual_complex2), intent(in) :: u(:)
! !       integer :: iloc(1)
! !       type(dual_complex2) :: res
! !
! !       iloc=maxloc(u%z)
! !       res=u(iloc(1))
! !
! !   end function maxval_dc
! !
! !
! !   !-----------------------------------------
! !   ! Obtain the min of 2 to 4 dual_complex2 numbers
! !   !----------------------------------------
! !   elemental function min_dc_d2(val1, val2, val3, val4) result(res)
! !       type(dual_complex2), intent(in) :: val1, val2
! !       type(dual_complex2), intent(in), optional :: val3, val4
! !       type(dual_complex2) :: res
! !
! !       if (val1%z < val2%z) then
! !           res = val1
! !       else
! !           res = val2
! !       endif
! !       if(present(val3))then
! !          if(res%z > val3%z) res = val3
! !       endif
! !       if(present(val4))then
! !          if(res%z > val4%z) res = val4
! !       endif
! !
! !   end function min_dc_d
! !
! !
! !   !-----------------------------------------
! !   ! Obtain the min of a dual_complex2 and a double
! !   !----------------------------------------
! !   elemental function min_dc_r2(u, r) result(res)
! !       type(dual_complex2), intent(in) :: u
! !       real(wp), intent(in) :: r
! !       type(dual_complex2) :: res
! !
! !       if (u%z < r) then
! !           res = u
! !       else
! !           res = r
! !       endif
! !
! !   end function min_dc_r2
! !
! !
! ! !-----------------------------------------
! !   ! Obtain the min value of vector u
! !   !----------------------------------------
! !   function minval_dc2(u) result(res)
! !       type(dual_complex2), intent(in) :: u(:)
! !       integer :: iloc(1)
! !       type(dual_complex2) :: res
! !
! !       iloc=minloc(u%z)
! !       res=u(iloc(1))
! !
! !   end function minval_dc2
! !
! !
! !   !------------------------------------------------------
! !   !Returns the nearest integer to u%z, ELEMENTAL
! !   !------------------------------------------------------
! !   elemental function nint_dc2(u) result(res)
! !       type(dual_complex2), intent(in) :: u
! !       integer :: res
! !
! !       res=nint(u%z)
! !
! !   end function nint_dc2
! !
! !
! !   !----------------------------------------------------------------
! !   ! SIGN(a,b) with two dual_complex2 numbers as inputs,
! !   ! the result will be |a| if b%z>=0, -|a| if b%z<0,ELEMENTAL
! !   !----------------------------------------------------------------
! !   elemental function sign_dc_d2(val1, val2) result(res)
! !       type(dual_complex2), intent(in) :: val1, val2
! !       type(dual_complex2) :: res
! !
! !       if (val2%z < 0.0d0) then
! !           res = -abs(val1)
! !       else
! !           res =  abs(val1)
! !       endif
! !
! !    end function sign_dc_d2
! !
! !
! !   !----------------------------------------------------------------
! !   ! SIGN(a,b) with one real and one dual_complex2 number as inputs,
! !   ! the result will be |a| if b%z>=0, -|a| if b%z<0,ELEMENTAL
! !   !----------------------------------------------------------------
! !   elemental function sign_rd2(val1, val2) result(res)
! !       real(wp), intent(in) :: val1
! !       type(dual_complex2), intent(in) :: val2
! !       type(dual_complex2) :: res
! !
! !       if (val2%z < 0.0d0) then
! !           res = -abs(val1)
! !       else
! !           res = abs(val1)
! !       endif
! !
! !    end function sign_rd2
! !
! !
! !   !-----------------------------------------
! !   ! SIN of dual_complex2 numbers
! !   ! <res, dres> = sin(<u, du>) = <sin(u), cos(u) * du>
! !   !----------------------------------------
!   elemental function sin_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = sin(u%z)
!       res%dz = cos(u%z) * u%dz

!   end function sin_dc2
! !
! !
!   !-----------------------------------------
!   ! TAN of dual_complex2 numbers
!   ! <res, dres> = tan(<u, du>) = <tan(u), du / cos(u)^2>
!   !----------------------------------------
!   elemental function tan_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = tan(u%z)
!       res%dz = u%dz / cos(u%z)**2

!   end function tan_dc2
! !
! !
!   !-----------------------------------------
!   ! SQRT of dual_complex2 numbers
!   ! <res, dres> = sqrt(<u, du>) = <sqrt(u), du / (2 * sqrt(u))>
!   !----------------------------------------
!   elemental function sqrt_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res
!       integer :: j

!       res%z = sqrt(u%z)

!       if (res%z /= complex(0.0d0, 0.0d0)) then
!           res%dz = 0.5 * u%dz / res%z
!       else
!           do j = 1, dual_size2
!               if (u%dz(j) == complex(0.0d0, 0.0d0)) then
!                   res%dz(j) = complex(0.0d0, 0.0d0)
!               else
!                   res%dz(j) = set_NaN()
!               end if
!           end do
!       end if

!   end function sqrt_dc2


!   !-----------------------------------------
!   ! Sum of a dual_complex2 array
!   !-----------------------------------------
!   function sum_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u(:)
!       type(dual_complex2) :: res
!       integer :: j

!       res%z = sum(u%z)
!       do j = 1, dual_size2
!           res%dz(j) = sum(u%dz(j))
!       end do

!   end function sum_dc2
! !
! !
! !   !-----------------------------------------
! !   ! Find the location of the max value in an
! !   ! array of dual_complex2 numbers
! !   !-----------------------------------------
! !   function maxloc_dc2(array) result(ind)
! !       type(dual_complex2), intent(in) :: array(:)
! !       integer :: ind(1)
! !
! !       ind = maxloc(array%z)
! !
! !   end function maxloc_dc2
! !
! !
! !   elemental function set_NaN() result(res)
! !       real(wp):: res
! !
! !       res = sqrt(negative_one)
! !
! !   end function set_NaN
! !
! !
!   !-----------------------------------------
!   ! Hyperbolic functions: sinh, cosh, tanh
!   ! and their inverses: asinh, acosh, atanh
!   !-----------------------------------------
!   !-----------------------------------------
!   ! SINH OF dual_complex2 numbers
!   ! <res, dres> = sinh(<u, du>) = <sinh(u), cosh(u) * du>
!   !-----------------------------------------
!   elemental function sinh_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = sinh(u%z)
!       res%dz = u%dz * cosh(u%z)

!   end function sinh_dc2

!   !-----------------------------------------
!   ! COSH OF dual_complex2 numbers
!   ! <res, dres> = cosh(<u, du>) = <cosh(u), sinh(u) * du>
!   !-----------------------------------------
!   elemental function cosh_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = cosh(u%z)
!       res%dz = u%dz * sinh(u%z)

!   end function cosh_dc2

!   !-----------------------------------------
!   ! TANH OF dual_complex2 numbers
!   ! <res, dres> = tanh(<u, du>) = <tanh(u), 1.0d0/cosh(u)**2 * du>
!   !-----------------------------------------
!   elemental function tanh_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = tanh(u%z)
!       res%dz = u%dz * 1.0d0/cosh(u%z)**2

!   end function tanh_dc2

!   !-----------------------------------------
!   ! ASINH OF dual_complex2 numbers
!   ! <res, dres> = asinh(<u, du>) = <asinh(u), 1/sqrt(u**2 + 1) * du>
!   !-----------------------------------------
!   elemental function asinh_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = asinh(u%z)
!       res%dz = u%dz * 1.0d0/sqrt(u%z**2 + 1.0d0)

!   end function asinh_dc2

!   !-----------------------------------------
!   ! ACOSH OF dual_complex2 numbers
!   ! <res, dres> = acosh(<u, du>) = <acosh(u), 1/sqrt(u**2 - 1) * du>
!   !-----------------------------------------
!   elemental function acosh_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = acosh(u%z)
!       if (u%z == complex(1.0d0, 0.0d0)) then
!           res%dz = set_Nan()  ! Undefined derivative (∞)
!       else
!           res%dz = u%dz * 1.0d0/sqrt(u%z**2 - 1.0d0)
!       end if


!   end function acosh_dc2

!   !-----------------------------------------
!   ! ATAHN OF dual_complex2 numbers
!   ! <res, dres> = atanh(<u, du>) = <atanh(u), 1/(1 - u**2) * du>
!   !-----------------------------------------
!   elemental function atanh_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z = atanh(u%z)
!       if (u%z == complex(1.0d0, 0.0d0) .or. u%z == complex(-1.0d0, 0.0d0)) then
!           res%dz = set_Nan()  ! Undefined derivative
!       else
!           res%dz = u%dz * 1.0d0/(1.0d0 - u%z**2)
!       end if

!   end function atanh_dc2


!     ! -----------------------------------------
!     ! ERF OF dual_complex2 numbers
!     ! <res, dres> = ERF(<u, du>) = <ERF(u), 2/sqrt(PI)*exp(-u**2) * du>
!     ! -----------------------------------------
!     elemental function erf_dc2(u) result(res)
!       type(dual_complex2), intent(in) :: u
!       type(dual_complex2) :: res

!       res%z  = erf(u%z)
!       res%dz = u%dz * 2.0/sqrt(PI) * exp(-u%z**2)

!     end function erf_dc2

end module dnad