! This module has been sourced from https://github.com/equipez/infnan on 25/10/2022. At the time of access the latest commit 
! [bc9213f] was on the 16/09/2021. The module has a GNU Lesser General Public License v3.0. 
! Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed. 
! No changes have been made to the module however the content of the file consts.f90 has been included in this file whilst 
! preserving it as a seperate module.

module consts_mod

implicit none
private
public :: SP, DP

integer, parameter :: SP = kind(0.0)
integer, parameter :: DP = kind(0.0D0)

end module consts_mod

module infnan_mod
    
! 1. This module provides functions for checking Inf/NaN. They aim to work even when compilers are
! invoked with aggressive optimization flags, such as `gfortran -Ofast`.
!
! 2. There are many ways to implement functions like `is_nan`. However, not all of them work with
! aggressive optimization flags. For example, for `gfortran 9.3.0`, the `ieee_is_nan` included in
! `ieee_arithmetic` does not work with `gfortran -Ofast`. See the following for discussions
! https://stackoverflow.com/questions/15944614
!
! 3. My choice of implementation is totally empirical, in the sense that I have not studied in-depth
! what the aggressive optimization flags really do, but only made some tests and found the
! implementations that worked correctly. In other words, I do not know why my implementation works
! but other implementations may not. The story may change when compilers are changed/updated.
!
! 4. Do NOT change the functions without thorough testing. Their implementation is delicate. For
! example, when compilers are invoked with aggressive optimization flags,
! (X <= HUGE(X) .AND. X >= -HUGE(X)) may differ from (ABS(X) <= HUGE(X)) ,
! (X > HUGE(X) .OR. X < -HUGE(X)) may differ from (ABS(X) > HUGE(X)) ,
! (ABS(X) > HUGE(X) .AND. X > 0) may differ from (X > HUGE(X)) .
!
! 5. Even though the functions involve invocation of ABS and HUGE, their performance (in terms of
! CPU time) turns out comparable to or even better than the functions in `ieee_arithmetic`.


use consts_mod, only : SP, DP
implicit none
private
public :: is_nan, is_finite, is_inf, is_posinf, is_neginf


interface is_nan
    module procedure is_nan_sp, is_nan_dp
end interface is_nan

interface is_finite
    module procedure is_finite_sp, is_finite_dp
end interface is_finite

interface is_posinf
    module procedure is_posinf_sp, is_posinf_dp
end interface is_posinf

interface is_neginf
    module procedure is_neginf_sp, is_neginf_dp
end interface is_neginf

interface is_inf
    module procedure is_inf_sp, is_inf_dp
end interface is_inf


contains


elemental pure function is_nan_sp(x) result(y)
implicit none
!$omp declare target
real(SP), intent(in) :: x
logical :: y
y = (.not. (x <= huge(x) .and. x >= -huge(x))) .and. (.not. abs(x) > huge(x))
end function is_nan_sp

elemental pure function is_nan_dp(x) result(y)
implicit none
!$omp declare target
real(DP), intent(in) :: x
logical :: y
y = (.not. (x <= huge(x) .and. x >= -huge(x))) .and. (.not. abs(x) > huge(x))
end function is_nan_dp


elemental pure function is_finite_sp(x) result(y)
implicit none
!$omp declare target
real(SP), intent(in) :: x
logical :: y
y = (x <= huge(x) .and. x >= -huge(x))
end function is_finite_sp

elemental pure function is_finite_dp(x) result(y)
implicit none
!$omp declare target
real(DP), intent(in) :: x
logical :: y
y = (x <= huge(x) .and. x >= -huge(x))
end function is_finite_dp


elemental pure function is_inf_sp(x) result(y)
implicit none
!$omp declare target
real(SP), intent(in) :: x
logical :: y
y = (abs(x) > huge(x))
end function is_inf_sp

elemental pure function is_inf_dp(x) result(y)
implicit none
!$omp declare target
real(DP), intent(in) :: x
logical :: y
y = (abs(x) > huge(x))
end function is_inf_dp


elemental pure function is_posinf_sp(x) result(y)
implicit none
!$omp declare target
real(SP), intent(in) :: x
logical :: y
y = (abs(x) > huge(x)) .and. (x > 0)
end function is_posinf_sp

elemental pure function is_posinf_dp(x) result(y)
implicit none
!$omp declare target
real(DP), intent(in) :: x
logical :: y
y = (abs(x) > huge(x)) .and. (x > 0)
end function is_posinf_dp


elemental pure function is_neginf_sp(x) result(y)
implicit none
!$omp declare target
real(SP), intent(in) :: x
logical :: y
y = (abs(x) > huge(x)) .and. (x < 0)
end function is_neginf_sp

elemental pure function is_neginf_dp(x) result(y)
implicit none
!$omp declare target
real(DP), intent(in) :: x
logical :: y
y = (abs(x) > huge(x)) .and. (x < 0)
end function is_neginf_dp


end module infnan_mod