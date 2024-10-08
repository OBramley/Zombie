module mod_types
    use iso_fortran_env, only:  int8, int16, int32, int64, real32, real64, input_unit, output_unit, error_unit
    implicit none
    integer, parameter :: sp        = real32
    integer, parameter :: dp        = real64
    integer, parameter :: int8_t    = int8
    integer, parameter :: int16_t   = int16
    integer, parameter :: int32_t   = int32
    integer, parameter :: stdin     = input_unit
    integer, parameter :: stdout    = output_unit
    integer, parameter :: stderr    = error_unit
end module