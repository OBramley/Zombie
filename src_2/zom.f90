MODULE zom 

    use globvars
    use alarrays

    contains

    subroutine genzs(zstore,zomtyp)
        
        implicit none
        type(zombiest), dimension(:), intent(inout)::zstore
        character(LEN=2), intent(in):: zomtyp

        if (errorflag .ne. 0) return

        select case(zomtyp)
            case('HF')
                call gen_hf_zs(zstore)
            case('RN')
                call gen_ran_zs(zstore)
            case('BB')
                call gen_biased_zs(zstore)
            case default
                write(0,"(a)") "Error! Initial zombie type method not recognised!"
                write(0,"(a)") "This should have been caught at the read stage!"
                errorflag = 1
                return
        end select

        return
    
    end subroutine genzs

END MODULE zom