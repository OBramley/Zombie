MODULE outputs
    use globvars

    contains

    subroutine matrixwriter(out,size,filenm)

        implicit none

        complex(kind=8),dimension(:,:),intent(in)::out
        integer,intent(in)::size
        character(LEN=8),intent(in)::filenm
        integer::ierr,j,k

        ierr=0

        open(unit=200,file=filenm,status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening matrix file. ierr had value ", ierr
            errorflag=1
            return
        end if

        do j=1, size
            write(200,*) (out(j,k),k=1,size)
        end do

        close(200)

        return

    end subroutine matrixwriter

    subroutine zombiewriter(zom,num)

        implicit none

        type(zombiest),intent(in)::zom 
        integer,intent(in)::num
        character(LEN=15)::filenm
        integer::ierr,zomnum,j
        character(LEN=4)::nums

        ierr=0

        write(nums,"(i4.4)")num

        filenm = "zombie_"//trim(nums)//".csv"

        zomnum=300+num

        open(unit=zomnum,file=trim(filenm),status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening zombie state file. ierr had value ", ierr
            errorflag=1
            return
        end if

        write(zomnum,*) (zom%dead(j),j=1,norb)
        write(zomnum,*) (zom%alive(j),j=1,norb)

        close(zomnum)

        return

    end subroutine zombiewriter



END MODULE outputs