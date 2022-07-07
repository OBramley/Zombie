MODULE electrons

use globvars
use alarrys

contains

subroutine electronintegrals(elecs)

    implicit none

    type(elecintrgl), intent(inout)::elecs

    integer: ierr

    if (errorflag .ne. 0) return
    ierr = 0

    call allocintgrl(elecs)


end subroutine electronintegrals

subroutine spattospin1(elecs)

    implicit none

    type(elecintrgl), intent(inout)::elecs
    real(kind=8), dimension(:,:),allocatable ::h1ea
    integer:: ierr, j, k,jj,kk, nspao 

    allocate (h1ea(norb,norb), stat=ierr)
    if (ierr/=0) then
        write(0,"(a,i0)") "Error in h1ea allocation. ierr had value ", ierr
        errorflag=1
        return
    end if

    open(unit=130, file='h1ea.csv',status='old',iostat=ierr)
    if (ierr.ne.0) then
        write(0,"(a)") 'Error in opening h1ea.csv file'
        errorflag = 1
        return
    end if

    read(130,*) ((h1ea(j,k),k=1,norb),j=1,norb)
    close(130)

    nspao = int(norb/2)

    do j=1,nspao
        do k=1, nspao
            jj=j*2
            kk=k*2
            elecs%h1ei(jj,kk)=h1ea(j,k)
            elecs%h1ei(jj+1,kk+1)=elecs%h1ei(jj,kk)
        end do
    end do
    
    deallocate(h1ea, stat=ierr)
    if (ierr/=0) then
        write(0,"(a,i0)") "Error in h1ea deallocation. ierr had value ", ierr
        errorflag=1
        return
    end if

    return

end subroutine spattospin1





END MODULE electrons