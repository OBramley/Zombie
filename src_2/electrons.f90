MODULE electrons

use globvars

contains

! Program to generate one and two electron integrals from PyScf
! Some of this code has been adapted from George Booth
subroutine electronintegrals(elecs)

    implicit none

    type(elecintrgl), intent(inout)::elecs

    integer:: ierr, nlines
    

    if (errorflag .ne. 0) return
    ierr = 0

    ! n=0
    nlines = lines(nlines)
    call spattospin1(elecs,nlines)
    call spattospin2(elecs,nlines)

    open(unit=128, file='/integrals/hnuc.csv',status='old',iostat=ierr)
    if (ierr.ne.0) then
        write(0,"(a)") 'Error in opening hnuc.csv file'
        errorflag = 1
        return
    end if
    read(128,*, iostat=ierr)elecs%hnuc
    if (ierr .ne. 0) then
        write(0,"(a)") 'Error reading Hnuc'
        errorflag = 1
        return
      end if
    close(128)

    return

end subroutine electronintegrals

integer function lines(nlines)
    implicit none

    integer, intent(INOUT):: nlines
    integer:: ierr
    
    ierr=0
    nlines=0
    open(unit=129, file='/integrals/h1ea.csv',status='old',iostat=ierr)
    if (ierr.ne.0) then
        write(0,"(a)") 'Error in opening h1ea.csv file'
        errorflag = 1
        return
    end if

    do 
        read(129,*, iostat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in counting h1ea rows. ierr had value ", ierr
            errorflag=1
            return
        end if
        nlines=nlines+1
    end do

    
    return 


end function lines


subroutine spattospin1(elecs,nlines)

    implicit none

    type(elecintrgl), intent(inout)::elecs
    integer, intent(in)::nlines
    real(kind=8), dimension(:,:),allocatable ::h1ea
    integer:: ierr, j, k,jj,kk, nspao

    ierr=0 
    allocate (h1ea(nlines,nlines), stat=ierr)
    if (ierr/=0) then
        write(0,"(a,i0)") "Error in h1ea allocation. ierr had value ", ierr
        errorflag=1
        return
    end if

    open(unit=130, file='/integrals/h1ea.csv',status='old',iostat=ierr)
    if (ierr.ne.0) then
        write(0,"(a)") 'Error in opening h1ea.csv file'
        errorflag = 1
        return
    end if

    
    read(130,*) ((h1ea(j,k),k=1,nlines),j=1,nlines)
    close(130)

    nspao = int(norb/2)

    do j=1,nspao
        do k=1, nspao
            jj=(j*2)-1
            kk=(k*2)-1
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

subroutine spattospin2(elecs,nlines)

    implicit none

    type(elecintrgl), intent(inout)::elecs
    integer, intent(in)::nlines
    real(kind=8), dimension(:,:,:,:),allocatable ::h2ea
    integer:: ierr, j, k, l, m, jj, kk, ll, mm, nspao
    character(len=4)::val1,val2

    ierr=0 
    allocate (h2ea(nlines,nlines,nlines,nlines), stat=ierr)
    if (ierr/=0) then
        write(0,"(a,i0)") "Error in h2ea allocation. ierr had value ", ierr
        errorflag=1
        return
    end if

    nspao = int(norb/2)
    do j=1,norb
        do k=1, norb
            write(val1,'(i0)')j
            write(val2,'(i0)')k
            open(unit=(131+j+k), file='/integrals/h2ea_'//trim(val1)//'_'//trim(val2)//'.csv', status='old',iostat=ierr)
            if (ierr.ne.0) then
                write(0,"(a)") 'Error in opening h1ea.csv file'
                errorflag = 1
                return
            end if
            read((131+j+k),*) ((h2ea(j,k,l,m),k=1,nlines),l=1,nlines)
            close(131+j+k)
        end do
    end do

    do j=1,nspao
        jj=(2*j)-1
        do k=1, nspao
            kk=(2*k)-1
            do l=1, nspao
                ll=(2*l)-1
                do m=1, nspao
                    mm=(2*m)-1
                    elecs%h2ei(jj,ll,kk,mm)=h2ea(j,k,l,m)
                    elecs%h2ei(jj+1,ll,kk+1,mm)=h2ea(j,k,l,m)
                    elecs%h2ei(jj,ll+1,kk,mm+1)=h2ea(j,k,l,m)
                    elecs%h2ei(jj+1,ll+1,kk+1,mm+1)=h2ea(j,k,l,m)
                end do
            end do
        end do
    end do

    deallocate(h2ea, stat=ierr)
    if (ierr/=0) then
        write(0,"(a,i0)") "Error in h1ea deallocation. ierr had value ", ierr
        errorflag=1
        return
    end if

    return
                    
end subroutine spattospin2
    





END MODULE electrons