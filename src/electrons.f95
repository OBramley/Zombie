subroutine spatospin(H1ea,norb,elecs,n)
    implicit none
    integer, intent(in) :: norb
    complex(kind=8), dimension(n,n) :: H1ea
    real(kind=8), dimension(norb,norb), intent(out) :: elecs
    !f2py intent(in) :: norb
    !f2py intent(hide), depend(H1ea) :: n=shape(H1ea,0)
    !f2py intent(out) output
    integer :: i, j, nspao, ii, jj, n

    elecs=0.00d0
    nspao=int(norb/2)
    
    do i=1, nspao
        do j=1, nspao
            ii=i*2
            jj=j*2
            elecs(ii,jj)=H1ea(i,j)
            elecs(ii+1,jj+1)=H1ea(i,j)
        end do
    end do
    return
end subroutine spatospin






