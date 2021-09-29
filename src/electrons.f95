subroutine spatospin(H1ea,norb,elecs)
    implicit none
    integer, intent(in) :: norb
    real(kind=8), dimension(28,28), intent(in) :: H1ea
    real(kind=8), dimension(norb,norb), intent(out) :: elecs
    integer :: i, j, nspao, ii, jj

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






