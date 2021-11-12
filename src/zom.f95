module zombie
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains the subroutines need to create a set of zombie states. The module also 
! contains the alogrithm for creating the Hamiltonian matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    implicit none

    contains

    subroutine zommake(norb, ndet, zombie, zomtyp)
    implicit none
        integer,intent(in) :: norb, ndet, zomtyp,
        complex, dimension(ndet,norb,2), intent(out) :: zombie
        integer::i,j
        

        ! 1 HF zombie states
        if(zomtyp.eq.1) then
            do i=1, ndet
                do j=1, norb
                
                end do
            end do
        ! 2 Random Zombie states
        else if (zomtyp.eq.2) then
            do i=1, ndet
                do j=1, norb
                
                end do
            end do
        ! 3 Biased Zombie states
        else if (zomtyp.eq.3) then
            do i=1, ndet
                do j=1, norb
                
                end do
            end do
        end if
        

    end subroutine zommake


    
end module zombie