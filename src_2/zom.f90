MODULE zom 

    use globvars
    use outputs

    contains


    subroutine gen_ran_zs(zstore,num)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        integer, intent(in)::num
        integer::j,k
        real(kind=8)::dummy,temp
        real(kind=8),dimension(norb)::rands

        if (errorflag .ne. 0) return

        dummy=1
        do j=1,num
            call random_number(rands)
            do k=1,norb
                dummy=(2*pirl*rands(k))
                temp=sin(dummy)
                zstore(j)%alive(k)=cmplx(temp,0.0d0,kind=8)
                temp=cos(dummy)
                zstore(j)%dead(K)=cmplx(temp,0.0d0,kind=8)
            end do
        end do

        return

    end subroutine gen_ran_zs

    subroutine gen_hf_zs(zstore)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        integer::j,k,total,ierr,count
        integer, allocatable, dimension(:,:)::combs

        if (errorflag .ne. 0) return
        ierr=0

        zstore(1)%dead(1:norb)=(1.0d0,0.0d0)
        count=2
        do j=1,norb
            total=int(gamma(dble(norb+1))/(gamma(dble(j+1))*gamma(dble(norb-j)+1)))
            allocate(combs(total,nel),stat=ierr)
            if(ierr/=0) then
                write(0,"(a,i0)") "Error in combination matrix allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
            call combinations(norb,j,combs)
            do k=1, total
                call zomhf(zstore(count),combs(k,:))
                count=count+1
            end do
            deallocate(combs,stat=ierr)
            if(ierr/=0) then
                write(0,"(a,i0)") "Error in combination matrix deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if
        end do

        return  
    
    end subroutine gen_hf_zs

    subroutine combinations(n_max,m_max,final)

        ! code adapted from https://rosettacode.org/wiki/Combinations#Fortran accessed 11:45 18/07/2022
        
        implicit none

        integer, intent(in)::n_max,m_max
        integer,dimension(:,:),intent(inout)::final
        integer, dimension (m_max) :: comb
        integer:: p


        p=1
        call gen (1,final,p)

        contains 

        recursive subroutine gen (m,final,p)
   
            implicit none
            integer, intent (in) :: m
            integer, intent(inout)::p
            integer, dimension(:,:), intent(inout)::final
            integer :: n
            
            if (m > m_max) then
            final(p,1:m_max)=comb(1:m_max)
            p=p+1
            else
            do n = 1, n_max
                if ((m == 1) .or. (n > comb (m - 1))) then
                comb (m) = n
                call gen (m + 1,final,p)
                end if
            end do
            end if
        
        end subroutine gen

    end subroutine combinations


    subroutine zomhf(zom,occ)

        implicit none
        type(zombiest),intent(inout)::zom 
        integer, dimension(:), intent(in)::occ
        integer::j

        if (errorflag .ne. 0) return

        zom%dead(1:norb)=(1.0d0,0.0d0)

        do j=1, size(occ)
            zom%alive(occ(j))=(1.0d0,0.0d0)
            zom%dead(occ(j))=(0.0d0,0.0d0)
        end do

        return

    end subroutine zomhf


    subroutine genzf(zstore,num)
        
        implicit none
        type(zombiest), dimension(:), intent(inout)::zstore
        integer, intent(in)::num
        integer::j

        if (errorflag .ne. 0) return

        select case(zst)
            case('HF')
                call gen_hf_zs(zstore)
            case('RN')
                call gen_ran_zs(zstore,num)
            case('BB')
                call gen_biased_zs(zstore)
            case default
                write(0,"(a)") "Error! Initial zombie type method not recognised!"
                write(0,"(a)") "This should have been caught at the read stage!"
                errorflag = 1
                return
        end select

        do j=1,num
            call zombiewriter(zstore(j),j)
        end do

        return
    
    end subroutine genzf

END MODULE zom