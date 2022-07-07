MODULE operators

    use globvars
    use alarrays

    contains


    ! Overlap 
    real function overlap(z1,z2)
        implicit none

        type(zombiest),intent(in)::z1,z2
        real(kind=8)::ovrl
        real(kind=8)::tt
        integer:: j

        if (errorflag .ne. 0) return

        ovrl=1.0d0
        do j=1, size(z1%alive)
            tt=dble(dconjg(z1%alive(j))*z2%alive(j))+(dconjg(z1%dead(j))*z2%dead(j))
            if(tt==0.0) then
                ovrl = 0.0d0
                EXIT 
            end if
            ovrl=ovrl*tt
        end do
        overlap=ovrl
        return

    end function overlap

    ! Creation operator 
    subroutine cr(zs, iorb)

        implicit none
        type(zombiest),intent(inout)::zs
        integer, intent(in)::iorb
        integer:: j

        if (errorflag .ne. 0) return

        zs%alive(iorb)=zs%dead(iorb)
        zs%dead(iorb)=(0.0d0,0.0d0)
        
        do j=1, iorb-1
            zs%alive(j)=zs%alive(j)*(-1.0)
        end do

        return

    end subroutine cr

    ! Annihilarion operator 
    subroutine an(zs, iorb)

        implicit none
        type(zombiest),intent(inout)::zs
        integer, intent(in)::iorb
        integer:: j

        if (errorflag .ne. 0) return

        zs%dead(iorb)=zs%alive(iorb)
        zs%alive(iorb)=(0.0d0,0.0d0)
        
        do j=1, iorb-1
            zs%alive(j)=zs%alive(j)*(-1.0)
        end do

        return

    end subroutine an

    ! Application of number operator
    real function numf(z1,z2)
        implicit none

        type(zombiest),intent(in)::z1,z2
        real(kind=8)::temp
        real(kind=8),dimension(:),allocatable::cc, dd, mult, multb
        integer:: j, ierr,iorb 

        if (errorflag .ne. 0) return
        ierr=0

        iorb=size(z1%alive)
        allocate(cc(size(z1%alive)),stat=ierr)
        if(ierr==0) allocate(dd(size(z1%alive)),stat=ierr)
        if(ierr==0) allocate(mult(size(z1%alive)),stat=ierr)
        if(ierr==0) allocate(multb(size(z1%alive)),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in allocation of arrays in number operator. ierr had value ", ierr
            errorflag=1
            return
        end if

        do j=1, iorb
            mult(j)=dble(dconjg(z1%alive(j))*z2%alive(j))
            multb(j)=mult(j) + dble(dconjg(z1%dead(j))*z2%dead(j))
        end do

        cc(1)=mult(1)
        dd(iorb)=multb(iorb)
        do j=2, iorb
            cc(j)=cc(j-1)*multb(j)
        end do
        do j=(iorb-1),1,-(1)
            dd(j)=dd(j+1)*multb(j)
        end do
        temp=mult(1)*dd(2)
        do j=2, (iorb-2)
            temp = temp+cc(j-1)*mult(j)*dd(j+1)
        end do
        temp=temp+cc(iorb-1)*mult(iorb)

        deallocate(cc,stat=ierr)
        if(ierr==0) deallocate(dd,stat=ierr)
        if(ierr==0) deallocate(mult,stat=ierr)
        if(ierr==0) deallocate(multb,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in deallocation of arrays in number operator. ierr had value ", ierr
            errorflag=1
            return
        end if

        numf=temp
        return

    end function numf

    ! Number operator on a specific orbital 
    subroutine num(zs,iorb)
        implicit none

        type(zombiest),intent(inout)::zs
        integer, intent(in)::iorb

        if (errorflag .ne. 0) return

        zs%dead(iorb)=(0.0d0,0.0d0)

        return
    end subroutine num

    ! N squared operator
    real function nsq(z1,z2)
        implicit none

        type(zombiest),intent(in)::z1,z2
        real(kind=8)::temp
        type(zombiest),dimension(:),allocatable::zt
        integer:: j, ierr

        if (errorflag .ne. 0) return
        ierr=0

        call alloczs(zt,norb)

        do j=1, norb
            zt(j)%alive(1:norb)=z2%alive(1:norb)
            zt(j)%dead(1:norb)=z2%dead(1:norb)
        end do
        temp=0.0d0

        do j=j, norb
            call num(zt(j),j)
            temp = temp + numf(z1,zt(j))
        end do

        call dealloczs(zt)
        nsq=temp
        return
    end function nsq

    ! Checks if a given zombie state is equal to a single determiant an the entreis for each
    ! spin orbital are either 0 and 1 or 0 and -1
    logical function isdet(zs)
        
        implicit none

        type(zombiest),intent(in)::zs
        integer:: j

        do j=1, norb
            if(zs%dead(j)==(0.0d0,0.0d0))then
                if((zs%alive(j)==(1.0d0,0.0d0)).or.(zs%alive(j)==(-1.0d0,0.0d0))) then
                    CYCLE
                else
                    isdet=.false.
                    return
                end if
            else if((zs%dead(j)==(1.0d0,0.0d0)).or.(zs%dead(j)==(-1.0d0,0.0d0))) then
                if(zs%alive(j)==(0.0d0,0.0d0)) then
                    CYCLE
                else
                    isdet=.false.
                    return
                end if
            else 
                isdet=.false.
                return
            end if
        end do
        
        isdet=.true.
        return
    end function isdet

        ! Determines if a given zombie state vanishes. Returns True if state vanishes
        ! False otherwise
    logical function iszero(zs)
        implicit none

        type(zombiest),intent(in)::zs
        integer:: j
        real(kind=8)::tt 

        do j=1, size(zs%alive)
            tt=dconjg(zs%dead(j))*zs%dead(j) + dconjg(zs%alive(j))*zs%alive(j)
            if(tt==0.0)then
                iszero=.true.
            end if
        end do

        iszero=.false.
        return
    end function iszero

    






    























END MODULE operators