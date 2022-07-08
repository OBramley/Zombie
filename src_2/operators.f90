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
    
    ! NEED TO CHECK NUM TO DET AND DET TO NUM WORK PROPERLY

    ! Turns an integer number 0 <= j <= 2**norb-1
    ! Into an arrange length norb with orbital occupancy
    ! by converting the integer into binary 
    function numtodet(j,orb) result(bini)
        
        implicit none
        integer, intent(in)::j,orb
        integer, dimension(orb)::bini
        integer::jt,k

        if(j>=2**orb) then
            write(0) "j too big"
            errorflag=1
            return
        end if

        jt=j
        do k=1, orb
            bini(k)=mod(jt,2)
            jt = jt-bini(k)
            jt = jt/2
        end do

        return

    end function numtodet

    ! Turns an integer array of 0s and 1s length norb
    ! into its corresponding binary number
    ! The reverse of numtodet

    integer function dettonum(bini,orb)

        implicit none
        integer, intent(in)::orb
        integer, dimension(orb),intent(in)::bini
        integer::temp, j 

        temp=0
        do j=1,orb
            temp = temp + ((2**j)*bini(j))
        end do

        dettonum=temp
        return
    end function dettonum

    ! Fast algorithm for application of sz operator O(N) steps

    real function szf(z1,z2)

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
            write(0,"(a,i0)") "Error in allocation of arrays in szf operator. ierr had value ", ierr
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
            temp = temp+((cc(j-1)*mult(j)*dd(j+1))*(-1**j))
        end do
        temp=temp-(cc(iorb-1)*mult(iorb))
        temp=temp*0.5

        deallocate(cc,stat=ierr)
        if(ierr==0) deallocate(dd,stat=ierr)
        if(ierr==0) deallocate(mult,stat=ierr)
        if(ierr==0) deallocate(multb,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in deallocation of arrays in sz. ierr had value ", ierr
            errorflag=1
            return
        end if

        szf=temp

        return

    end function szf

    ! Computing <zs1 | S_z^2 | zs2 > O(M^2) not O(M^3)
    real function sz2f(z1,z2)

        implicit none

        type(zombiest),intent(in)::z1,z2
        real(kind=8)::temp
        type(zombiest),dimension(:),allocatable::zt
        integer:: j, ierr,orb 

        if (errorflag .ne. 0) return
        ierr=0

        temp=0.0
        orb=size(z1%alive)
        call alloczs(zt,orb)

        do j=1, norb
            zt(j)%alive(1:orb)=z2%alive(1:orb)
            zt(j)%dead(1:orb)=z2%dead(1:orb)
        end do
        temp=0.0d0

        do j=1, orb
            call num(zt(j),j)
            temp = temp +(szf(z1,zt(j))*(-1**j))
        end do

        call dealloczs(zt)

        sz2f=temp*0.5
        return

    end function sz2f

    ! Fastest calcualtion of <zs1 |S_+S_- |zs2>

    real function spsmfast(z1,z2)

        implicit none

        type(zombiest),intent(in)::z1,z2
        real(kind=8)::tot,p1,p2,p3
        real(kind=8),dimension(:),allocatable::cc, dd, ss, tt
        integer:: j, ierr,orb,kmax,a,b,k,l  

        if (errorflag .ne. 0) return
        ierr=0

        orb=(size(z1%alive))
        kmax=orb/2
        allocate(cc(size(z1%alive)),stat=ierr)
        if(ierr==0) allocate(dd(size(z1%alive)),stat=ierr)
        if(ierr==0) allocate(ss(size(z1%alive)),stat=ierr)
        if(ierr==0) allocate(tt(size(z1%alive)),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in allocation of arrays in S- S+ operator. ierr had value ", ierr
            errorflag=1
            return
        end if

        do j=1, kmax
            a=(2*j)-1
            b=a+1
            cc(j)=(dble(dconjg(z1%alive(a))*z2%alive(a))+dble(dconjg(z1%dead(a))*z2%dead(a))) &
                    *(dble(dconjg(z1%alive(b))*z2%alive(b))+dble(dconjg(z1%dead(b))*z2%dead(b)))
            dd(j)=(dble(dconjg(z1%dead(a))*z2%dead(a)))*(dble(dconjg(z1%alive(b))*z2%alive(b)))
            ss(j)=(dble(dconjg(z1%alive(a))*z2%dead(a)))*(dble(dconjg(z1%dead(b))*z2%alive(b)))
            tt(j)=(dble(dconjg(z1%alive(a))*z2%alive(a)))*(dble(dconjg(z1%dead(b))*z2%dead(b)))
        end do
          
        tot=0.0

        do j=1, kmax
            do k=j, kmax
                p1=1
                p2=1
                if(j==1)then
                    p1=1
                else if (j==2) then
                    p1=cc(0)
                else
                    do l=1, j-1
                        p1=p1*cc(l)
                    end do
                end if

                if(k==kmax) then
                    p2=1
                else if(k==kmax-1)then
                    p2=cc(kmax)
                else 
                    do l=k+1, kmax-1
                        p2=p2*cc(l)
                    end do
                end if

                if(j==1)then
                    tot=tot+(p1*p2*tt(j))
                else if (k==j+1) then
                    tot=tot+(p1*p2*ss(j)*dd(k))+(p1*p2*ss(k)*dd(j))
                else if (k==j+2) then
                    tot = tot +(p1*p2*ss(j)*dd(k)*cc(j+1))+(p1*p2*ss(k)*dd(j)*cc(j+1))
                else 
                    p3=1
                    do l=j+1,k-1
                        p3=p3*cc(l)
                    end do
                    tot = tot +(p1*p2*ss(j)*dd(k)*p3)+(p1*p2*ss(k)*dd(j)*p3)
                end if
            end do
        end do

        deallocate(cc,stat=ierr)
        if(ierr==0) deallocate(dd,stat=ierr)
        if(ierr==0) deallocate(ss,stat=ierr)
        if(ierr==0) deallocate(tt,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in deallocation of arrays in S- S+ operator. ierr had value ", ierr
            errorflag=1
            return
        end if

        spsmfast=tot
        return
    end function spsmfast

    real function stotfast(z1,z2)
        
        implicit none
        type(zombiest),intent(in)::z1,z2
        stotfast=spsmfast(z1,z2)-szf(z1,z2)+sz2f(z1,z2)
        return
    end function stotfast

    




























END MODULE operators