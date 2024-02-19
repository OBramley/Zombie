MODULE operators

    use mod_types
    use globvars
    use alarrays

    contains

    
    

    ! Creation operator 
    subroutine cr(zs, iorb)

        implicit none
        type(zombiest),intent(inout)::zs
        integer, intent(in)::iorb
    
        if (errorflag .ne. 0) return

        zs%val(iorb)=zs%val(iorb+norb)
        zs%val(norb+iorb)=0
        zs%val(1:iorb-1)=(-1)*zs%val(1:iorb-1)
       
        return

    end subroutine cr

    ! Annihilarion operator 
    subroutine an(zs, iorb)

        implicit none
        type(zombiest),intent(inout)::zs
        integer, intent(in)::iorb
        if (errorflag .ne. 0) return

        zs%val(norb+iorb)=zs%val(iorb)
        zs%val(iorb)=0
        zs%val(1:iorb-1)=(-1)*zs%val(1:iorb-1)

        return

    end subroutine an

    ! Application of number operator
    real(wp) function numf(z1,z2)
        implicit none

        type(zombiest),intent(in)::z1,z2
        real(wp)::temp
        real(wp),dimension(norb)::cc, dd, mult, multb
        integer:: j

        do j=1, norb
            mult(j)=(z1%val(j))*z2%val(j)
            multb(j)=mult(j) + (z1%val(j+norb))*z2%val(j+norb)
        end do

        cc(1)=mult(1)
        dd(norb)=multb(norb)
        do j=2, norb
            cc(j)=cc(j-1)*multb(j)
        end do
        do j=(norb-1),1,-(1)
            dd(j)=dd(j+1)*multb(j)
        end do
        ! print*, "cc", cc
        ! print*, "dd", dd
        ! print*, "mult", mult
        ! print*, "multb", multb

        temp=mult(1)*dd(2)
        do j=2, (norb-2)
            temp = temp+(cc(j-1)*mult(j)*dd(j+1))
        end do
        temp=temp+cc(norb-1)*mult(norb)

        numf=temp
        return

    end function numf

    ! Number operator on a specific orbital 
    subroutine num(zs,iorb)
        implicit none

        type(zombiest),intent(inout)::zs
        integer, intent(in)::iorb

        if (errorflag .ne. 0) return

        zs%val(iorb+norb)=0.0d0

        return
    end subroutine num

    ! N squared operator
    real(wp) function nsq(z1,z2)
        implicit none

        type(zombiest),intent(in)::z1,z2
        real(wp)::temp
        type(zombiest),dimension(:),allocatable::zt
        integer:: j

        call alloczs(zt,norb)
       
        do j=1, norb
            zt(j)%val=z2%val
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
            if(zs%val(j+norb).lt.1.0d-15)then
                if((abs(zs%val(j)-1.0d0).lt.1.0d-15).or.(abs(zs%val(j)+1.0d0).lt.1.0d-15))then
                    CYCLE
                else
                    isdet=.false.
                    return
                end if
                
            else if((abs(zs%val(j+norb)-1.0d0).lt.1.0d-15).or.(abs(zs%val(j+norb)+1.0d0).lt.1.0d-15))then
                if(zs%val(j).lt.1.0d-15) then
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
        real(wp)::tt 
        do j=1, norb
            tt=zs%val(j+norb) + zs%val(j)
            if(tt.lt.1.0d-15)then
                iszero=.true.
                return
            end if
        end do

        iszero=.false.
        return
    end function iszero

    logical function occ_iszero(zs)

        implicit none

        type(zombiest),intent(in)::zs
        integer::j
        real(wp)::tt 

        occ_iszero=.false.
        
        do j=1, norb
            tt=(zs%val(j)*zs%val(j))+(zs%val(j+norb)*zs%val(j+norb))
            if(tt.lt.1.0d-15)then
                occ_iszero=.true.
                exit
            end if
        end do
        
        return
    end function occ_iszero
    
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

    real(wp) function szf(z1,z2)

        implicit none

        type(zombiest),intent(in)::z1,z2
        real(wp)::temp
        real(wp),dimension(norb)::cc, dd, mult, multb
        integer:: j

        do j=1, norb
            mult(j)=(z1%val(j))*z2%val(j)
            multb(j)=mult(j) + (z1%val(j+norb))*z2%val(j+norb)
        end do

        cc(1)=mult(1)
        dd(norb)=multb(norb)
        do j=2, norb
            cc(j)=cc(j-1)*multb(j)
        end do
        do j=(norb-1),1,-(1)
            dd(j)=dd(j+1)*multb(j)
        end do
        temp=mult(1)*dd(2)
        do j=2, (norb-2)
            temp = temp+((cc(j-1)*mult(j)*dd(j+1))*(-1**j))
        end do
        temp=temp-(cc(norb-1)*mult(norb))
        temp=temp*0.5
        szf=temp

        return

    end function szf

    ! Computing <zs1 | S_z^2 | zs2 > O(M^2) not O(M^3)
    real(wp) function sz2f(z1,z2)

        implicit none

        type(zombiest),intent(in)::z1,z2
        real(wp)::temp
        type(zombiest),dimension(:),allocatable::zt
        integer:: j,orb 

     
        temp=0.0d0
        call alloczs(zt,norb)

        do j=1, norb
            zt(j)%val=z2%val
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

    real(wp) function spsmfast(z1,z2)

        implicit none

        type(zombiest),intent(in)::z1,z2
        real(wp)::tot,p1,p2,p3
        real(wp),dimension(norb)::cc, dd, ss, tt
        integer:: j,kmax,a,b,k,l  

        
        kmax=norb/2
        
        do j=1, kmax
            a=(2*j)-1
            b=a+1
            cc(j)=((z1%val(a)*z2%val(a))+(z1%val(a+norb)*z2%val(a+norb))) &
                    *(z1%val(b)*z2%val(b))+(z1%val(b+norb)*z2%val(b+norb))
            dd(j)=(z1%val(a+norb)*z2%val(a+norb))*(z1%val(b)*z2%val(b))
            ss(j)=(z1%val(a)*z2%val(a+norb))*(z1%val(b+norb)*z2%val(b))
            tt(j)=(z1%val(a)*z2%val(a))*(z1%val(b+norb)*z2%val(b+norb))
        end do
          
        tot=0.0

        do j=1, kmax
            do k=j, kmax
                p1=1
                p2=1
                if(j==1)then
                    p1=1
                else if (j==2) then
                    p1=cc(1)
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

        spsmfast=tot
        return
    end function spsmfast

    real(wp) function stotfast(z1,z2)
        
        implicit none
        type(zombiest),intent(in)::z1,z2
        
        stotfast=spsmfast(z1,z2)-szf(z1,z2)+sz2f(z1,z2)
        return
    end function stotfast

    







END MODULE operators