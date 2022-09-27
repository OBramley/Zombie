MODULE operators

    use globvars
    use alarrays

    contains

    ! Overlap 
    complex(kind=8) function overlap2(z1,z2)
    ! This method has a lower scaling but is not faster that the brute force method in fortran
        implicit none

        type(zombiest),intent(in)::z1,z2
        complex(kind=8)::ovrl
        complex(kind=8)::tt
        integer:: j

        if (errorflag .ne. 0) return

        ovrl=1.0d0
        do j=1, size(z1%alive)
            tt=(conjg(z1%alive(j))*z2%alive(j))+(conjg(z1%dead(j))*z2%dead(j))
            if(tt==0.0) then
                ovrl = (0.0d0,0.0d0)
                EXIT 
            end if
            ovrl=ovrl*tt
        end do
        overlap2=ovrl
        return

    end function overlap2

    complex(kind=8) function overlap(z1,z2)
    ! This way of calcualting the overlap is faster because fortran is very efficient at matrix multiplicaiton
    implicit none

    type(zombiest),intent(in)::z1,z2

    if (errorflag .ne. 0) return

    overlap=product((conjg(z1%alive)*z2%alive)+((conjg(z1%dead))*z2%dead))

    return

    end function overlap

    ! computes the vector of values formed by the derivative of the overlap with respect to each orbital
    function diff_overlap(z1_alive,z1_dead,z2_alive,z2_dead)
        implicit none
        real(kind=8),dimension(:)::z1_alive,z1_dead,z2_alive,z2_dead
        real(kind=8),dimension(norb)::diff_overlap

        diff_overlap=(z1_alive*z2_alive)+(z1_dead*z2_dead)
        return
    end function diff_overlap

    ! Creation operator 
    subroutine cr(zs, iorb)

        implicit none
        type(zombiest),intent(inout)::zs
        integer, intent(in)::iorb
        integer:: j

        ! if (errorflag .ne. 0) return

        zs%alive(iorb)=zs%dead(iorb)
        zs%dead(iorb)=(0.0d0,0.0d0)
        zs%diffalive(iorb)=zs%diffdead(iorb)
        zs%diffdead(iorb)=0.0d0
        
        do j=1, iorb-1
            zs%alive(j)=zs%alive(j)*(-1.0)
            zs%diffalive(j)=zs%diffalive(j)*(-1.0)
        end do

        return

    end subroutine cr

    ! Annihilarion operator 
    subroutine an(zs, iorb)

        implicit none
        type(zombiest),intent(inout)::zs
        integer, intent(in)::iorb
        integer:: j

        ! if (errorflag .ne. 0) return

        zs%dead(iorb)=zs%alive(iorb)
        zs%alive(iorb)=(0.0d0,0.0d0)
        zs%diffdead(iorb)=zs%diffalive(iorb)
        zs%diffalive(iorb)=0.0d0
        
        do j=1, iorb-1
            zs%alive(j)=zs%alive(j)*(-1.0)
            zs%diffalive(j)=zs%diffalive(j)*(-1.0)
        end do

        return

    end subroutine an

    ! Application of number operator
    complex(kind=8) function numf(z1,z2)
        implicit none

        type(zombiest),intent(in)::z1,z2
        complex(kind=8)::temp
        complex(kind=8),dimension(:),allocatable::cc, dd, mult, multb
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
            mult(j)=(conjg(z1%alive(j))*z2%alive(j))
            multb(j)=mult(j) + (conjg(z1%dead(j))*z2%dead(j))
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
    complex(kind=8) function nsq(z1,z2)
        implicit none

        type(zombiest),intent(in)::z1,z2
        complex(kind=8)::temp
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

        if (errorflag .ne. 0) return

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
        complex(kind=8)::tt 

        if (errorflag .ne. 0) return

        do j=1, size(zs%alive)
            tt=conjg(zs%dead(j))*zs%dead(j) + conjg(zs%alive(j))*zs%alive(j)
            if(tt==(0.0,0.0))then
                iszero=.true.
                return
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

        if (errorflag .ne. 0) return

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

        if (errorflag .ne. 0) return

        temp=0
        do j=1,orb
            temp = temp + ((2**j)*bini(j))
        end do

        dettonum=temp
        return
    end function dettonum

    ! Fast algorithm for application of sz operator O(N) steps

    complex(kind=8) function szf(z1,z2)

        implicit none

        type(zombiest),intent(in)::z1,z2
        complex(kind=8)::temp
        complex(kind=8),dimension(:),allocatable::cc, dd, mult, multb
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
            mult(j)=(conjg(z1%alive(j))*z2%alive(j))
            multb(j)=mult(j) + (conjg(z1%dead(j))*z2%dead(j))
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
    complex(kind=8) function sz2f(z1,z2)

        implicit none

        type(zombiest),intent(in)::z1,z2
        complex(kind=8)::temp
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

    complex(kind=8) function spsmfast(z1,z2)

        implicit none

        type(zombiest),intent(in)::z1,z2
        complex(kind=8)::tot,p1,p2,p3
        complex(kind=8),dimension(:),allocatable::cc, dd, ss, tt
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
            cc(j)=((conjg(z1%alive(a))*z2%alive(a))+(conjg(z1%dead(a))*z2%dead(a))) &
                    *((conjg(z1%alive(b))*z2%alive(b))+(conjg(z1%dead(b))*z2%dead(b)))
            dd(j)=((conjg(z1%dead(a))*z2%dead(a)))*((conjg(z1%alive(b))*z2%alive(b)))
            ss(j)=((conjg(z1%alive(a))*z2%dead(a)))*((conjg(z1%dead(b))*z2%alive(b)))
            tt(j)=((conjg(z1%alive(a))*z2%alive(a)))*((conjg(z1%dead(b))*z2%dead(b)))
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

    complex(kind=8) function stotfast(z1,z2)
        
        implicit none
        type(zombiest),intent(in)::z1,z2

        if (errorflag .ne. 0) return
        
        stotfast=spsmfast(z1,z2)-szf(z1,z2)+sz2f(z1,z2)
        return
    end function stotfast

    ! Finding sum_k <zom1 | b_k | zom2> vec_k for all k in norb
    complex(kind=8) function z_an_z3(z1,z2,vec)

        implicit none
        type(zombiest),intent(in)::z1,z2
        type(zombiest)::vmult
        real(kind=8),dimension(norb),intent(in)::vec
        complex(kind=8),dimension(norb)::gg,hh
        complex(kind=8)::tot
        integer::j,gmax,hmin,ierr

        if (errorflag .ne. 0) return
        ierr=0

        call alloczf(vmult)

        do j=1, norb
            vmult%alive(j)=conjg(z1%alive(j))*z2%alive(j)
            vmult%dead(j)=conjg(z1%dead(j))*z2%dead(j)
        end do


        gg(1:norb)=(0.0,0.0)
        hh(1:norb)=(0.0,0.0)
        gmax=norb
        gg(1)=vmult%dead(1)-vmult%alive(1)

        do j=2, norb
            gg(j)=gg(j-1)*(vmult%dead(j)-vmult%alive(j))
            if(gg(j)==(0.0,0.0))then
                gmax=j
                EXIT 
            end if
        end do
        

        hmin=0
        hh(norb) = vmult%dead(norb)+vmult%alive(norb)
        do j=(norb-1),1,-(1)
            hh(j)=hh(j+1)*(vmult%dead(j)+vmult%alive(j))
            if(hh(j)==(0.0,0.0))then
                hmin=j
                EXIT 
            end if
        end do


        tot=(0.0,0.0)
        if (gmax < hmin) then
            z_an_z3=tot
            return
        end if

        if(vec(1).ne.0) then
            tot = tot+(conjg(z1%dead(1))*z2%alive(1)*hh(2)*vec(1))
        end if
        do j=2,norb-1
            if(vec(j).ne.0.0) then
                tot = tot+ (gg(j-1)*conjg(z1%dead(j))*z2%alive(j)*hh(j+1)*vec(j))
            end if

        end do

        if(vec(norb).ne.0) then
            tot = tot +(gg(norb-1)*conjg(z1%dead(norb))*z2%alive(norb)*vec(norb))
        end if
        
        call dealloczf(vmult)
    
        z_an_z3=tot
        return

    end function z_an_z3
        

    function cmplx_inv(mat,size)

        implicit none

        complex(kind=8),dimension(:,:),intent(in)::mat
        integer, intent(in):: size
        integer, allocatable,dimension(:)::IPIV1
        complex(kind=8),allocatable,dimension(:)::WORK1
        complex(kind=8),dimension(size,size)::cmplx_inv
        integer:: ierr
        allocate(IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        allocate(WORK1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if   

        call ZGETRF(size,size,mat,size,IPIV1,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if
        call ZGETRI(size,mat,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if

        deallocate(IPIV1,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if

        deallocate(WORK1,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if

        cmplx_inv=mat

        return

    end function cmplx_inv








END MODULE operators