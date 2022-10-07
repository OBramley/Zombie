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
            tt=(conjg(z1%sin(j))*z2%sin(j))+(conjg(z1%cos(j))*z2%cos(j))
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

     
  
        overlap=product(((conjg(z1%sin)*z2%sin)*(z1%alive*z2%alive))+((conjg(z1%cos)*z2%cos)*(z1%dead*z2%dead)))
        
        return

    end function overlap


    function diff_overlap_cran(z1,z2,dtype,zomt,annihilate2,create2,occupancy)

        implicit none

        type(zombiest),intent(in)::z1,z2
        real(kind=8),dimension(2,norb)::diff_overlap_cran
        complex(kind=8),dimension(:,:)::zomt
        integer,intent(in)::dtype,annihilate2,create2
        integer,dimension(:,:),intent(in)::occupancy
        real(kind=8),dimension(norb)::bra_prod,ket_prod,prod
        integer::j
        real(kind=8),dimension(norb)::temp,temp2

        if (errorflag .ne. 0) return

        prod=real(((conjg(z1%sin)*zomt(1,:)))+((conjg(z1%cos)*zomt(2,:))))
        diff_overlap_cran=0
        ! Differntial of overlap of the same ZS 0 unless an operator has acted on ZS
        if(dtype.eq.0)then
            if(annihilate2.eq.create2)then
                if((real(z2%cos(annihilate2)).eq.0).and.(real(z2%sin(annihilate2)).eq.1).or.(real(z2%sin(annihilate2)).eq.0))then
                    diff_overlap_cran(1,annihilate2)=0
                else
                    bra_prod=prod           !dead amplitude is zero
                    bra_prod(annihilate2)=sin(2*z1%phi(annihilate2))*occupancy(1,annihilate2)
                    diff_overlap_cran(1,annihilate2)=product(bra_prod)   
                end if 
            else if(annihilate2.ne.create2)then
                bra_prod=prod
                if((real(z2%cos(annihilate2)).eq.0).and.(real(z2%sin(annihilate2)).eq.1)) then
                    bra_prod(annihilate2)=-1*occupancy(2,annihilate2)
                else if((real(z2%sin(annihilate2)).eq.0).and.(real(z2%cos(annihilate2)).eq.1)) then
                    bra_prod(annihilate2)=occupancy(2,annihilate2)
                else
                    bra_prod(annihilate2)=cos(2*z1%phi(annihilate2))*occupancy(2,annihilate2)
                end if
                diff_overlap_cran(1,annihilate2)=product(bra_prod)   
                bra_prod=prod
                if((real(z2%cos(create2)).eq.0).and.(real(z2%sin(create2)).eq.1)) then
                    bra_prod(create2)=-1*occupancy(1,create2)
                else if((real(z2%sin(create2)).eq.0).and.(real(z2%cos(create2)).eq.1)) then
                    bra_prod(create2)=occupancy(1,create2)
                else
                    bra_prod(create2)=cos(2*z1%phi(create2))*occupancy(1,create2)
                end if               !dead amplitude is zero
                diff_overlap_cran(1,create2)=product(bra_prod)  
            end if
            diff_overlap_cran(2,:)=diff_overlap_cran(1,:) 
        else 
            bra_prod=real(z1%cos*z2%sin*occupancy(1,:)-z1%sin*z2%cos*occupancy(2,:))
            ket_prod=real(z1%sin*z2%cos*occupancy(1,:)-z1%cos*z2%sin*occupancy(2,:))                  
            if(annihilate2.eq.create2)then  !dead amplitude is zero
                bra_prod(annihilate2)=real(z1%cos(annihilate2)*z2%sin(annihilate2)*occupancy(1,annihilate2))
                ket_prod(annihilate2)=real(z1%sin(annihilate2)*z2%cos(annihilate2)*occupancy(1,annihilate2))
            else                   
                bra_prod(annihilate2)=-real(z1%sin(annihilate2)*z2%sin(annihilate2))*occupancy(2,annihilate2)
                ket_prod(annihilate2)=real(z1%cos(annihilate2)*z2%cos(annihilate2))*occupancy(2,annihilate2)!alive amplitude is zero
                bra_prod(create2)=real(z1%cos(create2)*z2%cos(create2))*occupancy(1,create2) 
                ket_prod(create2)=-real(z1%sin(create2)*z2%sin(create2))*occupancy(1,create2) !dead amplitude is zero
            end if
            do j=1,norb
                temp=prod
                temp2=prod
                temp(j)=bra_prod(j)
                temp2(j)=ket_prod(j)
                diff_overlap_cran(1,j)=product(temp)
                diff_overlap_cran(2,j)=product(temp2)
            end do   
        end if
        
        RETURN


    end function diff_overlap_cran


    ! computes the vector of values formed by the derivative of the overlap with respect to each orbital. 
    ! Does not have capability to deal with states where creation and annihilation operators have acted
    function diff_overlap(z1,z2)
    
        implicit none

        type(zombiest),intent(in)::z1,z2
        real(kind=8),dimension(2,norb)::diff_overlap
        real(kind=8),dimension(norb)::bra_prod,ket_prod,prod
        integer::j
        real(kind=8),dimension(norb)::temp,temp2

        if (errorflag .ne. 0) return

       
            prod=real((conjg(z1%sin)*z2%sin)+(conjg(z1%cos)*z2%cos))
            bra_prod=real(conjg(z1%cos)*z2%sin)-real(conjg(z1%sin)*z2%cos)
            ket_prod=real(conjg(z1%sin)*z2%cos)-real(conjg(z1%cos)*z2%sin)
            do j=1,norb
                temp=prod
                temp2=prod
                temp(j)=bra_prod(j)
                temp2(j)=ket_prod(j)
                diff_overlap(1,j)=product(temp)
                diff_overlap(2,j)=product(temp2)
            end do

        return
        
    end function diff_overlap


    ! Creation operator 
    subroutine cr(zs, iorb)

        implicit none
        type(zombiest),intent(inout)::zs
        integer, intent(in)::iorb
        ! integer:: j

        if (errorflag .ne. 0) return

        zs%alive(iorb)=zs%dead(iorb)
        zs%dead(iorb)=0
        zs%alive(1:iorb-1)=(-1)*zs%alive(1:iorb-1)
        zs%phi(iorb)=zs%phi(iorb)+0.5*pirl
        zs%sin(iorb)=sin(cmplx(zs%phi(iorb),0.0,kind=8))

        return

    end subroutine cr

    ! Annihilarion operator 
    subroutine an(zs, iorb)

        implicit none
        type(zombiest),intent(inout)::zs
        integer, intent(in)::iorb
        ! integer:: j

        if (errorflag .ne. 0) return

        zs%dead(iorb)=zs%alive(iorb)
        zs%alive(iorb)=0
        zs%alive(1:iorb-1)=(-1)*zs%alive(1:iorb-1)
        zs%phi(iorb)=zs%phi(iorb)-0.5*pirl
        zs%cos(iorb)=cos(cmplx(zs%phi(iorb),0.0,kind=8))


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
            mult(j)=(conjg(z1%sin(j))*z2%sin(j))
            multb(j)=mult(j) + (conjg(z1%cos(j))*z2%cos(j))
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
        integer::tt 

        if (errorflag .ne. 0) return

        do j=1, size(zs%alive)
            tt=zs%dead(j) + zs%alive(j)
            if(tt==0)then
                iszero=.true.
                return
            end if
        end do

        iszero=.false.
        return
    end function iszero

    logical function occ_iszero(mat)

        implicit none
        complex(kind=8),dimension(:,:),intent(in)::mat
        integer::j
        complex(kind=8)::tt 
    
        do j=1, size(mat(1,:))
            tt=conjg(mat(1,j))*mat(1,j) + conjg(mat(2,j))*mat(2,j)
            if(tt==(0.0,0.0))then
                occ_iszero=.true.
                return
            end if
        end do

        occ_iszero=.false.
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
            mult(j)=(conjg(z1%sin(j))*z2%sin(j))
            multb(j)=mult(j) + (conjg(z1%cos(j))*z2%cos(j))
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
            zt(j)%sin(1:orb)=z2%sin(1:orb)
            zt(j)%cos(1:orb)=z2%cos(1:orb)
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
            cc(j)=((conjg(z1%sin(a))*z2%sin(a))+(conjg(z1%cos(a))*z2%cos(a))) &
                    *((conjg(z1%sin(b))*z2%sin(b))+(conjg(z1%cos(b))*z2%cos(b)))
            dd(j)=((conjg(z1%cos(a))*z2%cos(a)))*((conjg(z1%sin(b))*z2%sin(b)))
            ss(j)=((conjg(z1%sin(a))*z2%cos(a)))*((conjg(z1%cos(b))*z2%sin(b)))
            tt(j)=((conjg(z1%sin(a))*z2%sin(a)))*((conjg(z1%cos(b))*z2%cos(b)))
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

        complex(kind=8),dimension(:,:),intent(in)::z1,z2
        ! type(zombiest),intent(in)::z1,z2
        complex(kind=8),dimension(:,:),allocatable::vmult
        ! type(zombiest)::vmult
        real(kind=8),dimension(norb),intent(in)::vec
        complex(kind=8),dimension(norb)::gg,hh
        complex(kind=8)::tot
        integer::j,gmax,hmin,ierr

        if (errorflag .ne. 0) return
        ierr=0

        ! call alloczf(vmult)
        allocate(vmult(2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in vmult allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        
        vmult=conjg(z1)*(z2)
        

        gg(1:norb)=(0.0,0.0)
        hh(1:norb)=(0.0,0.0)
        gmax=norb
        gg(1)=vmult(2,1)-vmult(1,1)

        do j=2, norb
            gg(j)=gg(j-1)*(vmult(2,j)-vmult(1,j))
            if(gg(j)==(0.0,0.0))then
                gmax=j
                EXIT 
            end if
        end do
        
        hmin=0
        hh(norb) = vmult(2,norb)+vmult(1,norb)
        do j=(norb-1),1,-(1)
            hh(j)=hh(j+1)*(vmult(2,j)+vmult(1,j))
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
            tot = tot+(conjg(z1(2,1))*z2(1,1)*hh(2)*vec(1))
        end if

        do j=2,norb-1
            if(vec(j).ne.0.0) then
                tot = tot+ (gg(j-1)*conjg(z1(2,j))*z2(1,j)*hh(j+1)*vec(j))
            end if
        end do

        if(vec(norb).ne.0) then
            tot = tot +(gg(norb-1)*conjg(z1(2,norb))*z2(1,norb)*vec(norb))
        end if

        deallocate(vmult)
    
        z_an_z3=tot
    
        return

    end function z_an_z3

    function z_an_z3_diff(z1,z2,vec,dtype,occupancy,annihilate1,annihilate1_2,annihilate2,zom1,zom2)

        implicit none
        complex(kind=8),dimension(:,:),intent(in)::z1,z2
        real(kind=8),dimension(:,:),allocatable::vmult_dd,vmult_1d, vmult_2d,vmult
        real(kind=8),dimension(norb),intent(in)::vec
        integer,dimension(:,:),intent(in)::occupancy
        integer, intent(in)::dtype,annihilate1,annihilate1_2,annihilate2
        type(zombiest),intent(in)::zom1,zom2
        real(kind=8),dimension(2,norb)::z_an_z3_diff
        real(kind=8),dimension(norb)::temp,temp2
        real(kind=8),dimension(norb)::gg_1,hh_1,gg_2,hh_2
        real(kind=8)::tot1, tot2
        integer::j,k,gmax1,hmin1,gmax2,hmin2,ierr,breakflag

        if (errorflag .ne. 0) return
        ierr=0

        allocate(vmult(2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in vmult allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        vmult=real(conjg(z1)*z2)
       
        if(dtype.eq.0) then !Differentiation when zombie states are the same
            allocate(vmult_dd(2,norb),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult_dd allocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
            temp=0
            do j=1, norb    !Differentiating w.r.t to orbital j
                breakflag=0
                vmult_dd=vmult
                if(annihilate1.eq.annihilate1_2)then
                    temp(j)=0
                    cycle
                else if((annihilate1.eq.j).or.(annihilate1_2.eq.j))then
                    if(annihilate2.eq.j)then    !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(1j)
                        vmult_dd(1,j)=0
                        vmult_dd(2,j)=2*real(zom1%sin(j)*zom1%cos(j))*occupancy(2,j) !(sin2x=2sinxcosx)
                    else
                        vmult_dd(1,j)=0        !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(0j)
                        vmult_dd(2,j)=(1.0-2*((real(zom1%sin(j)))**2))*occupancy(2,j)    !(cos2x = cos^x -sin^2x)
                    end if
                else if((annihilate1.ne.j).or.(annihilate1_2.ne.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(0j)*a^(a)_(1j)
                        vmult_dd(1,j)=0
                        vmult_dd(2,j)=(1.0-2*((real(zom1%sin(j)))**2))*occupancy(2,j)
                    else
                        breakflag=j !before diff alive:a^(a)_(1j)*a^(a)_(1j) dead:a^(a)_(0j)*a^(a)_(0j)
                    end if          !Unless an opeator acts at position j this evaluates to 0
                end if
            
                gg_1(1:norb)=(0.0,0.0)
                hh_1(1:norb)=(0.0,0.0)
                gmax1=norb
                gg_1(1)=vmult_dd(2,1)-vmult_dd(1,1)

                do k=2, norb
                    gg_1(k)=gg_1(k-1)*(vmult_dd(2,k)-vmult_dd(1,k))
                    if(gg_1(k)==(0.0))then
                        gmax1=k
                        EXIT 
                    end if
                end do
                
                hmin1=0
                hh_1(norb) = vmult_dd(2,norb)+vmult_dd(1,norb)
                do k=(norb-1),1,-(1)
                    hh_1(k)=hh_1(k+1)*(vmult_dd(2,k)+vmult_dd(1,k))
                    if(hh_1(k)==(0.0))then
                        hmin1=k
                        EXIT 
                    end if
                end do

                tot1=(0.0)
                if (gmax1 < hmin1) then
                    temp(j)=tot1
                    cycle
                end if

                if((breakflag.eq.0).or.(breakflag.eq.1))then
                    if(vec(1).ne.0) then
                        if(breakflag.eq.1)then
                            temp(j)=(1.0-2*((real(zom1%sin(1)))**2))*occupancy(2,1)*hh_1(2)*vec(1)
                            cycle
                        end if
                        if(j.eq.1)then
                            if((annihilate1.eq.1).or.(annihilate1_2.eq.1))then
                                if(annihilate2.ne.1)then
                                    tot1 = tot1+(2*real(zom1%sin(1)*zom1%cos(1))*occupancy(2,1)*hh_1(2)*vec(1))
                                end if
                            else
                                if(annihilate2.ne.1)then
                                    tot1 = tot1+((1.0-2*((real(zom1%sin(1)))**2))*occupancy(2,1)*hh_1(2)*vec(1))
                                end if 
                            end if
                        else
                            tot1 = tot1+(REAL(conjg(z1(2,1))*z2(1,1))*hh_1(2)*vec(1))
                        end if
                    end if
                else if((breakflag.lt.norb))then
                    if(breakflag.ne.0)then
                        temp(j)=(gg_1(j-1)*((1.0-2*((real(zom1%sin(j)))**2))*occupancy(2,j))*hh_1(j+1)*vec(j))
                        cycle
                    end if
                    do k=2,norb-1
                        if(vec(k).ne.0.0) then
                            if(k.eq.j)then
                                if((annihilate1.eq.k).or.(annihilate1_2.eq.k))then
                                    if(annihilate2.ne.k)then
                                        tot1 = tot1+(gg_1(k-1)*(2*real(zom1%sin(k)*zom1%cos(k))*occupancy(2,k))*hh_1(k+1)*vec(k))
                                    end if
                                else
                                    if(annihilate2.ne.k)then
                                        tot1 = tot1+(gg_1(k-1)*((1.0-2*((real(zom1%sin(k)))**2))*occupancy(2,k))*hh_1(k+1)*vec(k))
                                    end if 
                                end if
                            else
                                tot1 = tot1+ (gg_1(k-1)*REAL(conjg(z1(2,k))*z2(1,k))*hh_1(k+1)*vec(k))
                            end if
                        end if
                    end do
                else
                    if(breakflag.eq.norb)then
                        temp(j)=(gg_1(norb-1)*((1.0-2*((real(zom1%sin(norb)))**2))*occupancy(2,norb))*vec(norb))
                        cycle
                    end if
                    if(vec(norb).ne.0) then
                        if(norb.eq.j)then
                            if((annihilate1.eq.norb).or.(annihilate1_2.eq.norb))then
                                if(annihilate2.ne.norb)then
                                    tot1 = tot1 +(gg_1(norb-1)*(2*real(zom1%sin(norb)*zom1%cos(norb))*occupancy(2,norb))*vec(norb))
                                end if
                            else
                                if(annihilate2.ne.norb)then
                                    tot1 = tot1+ (gg_1(norb-1)*((1.0-2*((real(zom1%sin(norb)))**2))*occupancy(2,norb))*vec(norb))
                                end if    
                            end if
                        else
                            tot1 = tot1 +(gg_1(norb-1)*REAL(conjg(z1(2,norb))*z2(1,norb))*vec(norb))
                        end if
                    end if
                end if

                temp(j)=tot1
            end do
            z_an_z3_diff(1,:)=temp(:)
            z_an_z3_diff(2,:)=temp(:)
            deallocate(vmult,stat=ierr)
            if(ierr==0) deallocate(vmult_dd,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult deallocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
        
        else !Differentiation when zombie states are not the same
            allocate(vmult_1d(2,norb),stat=ierr)
            if(ierr==0)allocate(vmult_2d(2,norb),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult_1d/vmult_2d allocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
            temp=0
            temp2=0
            do j=1, norb
                vmult_1d=vmult
                vmult_2d=vmult
                if(annihilate1.eq.annihilate1_2)then
                    temp(j)=0
                    cycle
                else if((annihilate1.eq.j).or.(annihilate1_2.eq.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                        vmult_1d(1,j)=0
                        vmult_1d(2,j)=real(zom1%cos(j)*zom2%sin(j))*occupancy(2,j)
                        vmult_2d(1,j)=0
                        vmult_2d(2,j)=real(zom1%sin(j)*zom2%cos(j))*occupancy(2,j)
                    else !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(0j)
                        vmult_1d(1,j)=0
                        vmult_1d(2,j)=real(zom1%cos(j)*zom2%cos(j))*occupancy(2,j)
                        vmult_2d(1,j)=0
                        vmult_2d(2,j)=-real(zom1%sin(j)*zom2%sin(j))*occupancy(2,j)
                    end if
                else if((annihilate1.ne.j).or.(annihilate1_2.ne.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                        vmult_1d(1,j)=0
                        vmult_1d(2,j)=-real(zom1%sin(j)*zom2%sin(j))*occupancy(2,j)
                        vmult_2d(1,j)=0
                        vmult_2d(2,j)=real(zom1%cos(j)*zom2%cos(j))*occupancy(2,j)
                    else    !before diff alive:a^(a)_(1j)*a^(b)_(1j) dead:a^(a)_(0j)*a^(b)_(0j)
                        ! breakflag=j
                        vmult_1d(1,j)=real(zom1%cos(j)*zom2%sin(j))*occupancy(1,j)
                        vmult_1d(2,j)=-real(zom1%sin(j)*zom2%cos(j))*occupancy(2,j)
                        vmult_2d(1,j)=real(zom1%sin(j)*zom2%cos(j))*occupancy(1,j)
                        vmult_2d(2,j)=-real(zom1%cos(j)*zom2%sin(j))*occupancy(2,j)
                    end if
                end if

                gg_1(1:norb)=(0.0,0.0)
                hh_1(1:norb)=(0.0,0.0)
                gg_2(1:norb)=(0.0,0.0)
                hh_2(1:norb)=(0.0,0.0)
                gmax1=norb
                gmax2=norb
                gg_1(1)=(vmult_1d(2,1))-(vmult_1d(1,1))
                gg_2(1)=(vmult_2d(2,1))-(vmult_2d(1,1))

                do k=2, norb
                    gg_1(k)=gg_1(k-1)*((vmult_1d(2,k)-vmult_1d(1,k)))
                    if(gg_1(k)==(0.0,0.0))then
                        gmax1=k
                        EXIT 
                    end if
                end do

                do k=2, norb
                    gg_2(k)=gg_2(k-1)*((vmult_2d(2,k)-vmult_2d(1,k)))
                    if(gg_2(k)==(0.0,0.0))then
                        gmax2=k
                        EXIT 
                    end if
                end do
                
                hmin1=0
                hmin2=0
                hh_1(norb) = (vmult_1d(2,norb)+vmult_1d(1,norb))
                hh_2(norb) = (vmult_2d(2,norb)+vmult_2d(1,norb))
                do k=(norb-1),1,-(1)
                    hh_1(k)=hh_1(k+1)*(vmult_1d(2,k)+vmult_1d(1,k))
                    if(hh_1(k)==(0.0,0.0))then
                        hmin1=k
                        EXIT 
                    end if
                end do

                do k=(norb-1),1,-(1)
                    hh_2(k)=hh_2(k+1)*(vmult_2d(2,k)+vmult_2d(1,k))
                    if(hh_2(k)==(0.0,0.0))then
                        hmin2=k
                        EXIT 
                    end if
                end do

                tot1=(0.0)
                tot2=(0.0)
                if((gmax1 < hmin1).and.(gmax2 < hmin2))then
                    temp(j)=tot1
                    temp2(j)=tot2
                    cycle
                end if

                if(vec(1).ne.0) then
                    if(j.eq.1)then
                        if((annihilate1.eq.1).or.(annihilate1_2.eq.1))then
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot1 = tot1+((REAL(zom1%cos(1)*zom2%sin(1))*occupancy(2,1))*hh_1(2)*vec(1))
                                tot2 = tot2 +((REAL(zom2%sin(1)*zom1%cos(1))*occupancy(2,1))*hh_2(2)*vec(1))
                            end if
                        else
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot1 = tot1+((REAL(-zom1%sin(1)*zom2%sin(1))*occupancy(2,1))*hh_1(2)*vec(1))
                                tot2 = tot2 +((REAL(zom2%cos(1)*zom1%cos(1))*occupancy(2,1))*hh_2(2)*vec(1))
                            end if 
                        end if
                    else
                        tot1 = tot1+(REAL(conjg(z1(2,1))*z2(1,1))*hh_1(2)*vec(1))
                        tot2 = tot2+(REAL(conjg(z1(2,1))*z2(1,1))*hh_2(2)*vec(1))
                    end if
                end if
                do k=2,norb-1
                    if(vec(k).ne.0.0) then
                        if(k.eq.j)then
                            if((annihilate1.eq.k).or.(annihilate1_2.eq.k))then
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                    tot1 = tot1+(gg_1(k-1)*(REAL(zom1%cos(k)*zom2%sin(k))*occupancy(2,k))*hh_1(k+1)*vec(k))
                                    tot2 = tot2 +(gg_2(k-1)*(REAL(zom2%sin(k)*zom1%cos(k))*occupancy(2,k))*hh_2(k+1)*vec(k))
                                end if
                            else
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                    tot1 = tot1+(gg_1(k-1)*(REAL(-zom1%sin(k)*zom2%sin(k))*occupancy(2,k))*hh_1(k+1)*vec(k))
                                    tot2 = tot2 +(gg_2(k-1)*(REAL(zom2%cos(k)*zom1%cos(k))*occupancy(2,k))*hh_2(k+1)*vec(k))
                                end if 
                            end if
                        else
                            tot1 = tot1+ (gg_1(k-1)*REAL(conjg(z1(2,k))*z2(1,k))*hh_1(k+1)*vec(k))
                            tot2 = tot2+ (gg_2(k-1)*REAL(conjg(z1(2,k))*z2(1,k))*hh_2(k+1)*vec(k))
                        end if
                    end if
                end do

                if(vec(norb).ne.0) then
                    if(norb.eq.j)then
                        if((annihilate1.eq.norb).or.(annihilate1_2.eq.norb))then
                            if(annihilate2.ne.norb)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot1 = tot1+(gg_1(norb-1)*(REAL(zom1%cos(norb)*zom2%sin(norb))*occupancy(2,norb))*vec(norb))
                                tot2 = tot2 +(gg_2(norb-1)*(REAL(zom2%sin(norb)*zom1%cos(norb))*occupancy(2,norb))*vec(norb))
                            end if
                        else
                            if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot1 = tot1+(gg_1(norb-1)*(REAL(-zom1%sin(norb)*zom2%sin(norb))*occupancy(2,norb))*vec(norb))
                                tot2 = tot2 +(gg_2(norb-1)*(REAL(zom2%cos(norb)*zom1%cos(norb))*occupancy(2,norb))*vec(norb))
                            end if 
                        end if
                    else
                        tot1 = tot1 +(gg_1(norb-1)*REAL(conjg(z1(2,norb))*z2(1,norb))*vec(norb))
                        tot2 = tot2 +(gg_2(norb-1)*REAL(conjg(z1(2,norb))*z2(1,norb))*vec(norb))
                    end if
                end if

                temp(j)=tot1
                temp2(j)=tot2
                if((gmax1 < hmin1))then
                    temp(j)=0.0
                end if
                if((gmax2 < hmin2))then
                    temp2(j)=0.0
                end if
            end do
            z_an_z3_diff(1,:)=temp(:)
            z_an_z3_diff(2,:)=temp2(:)
            deallocate(vmult,stat=ierr)
            if(ierr==0) deallocate(vmult_1d,stat=ierr)
            if(ierr==0) deallocate(vmult_2d,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult deallocation . ierr had value ", ierr
                errorflag=1
                return
            end if 

        end if

        return

    end function z_an_z3_diff
        











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