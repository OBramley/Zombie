    MODULE vals 

        type zombiest
            real(kind=8), dimension(10)::sin
            real(kind=8), dimension(10)::cos
            real(kind=8),dimension(10)::phi
            real(kind=8),dimension(10)::img
            real(kind=8),dimension(0:20)::val
            integer(kind=1),dimension(10)::dead
            integer(kind=1),dimension(10)::alive
            integer::update_num
        end type zombiest


        type oprts
            integer(kind=2),dimension(:,:),allocatable::alive,dead
            integer(kind=1),dimension(:,:),allocatable::neg_alive,neg_dead

            integer(kind=2),dimension(:,:,:),allocatable::alive_diff,dead_diff
            integer(kind=1),dimension(:,:,:),allocatable::neg_alive_diff,neg_dead_diff
            integer(kind=2),dimension(:,:),allocatable::dcnt
        end type oprts


        ! Type defining the 1&2 electron integrals
        type elecintrgl
            integer::h1_num
            integer::h2_num
            real(kind=8),allocatable,  dimension(:,:)::h1ei_old
            real(kind=8), allocatable,dimension(:,:,:,:)::h2ei_old
            real(kind=8), allocatable,dimension(:)::h1ei
            real(kind=8), allocatable,dimension(:)::h2ei
            real(kind=8),allocatable, dimension(:)::h2ei_grad
            real(kind=8) :: hnuc
            
        end type elecintrgl

        integer::ndet       ! Number of Zombie states
        integer::norb       ! Number of spin orbitals
        integer::nel        ! Number of electrons in molecule
        integer::errorflag
        character(LEN=1)::GDflg      
        real(kind=8)::pirl      ! pi
        real(kind=8) :: sqrtpi  ! square root of pi      
        
        contains 
        subroutine init()
            implicit none
            ndet=10
            norb=10
            nel=6
            errorflag=0
            GDflg='y'
            sqrtpi = 1.7724538509055160272981674833411451827975494561223871d0
            pirl = sqrtpi**2.0d0
            
        end subroutine

    End Module vals 

    Module stuff 

        use vals
        contains
    
        ! calculates indvidual hamiltonian elements taking in two Zombie states and a set of 
        ! creation and annihilation operations
        real(kind=8) function haml_vals(z1d,z2d,ops,el,len)

        implicit none 
        real(kind=8),dimension(0:),intent(in)::z1d,z2d
        real(kind=8),dimension(:),intent(in)::el
        integer,intent(in)::len
        type(oprts),intent(in)::ops
        real(kind=8)::ov
        integer::j,k

        
        haml_vals=0.0
        !$omp parallel private(j,k,ov) shared(ops,z1d,z2d)
        !$omp do simd reduction(+:haml_vals) 
        do j=1,len
            ov=1.0
            !!$omp do simd reduction(*:ov)
            do k=1, norb
                ov=ov*((z1d(k)*z2d(ops%alive(k,j))*ops%neg_alive(k,j))+(z1d(k+norb)*z2d(ops%dead(k,j))*ops%neg_dead(k,j)))
            end do
            ! if(j==21)then
            !     print*,ops%alive(:,j)
            !     print*,ops%dead(:,j)
            !     print*,ops%neg_alive(:,j)
            !     print*,ops%neg_dead(:,j)
            ! end if
            ! print*,ov
            print*,el(j)*ov
            !!$omp end do simd
            haml_vals=haml_vals+(ov*el(j))
            ! haml_vals=haml_vals+(el(j))
        end do
        !$omp end do simd
        !$omp end parallel 
        
        return 
    
        end function haml_vals
        
        
        ! calculates indvidual hamliltonian elements taking in two Zombie states and a set of 
        ! creation and annihilation operations
        real(kind=8) function haml_val_grad(z1d,z2d,ops,el,orb)

        implicit none 
        real(kind=8),dimension(0:),intent(in)::z1d,z2d
        real(kind=8),dimension(:),intent(in)::el
        integer,intent(in)::orb
        type(oprts),intent(in)::ops
        real(kind=8)::ov
        integer::j,k,len

        len=ops%dcnt(0,orb)
        haml_val_grad=0.0
        !$omp parallel private(j,k,ov) shared(ops,z1d,z2d)
        !$omp do simd reduction(+:haml_val_grad) 
        do j=1,len
            ov=1.0
            !!$omp do simd reduction(*:ov)
            do k=1, norb
                ov=ov*((z1d(k)*z2d(ops%alive_diff(orb,k,(ops%dcnt(j,orb))))*ops%neg_alive_diff(orb,k,(ops%dcnt(j,orb))))&
                +(z1d(k+norb)*z2d(ops%dead_diff(orb,k,(ops%dcnt(j,orb))))*ops%neg_dead_diff(orb,k,(ops%dcnt(j,orb))))) 
            end do
            !!$omp end do simd
            haml_val_grad=haml_val_grad+(ov*el(ops%dcnt(j,orb)))
        end do
        !$omp end do simd
        !$omp end parallel 
        
        return 
    
    end function haml_val_grad

    subroutine haml_grad_rc(hcol,z1d,zstore,an_cr,an2_cr2,an2_cr2_diff,elecs,state,orb)

        implicit none
        real(kind=8),dimension(:),intent(inout)::hcol 
        type(zombiest),dimension(:),intent(in)::zstore
        real(kind=8),dimension(0:),intent(in)::z1d
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2,an2_cr2_diff
        integer,intent(in)::state,orb
        real(kind=8)::h1etot,h2etot
        integer::j
        

        !$omp parallel 
        !$omp single
        do j=1,1!ndet
            if(j.ne.state)then
                ! ! Differentiating the bra 1 el
                ! !$omp task firstprivate(h1etot,j) shared(hcol,zstore,an_cr,elecs,z1d)
                ! h1etot = haml_vals(z1d,zstore(j)%val,an_cr,elecs%h1ei,elecs%h1_num)
                ! !$omp atomic
                ! hcol(j)=hcol(j)+h1etot
                ! !$omp end atomic
                ! !$omp end task
                !Differentiating the bra 2 el
                !$omp task firstprivate(h2etot,j) shared(hcol,zstore,an2_cr2,elecs,z1d)
                h2etot = haml_vals(z1d,zstore(j)%val,an2_cr2,elecs%h2ei,elecs%h2_num)
                !$omp atomic
                hcol(j)=hcol(j)+(0.5*h2etot)
                !$omp end atomic
                !$omp end task

            else
                ! !Differentiaitn hamiltonian element (a,a) only placed in hamiltonian column
                ! !$omp task firstprivate(h1etot,j) shared(hcol,zstore,an_cr,elecs)
                ! h1etot = haml_val_grad(zstore(state)%val,zstore(state)%val,an_cr,elecs%h1ei,orb)
                ! !$omp atomic
                ! hcol(j)=hcol(j)+h1etot
                ! !$omp end atomic
                ! !$omp end task
                !$omp task firstprivate(h2etot,j) shared(hcol,zstore,an2_cr2,elecs)
                h2etot = haml_val_grad(zstore(state)%val,zstore(state)%val,an2_cr2,elecs%h2ei,orb)
                !$omp atomic
                hcol(j)=hcol(j)+(0.5*h2etot)
                !$omp end atomic
                !$omp end task
            end if
        end do 
        !$omp end single
        !$omp end parallel  
        return

    end subroutine haml_grad_rc

! Hamiltonian calcualtion - calcualtes the gradient of ther hamliltonian w.r.t one zombie state 
    subroutine haml_grad(haml_diff,zstore,elecs,an_cr,an2_cr2,an2_cr2_diff,state) 

        implicit none
        
        real(kind=8),dimension(:,:),intent(inout)::haml_diff 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2,an2_cr2_diff
        integer,intent(in)::state
        real(kind=8),dimension(0:2*norb)::z1d
  
        integer::j,ierr
    
        if (errorflag .ne. 0) return 
        ierr=0

        !$omp parallel do &
        !$omp & private(j) &
        !$omp & shared(elecs,zstore,an_cr,an2_cr2,haml_diff)
        do j=2,2!norb
            z1d(0:2*norb)=zstore(state)%val(0:2*norb)
            z1d(j)=zstore(state)%cos(j)
            z1d(j+norb)=zstore(state)%sin(j)*(-1)
            call haml_grad_rc(haml_diff(:,j),z1d,zstore,an_cr,an2_cr2,an2_cr2_diff,elecs,state,j)
        end do
        !$omp end parallel do

        print*,haml_diff(1,:)
        
    end subroutine haml_grad

    subroutine two_elec_part_grad(zs1sin,zs1cos,zs2sin,zs2cos,z2l,z1jk,h2ei,occupancy_2an,occupancy_an,elecs)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        type(elecintrgl)::elecs
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:,:),intent(in)::z1jk
        real(kind=8), dimension(:,:,:,:), intent(in)::h2ei
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an
        integer,dimension(:,:,:),intent(in)::occupancy_an
        real(kind=8),dimension(norb,norb,norb,norb)::h2etot_diff
        integer,allocatable,dimension(:,:)::occupancy
        real(kind=8),dimension(norb)::tem
        integer::k,l,j,jspin,equal,p,cnt

        allocate(occupancy(2,norb))
        tem=1
        h2etot_diff=0.0
        equal=2
        cnt=1
        ! print*,zs1sin
        ! print*,zs1cos
        ! print*,real(zs2sin)
        ! print*,real(zs2cos)
        do j=1, norb

            if(modulo(j,2)==0)then
                jspin=2
            else
                jspin=1
            end if
            do k=1, norb
                if(j.eq.k)then
                    cycle
                end if
                do l=jspin, norb, 2
                    ! do p=1,norb 
                        ! if(elecs%h2ei_old(j,k,l,p).ne.0)then
                        !     print*,j,k,l,p
                        !     ! print*,elecs%h2ei_old(j,k,l,p),elecs%h2ei(cnt)
                        !     cnt=cnt+1
                        ! end if 
                    ! end do 
                    ! print*,j,k,l
                    occupancy=occupancy_2an(j,k,:,:)*occupancy_an(l,:,:)
                    
                    h2etot_diff(j,k,l,:) = z_an_z3_diff_2(norb,z1jk(j,k,:,:),z2l(l,:,:),h2ei(j,k,l,:),&     
                    2,occupancy,j,k,l,zs1sin,zs1cos,zs2sin,zs2cos)
                end do
            end do
        end do
        ! stop
        ! print*,cnt

        tem=0.0
        do j=1,norb
            do k=1, norb 
                do l=1,norb
                    tem=tem+h2etot_diff(j,k,l,:) 
                end do 
            end do 
        end do
        tem=tem*0.5
        print*,tem
        deallocate(occupancy)
        
    

    end subroutine two_elec_part_grad

    function z_an_z3_diff_2(len,z1,z2,vec,dtype,occupancy,annihilate1,annihilate1_2,annihilate2,zs1sin,zs1cos,zs2sin,zs2cos)

        implicit none
        integer, intent(in)::len
        complex(kind=8),dimension(:,:),intent(in)::z1,z2
        real(kind=8),dimension(:,:),allocatable::vmult_dd,vmult_1d, vmult_2d,vmult
        real(kind=8),dimension(len),intent(in)::vec
        integer,dimension(:,:),intent(in)::occupancy
        integer, intent(in)::dtype,annihilate1,annihilate1_2,annihilate2
        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        real(kind=8),dimension(len)::z_an_z3_diff_2
        real(kind=8),allocatable,dimension(:)::temp,temp2,temp0
        real(kind=8),allocatable,dimension(:)::gg_1,hh_1,gg_2,hh_2,gg_0,hh_0
        real(kind=8)::tot1, tot2,tot0
        integer::j,k,gmax1,hmin1,gmax2,hmin2,ierr,breakflag,gmax0,hmin0

        if (errorflag .ne. 0) return
        ierr=0

        allocate(vmult(2,len),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in vmult allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        vmult=real(conjg(z1)*z2)
        
        z_an_z3_diff_2=0
        
        if(dtype.eq.1) then !Differentiation when zombie states are the same
            allocate(temp0(len))
            allocate(gg_0(len),stat=ierr)
            allocate(hh_0(len),stat=ierr)
            allocate(vmult_dd(2,len),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult_dd allocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
            temp0=0
            !$omp parallel do &
            !$omp & private(vmult_dd,gmax0,hmin0,breakflag,tot0,j,k) &
            !$omp shared(annihilate1,annihilate1_2,annihilate2,temp0,zs1sin,zs1cos,vmult,occupancy,vec,z1,z2)
            do j=1, len    !Differentiating w.r.t to orbital j
                
                breakflag=0
                vmult_dd=vmult
                if(annihilate1.eq.annihilate1_2)then
                    temp0(j)=0
                    cycle
                else if((annihilate1.eq.j).or.(annihilate1_2.eq.j))then
                    if(annihilate2.eq.j)then    !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(1j)
                        vmult_dd(1,j)=0
                        vmult_dd(2,j)=2*real(zs1sin(j)*zs1cos(j))*occupancy(2,j) !(sin2x=2sinxcosx)
                    else
                        vmult_dd(1,j)=0        !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(0j)
                        vmult_dd(2,j)=(1.0-2*((real(zs1sin(j)))**2))*occupancy(2,j)    !(cos2x = cos^x -sin^2x)
                    end if
                else if((annihilate1.ne.j).or.(annihilate1_2.ne.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(0j)*a^(a)_(1j)
                        vmult_dd(1,j)=0
                        vmult_dd(2,j)=(1.0-2*((real(zs1sin(j)))**2))*occupancy(2,j)
                    else
                        breakflag=j !before diff alive:a^(a)_(1j)*a^(a)_(1j) dead:a^(a)_(0j)*a^(a)_(0j)
                    end if          !Unless an opeator acts at position j this evaluates to 0
                end if
            
                gg_0(1:len)=(0.0,0.0)
                hh_0(1:len)=(0.0,0.0)
                gmax0=len
                gg_0(1)=vmult_dd(2,1)-vmult_dd(1,1)

                do k=2, len
                    gg_0(k)=gg_0(k-1)*(vmult_dd(2,k)-vmult_dd(1,k))
                    if(gg_0(k)==(0.0))then
                        gmax0=k
                        EXIT 
                    end if
                end do
                
                hmin0=0
                hh_0(len) = vmult_dd(2,len)+vmult_dd(1,len)
                do k=(len-1),1,-(1)
                    hh_0(k)=hh_0(k+1)*(vmult_dd(2,k)+vmult_dd(1,k))
                    if(hh_0(k)==(0.0))then
                        hmin0=k
                        EXIT 
                    end if
                end do

                tot0=(0.0)
                if (gmax0 < hmin0) then
                    temp0(j)=tot0
                    cycle
                end if

                if((breakflag.eq.0).or.(breakflag.eq.1))then
                    if(vec(1).ne.0) then
                        if(breakflag.eq.1)then
                            temp0(j)=(1.0-2*((real(zs1sin(1)))**2))*occupancy(2,1)*hh_0(2)*vec(1)
                            
                            cycle
                        end if
                        if(j.eq.1)then
                            if((annihilate1.eq.1).or.(annihilate1_2.eq.1))then
                                if(annihilate2.ne.1)then
                                    tot0 = tot0+(2*real(zs1sin(1)*zs1cos(1))*occupancy(2,1)*hh_0(2)*vec(1))
                                end if
                            else
                                if(annihilate2.ne.1)then
                                    tot0 = tot0+((1.0-2*((real(zs1sin(1)))**2))*occupancy(2,1)*hh_0(2)*vec(1))
                                end if 
                            end if
                        else
                            tot0 = tot0+(REAL(conjg(z1(2,1))*z2(1,1))*hh_0(2)*vec(1))
                        end if
                    end if
                else if((breakflag.lt.len))then
                    if(breakflag.ne.0)then
                        temp0(j)=(gg_0(j-1)*((1.0-2*((real(zs1sin(j)))**2))*occupancy(2,j))*hh_0(j+1)*vec(j))
                       
                        cycle
                    end if
                    do k=2,len-1
                        if(vec(k).ne.0.0) then
                            if(k.eq.j)then
                                if((annihilate1.eq.k).or.(annihilate1_2.eq.k))then
                                    if(annihilate2.ne.k)then
                                        tot0 = tot0+(gg_0(k-1)*(2*real(zs1sin(k)*zs1cos(k))*occupancy(2,k))*hh_0(k+1)*vec(k))
                                    end if
                                else
                                    if(annihilate2.ne.k)then
                                        tot0 = tot0+(gg_0(k-1)*((1.0-2*((real(zs1sin(k)))**2))*occupancy(2,k))*hh_0(k+1)*vec(k))
                                    end if 
                                end if
                            else
                                tot0 = tot0+ (gg_0(k-1)*REAL(conjg(z1(2,k))*z2(1,k))*hh_0(k+1)*vec(k))
                            end if
                        end if
                    end do
                else
                    if(breakflag.eq.len)then
                        temp0(j)=(gg_0(len-1)*((1.0-2*((real(zs1sin(len)))**2))*occupancy(2,len))*vec(len))
                      
                        cycle
                    end if
                    if(vec(len).ne.0) then
                        if(len.eq.j)then
                            if((annihilate1.eq.len).or.(annihilate1_2.eq.len))then
                                if(annihilate2.ne.len)then
                                    tot0 = tot0 +(gg_0(len-1)*(2*real(zs1sin(len)*zs1cos(len))*occupancy(2,len))*vec(len))
                                end if
                            else
                                if(annihilate2.ne.len)then
                                    tot0 = tot0+ (gg_0(len-1)*((1.0-2*((real(zs1sin(len)))**2))*occupancy(2,len))*vec(len))
                                end if    
                            end if
                        else
                            tot0 = tot0 +(gg_0(len-1)*REAL(conjg(z1(2,len))*z2(1,len))*vec(len))
                        end if
                    end if
                end if

                temp0(j)=tot0
            end do
            !$omp end parallel do

            z_an_z3_diff_2=temp0(:)
            z_an_z3_diff_2=temp0(:)
            deallocate(temp0,stat=ierr)
            deallocate(gg_0,stat=ierr)
            deallocate(hh_0,stat=ierr) 
            deallocate(vmult,stat=ierr)
            if(ierr==0) deallocate(vmult_dd,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult deallocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
        
        else if(dtype.eq.2)then !Differentiation when zombie states are not the same
            allocate(temp(len),stat=ierr)
            allocate(gg_1(len),stat=ierr)
            allocate(hh_1(len),stat=ierr)
            allocate(vmult_1d(2,len),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult_1d/vmult_2d allocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
            temp=0
            !$omp  parallel do &
            !$omp & private(vmult_1d,gmax1,hmin1,tot1,j,k) shared(annihilate1,&
            !$omp annihilate1_2,annihilate2,temp,zs1sin,zs1cos,zs2sin,zs2cos,vmult,occupancy,vec,z1,z2,len)
           
            do j=2,2!len
                
                vmult_1d=vmult
                if(annihilate1.eq.annihilate1_2)then
                    temp(j)=0
                    cycle
                else if((annihilate1.eq.j).or.(annihilate1_2.eq.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                        vmult_1d(1,j)=0
                        vmult_1d(2,j)=real(zs1cos(j)*zs2sin(j))*occupancy(2,j)
                        
                    else !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(0j)
                        vmult_1d(1,j)=0
                        vmult_1d(2,j)=real(zs1cos(j)*zs2cos(j))*occupancy(2,j)
                    end if
                else if((annihilate1.ne.j).or.(annihilate1_2.ne.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                        vmult_1d(1,j)=0
                        vmult_1d(2,j)=-real(zs1sin(j)*zs2sin(j))*occupancy(2,j)
                    else    !before diff alive:a^(a)_(1j)*a^(b)_(1j) dead:a^(a)_(0j)*a^(b)_(0j)
                        vmult_1d(1,j)=real(zs1cos(j)*zs2sin(j))*occupancy(1,j)
                        vmult_1d(2,j)=-real(zs1sin(j)*zs2cos(j))*occupancy(2,j)
                    end if
                end if
                
                gg_1(1:len)=(0.0,0.0)
                hh_1(1:len)=(0.0,0.0)
                gmax1=len
                gg_1(1)=(vmult_1d(2,1))-(vmult_1d(1,1))

                do k=2, len
                    gg_1(k)=gg_1(k-1)*((vmult_1d(2,k)-vmult_1d(1,k)))
                    if(gg_1(k)==(0.0,0.0))then
                        gmax1=k
                        EXIT 
                    end if
                end do
                
                hmin1=0
                
                hh_1(len) = (vmult_1d(2,len)+vmult_1d(1,len))
                
                do k=(len-1),1,-(1)
                    hh_1(k)=hh_1(k+1)*(vmult_1d(2,k)+vmult_1d(1,k))
                    if(hh_1(k)==(0.0,0.0))then
                        hmin1=k
                        EXIT 
                    end if
                end do

                tot1=(0.0)
                
                if((gmax1 < hmin1))then
                    temp(j)=tot1
                    cycle
                end if
                ! if((annihilate1.eq.1).and.(annihilate1_2.eq.3).and.(annihilate2.eq.3))then
                !     print*,occupancy(1,:)
                !     print*,occupancy(2,:)
                !     print*,vmult_1d(1,:)
                !     print*,vmult_1d(2,:)
                ! end if
                if(vec(1).ne.0) then
                !    print*,annihilate1,annihilate1_2,annihilate2,1
                    if(j.eq.1)then
                        if((annihilate1.eq.1).or.(annihilate1_2.eq.1))then
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot1 = tot1+((REAL(zs1cos(1)*zs2sin(1))*occupancy(2,1))*hh_1(2)*vec(1))
                        print*, annihilate1,annihilate1_2,annihilate2,1,((REAL(zs1cos(1)*zs2sin(1))*occupancy(2,1))*hh_1(2)*vec(1))
                            end if
                        else
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot1 = tot1+((REAL(-zs1sin(1)*zs2sin(1))*occupancy(2,1))*hh_1(2)*vec(1))
                        print*,annihilate1,annihilate1_2,annihilate2,1, ((REAL(-zs1sin(1)*zs2sin(1))*occupancy(2,1))*hh_1(2)*vec(1))
                            end if 
                        end if
                    else
                        tot1 = tot1+(REAL(conjg(z1(2,1))*z2(1,1))*hh_1(2)*vec(1))
                        print*,annihilate1,annihilate1_2,annihilate2,1,(REAL(conjg(z1(2,1))*z2(1,1))*hh_1(2)*vec(1))
                    end if
                    ! tot1=tot1+vec(1)
                   
                end if
                do k=2,len-1
                    if(vec(k).ne.0.0) then
                        ! print*,annihilate1,annihilate1_2,annihilate2,k
                        if(k.eq.j)then
                            if((annihilate1.eq.k).or.(annihilate1_2.eq.k))then
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                    tot1 = tot1+(gg_1(k-1)*(REAL(zs1cos(k)*zs2sin(k))*occupancy(2,k))*hh_1(k+1)*vec(k))
            print*, annihilate1,annihilate1_2,annihilate2,k,(gg_1(k-1)*(REAL(zs1cos(k)*zs2sin(k))*occupancy(2,k))*hh_1(k+1)*vec(k))
                                end if
                            else
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                    tot1 = tot1+(gg_1(k-1)*(REAL(-zs1sin(k)*zs2sin(k))*occupancy(2,k))*hh_1(k+1)*vec(k))
        print*,  annihilate1,annihilate1_2,annihilate2,k,(gg_1(k-1)*(REAL(-zs1sin(k)*zs2sin(k))*occupancy(2,k))*hh_1(k+1)*vec(k))
                                end if 
                            end if
                        else
                            tot1 = tot1+ (gg_1(k-1)*REAL(conjg(z1(2,k))*z2(1,k))*hh_1(k+1)*vec(k))
                    print*,  annihilate1,annihilate1_2,annihilate2,k,(gg_1(k-1)*REAL(conjg(z1(2,k))*z2(1,k))*hh_1(k+1)*vec(k))
                        end if
                        ! tot1=tot1+vec(k)
                    end if
                    
                end do

                if(vec(len).ne.0) then
                    ! print*,annihilate1,annihilate1_2,annihilate2,len
                    if(len.eq.j)then
                        if((annihilate1.eq.len).or.(annihilate1_2.eq.len))then
                            if(annihilate2.ne.len)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot1 = tot1+(gg_1(len-1)*(REAL(zs1cos(len)*zs2sin(len))*occupancy(2,len))*vec(len))
        print*, annihilate1,annihilate1_2,annihilate2,len,(gg_1(len-1)*(REAL(zs1cos(len)*zs2sin(len))*occupancy(2,len))*vec(len))
                            end if
                        else
                            if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot1 = tot1+(gg_1(len-1)*(REAL(-zs1sin(len)*zs2sin(len))*occupancy(2,len))*vec(len))
        print*, annihilate1,annihilate1_2,annihilate2,len, (gg_1(len-1)*(REAL(-zs1sin(len)*zs2sin(len))*occupancy(2,len))*vec(len))
                            end if 
                        end if
                    else
                        tot1 = tot1 +(gg_1(len-1)*REAL(conjg(z1(2,len))*z2(1,len))*vec(len))
                        print*, annihilate1,annihilate1_2,annihilate2,len,(gg_1(len-1)*REAL(conjg(z1(2,len))*z2(1,len))*vec(len))
                    end if
                !    tot1=tot1+vec(len)
                end if
                
                temp(j)=tot1
            end do
            ! print*,temp
            !$omp end parallel do
            z_an_z3_diff_2=temp(:)
            deallocate(temp,stat=ierr)
            deallocate(gg_1,stat=ierr)
            deallocate(hh_1,stat=ierr) 
            deallocate(vmult,stat=ierr)
            if(ierr==0) deallocate(vmult_1d,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult deallocation . ierr had value ", ierr
                errorflag=1
                return
            end if 

            
        else if(dtype.eq.3)then
            allocate(temp2(len),stat=ierr)
            allocate(gg_2(len),stat=ierr)
            allocate(hh_2(len),stat=ierr) 
            allocate(vmult_2d(2,len),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult_1d/vmult_2d allocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
            temp2=0
            !$omp parallel do &
            !$omp & private(vmult_2d,gmax2,hmin2,tot2,j,k) shared(annihilate1,&
            !$omp annihilate1_2,annihilate2,temp2,zs1sin,zs1cos,zs2sin,zs2cos,vmult,occupancy,vec,z1,z2)
            do j=1,len
                
                vmult_2d=vmult
                if(annihilate1.eq.annihilate1_2)then
                    temp2(j)=0
                    cycle
                else if((annihilate1.eq.j).or.(annihilate1_2.eq.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                        vmult_2d(1,j)=0
                        vmult_2d(2,j)=real(zs1sin(j)*zs2cos(j))*occupancy(2,j)
                    else !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(0j)
                        vmult_2d(1,j)=0
                        vmult_2d(2,j)=-real(zs1sin(j)*zs2sin(j))*occupancy(2,j)
                    end if
                else if((annihilate1.ne.j).or.(annihilate1_2.ne.j))then
                    if(annihilate2.eq.j)then !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                        vmult_2d(1,j)=0
                        vmult_2d(2,j)=real(zs1cos(j)*zs2cos(j))*occupancy(2,j)
                    else    !before diff alive:a^(a)_(1j)*a^(b)_(1j) dead:a^(a)_(0j)*a^(b)_(0j)
                        vmult_2d(1,j)=real(zs1sin(j)*zs2cos(j))*occupancy(1,j)
                        vmult_2d(2,j)=-real(zs1cos(j)*zs2sin(j))*occupancy(2,j)
                    end if
                end if

                gg_2(1:len)=(0.0,0.0)
                hh_2(1:len)=(0.0,0.0)
                gmax2=len
                gg_2(1)=(vmult_2d(2,1))-(vmult_2d(1,1))

                do k=2, len
                    gg_2(k)=gg_2(k-1)*((vmult_2d(2,k)-vmult_2d(1,k)))
                    if(gg_2(k)==(0.0,0.0))then
                        gmax2=k
                        EXIT 
                    end if
                end do
                

                hmin2=0
                hh_2(len) = (vmult_2d(2,len)+vmult_2d(1,len))
        
                do k=(len-1),1,-(1)
                    hh_2(k)=hh_2(k+1)*(vmult_2d(2,k)+vmult_2d(1,k))
                    if(hh_2(k)==(0.0,0.0))then
                        hmin2=k
                        EXIT 
                    end if
                end do

    
                tot2=(0.0)
                if((gmax2 < hmin2))then
                    temp2(j)=tot2
                    cycle
                end if

                if(vec(1).ne.0) then
                   
                    if(j.eq.1)then
                        if((annihilate1.eq.1).or.(annihilate1_2.eq.1))then
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot2 = tot2 +((REAL(zs2sin(1)*zs1cos(1))*occupancy(2,1))*hh_2(2)*vec(1))
                            end if
                        else
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot2 = tot2 +((REAL(zs2cos(1)*zs1cos(1))*occupancy(2,1))*hh_2(2)*vec(1))
                            end if 
                        end if
                    else
                        tot2 = tot2+(REAL(conjg(z1(2,1))*z2(1,1))*hh_2(2)*vec(1))
                    end if
                end if
                do k=2,len-1
                    if(vec(k).ne.0.0) then
                       
                        if(k.eq.j)then
                            if((annihilate1.eq.k).or.(annihilate1_2.eq.k))then
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                    tot2 = tot2 +(gg_2(k-1)*(REAL(zs2sin(k)*zs1cos(k))*occupancy(2,k))*hh_2(k+1)*vec(k))
                                end if
                            else
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                    tot2 = tot2 +(gg_2(k-1)*(REAL(zs2cos(k)*zs1cos(k))*occupancy(2,k))*hh_2(k+1)*vec(k))
                                end if 
                            end if
                        else
                            tot2 = tot2+ (gg_2(k-1)*REAL(conjg(z1(2,k))*z2(1,k))*hh_2(k+1)*vec(k))
                        end if
                    end if
                end do

                if(vec(len).ne.0) then
                    
                    if(len.eq.j)then
                        if((annihilate1.eq.len).or.(annihilate1_2.eq.len))then
                            if(annihilate2.ne.len)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot2 = tot2 +(gg_2(len-1)*(REAL(zs2sin(len)*zs1cos(len))*occupancy(2,len))*vec(len))
                            end if
                        else
                            if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot2 = tot2 +(gg_2(len-1)*(REAL(zs2cos(len)*zs1cos(len))*occupancy(2,len))*vec(len))
                            end if 
                        end if
                    else
                        tot2 = tot2 +(gg_2(len-1)*REAL(conjg(z1(2,len))*z2(1,len))*vec(len))
                    end if
                end if

                temp2(j)=tot2
                
            end do
            !$omp end parallel do
        
            z_an_z3_diff_2=temp2(:)
            deallocate(temp2,stat=ierr)
            deallocate(gg_2,stat=ierr)
            deallocate(hh_2,stat=ierr) 
            deallocate(vmult,stat=ierr)
            if(ierr==0) deallocate(vmult_2d,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in vmult deallocation . ierr had value ", ierr
                errorflag=1
                return
            end if 

        end if

        return

    end function z_an_z3_diff_2

    subroutine electronintegrals_old(elecs)

        implicit none
    
        type(elecintrgl), intent(inout)::elecs
    
        integer:: ierr, nlines
        
    
        if (errorflag .ne. 0) return
        ierr = 0
    
        nlines = lines(nlines)
        call spattospin1(elecs,nlines)
        call spattospin2(elecs,nlines)
    
        open(unit=128, file='integrals/hnuc.csv',status='old',iostat=ierr)
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
        
        write(6,"(a)") "1 & 2 electron integrals successfully generated"
    
        return
    
    end subroutine electronintegrals_old
    
    integer function lines(nlines)
        implicit none
    
        integer, intent(INOUT):: nlines
        integer:: ierr
        
        ierr=0
        nlines=0
        open(unit=129, file='integrals/h1ea.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a,i0)") 'Error in opening h1ea.csv file',ierr
            errorflag = 1
            return
        end if
    
        do 
            read(129,*, iostat=ierr)
            if(ierr<0)then
                ! write(0,"(a,i0)") "nlines has value ", nlines
                lines=nlines
                close(129)
                return
            else if (ierr/=0) then
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
    
    
        open(unit=130, file='integrals/h1ea.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a)") 'Error in opening h1ea.csv file'
            errorflag = 1
            return
        end if
      
        do j=1, nlines
            read(130,*,iostat=ierr) (h1ea(j,k),k=1,nlines)
        end do
    
        close(130)
    
    
        nspao = int(norb/2)
    
        !$omp parallel shared(elecs,h1ea) private(j,k,jj,kk)
        !$omp do 
        do j=1,nspao
            do k=1, nspao
                jj=(j*2)-1
                kk=(k*2)-1
                ! if(h1ea(j,k).ne.0) h1ea(j,k)=1
                elecs%h1ei_old(jj,kk)=h1ea(j,k)
                elecs%h1ei_old(jj+1,kk+1)=elecs%h1ei_old(jj,kk)
            end do
        end do
        !$omp end do
        !$omp end parallel
        
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
        do j=1,nlines
            do k=1, nlines
                write(val1,'(i0)')j
                write(val2,'(i0)')k
                open(unit=(131+j+k), file='integrals/h2ea_'//trim(val1)//'_'//trim(val2)//'.csv', status='old',iostat=ierr)
                if (ierr.ne.0) then
                    write(0,"(a)") 'Error in opening h2ea_'//trim(val1)//'_'//trim(val2)//'.csv file'
                    errorflag = 1
                    return
                end if
                do l=1, nlines
                    read((131+j+k),*) (h2ea(j,k,l,m),m=1,nlines)
                end do
                close(131+j+k)
            end do
        end do
    
        !$omp parallel shared(elecs,h2ea) private(j,k,l,m,jj,kk,ll,mm) 
        !$omp do
        do j=1,nspao
            jj=(2*j)-1
            do k=1, nspao
                kk=(2*k)-1
                do l=1, nspao
                    ll=(2*l)-1
                    do m=1, nspao
                        mm=(2*m)-1
                        ! if(h2ea(j,k,l,m).ne.0) h2ea(j,k,l,m)=1
                        elecs%h2ei_old(jj,ll,kk,mm)=h2ea(j,k,l,m)
                        elecs%h2ei_old(jj+1,ll,kk+1,mm)=h2ea(j,k,l,m)
                        elecs%h2ei_old(jj,ll+1,kk,mm+1)=h2ea(j,k,l,m)
                        elecs%h2ei_old(jj+1,ll+1,kk+1,mm+1)=h2ea(j,k,l,m)
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
       
    
        deallocate(h2ea, stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in h2ea deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if
    
        return
                        
    end subroutine spattospin2


    subroutine alloc_oprts(oper,n)

        implicit none 
        type(oprts),intent(inout)::oper 
        integer,intent(in)::n 
        integer::ierr,k
        integer(kind=2)::j,norbs
        if (errorflag .ne. 0) return

        ierr=0
        norbs=int(norb,kind=2)
        allocate(oper%alive(norb,n),oper%dead(norb,n),oper%neg_alive(norb,n),oper%neg_dead(norb,n),stat=ierr)
        if(GDflg.eq.'y')then 
            if(ierr==0) allocate(oper%alive_diff(norb,norb,n),oper%dead_diff(norb,norb,n),stat=ierr)
            if(ierr==0) allocate(oper%neg_alive_diff(norb,norb,n),oper%neg_dead_diff(norb,norb,n),stat=ierr)
            if(ierr==0) allocate(oper%dcnt(0:n,norb),stat=ierr)
        end if
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in operators allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        oper%neg_alive=1
        oper%neg_dead=1
        oper%alive=0
        oper%dead=0
        do k=1,n
            oper%alive(:,k)=[integer(kind=2)::(j,j=1,norbs)]
            oper%dead(:,k)=[integer(kind=2)::((j+norbs),j=1,norbs)]
        end do
        if(GDflg.eq.'y')then
            oper%dcnt=0
            oper%neg_alive_diff=1
            oper%neg_dead_diff=1
            do j=1,norbs
                oper%alive_diff(j,:,:)=oper%alive
                oper%dead_diff(j,:,:)=oper%dead
            end do
        end if 



    end subroutine alloc_oprts


    subroutine electronintegrals(elecs,an_cr,an2_cr2,an2_cr2_diff)

        implicit none
    
        type(elecintrgl), intent(inout)::elecs
        type(oprts),intent(inout)::an_cr,an2_cr2,an2_cr2_diff
        real::e1in,e2in
        integer:: ierr,e1,e2
        
    
        if (errorflag .ne. 0) return
    
    
        ierr = 0
        open(unit=129, file='integrals/h1e.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a)") 'Error in opening h1ec.csv file'
            errorflag = 1
            return
        end if
        read(129,*, iostat=ierr)e1in
        if (ierr .ne. 0) then
            write(0,"(a)") 'Error reading h1e'
            errorflag = 1
            return
          end if
        close(129)
    
        e1=int(e1in)
    
        open(unit=131, file='integrals/h2e.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a)") 'Error in opening h2e.csv file'
            errorflag = 1
            return
        end if
        read(131,*, iostat=ierr)e2in
        if (ierr .ne. 0) then
            write(0,"(a)") 'Error reading h1e'
            errorflag = 1
            return
          end if
        close(131)
    
        e2=int(e2in)
        ! e2=1243 !int(e2in)
        
        elecs%h1_num=e1
        elecs%h2_num=e2
        allocate(elecs%h1ei(e1))
        allocate(elecs%h2ei(e2))
        allocate(elecs%h2ei_grad(e2))
        
        call  one_electrons(elecs,an_cr,e1)
        call  two_electrons(elecs,an2_cr2,an2_cr2_diff,e2)
    
    
        open(unit=128, file='integrals/hnuc.csv',status='old',iostat=ierr)
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
        
        write(6,"(a)") "1 & 2 electron integrals successfully generated"
    
        return
    
    end subroutine electronintegrals
    
    ! 
    subroutine two_electrons(elecs,an2_cr2,an2_cr2_diff,e2)
    
        implicit none 
        type(elecintrgl), intent(inout)::elecs
        type(oprts),intent(inout)::an2_cr2,an2_cr2_diff
        integer,intent(in)::e2
        real(kind=8),dimension(e2+1,5)::read_in
        integer::ierr,l,k,an1,an2,cr1,cr2,j,jspin,cnt,p
        
        
        
        open(unit=131, file='integrals/h2e.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a)") 'Error in opening hnuc.csv file'
            errorflag = 1
            return
        end if
       
        do j=1,e2+1
            read(131,*) (read_in(j,k),k=1,5)
        end do 
        
    
        close(131)
        elecs%h2ei=read_in(2:,1)
        
        call alloc_oprts(an2_cr2,e2)
        ! call alloc_oprts(an2_cr2,1243)

        ! cnt=1
        ! do j=1, norb
        !     if(modulo(j,2)==0)then
        !         jspin=2
        !     else
        !         jspin=1
        !     end if
        !     do k=1, norb
        !         if(j.eq.k)then
        !             cycle
        !         end if
        !         do l=jspin, norb, 2
        !             do p=1,norb
        !                 if(elecs%h2ei_old(j,k,l,p).ne.0)then
        !                     elecs%h2ei(cnt)=elecs%h2ei_old(j,k,l,p)
        !                     cr1=j
        !                     cr2=k
        !                     an2=p
        !                     an1=l
        !                     !annihilation
        !                     an2_cr2%dead(an1,cnt)=an2_cr2%alive(an1,cnt)
        !                     an2_cr2%neg_dead(an1,cnt)=an2_cr2%neg_alive(an1,cnt)
        !                     an2_cr2%alive(an1,cnt)=0
        !                     an2_cr2%neg_alive(:an1-1,cnt)=int(an2_cr2%neg_alive(:an1-1,cnt)*(-1),kind=1)
        !                     !annihilation
        !                     an2_cr2%dead(an2,cnt)=an2_cr2%alive(an2,cnt)
        !                     an2_cr2%neg_dead(an2,cnt)=an2_cr2%neg_alive(an2,cnt)
        !                     an2_cr2%alive(an2,cnt)=0
        !                     an2_cr2%neg_alive(:an2-1,cnt)=int(an2_cr2%neg_alive(:an2-1,cnt)*(-1),kind=1) 
        !                     !creation 
        !                     an2_cr2%alive(cr1,cnt)=an2_cr2%dead(cr1,cnt)
        !                     an2_cr2%neg_alive(cr1,cnt)=an2_cr2%neg_dead(cr1,cnt)
        !                     an2_cr2%dead(cr1,cnt)=0
        !                     an2_cr2%neg_alive(:cr1-1,cnt)=int(an2_cr2%neg_alive(:cr1-1,cnt)*(-1),kind=1)
        !                     !creation 
        !                     an2_cr2%alive(cr2,cnt)=an2_cr2%dead(cr2,cnt)
        !                     an2_cr2%neg_alive(cr2,cnt)=an2_cr2%neg_dead(cr2,cnt)
        !                     an2_cr2%dead(cr2,cnt)=0
        !                     an2_cr2%neg_alive(:cr2-1,cnt)=int(an2_cr2%neg_alive(:cr2-1,cnt)*(-1),kind=1) 
        !                     cnt=cnt+1
        !                 end if
        !             end do 
        !         end do 
        !     end do
        ! end do 
        !
    
        do l=1,e2 
            cr1=int(read_in(l+1,2)) 
            cr2=int(read_in(l+1,3)) 
            an2=int(read_in(l+1,4)) 
            an1=int(read_in(l+1,5))
            if((an1==2).and.(an2>2))then
                an1=an2
                an2=2
            end if
            ! print*,cr1,cr2,an2,an1,elecs%h2ei(l)
            !annihilation
            an2_cr2%dead(an1,l)=an2_cr2%alive(an1,l)
            an2_cr2%neg_dead(an1,l)=an2_cr2%neg_alive(an1,l)
            an2_cr2%alive(an1,l)=0
            an2_cr2%neg_alive(:an1-1,l)=int(an2_cr2%neg_alive(:an1-1,l)*(-1),kind=1)
            !annihilation
            an2_cr2%dead(an2,l)=an2_cr2%alive(an2,l)
            an2_cr2%neg_dead(an2,l)=an2_cr2%neg_alive(an2,l)
            an2_cr2%alive(an2,l)=0
            an2_cr2%neg_alive(:an2-1,l)=int(an2_cr2%neg_alive(:an2-1,l)*(-1),kind=1) 
            !creation 
            an2_cr2%alive(cr1,l)=an2_cr2%dead(cr1,l)
            an2_cr2%neg_alive(cr1,l)=an2_cr2%neg_dead(cr1,l)
            an2_cr2%dead(cr1,l)=0
            an2_cr2%neg_alive(:cr1-1,l)=int(an2_cr2%neg_alive(:cr1-1,l)*(-1),kind=1)
            !creation 
            an2_cr2%alive(cr2,l)=an2_cr2%dead(cr2,l)
            an2_cr2%neg_alive(cr2,l)=an2_cr2%neg_dead(cr2,l)
            an2_cr2%dead(cr2,l)=0
            an2_cr2%neg_alive(:cr2-1,l)=int(an2_cr2%neg_alive(:cr2-1,l)*(-1),kind=1) 
            
        end do
        ! read_in=0
        if(GDflg.eq.'y')then
            
            ! do l=1,e2
            !     cr1=int(read_in(l+1,2))
            !     cr2=int(read_in(l+1,3))
            !     an2=int(read_in(l+1,4))
            !     an1=int(read_in(l+1,5))
               
            !     do j=1,norb 
            !         if((an2.eq.j).or.(cr2.eq.j).or.(an1.eq.j).or.(cr1.eq.j))then
            !             an2_cr2%dcnt(0,j)=int(an2_cr2%dcnt(0,j)+1,kind=2)
            !             an2_cr2%dcnt(an2_cr2%dcnt(0,j),j)=int(l,kind=2)
            !             if((an2.eq.j).or.(an1.eq.j))then
            !                 if((cr1.eq.j).or.(cr2.eq.j))then
            !                     ! print*,j,' equals ', cr2,' or ',cr1, ' and one equals ', an2, an1
            !                     an2_cr2%alive_diff(j,:,l)=an2_cr2%alive(:,l)
            !                     an2_cr2%dead_diff(j,:,l)=an2_cr2%dead(:,l)
            !                     an2_cr2%neg_alive_diff(j,:,l)=an2_cr2%neg_alive(:,l)
            !                     an2_cr2%neg_dead_diff(j,:,l)=an2_cr2%neg_dead(:,l)
            !                     an2_cr2%alive_diff(j,j,l)=int(j+norb,kind=2)
            !                     an2_cr2%neg_alive_diff(j,j,l)=int(2*an2_cr2%neg_alive_diff(j,j,l),kind=1)
    
            !                 else 
            !                     ! print*,j,' equals ', an2,' or ',an1, ' but not ', cr2, cr1
            !                     an2_cr2%alive_diff(j,:,l)=an2_cr2%alive(:,l)
            !                     an2_cr2%dead_diff(j,:,l)=an2_cr2%dead(:,l)
            !                     an2_cr2%neg_alive_diff(j,:,l)=an2_cr2%neg_alive(:,l)
            !                     an2_cr2%neg_dead_diff(j,:,l)=an2_cr2%neg_dead(:,l)
            !                     an2_cr2%alive_diff(j,j,l)=int(j,kind=2)
            !                     an2_cr2%dead_diff(j,j,l)=int(j+norb,kind=2)
            !                     an2_cr2%neg_alive_diff(j,j,l)=int(an2_cr2%neg_dead_diff(j,j,l)*(-1),kind=1)
            !                 end if
            !             else
            !                 an2_cr2%dcnt(0,j)=int(an2_cr2%dcnt(0,j)+1,kind=2)
            !                 ! print*,j,' equals ', cr2,' or ',cr1, ' but not ', an2, an1
            !                 an2_cr2%alive_diff(j,:,l)=an2_cr2%alive(:,l)
            !                 an2_cr2%dead_diff(j,:,l)=an2_cr2%dead(:,l)
            !                 an2_cr2%neg_alive_diff(j,:,l)=an2_cr2%neg_alive(:,l)
            !                 an2_cr2%neg_dead_diff(j,:,l)=an2_cr2%neg_dead(:,l)
            !                 an2_cr2%alive_diff(j,j,l)=int(j,kind=2)
            !                 an2_cr2%dead_diff(j,j,l)=int(j+norb,kind=2)
            !                 an2_cr2%neg_dead_diff(j,j,l)=an2_cr2%neg_alive_diff(j,j,l)
            !                 an2_cr2%neg_alive_diff(j,j,l)=int(an2_cr2%neg_alive_diff(j,j,l)*(-1),kind=1)
            !             end if
            !         end if 
                    
            !     end do
            ! end do
          
        end if
    
        return
    
    end subroutine two_electrons
    
    
    subroutine one_electrons(elecs,an_cr,e1)
    
        implicit none 
        type(elecintrgl), intent(inout)::elecs
        type(oprts),intent(inout)::an_cr
        integer,intent(in)::e1
        real(kind=8),dimension(e1+1,3)::read_in
        integer::ierr,l,k,an,cr,j
    
        
        open(unit=129, file='integrals/h1e.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a)") 'Error in opening hnuc.csv file'
            errorflag = 1
            return
        end if
    
        do j=1,e1+1
            read(129,*) (read_in(j,k),k=1,3)
        end do 
        close(129)
    
        
        elecs%h1ei=read_in(2:,1)
        
        
        call alloc_oprts(an_cr,e1)
    
        do l=1,e1
            an=int(read_in(l+1,2))
            cr=int(read_in(l+1,3))
            !annihilation
            an_cr%dead(an,l)=an_cr%alive(an,l)
            an_cr%neg_dead(an,l)=an_cr%neg_alive(an,l)
            an_cr%alive(an,l)=0
            an_cr%neg_alive(:an-1,l)=int(an_cr%neg_alive(:an-1,l)*(-1),kind=1)
            !creation 
            an_cr%alive(cr,l)=an_cr%dead(cr,l)
            an_cr%neg_alive(cr,l)=an_cr%neg_dead(cr,l)
            an_cr%dead(cr,l)=0
            an_cr%neg_alive(:cr-1,l)=int(an_cr%neg_alive(:cr-1,l)*(-1),kind=1)
        end do 
        if(GDflg.eq.'y')then
            do l=1,e1
                an=int(read_in(l+1,2))
                cr=int(read_in(l+1,3))
                do j=1,norb 
                    if((an.eq.j).or.(cr.eq.j))then
                        an_cr%dcnt(0,j)=int(an_cr%dcnt(0,j)+1,kind=2)
                        an_cr%dcnt(an_cr%dcnt(0,j),j)=int(l,kind=2)
                        if(an.eq.cr)then 
                            an_cr%alive_diff(j,:,l)=an_cr%alive(:,l)
                            an_cr%dead_diff(j,:,l)=an_cr%dead(:,l)
                            an_cr%neg_alive_diff(j,:,l)=an_cr%neg_alive(:,l)
                            an_cr%neg_dead_diff(j,:,l)=an_cr%neg_dead(:,l)
                            an_cr%alive_diff(j,j,l)=int(j+norb,kind=2)
                            an_cr%neg_alive_diff(j,j,l)=int(2*an_cr%neg_alive_diff(j,j,l),kind=1)
                        else if(j.eq.an)then
                            an_cr%alive_diff(j,:,l)=an_cr%alive(:,l)
                            an_cr%dead_diff(j,:,l)=an_cr%dead(:,l)
                            an_cr%neg_alive_diff(j,:,l)=an_cr%neg_alive(:,l)
                            an_cr%neg_dead_diff(j,:,l)=an_cr%neg_dead(:,l)
                            an_cr%alive_diff(j,j,l)=int(j,kind=2)
                            an_cr%dead_diff(j,j,l)=int(j+norb,kind=2)
                            an_cr%neg_alive_diff(j,j,l)=int(an_cr%neg_dead_diff(j,j,l)*(-1),kind=1)
                        else if(j.eq.cr)then
                            an_cr%alive_diff(j,:,l)=an_cr%alive(:,l)
                            an_cr%dead_diff(j,:,l)=an_cr%dead(:,l)
                            an_cr%neg_alive_diff(j,:,l)=an_cr%neg_alive(:,l)
                            an_cr%neg_dead_diff(j,:,l)=an_cr%neg_dead(:,l)
                            an_cr%alive_diff(j,j,l)=int(j,kind=2)
                            an_cr%dead_diff(j,j,l)=int(j+norb,kind=2)
                            an_cr%neg_dead_diff(j,j,l)=an_cr%neg_alive_diff(j,j,l)
                            an_cr%neg_alive_diff(j,j,l)=int(an_cr%neg_alive_diff(j,j,l)*(-1),kind=1)
                        end if
                    end if 
                end do
            end do 
        
        end if
        
        
        return
    
    
    end subroutine one_electrons

End Module stuff

program test 

    use stuff 
    use vals 


    implicit none 

    type(zombiest), dimension(3):: zstore
    type(elecintrgl)::elecs
    type(oprts)::an_cr,an2_cr2,an2_cr2_diff
    integer::j,k,ierr,l,state_size
    real(kind=8)::r
    integer,allocatable,dimension(:,:,:)::occupancy_an
    integer,allocatable,dimension(:,:,:,:)::occupancy_2an
    complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
    complex(kind=8),allocatable, dimension(:,:,:)::z2l
    real(kind=8),dimension(10,10)::haml_diff 
    complex(kind=8),dimension(10)::zs1sin,zs1cos,zs2sin,zs2cos
    integer, allocatable, dimension(:) :: state

    call init()

    allocate(occupancy_an(norb,2,norb),stat=ierr)
    allocate(occupancy_2an(norb,norb,2,norb),stat=ierr)
    allocate(z1jk(norb,norb,2,norb),stat=ierr)
    allocate(z2l(norb,2,norb),stat=ierr)

    

    occupancy_an=1
    occupancy_2an=1

    call random_seed( size=state_size )
    allocate(state(state_size))
    state = 20180815
    call random_seed(put=state)
       
    do j=1,norb
        occupancy_an(j,1,j)=0
        do l=j-1, 1, -1
            occupancy_an(j,1,l)=-1
        end do
    end do
    
    do j=1,norb
        do k=1, norb
            occupancy_2an(j,k,:,:)=occupancy_an(j,:,:)
            occupancy_2an(j,k,2,k)=occupancy_2an(j,k,1,k)
            occupancy_2an(j,k,1,k)=0
            do l=k-1,1,-1
                occupancy_2an(j,k,1,l)=occupancy_2an(j,k,1,l)*(-1)
            end do
        end do
    end do
      


    allocate(elecs%h1ei_old(10,10))
    allocate(elecs%h2ei_old(10,10,10,10))
        
    call electronintegrals_old(elecs)
    call electronintegrals(elecs,an_cr,an2_cr2,an2_cr2_diff)
    

    ! zstore(1)%phi(1:nel)=0.5*pirl
    ! zstore(1)%phi(nel+1:)=0
    ! zstore(1)%sin=0
    ! zstore(1)%cos=1
    ! zstore(1)%sin(1:nel)=1
    ! zstore(1)%cos(1:nel)=0
    ! zstore(1)%val(nel:norb)=zstore(1)%sin
    ! zstore(1)%val(norb+1:)=zstore(1)%cos
    
    do j=1,3
        do k=1,norb
            call random_number(r)
            zstore(j)%phi(k)=2*pirl*r 
        end do 
        zstore(j)%val(0)=0
        zstore(j)%sin=sin(zstore(j)%phi)
        zstore(j)%cos=cos(zstore(j)%phi)
        zstore(j)%val(1:norb)=zstore(j)%sin
        zstore(j)%val(norb+1:)=zstore(j)%cos
    end do
   
    do l=1, norb
        z2l(l,1,:)=zstore(2)%sin(:)
        z2l(l,2,:)=zstore(2)%cos(:)
        z2l(l,2,l)=z2l(l,1,l)
        z2l(l,1,l)=cmplx(0.0,0.0)
    end do

    do j=1, norb
        do k=1, norb
            z1jk(j,k,:,:)=z2l(j,:,:)
            z1jk(j,k,2,k)=z1jk(j,k,1,k)
            z1jk(j,k,1,k)=cmplx(0.0,0.0)
        end do
    end do

    do l=1, norb
        z2l(l,1,:)=zstore(1)%sin(:)
        z2l(l,2,:)=zstore(1)%cos(:)
        z2l(l,2,l)=z2l(l,1,l)
        z2l(l,1,l)=cmplx(0.0,0.0)
    end do
    
    z1jk=z1jk*occupancy_2an
    z2l=z2l*occupancy_an
    haml_diff=0
    call haml_grad(haml_diff,zstore,elecs,an_cr,an2_cr2,an2_cr2_diff,2)

    zs1sin=cmplx(zstore(2)%sin,0.0d0,kind=8)
    zs1cos=cmplx(zstore(2)%cos,0.0d0,kind=8)
    zs2sin=cmplx(zstore(1)%sin,0.0d0,kind=8)
    zs2cos=cmplx(zstore(1)%cos,0.0d0,kind=8)
    
    ! call two_elec_part_grad(zs1sin,zs1cos,zs2sin,zs2cos,z2l,z1jk,elecs%h2ei_old,occupancy_2an,occupancy_an,elecs) 
    


end program test