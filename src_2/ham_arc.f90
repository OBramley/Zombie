MODULE ham 

    use globvars
    use alarrays
    use operators
    use omp_lib
    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators
    subroutine he_row(ham,zstore,elecs,row,size,occupancy_2an,occupancy_an_cr,occupancy_an,passback)
        
        implicit none
        type(hamiltonian), intent(inout)::ham
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::row,size
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an,occupancy_an_cr
        integer,dimension(:,:,:),intent(in)::occupancy_an
        complex(kind=8),dimension(:,:,:),intent(inout)::passback
        complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
        complex(kind=8),allocatable, dimension(:,:,:)::z2l
        integer::j,k,l,m,ierr,equal
        complex(kind=8)::h1etot, h2etot
        real(kind=8),allocatable,dimension(:)::h1etot_diff_bra,h2etot_diff_bra
        real(kind=8),allocatable,dimension(:)::h1etot_diff_ket,h2etot_diff_ket
        real(kind=8),allocatable,dimension(:)::overlap_diff
        
        if (errorflag .ne. 0) return
        ierr = 0
        
        
        allocate(z1jk(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        if(GDflg.eq.'y')then
            allocate(h1etot_diff_bra(norb),stat=ierr)
            if(ierr==0) allocate(h2etot_diff_bra(norb),stat=ierr)
            if(ierr==0) allocate(h1etot_diff_ket(norb),stat=ierr)
            if(ierr==0) allocate(h2etot_diff_ket(norb),stat=ierr)
            if(ierr==0) allocate(overlap_diff(norb),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in gradient array vector allocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
            h1etot_diff_bra=0.0
            h1etot_diff_ket=0.0
            h2etot_diff_bra=0.0
            h2etot_diff_ket=0.0
        end if
    
        
        h1etot=cmplx(0.0,0.0)
        h2etot=cmplx(0.0,0.0)
        
        
        !$omp parallel shared(passback,zstore,z1jk, z2l) private(j,k,l)
        if(row.eq.1)then
            
            !$omp do simd
            do l=1, norb
                z2l(l,1,:)=zstore(1)%sin(:)
                z2l(l,2,:)=zstore(1)%cos(:)
                z2l(l,2,l)=z2l(l,1,l)
                z2l(l,1,l)=cmplx(0.0,0.0)
            end do
            !!$omp end do simd
        
        else
            z2l=passback
            
        end if
        !$omp do simd collapse(2)
        do j=1, norb
            do k=1, norb
                z1jk(j,k,:,:)=z2l(j,:,:)
                z1jk(j,k,2,k)=z1jk(j,k,1,k)
                z1jk(j,k,1,k)=cmplx(0.0,0.0)
            end do
        end do
    
        !!$omp end do simd
        !$omp end parallel
        z1jk=z1jk*occupancy_2an
        
        !!$omp flush(z1jk)
        !!$omp parallel private(j,k,l,z2l,h1etot,h2etot, &
        !!$omp h1etot_diff_bra,h2etot_diff_bra, h1etot_diff_ket,h2etot_diff_ket,&
        !!$omp  h1etot_diff,h2etot_diff,m) shared(z1jk,zstore,row,size,occupancy_2an,occupancy_an_cr,occupancy_an)
        !!$omp do schedule(dynamic)
        do m=row,size
        
            h1etot=cmplx(0.0,0.0)
            h2etot=cmplx(0.0,0.0)
            
            if(m.gt.row)then
                !$omp parallel shared(z2l)
                !$omp do simd
                do l=1, norb
                    z2l(l,1,:)=zstore(m)%sin(:)
                    z2l(l,2,:)=zstore(m)%cos(:)
                    z2l(l,2,l)=z2l(l,1,l)
                    z2l(l,1,l)=cmplx(0.0,0.0)
                end do
                !!$omp end do simd
                !$omp end parallel
                if(m.eq.row+1)then 
                    passback=z2l
                end if
            end if
            
            h1etot=(0.0,0.0)
            h2etot=(0.0,0.0)
            
            equal=9
            if(GDflg.eq.'y')then
                h1etot_diff_bra=0.0
                h1etot_diff_ket=0.0
                h2etot_diff_bra=0.0
                h2etot_diff_ket=0.0
                if(2.eq.row)then 
                    if(row.eq.m)then 
                        equal=1
                    else if(row.ne.m)then 
                        equal=2
                    end if
                else if(2.ne.row)then 
                    if(m.eq.2)then 
                        equal=3
                    end if
                else 
                    equal = 9
                end if
            end if
            
            
            call one_elec_part(zstore(row)%sin,zstore(row)%cos,z2l,h1etot,occupancy_an_cr,&
                elecs%h1ei,h1etot_diff_bra,h1etot_diff_ket,zstore(m)%sin,zstore(m)%cos,equal)
            
            
            z2l=z2l*occupancy_an
            !!$omp flush(z2l)
            
            call two_elec_part(zstore(row)%sin,zstore(row)%cos,z1jk,z2l,h2etot,occupancy_2an,occupancy_an,&
                    elecs%h2ei,h2etot_diff_bra,h2etot_diff_ket,zstore(m)%sin,zstore(m)%cos,equal)
            
            ham%ovrlp(row,m)=overlap(zstore(row),zstore(m))
            ham%ovrlp(m,row)= ham%ovrlp(row,m)
            ham%hjk(row,m)=h1etot+h2etot+(elecs%hnuc*ham%ovrlp(row,m))
            ham%hjk(m,row)=ham%hjk(row,m)
            
            if(GDflg.eq.'y') then
                if(row.eq.2)then
                    ham%diff_hjk(row,m,:)=h1etot_diff_bra+h2etot_diff_bra
                    if(m.eq.row)then
                        ham%diff_ovrlp(row,m,:) = 0
                    else
                        overlap_diff = diff_overlap_2(zstore(row)%sin,zstore(row)%cos,zstore(m)%sin,zstore(m)%cos,equal)
                        ham%diff_ovrlp(row,m,:) = overlap_diff
                    end if
                else if(m.eq.2)then
                    ham%diff_hjk(m,row,:)=h1etot_diff_ket+h2etot_diff_ket
                    overlap_diff = diff_overlap_2(zstore(row)%sin,zstore(row)%cos,zstore(m)%sin,zstore(m)%cos,equal)
                    ham%diff_ovrlp(m,row,:) = overlap_diff
                end if
            end if
            ! write(6,"(a,i0,a,i0,a)") "Hamiltonian row/comun ",row,'/', m, " completed"
        end do
        !!$omp end do
        !!$omp end parallel
        
        
        deallocate(z1jk,stat=ierr)
        if(ierr==0) deallocate(z2l,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        if(GDflg.eq.'y')then
            deallocate(h1etot_diff_bra,stat=ierr)
            if(ierr==0) deallocate(h2etot_diff_bra,stat=ierr)
            if(ierr==0) deallocate(h1etot_diff_ket,stat=ierr)
            if(ierr==0) deallocate(h2etot_diff_ket,stat=ierr)
            if(ierr==0) deallocate(overlap_diff,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in gradient array vector deallocation . ierr had value ", ierr
                errorflag=1
                return
            end if 
        end if

        
        return

    end subroutine he_row

    subroutine one_elec_part(zs1sin,zs1cos,z2l,h1etot,occupancy,h1ei,h1etot_diff_bra,h1etot_diff_ket,zs2sin,zs2cos,equal)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),intent(inout)::h1etot
        integer,dimension(:,:,:,:),intent(in)::occupancy
        real(kind=8),dimension(:),intent(inout)::h1etot_diff_bra,h1etot_diff_ket
        real(kind=8), dimension(:,:), intent(in)::h1ei
        integer,intent(in)::equal
        real(kind=8),allocatable,dimension(:,:,:)::h1etot_diff
        complex(kind=8),allocatable,dimension(:,:)::temp
        integer::j,k,l,len,ierr
        complex(kind=8),allocatable,dimension(:,:)::zomt
        real(kind=8),allocatable,dimension(:)::chng_prod,temp_prod,prod

        if (errorflag .ne. 0) return

        len=norb
        
        allocate(temp(len,len),stat=ierr)
        allocate(zomt(2,len),stat=ierr)
        if(equal.lt.4)then
            allocate(h1etot_diff(len,len,len),stat=ierr)
            allocate(prod(len))
            allocate(chng_prod(len))
            allocate(temp_prod(len))
            h1etot_diff_bra=0.0
            h1etot_diff_ket=0.0
            h1etot_diff=0.0
        end if
        h1etot=(0.0,0.0)
    
        temp=cmplx(0.0,0.0)
        
        if(equal.lt.4)then
            !$omp parallel do  collapse(2) &
            !$omp & private(j,k,l,prod,chng_prod,temp_prod,zomt) &
            !$omp & shared(h1ei,occupancy,z2l,zs1sin,zs1cos,temp,zs2sin,zs2cos,equal,h1etot_diff)
            do j=1, len
                do k=1, len
                    if(h1ei(j,k).ne.(0.0)) then 
                        zomt=z2l(j,:,:)
                        zomt(1,k)=zomt(2,k)
                        zomt(2,k)=cmplx(0.0,0.0)
                        zomt=zomt*occupancy(j,k,:,:)
                        temp(j,k)=product((conjg(zs1sin)*zomt(1,:))+(conjg(zs1cos)*zomt(2,:)))*h1ei(j,k)
                        prod=real(((conjg(zs1sin)*zomt(1,:)))+((conjg(zs1cos)*zomt(2,:))))
                        if(equal.eq.1)then
                            if(j.eq.k)then
                                if((real(zs2cos(j)).eq.0).and.(real(zs2sin(j)).eq.1).or.(real(zs2sin(j)).eq.0))then
                                    h1etot_diff(j,k,j)=0.0
                                else
                                    chng_prod=prod           !dead amplitude is zero
                                    chng_prod(j)=real(2*zs1sin(j)*zs1cos(j)*occupancy(j,k,1,j))
                                    h1etot_diff(j,k,j)=product(chng_prod)*h1ei(j,k)   
                                end if 
                            else if(j.ne.k)then
                                chng_prod=prod
                                if((real(zs2cos(j)).eq.0).and.(real(zs2sin(j)).eq.1)) then
                                    chng_prod(j)=-1*occupancy(j,k,2,j)
                                else if((real(zs2sin(j)).eq.0).and.(real(zs2cos(j)).eq.1)) then
                                    chng_prod(j)=occupancy(j,k,2,j)
                                else
                                    chng_prod(j)=real(((zs1cos(j)**2)-(zs1sin(j)**2))*occupancy(j,k,2,j))
                                end if
                                h1etot_diff(j,k,j)=product(chng_prod)*h1ei(j,k)
                                chng_prod=prod
                                if((real(zs2cos(k)).eq.0).and.(real(zs2sin(k)).eq.1)) then
                                    chng_prod(k)=-1*occupancy(j,k,1,k)
                                else if((real(zs2sin(k)).eq.0).and.(real(zs2cos(k)).eq.1)) then
                                    chng_prod(k)=occupancy(j,k,1,k)
                                else
                                    chng_prod(k)=real(((zs1cos(k)**2)-(zs1sin(k)**2))*occupancy(j,k,1,k))
                                end if               !dead amplitude is zero
                                h1etot_diff(j,k,k)=product(chng_prod)*h1ei(j,k)
                            end if
                        else if(equal.eq.2)then 
                            chng_prod=real(zs1cos*zs2sin*occupancy(j,k,1,:)-zs1sin*zs2cos*occupancy(j,k,2,:))            
                            if(j.eq.k)then  !dead amplitude is zero
                                chng_prod(j)=real(zs1cos(j)*zs2sin(j)*occupancy(j,k,1,j))
                            else                   
                                chng_prod(j)=-real(zs1sin(j)*zs2sin(j))*occupancy(j,k,2,j)
                                chng_prod(k)=real(zs1cos(k)*zs2cos(k))*occupancy(j,k,1,k) 
                            end if
                            !!$omp parallel do simd shared(prod,chng_prod) private(temp,j)
                            do l=1,len
                                temp_prod=prod
                                temp_prod(l)=chng_prod(l)
                                h1etot_diff(j,k,l)=product(temp_prod)*h1ei(j,k)
                            end do
                            !!$omp end parallel do simd 
                        else if(equal.eq.3)then
                            chng_prod=real(zs1sin*zs2cos*occupancy(j,k,1,:)-zs1cos*zs2sin*occupancy(j,k,2,:))                  
                            if(j.eq.k)then  !dead amplitude is zero
                                chng_prod(j)=real(zs1sin(j)*zs2cos(j)*occupancy(j,k,1,j))
                            else                   
                                chng_prod(j)=real(zs1cos(j)*zs2cos(j))*occupancy(j,k,2,j)!alive amplitude is zero
                                chng_prod(k)=-real(zs1sin(k)*zs2sin(k))*occupancy(j,k,1,k) !dead amplitude is zero
                            end if
                            !!$omp parallel do simd  shared(prod,ket_prod) private(temp_prod,j)
                            do l=1,len
                                temp_prod=prod
                                temp_prod(l)=chng_prod(l)
                                h1etot_diff(j,k,l)=product(temp_prod)*h1ei(j,k)
                            end do
                            !!$omp end parallel do simd   
                        end if
                    end if
                end do
            end do
            !$omp end parallel do          
        else
            !$omp parallel do  collapse(2) &
            !$omp & private(j,k,zomt) shared(h1ei,occupancy,z2l,zs1sin,zs1cos,temp)
            do j=1, len
                do k=1, len
                    if(h1ei(j,k).ne.(0.0)) then 
                        zomt=z2l(j,:,:)
                        zomt(1,k)=zomt(2,k)
                        zomt(2,k)=cmplx(0.0,0.0)
                        zomt=zomt*occupancy(j,k,:,:)
                        temp(j,k)=product((conjg(zs1sin)*zomt(1,:))+(conjg(zs1cos)*zomt(2,:)))*h1ei(j,k)
                    end if
                end do
            end do
            !$omp end parallel do 
        end 

        do j=1, norb
            do k=1, norb
                h1etot=h1etot+temp(j,k)
                if(equal.lt.4)then 
                    if(equal.eq.1)then
                        h1etot_diff_bra= h1etot_diff_bra+h1etot_diff(j,k,:)
                    else if(equal.eq.2)then
                        h1etot_diff_bra= h1etot_diff_bra+h1etot_diff(j,k,:)
                    else if(equal.eq.3)then
                        h1etot_diff_ket= h1etot_diff_ket+h1etot_diff(j,k,:)
                    end if
                end if
            end do
        end do

        if(equal.lt.4)then
            if(equal.eq.1)then
                h1etot_diff_ket=h1etot_diff_bra
            end if 
            deallocate(prod)
            deallocate(chng_prod)
            deallocate(temp_prod)
            deallocate(h1etot_diff,stat=ierr)
        end if
        
        
        deallocate(temp,stat=ierr)
        return

    end subroutine one_elec_part


    subroutine two_elec_part(zs1sin,zs1cos,z1jk,z2l,h2etot,occupancy_2an,occupancy_an,h2ei,&
        h2etot_diff_bra,h2etot_diff_ket,zs2sin,zs2cos,equal)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:,:),intent(in)::z1jk
        complex(kind=8),intent(inout)::h2etot
        integer,dimension(:,:,:,:),intent(in)::occupancy_2an
        integer,dimension(:,:,:),intent(in)::occupancy_an
        real(kind=8),dimension(norb),intent(inout)::h2etot_diff_bra,h2etot_diff_ket
        real(kind=8), dimension(:,:,:,:), intent(in)::h2ei
        integer,intent(in)::equal
        real(kind=8),allocatable,dimension(:,:,:,:)::h2etot_diff
        integer::j,k,l,ierr,len
        complex(kind=8),allocatable,dimension(:,:,:)::tot

        if (errorflag .ne. 0) return
        len=norb
        
        allocate(tot(len,len,len),stat=ierr)
        tot=cmplx(0.0,0.0)
        h2etot=(0.0,0.0)
        if(equal.lt.4)then
            allocate(h2etot_diff(len,len,len,len),stat=ierr)
            h2etot_diff=0.0
        end if
        
        
        
    
        !$omp parallel shared(z1jk,z2l,zs1sin,zs1cos,zs2sin,zs2cos,tot,occupancy_2an,occupancy_an,h2ei,equal,h2etot_diff) private(j)
        !$omp do
        do j=1, norb
            tot(j,:,:) = two_elec_part_body(norb,zs1sin,zs2sin,z2l,z1jk(j,:,:,:),h2ei(j,:,:,:),j)
            if(equal.lt.4)then
                h2etot_diff(j,:,:,:) = two_elec_part_grad(norb,zs1sin,zs1cos,zs2sin,zs2cos,z2l,z1jk(j,:,:,:),&
                h2ei(j,:,:,:),occupancy_2an(j,:,:,:),occupancy_an(:,:,:),j,equal)
            end if
        end do
        !$omp end  do
        !$omp end parallel
        
        do j=1, norb
            do k=1,norb
                do l=1, norb
                    h2etot=h2etot+tot(j,k,l)
                    if(equal.lt.4)then 
                        if(equal.eq.1)then
                            h2etot_diff_bra= h2etot_diff_bra+h2etot_diff(j,k,l,:)
                        else if(equal.eq.2)then
                            h2etot_diff_bra= h2etot_diff_bra+h2etot_diff(j,k,l,:)
                        else if(equal.eq.3)then
                            h2etot_diff_ket= h2etot_diff_ket+h2etot_diff(j,k,l,:)
                        end if
                    end if
                end do 
            end do
        end do


        h2etot=h2etot*0.5
        if(equal.lt.4)then
            if(equal.eq.1)then
                h2etot_diff_ket=h2etot_diff_bra
            end if 
            h2etot_diff_bra = h2etot_diff_bra*0.5
            h2etot_diff_ket = h2etot_diff_ket*0.5
            deallocate(h2etot_diff,stat=ierr)
        end if
        deallocate(tot,stat=ierr)
        
        return

    end subroutine two_elec_part

    function two_elec_part_body(len,zs1sin,zs2sin,z2l,z1jk,h2ei,j)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs2sin
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:),intent(in)::z1jk
        real(kind=8), dimension(:,:,:), intent(in)::h2ei
        complex(kind=8),dimension(:,:),allocatable::vmult,tot
        complex(kind=8),dimension(len,len)::two_elec_part_body
        complex(kind=8),allocatable,dimension(:)::gg,hh
        integer,intent(in)::len
        integer::j,k,l,p,jspin,gmax,hmin,ierr

        ierr=0

        allocate(tot(len,len),stat=ierr)
        ! allocate(two_elec_part_body(len,len),stat=ierr)
        allocate(gg(len),stat=ierr)
        allocate(hh(len),stat=ierr)
        allocate(vmult(2,len),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in vmult allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        
        tot=cmplx(0.0,0.0)

        if(zs1sin(j)==(0.0,0.0))then
            two_elec_part_body=(0.0,0.0)
            return
        end if
        if(modulo(j,2)==0)then
            jspin=2
        else
            jspin=1
        end if
        
        do k=1, len
            if(j.eq.k) cycle
            if(occ_iszero(z1jk(k,:,:)).eqv..true.)then
                CYCLE 
            end if
            
            !$omp parallel do  &
            !$omp & private(gmax,hmin,gg,hh,p,vmult) shared(z1jk,z2l,tot)
            do l=jspin, len, 2
                if(zs2sin(l)==(0.0,0.0))then

                    CYCLE
                else
                    vmult=conjg(z1jk(k,:,:))*(z2l(l,:,:))
                    
                    gg(1:len)=(0.0,0.0)
                    hh(1:len)=(0.0,0.0)
                    gmax=len
                    gg(1)=vmult(2,1)-vmult(1,1)

                    do p=2, len
                        gg(p)=gg(p-1)*(vmult(2,p)-vmult(1,p))
                        if(gg(p)==(0.0,0.0))then
                            gmax=p
                            EXIT 
                        end if
                    end do
                    
                    hmin=0
                    hh(len) = vmult(2,len)+vmult(1,len)
                    do p=(len-1),1,-(1)
                        hh(p)=hh(p+1)*(vmult(2,p)+vmult(1,p))
                        if(hh(p)==(0.0,0.0))then
                            hmin=p
                            EXIT 
                        end if
                    end do

                    tot(k,l)=(0.0,0.0)
                    if (gmax < hmin) then
                        tot(k,l)=(0.0,0.0)
                        cycle
                    end if

                    if(h2ei(k,l,1).ne.0) then
                        tot(k,l) = tot(k,l)+(conjg(z1jk(k,2,1))*z2l(l,1,1)*hh(2)*h2ei(k,l,1))
                    end if

                    do p=2,len-1
                        if(h2ei(k,l,p).ne.0.0) then
                            tot(k,l) = tot(k,l)+ (gg(p-1)*conjg(z1jk(k,2,p))*z2l(l,1,p)*hh(p+1)*h2ei(k,l,p))
                        end if
                    end do

                    if(h2ei(k,l,len).ne.0) then
                        tot(k,l) = tot(k,l) +(gg(len-1)*conjg(z1jk(k,2,len))*z2l(l,1,len)*h2ei(k,l,len))
                    end if
                end if
            end do
            !$omp end parallel do
        end do


        

        two_elec_part_body=tot

        deallocate(tot,stat=ierr)
        deallocate(vmult)
        deallocate(gg,stat=ierr)
        deallocate(hh,stat=ierr)

        return

    end function two_elec_part_body

    function two_elec_part_grad(len,zs1sin,zs1cos,zs2sin,zs2cos,z2l,z1jk,h2ei,occupancy_2an,occupancy_an,j,equal)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:),intent(in)::z1jk
        real(kind=8), dimension(:,:,:), intent(in)::h2ei
        integer,dimension(:,:,:),intent(in)::occupancy_2an
        integer,dimension(:,:,:),intent(in)::occupancy_an
        integer,intent(in)::j,equal,len
        real(kind=8),dimension(len,len,len)::two_elec_part_grad
        integer,allocatable,dimension(:,:)::occupancy
        integer::k,l,jspin

        allocate(occupancy(2,len))
        
        two_elec_part_grad=0.0

        if(modulo(j,2)==0)then
            jspin=2
        else
            jspin=1
        end if
        do k=1, len
            if(j.eq.k)then
                cycle
            end if
            do l=jspin, len, 2
                occupancy=occupancy_2an(k,:,:)*occupancy_an(l,:,:)
                two_elec_part_grad(k,l,:) = z_an_z3_diff_2(len,z1jk(k,:,:),z2l(l,:,:),h2ei(k,l,:),&
                equal,occupancy,j,k,l,zs1sin,zs1cos,zs2sin,zs2cos)
            end do
        end do

        deallocate(occupancy)
        
        return

    end function two_elec_part_grad

    ! Top level routine to allocate hamiltonian and overlap matrices 
    subroutine hamgen(ham,zstore,elecs,size,verb)

        implicit none

        type(hamiltonian), intent(inout)::ham 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,allocatable,dimension(:,:,:)::occupancy_an
        integer,allocatable,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an
        complex(kind=8),allocatable, dimension(:,:,:)::passback
        integer, allocatable,dimension(:)::IPIV1
        integer,intent(in)::size,verb
        complex(kind=8),allocatable,dimension(:)::WORK1
        real(kind=8),allocatable,dimension(:,:,:)::temp2
        integer:: j,k,l,ierr
        

        if (errorflag .ne. 0) return
        write(0,"(a)") "Begining of ham"
        allocate(occupancy_an(norb,2,norb),stat=ierr)
        if(ierr==0) allocate(occupancy_2an(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(occupancy_an_cr(norb,norb,2,norb),stat=ierr)
        if(ierr==0) allocate(passback(norb,2,norb),stat=ierr)
        if(ierr==0) allocate(temp2(ndet,ndet,norb),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in occupancy vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 

        occupancy_an=1
        occupancy_2an=1
        occupancy_an_cr=1

        
        !call omp_set_nested(.TRUE.)
        !$omp parallel private(j,k,l) shared(occupancy_an,occupancy_2an,occupancy_an_cr)
        !$omp do
        do j=1,norb
            occupancy_an(j,1,j)=0
            do l=j-1, 1, -1
                occupancy_an(j,1,l)=-1
            end do
        end do
        !!$omp end do
        !$omp barrier
        !$omp do collapse(2)
        do j=1,norb
            do k=1, norb
                occupancy_2an(j,k,:,:)=occupancy_an(j,:,:)
                occupancy_2an(j,k,2,k)=occupancy_2an(j,k,1,k)
                occupancy_2an(j,k,1,k)=0
                occupancy_an_cr(j,k,:,:)=occupancy_an(j,:,:)
                occupancy_an_cr(j,k,1,k)=occupancy_an_cr(j,k,2,k)
                occupancy_an_cr(j,k,2,k)=0
                do l=k-1,1,-1
                    occupancy_2an(j,k,1,l)=occupancy_2an(j,k,1,l)*(-1)
                    occupancy_an_cr(j,k,1,l)=occupancy_an_cr(j,k,1,l)*(-1)
                end do
            end do
        end do
        !!$omp end do simd
        !$omp end parallel
        write(0,"(a)") "occupancy done"
        
        do j=1, size
            call he_row(ham,zstore,elecs,j,size,occupancy_2an,occupancy_an_cr,occupancy_an,passback)
            if(verb.eq.1)then
                write(6,"(a,i0,a)") "Hamiltonian row ",j, " completed"
            end if  
        end do

        
        
        if(ierr==0) deallocate(passback,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in passback vector deallocation . ierr had value ", ierr
            errorflag=1
        end if

        ham%inv=ham%ovrlp
        allocate(IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
            errorflag=1
        end if 
        
        if (ierr==0) allocate(WORK1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
            errorflag=1
        end if   

        if (ierr==0) call ZGETRF(size,size,ham%inv,size,IPIV1,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if
        if (ierr==0) call ZGETRI(size,ham%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if

        if (ierr==0) deallocate(IPIV1,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV vector deallocation . ierr had value ", ierr
            errorflag=1
        end if

        if (ierr==0) deallocate(WORK1,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in WORK vector deallocation . ierr had value ", ierr
            errorflag=1
        end if
        
        !$omp parallel
        !$omp workshare
        ham%kinvh=matmul(ham%inv,ham%hjk)
        !$omp end workshare
        !$omp end parallel
        
        if(GDflg.eq.'y')then
            do k=1, ndet
                do l=1, ndet
                    if(l.eq.2)then
                        temp2(k,2,:)=matmul(REAL(ham%inv(k,:)),ham%diff_ovrlp(2,:,:))
                    else
                        temp2(k,l,:)=real(ham%inv(k,l))*ham%diff_ovrlp(2,l,:)
                    end if
                end do
            end do
            do k=1, ndet
                do l=1, ndet
                    ham%diff_invh(2,k,l,:)=matmul(transpose(temp2(k,:,:)),real(ham%kinvh(:,l)))*(-1)
                end do
            end do
        end if

        deallocate(temp2)
        
        return
        
    end subroutine hamgen


    !

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
            do j=1, len
                
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

                if(vec(1).ne.0) then
                    if(j.eq.1)then
                        if((annihilate1.eq.1).or.(annihilate1_2.eq.1))then
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot1 = tot1+((REAL(zs1cos(1)*zs2sin(1))*occupancy(2,1))*hh_1(2)*vec(1))
                            end if
                        else
                            if(annihilate2.ne.1)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot1 = tot1+((REAL(-zs1sin(1)*zs2sin(1))*occupancy(2,1))*hh_1(2)*vec(1))
                            end if 
                        end if
                    else
                        tot1 = tot1+(REAL(conjg(z1(2,1))*z2(1,1))*hh_1(2)*vec(1))
                    end if
                end if
                do k=2,len-1
                    if(vec(k).ne.0.0) then
                        if(k.eq.j)then
                            if((annihilate1.eq.k).or.(annihilate1_2.eq.k))then
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                    tot1 = tot1+(gg_1(k-1)*(REAL(zs1cos(k)*zs2sin(k))*occupancy(2,k))*hh_1(k+1)*vec(k))
                                end if
                            else
                                if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                    tot1 = tot1+(gg_1(k-1)*(REAL(-zs1sin(k)*zs2sin(k))*occupancy(2,k))*hh_1(k+1)*vec(k))
                                end if 
                            end if
                        else
                            tot1 = tot1+ (gg_1(k-1)*REAL(conjg(z1(2,k))*z2(1,k))*hh_1(k+1)*vec(k))
                        end if
                    end if
                end do

                if(vec(len).ne.0) then
                    if(len.eq.j)then
                        if((annihilate1.eq.len).or.(annihilate1_2.eq.len))then
                            if(annihilate2.ne.len)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
                                tot1 = tot1+(gg_1(len-1)*(REAL(zs1cos(len)*zs2sin(len))*occupancy(2,len))*vec(len))
                            end if
                        else
                            if(annihilate2.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
                                tot1 = tot1+(gg_1(len-1)*(REAL(-zs1sin(len)*zs2sin(len))*occupancy(2,len))*vec(len))
                            end if 
                        end if
                    else
                        tot1 = tot1 +(gg_1(len-1)*REAL(conjg(z1(2,len))*z2(1,len))*vec(len))
                    end if
                end if

                temp(j)=tot1
            end do
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
            do j=1, len
                
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

    function diff_overlap_2(zs1sin,zs1cos,zs2sin,zs2cos,dtype)
    
        implicit none
      
        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        integer,intent(in)::dtype
        real(kind=8),dimension(norb)::diff_overlap_2
        real(kind=8),allocatable,dimension(:)::bra_prod,ket_prod,prod
        integer::j
        real(kind=8),allocatable,dimension(:)::temp,temp2

        allocate(prod(norb))

        ! if (errorflag .ne. 0) return
        prod=real((conjg(zs1sin)*zs2sin)+(conjg(zs1cos)*zs2cos))
        if(dtype.eq.2)then
            allocate(temp(norb))
            allocate(bra_prod(norb))
            bra_prod=real(conjg(zs1cos)*zs2sin)-real(conjg(zs1sin)*zs2cos)
         
            !$omp parallel do  shared(prod,bra_prod) private(temp)
            do j=1,norb
                temp=prod
                temp(j)=bra_prod(j)
                diff_overlap_2(j)=product(temp)
            end do
            !$omp end parallel do   
         
            deallocate(bra_prod)
            deallocate(temp)
        else if(dtype.eq.3)then
            allocate(ket_prod(norb))
            allocate(temp2(norb))
            ket_prod=real(conjg(zs1sin)*zs2cos)-real(conjg(zs1cos)*zs2sin)
          
            !$omp parallel do  shared(prod,ket_prod) private(temp2)
            do j=1,norb
                temp2=prod
                temp2(j)=ket_prod(j)
                diff_overlap_2(j)=product(temp2)
            end do
            !$omp end parallel do   
           
            deallocate(ket_prod)
            deallocate(temp2)
        end if
        deallocate(prod)
        return
        
    end function diff_overlap_2

END MODULE ham



    ! ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators
    ! subroutine he_row(ham,zstore,elecs,row,size,occupancy_2an,occupancy_an_cr,occupancy_an,passback)
        
    !     implicit none
    !     type(hamiltonian), intent(inout)::ham
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     type(elecintrgl),intent(in)::elecs
    !     integer,intent(in)::row,size
    !     integer,dimension(:,:,:,:),intent(in)::occupancy_2an,occupancy_an_cr
    !     integer,dimension(:,:,:),intent(in)::occupancy_an
    !     complex(kind=8),dimension(:,:,:),intent(inout)::passback
    !     complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
    !     complex(kind=8),allocatable, dimension(:,:,:)::z2l
    !     integer::j,k,l,m,ierr,equal
    !     complex(kind=8)::h1etot, h2etot
    !     real(kind=8),dimension(norb)::h1etot_diff_bra,h2etot_diff_bra
    !     real(kind=8),dimension(norb)::h1etot_diff_ket,h2etot_diff_ket
    !     real(kind=8),dimension(norb)::overlap_diff
        
    !     if (errorflag .ne. 0) return
    !     ierr = 0

    !     allocate(z1jk(norb,norb,2,norb),stat=ierr)
    !     if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
    !         errorflag=1
    !         return
    !     end if 

    !     h1etot=cmplx(0.0,0.0)
    !     h2etot=cmplx(0.0,0.0)
    !     h1etot_diff_bra=0.0
    !     h1etot_diff_ket=0.0
    !     h2etot_diff_bra=0.0
    !     h2etot_diff_ket=0.0
    
   
    !     !$omp parallel shared(passback,zstore,z1jk, z2l) private(j,k,l)
    !     if(row.eq.1)then
            
    !         !$omp do simd
    !         do l=1, norb
    !             z2l(l,1,:)=zstore(1)%sin(:)
    !             z2l(l,2,:)=zstore(1)%cos(:)
    !             z2l(l,2,l)=z2l(l,1,l)
    !             z2l(l,1,l)=cmplx(0.0,0.0)
    !         end do
    !         !$omp end do simd
           
    !     else
    !         z2l=passback
          
    !     end if
    !     !$omp do simd collapse(2)
    !     do j=1, norb
    !         do k=1, norb
    !             z1jk(j,k,:,:)=z2l(j,:,:)
    !             z1jk(j,k,2,k)=z1jk(j,k,1,k)
    !             z1jk(j,k,1,k)=cmplx(0.0,0.0)
    !         end do
    !     end do
    !     !$omp end do simd
    !     !$omp end parallel
    !     z1jk=z1jk*occupancy_2an

    !     !!$omp flush(z1jk)
    !     !!$omp parallel private(j,k,l,z2l,h1etot,h2etot, &
    !     !!$omp h1etot_diff_bra,h2etot_diff_bra, h1etot_diff_ket,h2etot_diff_ket,&
    !     !!$omp  h1etot_diff,h2etot_diff,m) shared(z1jk,zstore,row,size,occupancy_2an,occupancy_an_cr,occupancy_an)
    !     !!$omp do schedule(dynamic)
    !     do m=row,size
      
    !         h1etot=cmplx(0.0,0.0)
    !         h2etot=cmplx(0.0,0.0)
          
            
    !         if(m.gt.row)then
               
    !             !$omp parallel shared(z2l)
    !             !$omp do simd
    !             do l=1, norb
    !                 z2l(l,1,:)=zstore(m)%sin(:)
    !                 z2l(l,2,:)=zstore(m)%cos(:)
    !                 z2l(l,2,l)=z2l(l,1,l)
    !                 z2l(l,1,l)=cmplx(0.0,0.0)
    !             end do
    !             !$omp end do simd
    !             !$omp end parallel
    !             if(m.eq.row+1)then 
    !                 passback=z2l
    !             end if
    !         end if
          
    !         h1etot=(0.0,0.0)
    !         h2etot=(0.0,0.0)
            
    !         equal=9
    !         if(GDflg.eq.'y')then
    !             h1etot_diff_bra=0.0
    !             h1etot_diff_ket=0.0
    !             h2etot_diff_bra=0.0
    !             h2etot_diff_ket=0.0
    !             if(2.eq.row)then 
    !                 if(row.eq.m)then 
    !                     equal=1
    !                 else if(row.ne.m)then 
    !                     equal=2
    !                 end if
    !             else if(2.ne.row)then 
    !                 if(m.eq.2)then 
    !                     equal=3
    !                 end if
    !             else 
    !                 equal = 9
    !             end if
    !         end if
        
            
    !         call one_elec_part(zstore(row),z2l,h1etot,occupancy_an_cr,&
    !             elecs%h1ei,h1etot_diff_bra,h1etot_diff_ket,zstore(m),equal)
                       
    !         z2l=z2l*occupancy_an
    !         !!$omp flush(z2l)
         
    !         call two_elec_part(zstore(row),z1jk,z2l,h2etot,occupancy_2an,occupancy_an,&
    !                 elecs%h2ei,h2etot_diff_bra,h2etot_diff_ket,zstore(m),equal)
            
    !         ham%ovrlp(row,m)=overlap(zstore(row),zstore(m))
    !         ham%ovrlp(m,row)= ham%ovrlp(row,m)
    !         ham%hjk(row,m)=h1etot+h2etot+(elecs%hnuc*ham%ovrlp(row,m))
    !         ham%hjk(m,row)=ham%hjk(row,m)
            
    !         if(GDflg.eq.'y') then
    !             if(row.eq.2)then
    !                 ham%diff_hjk(row,m,:)=h1etot_diff_bra+h2etot_diff_bra
    !                 if(m.eq.row)then
    !                     ham%diff_ovrlp(row,m,:) = 0
    !                 else
    !                     overlap_diff = diff_overlap(zstore(row),zstore(m),equal)
    !                     ham%diff_ovrlp(row,m,:) = overlap_diff
    !                 end if
    !             else if(m.eq.2)then
    !                 ham%diff_hjk(m,row,:)=h1etot_diff_ket+h2etot_diff_ket
    !                 overlap_diff = diff_overlap(zstore(row),zstore(m),equal)
    !                 ham%diff_ovrlp(m,row,:) = overlap_diff
    !             end if
    !         end if
    !         ! write(6,"(a,i0,a,i0,a)") "Hamiltonian row/comun ",row,'/', m, " completed"
    !     end do
    !     !!$omp end do
    !     !!$omp end parallel
     

    !     deallocate(z1jk,stat=ierr)
    !     if(ierr==0) deallocate(z2l,stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
    !         errorflag=1
    !         return
    !     end if 

        
    !     return

    ! end subroutine he_row

    ! subroutine one_elec_part(zs1,z2l,h1etot,occupancy,h1ei,h1etot_diff_bra,h1etot_diff_ket,zs2,equal)

    !     implicit none

    !     type(zombiest),intent(in)::zs1,zs2
    !     complex(kind=8),dimension(:,:,:),intent(in)::z2l
    !     complex(kind=8),intent(inout)::h1etot
    !     integer,dimension(:,:,:,:),intent(in)::occupancy
    !     real(kind=8),dimension(norb),intent(inout)::h1etot_diff_bra,h1etot_diff_ket
    !     real(kind=8), dimension(:,:), intent(in)::h1ei
    !     integer,intent(in)::equal
    !     real(kind=8),dimension(norb,norb,norb)::h1etot_diff
    !     complex(kind=8),dimension(norb,norb)::temp
    !     integer::j,k

    !     if (errorflag .ne. 0) return
    !     h1etot=(0.0,0.0)
    !     h1etot_diff_bra=0.0
    !     h1etot_diff_ket=0.0
    !     h1etot_diff=0.0
    !     temp=cmplx(0.0,0.0)
    !     !$omp parallel private(j,k) shared(h1ei,occupancy,z2l,zs1,zs2,equal,temp,h1etot_diff)
    !     !$omp do simd collapse(2)
    !     do j=1, norb
    !         do k=1, norb
    !             if(h1ei(j,k).ne.(0.0)) then 
    !                 temp(j,k)=one_elec_body(zs1,z2l(j,:,:),occupancy(j,k,:,:),h1ei(j,k),k)
    !                 if(equal.lt.4)then 
    !                     h1etot_diff(j,k,:)=one_elec_body_grad(zs1,zs2,z2l(j,:,:),occupancy(j,k,:,:),h1ei(j,k),k,j,equal)
    !                 end if
    !             end if
    !         end do
    !     end do
    !     !$omp end do simd
    !     !$omp end parallel
       
    !     do j=1, norb
    !         do k=1, norb
    !             h1etot=h1etot+temp(j,k)
    !             if(equal.lt.4)then 
    !                 if(equal.eq.1)then
    !                     h1etot_diff_bra= h1etot_diff_bra+h1etot_diff(j,k,:)
    !                 else if(equal.eq.2)then
    !                     h1etot_diff_bra= h1etot_diff_bra+h1etot_diff(j,k,:)
    !                 else if(equal.eq.3)then
    !                     h1etot_diff_ket= h1etot_diff_ket+h1etot_diff(j,k,:)
    !                 end if
    !             end if
    !         end do
    !     end do

    !     if(equal.eq.1)then
    !         h1etot_diff_ket=h1etot_diff_bra
    !     end if 
     
    !     return

    ! end subroutine one_elec_part

    ! complex(kind=8) function one_elec_body(zs1,z2l,occupancy,h1ei,k)

    !     implicit none
    !     type(zombiest),intent(in)::zs1
    !     integer,dimension(:,:),intent(in)::occupancy
    !     complex(kind=8),dimension(:,:),intent(in)::z2l
    !     real(kind=8), intent(in)::h1ei
    !     integer,intent(in)::k
    !     complex(kind=8),dimension(2,norb)::zomt

    !     zomt=z2l
    !     zomt(1,k)=zomt(2,k)
    !     zomt(2,k)=cmplx(0.0,0.0)
    !     zomt=zomt*occupancy
    !     one_elec_body=product((conjg(zs1%sin)*zomt(1,:))+(conjg(zs1%cos)*zomt(2,:)))*h1ei
     
    !     return

    ! end function one_elec_body

    ! function one_elec_body_grad(zs1,zs2,z2l,occupancy,h1ei,k,j,equal)

    !     implicit none

    !     type(zombiest),intent(in)::zs1,zs2
    !     real(kind=8),dimension(norb)::one_elec_body_grad
    !     integer,dimension(:,:),intent(in)::occupancy
    !     complex(kind=8),dimension(:,:),intent(in)::z2l
    !     real(kind=8), intent(in)::h1ei
    !     integer,intent(in)::k,j,equal
    !     complex(kind=8),dimension(2,norb)::zomt

    !     zomt=z2l
    !     zomt(1,k)=zomt(2,k)
    !     zomt(2,k)=cmplx(0.0,0.0)
    !     zomt=zomt*occupancy
    !     one_elec_body_grad = diff_overlap_cran(zs1,zs2,equal,zomt,j,k,occupancy)*h1ei
        
    !     return

    ! end function one_elec_body_grad
        
    ! subroutine two_elec_part(zs1,z1jk,z2l,h2etot,occupancy_2an,occupancy_an,h2ei,h2etot_diff_bra,h2etot_diff_ket,zs2,equal)

    !     implicit none

    !     type(zombiest),intent(in)::zs1,zs2
    !     complex(kind=8),dimension(:,:,:),intent(in)::z2l
    !     complex(kind=8),dimension(:,:,:,:),intent(in)::z1jk
    !     complex(kind=8),intent(inout)::h2etot
    !     integer,dimension(:,:,:,:),intent(in)::occupancy_2an
    !     integer,dimension(:,:,:),intent(in)::occupancy_an
    !     real(kind=8),dimension(norb),intent(inout)::h2etot_diff_bra,h2etot_diff_ket
    !     real(kind=8), dimension(:,:,:,:), intent(in)::h2ei
    !     integer,intent(in)::equal
    !     real(kind=8),dimension(norb,norb,norb,norb)::h2etot_diff
    !     integer::j,k,l
    !     complex(kind=8),dimension(norb,norb,norb)::tot

    !     if (errorflag .ne. 0) return
        
    !     h2etot=(0.0,0.0)
    !     h2etot_diff_bra=0.0
    !     h2etot_diff_ket=0.0
    !     h2etot_diff=0.0
    !     tot=cmplx(0.0,0.0)
       
    !     !$omp parallel shared(z1jk,z2l,zs1,zs2,tot,occupancy_2an,occupancy_an,h2ei,equal,h2etot_diff) private(j)
    !     !$omp do
    !     do j=1, norb
    !         tot(j,:,:) = two_elec_part_body(zs1,zs2,z2l,z1jk(j,:,:,:),h2ei(j,:,:,:),j)
    !         if(equal.lt.4)then
    !             h2etot_diff(j,:,:,:) = two_elec_part_grad(zs1,zs2,z2l,z1jk(j,:,:,:),&
    !             h2ei(j,:,:,:),occupancy_2an(j,:,:,:),occupancy_an(:,:,:),j,equal)
    !         end if
    !     end do
    !     !$omp end  do
    !     !$omp end parallel

    !     do j=1, norb
    !         do k=1,norb
    !             do l=1, norb
    !                 h2etot=h2etot+tot(j,k,l)
    !                 if(equal.lt.4)then 
    !                     if(equal.eq.1)then
    !                         h2etot_diff_bra= h2etot_diff_bra+h2etot_diff(j,k,l,:)
    !                     else if(equal.eq.2)then
    !                         h2etot_diff_bra= h2etot_diff_bra+h2etot_diff(j,k,l,:)
    !                     else if(equal.eq.3)then
    !                         h2etot_diff_ket= h2etot_diff_ket+h2etot_diff(j,k,l,:)
    !                     end if
    !                 end if
    !             end do 
    !         end do
    !     end do

    !     if(equal.eq.1)then
    !         h2etot_diff_ket=h2etot_diff_bra
    !     end if 

    !     h2etot=h2etot*0.5
    !     h2etot_diff_bra = h2etot_diff_bra*0.5
    !     h2etot_diff_ket = h2etot_diff_ket*0.5
     
    !     return

    ! end subroutine two_elec_part

    ! function two_elec_part_body(zs1,zs2,z2l,z1jk,h2ei,j)

    !     implicit none

    !     type(zombiest),intent(in)::zs1,zs2
    !     complex(kind=8),dimension(:,:,:),intent(in)::z2l
    !     complex(kind=8),dimension(:,:,:),intent(in)::z1jk
    !     real(kind=8), dimension(:,:,:), intent(in)::h2ei
    !     integer,intent(in)::j
    !     complex(kind=8),dimension(norb,norb)::tot,two_elec_part_body
    !     integer::k,l,jspin

    !     tot=cmplx(0.0,0.0)

    !     if(zs1%sin(j)==(0.0,0.0))then
    !         two_elec_part_body=(0.0,0.0)
    !         return
    !     end if
    !     if(modulo(j,2)==0)then
    !         jspin=2
    !     else
    !         jspin=1
    !     end if
    !     do k=1, norb
    !         if(j.eq.k) cycle
    !         if(occ_iszero(z1jk(k,:,:)).eqv..true.)then
    !             CYCLE
    !         end if
    !         do l=jspin, norb, 2
    !             if(zs2%sin(l)==(0.0,0.0))then
    !                 CYCLE
    !             end if
    !             tot(k,l) = z_an_z3(z1jk(k,:,:),z2l(l,:,:),h2ei(k,l,:))
    !         end do
    !     end do

    !     two_elec_part_body=tot

    !     return

    ! end function two_elec_part_body

    ! function two_elec_part_grad(zs1,zs2,z2l,z1jk,h2ei,occupancy_2an,occupancy_an,j,equal)

    !     implicit none

    !     type(zombiest),intent(in)::zs1,zs2
    !     complex(kind=8),dimension(:,:,:),intent(in)::z2l
    !     complex(kind=8),dimension(:,:,:),intent(in)::z1jk
    !     real(kind=8), dimension(:,:,:), intent(in)::h2ei
    !     integer,dimension(:,:,:),intent(in)::occupancy_2an
    !     integer,dimension(:,:,:),intent(in)::occupancy_an
    !     integer,intent(in)::j,equal
    !     real(kind=8),dimension(norb,norb,norb)::two_elec_part_grad
    !     integer,dimension(2,norb)::occupancy
    !     integer::k,l,jspin

    !     two_elec_part_grad=0.0

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
    !             occupancy=occupancy_2an(k,:,:)*occupancy_an(l,:,:)
    !             two_elec_part_grad(k,l,:) = z_an_z3_diff(z1jk(k,:,:),z2l(l,:,:),h2ei(k,l,:),equal,occupancy,j,k,l,zs1,zs2)
    !         end do
    !     end do

    !     return

    ! end function two_elec_part_grad

    ! ! Top level routine to allocate hamiltonian and overlap matrices 
    ! subroutine hamgen(ham,zstore,elecs,size,verb)

    !     implicit none

    !     type(hamiltonian), intent(inout)::ham 
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     type(elecintrgl),intent(in)::elecs
    !     integer,allocatable,dimension(:,:,:)::occupancy_an
    !     integer,allocatable,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an
    !     complex(kind=8),allocatable, dimension(:,:,:)::passback
    !     integer, allocatable,dimension(:)::IPIV1
    !     integer,intent(in)::size,verb
    !     complex(kind=8),allocatable,dimension(:)::WORK1
    !     real(kind=8),dimension(ndet,ndet,norb)::temp2
    !     integer:: j,k,l,ierr

    !     if (errorflag .ne. 0) return
        
    !     allocate(occupancy_an(norb,2,norb),stat=ierr)
    !     if(ierr==0) allocate(occupancy_2an(norb,norb,2,norb),stat=ierr)
    !     if(ierr==0) allocate(occupancy_an_cr(norb,norb,2,norb),stat=ierr)
    !     if(ierr==0) allocate(passback(norb,2,norb),stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in occupancy vector allocation . ierr had value ", ierr
    !         errorflag=1
    !         return
    !     end if 

    !     occupancy_an=1
    !     occupancy_2an=1
    !     occupancy_an_cr=1
       
       
    !     !call omp_set_nested(.TRUE.)
    !     !$omp parallel private(j,k,l) shared(occupancy_an,occupancy_2an,occupancy_an_cr)
    !     !$omp do
    !     do j=1,norb
    !         occupancy_an(j,1,j)=0
    !         do l=j-1, 1, -1
    !             occupancy_an(j,1,l)=-1
    !         end do
    !     end do
    !     !$omp end do
    !     !$omp barrier
    !     !$omp do simd collapse(2)
    !     do j=1,norb
    !         do k=1, norb
    !             occupancy_2an(j,k,:,:)=occupancy_an(j,:,:)
    !             occupancy_2an(j,k,2,k)=occupancy_2an(j,k,1,k)
    !             occupancy_2an(j,k,1,k)=0
    !             occupancy_an_cr(j,k,:,:)=occupancy_an(j,:,:)
    !             occupancy_an_cr(j,k,1,k)=occupancy_an_cr(j,k,2,k)
    !             occupancy_an_cr(j,k,2,k)=0
    !             do l=k-1,1,-1
    !                 occupancy_2an(j,k,1,l)=occupancy_2an(j,k,1,l)*(-1)
    !                 occupancy_an_cr(j,k,1,l)=occupancy_an_cr(j,k,1,l)*(-1)
    !             end do
    !         end do
    !     end do
    !     !$omp end do simd
    !     !$omp end parallel
    
       
    !     do j=1, size
    !         call he_row(ham,zstore,elecs,j,size,occupancy_2an,occupancy_an_cr,occupancy_an,passback)
    !         if(verb.eq.1)then
    !             write(6,"(a,i0,a)") "Hamiltonian row ",j, " completed"
    !         end if  
    !     end do

        
     
    !    if(ierr==0) deallocate(passback,stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in passback vector deallocation . ierr had value ", ierr
    !         errorflag=1
    !     end if

    !     ham%inv=ham%ovrlp
    !     allocate(IPIV1(size),stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
    !         errorflag=1
    !     end if 
        
    !     if (ierr==0) allocate(WORK1(size),stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
    !         errorflag=1
    !     end if   

    !     if (ierr==0) call ZGETRF(size,size,ham%inv,size,IPIV1,ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)")"Error in ZGETRF",ierr
    !     end if
    !     if (ierr==0) call ZGETRI(size,ham%inv,size,IPIV1,WORK1,size,ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)")"Error in ZGETRF",ierr
    !     end if

    !     if (ierr==0) deallocate(IPIV1,stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in IPIV vector deallocation . ierr had value ", ierr
    !         errorflag=1
    !     end if

    !     if (ierr==0) deallocate(WORK1,stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in WORK vector deallocation . ierr had value ", ierr
    !         errorflag=1
    !     end if
        
    !     !$omp parallel
    !     !$omp workshare
    !     ham%kinvh=matmul(ham%inv,ham%hjk)
    !     !$omp end workshare
    !     !$omp end parallel
       
    !     if(GDflg.eq.'y')then
    !         do k=1, ndet
    !             do l=1, ndet
    !                 if(l.eq.2)then
    !                     temp2(k,2,:)=matmul(REAL(ham%inv(k,:)),ham%diff_ovrlp(2,:,:))
    !                 else
    !                     temp2(k,l,:)=real(ham%inv(k,l))*ham%diff_ovrlp(2,l,:)
    !                 end if
    !             end do
    !         end do
    !         do k=1, ndet
    !             do l=1, ndet
    !                 ham%diff_invh(2,k,l,:)=matmul(transpose(temp2(k,:,:)),real(ham%kinvh(:,l)))*(-1)
    !             end do
    !         end do
    !     end if
      
    !     return
        
    ! end subroutine hamgen

    
    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators
    ! subroutine he_row_gpu(ph1ei,ph2ei,phnuc,psin,pcos,pphi,occupancy_2an,&
    !     occupancy_an_cr,occupancy_an,passback,phjk,povrlp,pdiff_hjk,pdiff_ovrlp,row,size)
        
    !     implicit none
        
    !     real(kind=8),dimension(:,:),intent(in)::ph1ei
    !     real(kind=8), dimension(:,:,:,:),intent(in)::ph2ei
    !     real(kind=8),intent(in) :: phnuc
    !     complex(kind=8), dimension(:,:),intent(in)::psin
    !     complex(kind=8), dimension(:,:),intent(in)::pcos
    !     real(kind=8),dimension(:,:),intent(in)::pphi
    !     complex(kind=8), dimension(:,:), intent(inout)::phjk
    !     complex(kind=8), dimension(:,:), intent(inout)::povrlp
    !     real(kind=8), dimension(:,:,:), intent(inout)::pdiff_hjk
    !     real(kind=8), dimension(:,:,:), intent(inout)::pdiff_ovrlp
    !     integer,intent(in)::row,size
    !     integer,dimension(:,:,:,:),intent(in)::occupancy_2an,occupancy_an_cr
    !     integer,dimension(:,:,:),intent(in)::occupancy_an
    !     complex(kind=8),dimension(:,:,:),intent(inout)::passback
    !     complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
    !     complex(kind=8),allocatable, dimension(:,:,:)::z2l
    !     integer::j,k,l,m,ierr,equal
    !     complex(kind=8)::h1etot, h2etot
    !     real(kind=8),dimension(norb)::h1etot_diff_bra,h2etot_diff_bra
    !     real(kind=8),dimension(norb)::h1etot_diff_ket,h2etot_diff_ket
    !     real(kind=8),dimension(norb)::bra_prod,ket_prod,prod,temp
        
    !     !!$omp declare target 
    !     if (errorflag .ne. 0) return
    !     ierr = 0

        
    !     allocate(z1jk(norb,norb,2,norb),stat=ierr)
    !     if(ierr==0) allocate(z2l(norb,2,norb),stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in annihilation and creation array vector allocation . ierr had value ", ierr
    !         errorflag=1
    !         return
    !     end if
  

    !     !$omp target teams map(alloc:z1jk(norb,norb,2,norb),z2l(norb,2,norb),prod(norb),bra_prod(norb),temp(norb),ket_prod(norb),&
    !     !$omp h1etot_diff_bra(norb),h2etot_diff_bra(norb),h1etot_diff_ket(norb),h2etot_diff_ket(norb)) &
    !     !$omp & map(to:equal,row,m) map(from:h1etot,h2etot)&
    !     !$omp & num_teams(max_teams) thread_limit(threadpteam)
    !     h1etot=cmplx(0.0,0.0)
    !     h2etot=cmplx(0.0,0.0)
    !     h1etot_diff_bra(:)=0.0
    !     h1etot_diff_ket(:)=0.0
    !     h2etot_diff_bra(:)=0.0
    !     h2etot_diff_ket(:)=0.0
    
    
    !     if(row.eq.1)then
    !         !$omp distribute parallel do simd
    !         do l=1, norb
    !             z2l(l,1,:)=psin(1,:) 
    !             z2l(l,2,:)=pcos(1,:)
    !             z2l(l,2,l)=z2l(l,1,l)
    !             z2l(l,1,l)=cmplx(0.0,0.0)
    !         end do
    !         !$omp end distribute parallel do simd
    !     else
    !         z2l=passback
    !     end if
        
    !     !$omp distribute parallel do simd collapse(2)
    !     do j=1, norb
    !         do k=1, norb
    !             z1jk(j,k,:,:)=z2l(j,:,:)
    !             z1jk(j,k,2,k)=z1jk(j,k,1,k)
    !             z1jk(j,k,1,k)=cmplx(0.0,0.0)
    !         end do
    !     end do
    !     !$omp end distribute parallel do simd
    !     z1jk=z1jk*occupancy_2an
  

    !     do m=row,size
        
    !         h1etot=cmplx(0.0,0.0)
    !         h2etot=cmplx(0.0,0.0)
    !         h1etot_diff_bra(:)=0.0
    !         h1etot_diff_ket(:)=0.0
    !         h2etot_diff_bra(:)=0.0
    !         h2etot_diff_ket(:)=0.0
    
    !         if(m.gt.row)then
    !             !$omp distribute parallel do simd
    !             do l=1, norb
    !                 z2l(l,1,:)=psin(m,:) 
    !                 z2l(l,2,:)=pcos(m,:)
    !                 z2l(l,2,l)=z2l(l,1,l)
    !                 z2l(l,1,l)=cmplx(0.0,0.0)
    !             end do
    !             !$omp end distribute parallel do simd
    !             if(m.eq.row+1)then 
    !                 passback=z2l
    !             end if
    !         end if
            
            
    !         equal=9
    !         if(2.eq.row)then 
    !             if(row.eq.m)then 
    !                 equal=1
    !             else if(row.ne.m)then 
    !                 equal=2
    !             end if
    !         else if(2.ne.row)then 
    !             if(m.eq.2)then 
    !                 equal=3
    !             end if
    !         else 
    !             equal = 9
    !         end if
          
         
    !         call one_elec_part_gpu(psin(row,:),psin(m,:),pcos(row,:),pcos(m,:),pphi(row,:),&
    !             z2l,h1etot,occupancy_an_cr,ph1ei,h1etot_diff_bra,h1etot_diff_ket,equal) 

    !         z2l=z2l*occupancy_an
          
    !         call two_elec_part_gpu(psin(row,:),psin(m,:),pcos(row,:),pcos(m,:),z1jk,z2l,h2etot,occupancy_2an,&
    !         occupancy_an,ph2ei,h2etot_diff_bra,h2etot_diff_ket,equal)
         
    !         povrlp(row,m)=product(((conjg(psin(row,:))*psin(m,:)))+((conjg(pcos(row,:))*pcos(m,:))))
    !         povrlp(m,row)= povrlp(row,m)
    !         phjk(row,m)=h1etot+h2etot+(phnuc*povrlp(row,m))
    !         phjk(m,row)=phjk(row,m)
            
    !         if(GDflg.eq.'y') then
    !             if(row.eq.2)then
    !                 temp=h1etot_diff_bra(:)+h2etot_diff_bra(:)
    !                 pdiff_hjk(row,m,:)= temp 
    !                 if(m.eq.row)then
    !                     pdiff_ovrlp(row,m,:) = 0
    !                 else
    !                     prod=real((conjg(psin(row,:))*psin(m,:))+(conjg(pcos(row,:))*pcos(m,:)))
    !                     bra_prod=real(conjg(pcos(row,:))*psin(m,:))-real(conjg(psin(row,:))*pcos(m,:))
    !                     !$omp distribute parallel do
    !                     do j=1,norb
    !                         temp=prod
    !                         temp(j)=bra_prod(j)
    !                         pdiff_ovrlp(row,m,j)=product(temp)
    !                     end do
                 
    !                 end if
    !             else if(m.eq.2)then
    !                 temp=h1etot_diff_ket(:)+h2etot_diff_ket(:)
    !                 pdiff_hjk(m,row,:)=temp!h1etot_diff_ket(:)+h2etot_diff_ket(:)
    !                 prod=real((conjg(psin(row,:))*psin(m,:))+(conjg(pcos(row,:))*pcos(m,:)))
    !                 ket_prod=real(conjg(psin(row,:))*pcos(m,:))-real(conjg(pcos(row,:))*psin(m,:))
    !                 !$omp distribute parallel do
    !                 do j=1,norb
    !                     temp=prod
    !                     temp(j)=ket_prod(j)
    !                     pdiff_ovrlp(m,row,j)=product(temp)
    !                 end do
                    
    !             end if
             
    !         end if
          
    !     end do
    !     !$omp end target teams
       
    !     deallocate(z1jk,stat=ierr)
    !     if(ierr==0) deallocate(z2l,stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in annihilation and creation array vector deallocation . ierr had value ", ierr
    !         errorflag=1
    !         return
    !     end if 

       
    !     return

    ! end subroutine he_row_gpu

    ! subroutine one_elec_part_gpu(psin1,psin2,pcos1,pcos2,pphi1,z2l,h1etot,occupancy,&
    !     ph1ei,h1etot_diff_bra,h1etot_diff_ket,equal)
    !     !!$omp declare target 
    !     implicit none

    !     complex(kind=8),dimension(:,:,:),intent(in)::z2l
    !     complex(kind=8),intent(inout)::h1etot
    !     integer,dimension(:,:,:,:),intent(in)::occupancy
    !     real(kind=8),dimension(:),intent(inout)::h1etot_diff_bra,h1etot_diff_ket
    !     real(kind=8),dimension(:,:),intent(in)::ph1ei
    !     integer,intent(in)::equal
    !     complex(kind=8), dimension(:),intent(in)::psin1,psin2
    !     complex(kind=8), dimension(:),intent(in)::pcos1,pcos2
    !     real(kind=8),dimension(:),intent(in)::pphi1
    !     complex(kind=8),dimension(2,norb)::zomt
    !     real(kind=8),dimension(norb)::temp1,temp2,bra_prod,ket_prod,prod,total
    !     integer::j,k,l
    !     complex(kind=8)::totc
        
       
  
    !     if (errorflag .ne. 0) return
    !     h1etot=(0.0,0.0)
    !     h1etot_diff_bra=0.0
    !     h1etot_diff_ket=0.0
    !     totc=cmplx(0.0,0.0)
        
       
 
    !     !$omp target teams distribute parallel do simd collapse(2) reduction(+:totc,total)&
    !     !$omp map(alloc:zomt(2,norb),bra_prod(norb),ket_prod(norb),prod(norb),temp1(norb),temp2(norb))&
    !     !$omp & map(tofrom:totc,total)&
    !     !$omp & num_teams(max_teams) thread_limit(threadpteam) &
    !     !$omp & private(temp1,temp2,bra_prod,ket_prod,prod) &
    !     !$omp & shared(ph1ei,psin1,psin2,pcos1,pcos2,pphi1,occupancy,equal,z2l)
    !     do j=1, norb
    !         do k=1, norb
    !             if(ph1ei(j,k).eq.0.0)then
    !                 cycle
    !             end if
    !             zomt(:,:)=z2l(j,:,:)
    !             zomt(1,k)=zomt(2,k)
    !             zomt(2,k)=cmplx(0.0,0.0)
    !             zomt=zomt*occupancy(j,k,:,:)
    !             totc=totc+product((conjg(psin1)*zomt(1,:))+(conjg(pcos1)*zomt(2,:)))*ph1ei(j,k)
    !             if(GDflg.eq.'y')then
    !                 if(equal.lt.4)then
    !                     prod=real(((conjg(psin1)*zomt(1,:)))+((conjg(pcos1)*zomt(2,:))))
    !                     ! Differntial of overlap of the same ZS 0 unless an operator has acted on ZS
    !                     if(equal.eq.1)then
    !                         if(j.eq.k)then
    !                             if((real(pcos2(j)).eq.0).and.(real(psin2(j)).eq.1).or.(real(psin2(j)).eq.0))then
    !                                 total(j)=total(j)+0
    !                             else
    !                                 bra_prod=prod           !dead amplitude is zero
    !                                 bra_prod(j)=sin(2*pphi1(j))*occupancy(j,k,1,j)
    !                                 total(j)= total(j)+product(bra_prod)*ph1ei(j,k)
    !                             end if 
    !                         else if(j.ne.k)then
    !                             bra_prod=prod
    !                             if((real(pcos2(j)).eq.0).and.(real(psin2(j)).eq.1)) then
    !                                 bra_prod(j)=-1*occupancy(j,k,2,j)
    !                             else if((real(psin2(j)).eq.0).and.(real(pcos2(j)).eq.1)) then
    !                                 bra_prod(j)=occupancy(j,k,2,j)
    !                             else
    !                                 bra_prod(j)=cos(2*pphi1(j))*occupancy(j,k,2,j)
    !                             end if
    !                             total(j)= total(j) + product(bra_prod)*ph1ei(j,k)      
    !                             bra_prod=prod
    !                             if((real(pcos2(k)).eq.0).and.(real(psin2(k)).eq.1)) then
    !                                 bra_prod(k)=-1*occupancy(j,k,1,k)
    !                             else if((real(psin2(k)).eq.0).and.(real(pcos2(k)).eq.1)) then
    !                                 bra_prod(k)=occupancy(j,k,1,k)
    !                             else
    !                                 bra_prod(k)=cos(2*pphi1(k))*occupancy(j,k,1,k)
    !                             end if               !dead amplitude is zero
    !                             total(k)= total(k) + product(bra_prod)*ph1ei(j,k)       
    !                         end if
    !                     else if(equal.eq.2)then 
    !                         bra_prod=real(pcos1*psin2*occupancy(j,k,1,:)-psin1*pcos2*occupancy(j,k,2,:))            
    !                         if(j.eq.k)then  !dead amplitude is zero
    !                             bra_prod(j)=real(pcos1(j)*psin2(j)*occupancy(j,k,1,j))
    !                         else                   
    !                             bra_prod(j)=-real(psin1(j)*psin2(j))*occupancy(j,k,2,j)
    !                             bra_prod(k)=real(pcos1(k)*pcos2(k))*occupancy(j,k,1,k) 
    !                         end if
    !                         do l=1,norb
    !                             temp1=prod
    !                             temp1(l)=bra_prod(l)
    !                             total(l)= total(l)+product(temp1)*ph1ei(j,k)   
    !                         end do 
    !                     else if(equal.eq.3)then
    !                         ket_prod=real(psin1*pcos2*occupancy(j,k,1,:)-pcos1*psin2*occupancy(j,k,2,:))                  
    !                         if(j.eq.k)then  !dead amplitude is zero
    !                             ket_prod(j)=real(psin1(j)*pcos2(j)*occupancy(j,k,1,j))
    !                         else                   
    !                             ket_prod(j)=real(pcos1(j)*pcos2(j))*occupancy(j,k,2,j)!alive amplitude is zero
    !                             ket_prod(k)=-real(psin1(k)*psin2(k))*occupancy(j,k,1,k) !dead amplitude is zero
    !                         end if
    !                         do l=1,norb
    !                             temp2=prod
    !                             temp2(l)=ket_prod(l)
    !                             total(l)= total(l) + product(temp2)*ph1ei(j,k)   
    !                         end do   
    !                     end if  
    !                 end if
    !             end if
    !         end do
    !     end do
    !     !$omp end target teams  distribute parallel do simd
 
       
    !     h1etot=totc
    !     if(equal.eq.1)then
    !         h1etot_diff_bra=total
    !         h1etot_diff_ket=total
    !     else if(equal.eq.2)then
    !         h1etot_diff_bra=total
    !     else if(equal.eq.3)then
    !         h1etot_diff_ket=total
    !     end if
      
    !     return

    ! end subroutine one_elec_part_gpu
        
    ! subroutine two_elec_part_gpu(psin1,psin2,pcos1,pcos2,z1jk,z2l,h2etot,occupancy_2an,&
    !     occupancy_an,ph2ei,h2etot_diff_bra,h2etot_diff_ket,equal)
    !     !!$omp declare target 
    !     implicit none

    !     complex(kind=8),dimension(:,:,:),intent(in)::z2l
    !     complex(kind=8),dimension(:,:,:,:),intent(in)::z1jk
    !     complex(kind=8),intent(inout)::h2etot
    !     integer,dimension(:,:,:,:),intent(in)::occupancy_2an
    !     integer,dimension(:,:,:),intent(in)::occupancy_an
    !     real(kind=8),dimension(norb),intent(inout)::h2etot_diff_bra,h2etot_diff_ket
    !     real(kind=8), dimension(:,:,:,:), intent(in)::ph2ei
    !     integer,intent(in)::equal
    !     complex(kind=8), dimension(:),intent(in)::psin1,psin2
    !     complex(kind=8), dimension(:),intent(in)::pcos1,pcos2
    !     integer,dimension(2,norb)::occupancy
    !     integer::j,k,l,jspin
    !     real(kind=8)::tot
    !     complex(kind=8)::totc,totc2
    !     complex(kind=8),dimension(2,norb)::vmult
    !     real(kind=8),dimension(2,norb)::vmult_dd,vmultr
    !     complex(kind=8),dimension(norb)::gg,hh
    !     real(kind=8),dimension(norb)::gg_1,hh_1,gg_2,hh_2
    !     integer::n,p,gmax,hmin,gmax1,hmin1,gmax2,hmin2,breakflag

    !     real(kind=8),dimension(norb)::temp,total
    !     if (errorflag .ne. 0) return

        
    !     h2etot=(0.0,0.0)
    !     h2etot_diff_bra=0.0
    !     h2etot_diff_ket=0.0
    !     tot=0.0
    !     totc=(0.0,0.0)
    !     totc2=(0.0,0.0)
    !     !$omp target teams distribute parallel do simd reduction(+:totc,totc2) map(to:gmax,hmin,jspin,totc) &
    !     !$omp & map(alloc:vmult(2,norb),gg(norb),hh(norb),occupancy(2,norb)) map(tofrom:totc2) &
    !     !$omp & private(j,k,l,jspin,occupancy,gg,hh,gmax,hmin,vmult,temp) &
    !     !$omp & shared(z2l,z1jk,occupancy_2an,occupancy_an,ph2ei)&
    !     !$omp & num_teams(max_teams) thread_limit(threadpteam)
    !     do j=1, norb
    !         if(psin1(j)==(0.0,0.0))then
    !             CYCLE
    !         end if
    !         if(modulo(j,2)==0)then
    !             jspin=2
    !         else
    !             jspin=1
    !         end if
    !         do k=1, norb
    !             if(j.eq.k) cycle
    !             if(occ_iszero(z1jk(j,k,:,:)).eqv..true.)then
    !                 CYCLE
    !             end if
    !             do l=jspin, norb, 2
    !                 if(psin2(l)==(0.0,0.0))then
    !                     CYCLE
    !                 end if
                
    !                 vmult=conjg(z1jk(j,k,:,:))*(z2l(l,:,:))
    
    !                 gg(1:norb)=(0.0,0.0)
    !                 hh(1:norb)=(0.0,0.0)
    !                 gmax=norb
    !                 gg(1)=vmult(2,1)-vmult(1,1)

    !                 do n=2, norb
    !                     gg(n)=gg(n-1)*(vmult(2,n)-vmult(1,n))
    !                     if(gg(n)==(0.0,0.0))then
    !                         gmax=n
    !                         EXIT 
    !                     end if
    !                 end do
                    
    !                 hmin=0
    !                 hh(norb) = vmult(2,norb)+vmult(1,norb)
    !                 do n=(norb-1),1,-(1)
    !                     hh(n)=hh(n+1)*(vmult(2,n)+vmult(1,n))
    !                     if(hh(n)==(0.0,0.0))then
    !                         hmin=n
    !                         EXIT 
    !                     end if
    !                 end do
    !                 totc=(0.0,0.0)
    !                 if (gmax < hmin) then
    !                     totc2=totc2+totc
    !                     cycle
    !                 end if

    !                 if(ph2ei(j,k,l,1).ne.0) then
    !                     totc = totc+(conjg(z1jk(j,k,2,1))*z2l(l,1,1)*hh(2)*ph2ei(j,k,l,1))
    !                 end if

    !                 do n=2,norb-1
    !                     if(ph2ei(j,k,l,n).ne.0.0) then
    !                         totc = totc+ (gg(n-1)*conjg(z1jk(j,k,2,n))*z2l(l,1,n)*hh(n+1)*ph2ei(j,k,l,n))
    !                     end if
    !                 end do

    !                 if(ph2ei(j,k,l,norb).ne.0) then
    !                     totc = totc +(gg(norb-1)*conjg(z1jk(j,k,2,norb))*z2l(l,1,norb)*ph2ei(j,k,l,norb))
    !                 end if
    !                 totc2=totc2+totc
    !             end do
    !         end do
    !     end do
    !    !$omp end target teams distribute parallel do simd
    !     h2etot=totc2
       
    !     if(GDflg.eq.'y')then
    !         if(equal.lt.4)then 
    !             h2etot_diff_ket=0.0
    !             h2etot_diff_bra=0.0
    !             !$omp target teams distribute parallel do reduction(+:total,tot) &
    !             !$omp map(to:gg_1,hh_1,gg_2,hh_2,gmax1,hmin1,gmax2,hmin2,jspin,breakflag) &
    !             !$omp & map(alloc:vmultr(2,norb),vmult_dd(2,norb),occupancy(2,norb),temp(norb)) map(tofrom: total) &
    !             !$omp & private(gg_1,hh_1,gg_2,hh_2,gmax1,gmax2,hmin1,hmin2,jspin,occupancy,breakflag,vmultr,vmult_dd,temp)&
    !             !$omp & shared(ph2ei,z1jk,z2l) &
    !             !$omp & num_teams(max_teams) thread_limit(threadpteam)
    !             do j=1, norb
    !                 if(modulo(j,2)==0)then
    !                     jspin=2
    !                 else
    !                     jspin=1
    !                 end if
    !                 do k=1, norb
    !                     if(j.eq.k)then
    !                         cycle
    !                     end if
    !                     do l=jspin, norb, 2
                            
    !                         occupancy=occupancy_2an(j,k,:,:)*occupancy_an(l,:,:)
    !                          !j,k,l=annihilate1,annihilate1_2,annihilate2
    !                         if(equal.eq.1) then !Differentiation when zombie states are the same
    !                             temp=0
    !                             vmultr=real(z1jk(j,k,:,:))*real(z2l(l,:,:))
    !                             do n=1, norb    !Differentiating w.r.t to orbital n
    !                                 breakflag=0
    !                                 vmult_dd=vmultr
    !                                 if(j.eq.k)then
    !                                     temp(n)=0
    !                                     cycle
    !                                 else if((j.eq.n).or.(k.eq.n))then
    !                                     if(l.eq.n)then    !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(1j)
    !                                         vmult_dd(1,n)=0
    !                                         vmult_dd(2,n)=2*real(psin1(n)*pcos1(n))*occupancy(2,n) !(sin2x=2sinxcosx)
    !                                     else
    !                                         vmult_dd(1,n)=0        !before diff alive:0 dead:a^(a)_(1j)*a^(a)_(0j)
    !                                         vmult_dd(2,n)=(1.0-2*((real(psin1(n)))**2))*occupancy(2,n)    !(cos2x = cos^x -sin^2x)
    !                                     end if
    !                                 else if((j.ne.n).or.(k.ne.n))then
    !                                     if(l.eq.n)then !before diff alive:0 dead:a^(a)_(0j)*a^(a)_(1j)
    !                                         vmult_dd(1,n)=0
    !                                         vmult_dd(2,n)=(1.0-2*((real(psin1(n)))**2))*occupancy(2,n)
    !                                     else
    !                                         breakflag=n !before diff alive:a^(a)_(1j)*a^(a)_(1j) dead:a^(a)_(0j)*a^(a)_(0j)
    !                                     end if          !Unless an opeator acts at position j this evaluates to 0
    !                                 end if
                                
    !                                 gg_1(1:norb)=(0.0,0.0)
    !                                 hh_1(1:norb)=(0.0,0.0)
    !                                 gmax1=norb
    !                                 gg_1(1)=vmult_dd(2,1)-vmult_dd(1,1)
                    
    !                                 do p=2, norb
    !                                     gg_1(p)=gg_1(p-1)*(vmult_dd(2,p)-vmult_dd(1,p))
    !                                     if(gg_1(p)==(0.0))then
    !                                         gmax1=p
    !                                         EXIT 
    !                                     end if
    !                                 end do
                                    
    !                                 hmin1=0
    !                                 hh_1(norb) = vmult_dd(2,norb)+vmult_dd(1,norb)
    !                                 do p=(norb-1),1,-(1)
    !                                     hh_1(p)=hh_1(p+1)*(vmult_dd(2,p)+vmult_dd(1,p))
    !                                     if(hh_1(p)==(0.0))then
    !                                         hmin1=p
    !                                         EXIT 
    !                                     end if
    !                                 end do
                    
    !                                 tot=0.0
    !                                 if (gmax1 < hmin1) then
    !                                     temp(n)=tot
    !                                     cycle
    !                                 end if
                    
    !                                 if((breakflag.eq.0).or.(breakflag.eq.1))then
    !                                     if(ph2ei(j,k,l,1).ne.0) then
    !                                         if(breakflag.eq.1)then
    !                                             temp(n)=(1.0-2*((real(psin1(1)))**2))*&
    !                                             occupancy(2,1)*hh_1(2)*ph2ei(j,k,l,1)
    !                                             cycle
    !                                         end if
    !                                         if(n.eq.1)then
    !                                             if((j.eq.1).or.(k.eq.1))then
    !                                                 if(l.ne.1)then
    !                                                     tot = tot+&
    !                                                     (2*real(psin1(1)*pcos1(1))*occupancy(2,1)*hh_1(2)*ph2ei(j,k,l,1))
    !                                                 end if
    !                                             else
    !                                                 if(l.ne.1)then
    !                                                     tot = tot+&
    !                                                     ((1.0-2*((real(psin1(1)))**2))*occupancy(2,1)*hh_1(2)*ph2ei(j,k,l,1))
    !                                                 end if 
    !                                             end if
    !                                         else
    !                                             tot = tot+&
    !                                             (REAL(conjg(z1jk(j,k,2,1))*z2l(l,1,1))*hh_1(2)*ph2ei(j,k,l,1))
    !                                         end if
    !                                     end if
    !                                 else if((breakflag.lt.norb))then
    !                                     if(breakflag.ne.0)then
    !                                         temp(n)=(gg_1(n-1)*((1.0-2*((real(psin1(n)))**2))*&
    !                                         occupancy(2,n))*hh_1(n+1)*ph2ei(j,k,l,n))
    !                                         cycle
    !                                     end if
    !                                     do p=2,norb-1
    !                                         if(ph2ei(j,k,l,p).ne.0.0) then
    !                                             if(p.eq.n)then
    !                                                 if((j.eq.p).or.(k.eq.p))then
    !                                                     if(l.ne.p)then
    !                                                         tot = tot+&
    !                                                         (gg_1(p-1)*(2*real(psin1(p)*pcos1(p))*&
    !                                                         occupancy(2,p))*hh_1(p+1)*ph2ei(j,k,l,p))
    !                                                     end if
    !                                                 else
    !                                                     if(l.ne.p)then
    !                                                         tot = tot+&
    !                                                         (gg_1(p-1)*((1.0-2*((real(psin1(p)))**2))*&
    !                                                         occupancy(2,p))*hh_1(p+1)*ph2ei(j,k,l,p))
    !                                                     end if 
    !                                                 end if
    !                                             else
    !                                                 tot = tot+&
    !                                                 (gg_1(p-1)*REAL(conjg(z1jk(j,k,2,p))*z2l(l,1,p))*hh_1(p+1)*ph2ei(j,k,l,p))
    !                                             end if
    !                                         end if
    !                                     end do
    !                                 else
    !                                     if(breakflag.eq.norb)then
    !                                         temp(n)=(gg_1(norb-1)*((1.0-2*((real(psin1(norb)))**2))&
    !                                         *occupancy(2,norb))*ph2ei(j,k,l,norb))
    !                                         cycle
    !                                     end if
    !                                     if(ph2ei(j,k,l,norb).ne.0) then
    !                                         if(norb.eq.n)then
    !                                             if((j.eq.norb).or.(k.eq.norb))then
    !                                                 if(l.ne.norb)then
    !                                                     tot = tot+&
    !                                                     (gg_1(norb-1)*(2*real(psin1(norb)*pcos1(norb))*&
    !                                                     occupancy(2,norb))*ph2ei(j,k,l,norb))
    !                                                 end if
    !                                             else
    !                                                 if(l.ne.norb)then
    !                                                     tot = tot+&
    !                                                      (gg_1(norb-1)*((1.0-2*((real(psin1(norb)))**2))*&
    !                                                      occupancy(2,norb))*ph2ei(j,k,l,norb))
    !                                                 end if    
    !                                             end if
    !                                         else
    !                                             tot = tot+&
    !                                             (gg_1(norb-1)*REAL(conjg(z1jk(j,k,2,norb))*z2l(l,1,norb))*ph2ei(j,k,l,norb))
    !                                         end if
    !                                     end if
    !                                 end if
    !                                 temp(n)=tot
                                    
    !                             end do
    !                             total=total+temp(:)
                                
                               
                                
                        
    !                         else if(equal.eq.2)then !Differentiation when zombie states are not the same
    !                             temp=0
    !                             vmultr=real(z1jk(j,k,:,:))*real(z2l(l,:,:))
    !                             do n=1, norb
    !                                 vmult_dd=vmultr
    !                                 if(j.eq.k)then
    !                                     temp(n)=0
    !                                     cycle
    !                                 else if((j.eq.n).or.(k.eq.n))then
    !                                     if(l.eq.n)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
    !                                         vmult_dd(1,n)=0
    !                                         vmult_dd(2,n)=real(pcos1(n)*psin2(n))*occupancy(2,n)
    !                                     else !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(0j)
    !                                         vmult_dd(1,n)=0
    !                                         vmult_dd(2,n)=real(pcos1(n)*pcos2(n))*occupancy(2,n)
    !                                     end if
    !                                 else if((j.ne.n).or.(k.ne.n))then
    !                                     if(l.eq.n)then !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
    !                                         vmult_dd(1,n)=0
    !                                         vmult_dd(2,n)=-real(psin1(n)*psin2(n))*occupancy(2,n)
    !                                     else    !before diff alive:a^(a)_(1j)*a^(b)_(1j) dead:a^(a)_(0j)*a^(b)_(0j)
    !                                         vmult_dd(1,n)=real(pcos1(n)*psin2(n))*occupancy(1,n)
    !                                         vmult_dd(2,n)=-real(psin1(n)*pcos2(n))*occupancy(2,n)
    !                                     end if
    !                                 end if
                    
    !                                 gg_1(1:norb)=(0.0,0.0)
    !                                 hh_1(1:norb)=(0.0,0.0)
    !                                 gmax1=norb
    !                                 gg_1(1)=(vmult_dd(2,1))-(vmult_dd(1,1))
                    
    !                                 do p=2, norb
    !                                     gg_1(p)=gg_1(p-1)*((vmult_dd(2,p)-vmult_dd(1,p)))
    !                                     if(gg_1(p)==(0.0,0.0))then
    !                                         gmax1=p
    !                                         EXIT 
    !                                     end if
    !                                 end do
                    
                                    
    !                                 hmin1=0
                                    
    !                                 hh_1(norb) = (vmult_dd(2,norb)+vmult_dd(1,norb))
                                    
    !                                 do p=(norb-1),1,-(1)
    !                                     hh_1(p)=hh_1(p+1)*(vmult_dd(2,p)+vmult_dd(1,p))
    !                                     if(hh_1(p)==(0.0,0.0))then
    !                                         hmin1=p
    !                                         EXIT 
    !                                     end if
    !                                 end do
                    
    !                                 tot=0.0
    !                                 if((gmax1 < hmin1))then
    !                                     temp(n)=tot
    !                                     cycle
    !                                 end if
                    
    !                                 if(ph2ei(j,k,l,1).ne.0) then
    !                                     if(n.eq.1)then
    !                                         if((j.eq.1).or.(k.eq.1))then
    !                                             if(l.ne.1)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
    !                                                 tot = tot+&
    !                                                 ((REAL(pcos1(1)*psin2(1))*occupancy(2,1))*hh_1(2)*ph2ei(j,k,l,1))
    !                                             end if
    !                                         else
    !                                             if(l.ne.1)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
    !                                                 tot = tot+&
    !                                                 ((REAL(-psin1(1)*psin2(1))*occupancy(2,1))*hh_1(2)*ph2ei(j,k,l,1))
    !                                             end if 
    !                                         end if
    !                                     else
    !                                         tot = tot+&
    !                                         (REAL(conjg(z1jk(j,k,2,1))*z2l(l,1,1))*hh_1(2)*ph2ei(j,k,l,1))
    !                                     end if
    !                                 end if
    !                                 do p=2,norb-1
    !                                     if(ph2ei(j,k,l,p).ne.0.0) then
    !                                         if(p.eq.n)then
    !                                             if((j.eq.p).or.(k.eq.p))then
    !                                                 if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
    !                                                     tot = tot+&
    !                                                     (gg_1(p-1)*(REAL(pcos1(p)*psin2(p))*&
    !                                                     occupancy(2,p))*hh_1(p+1)*ph2ei(j,k,l,p))
    !                                                 end if
    !                                             else
    !                                                 if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
    !                                                     tot = tot+&
    !                                                     (gg_1(p-1)*(REAL(-psin1(p)*psin2(p))*occupancy(2,p))*&
    !                                                     hh_1(p+1)*ph2ei(j,k,l,p))
    !                                                 end if 
    !                                             end if
    !                                         else
    !                                             tot= tot+&
    !                                              (gg_1(p-1)*REAL(conjg(z1jk(j,k,2,p))*z2l(l,1,p))*hh_1(p+1)*ph2ei(j,k,l,p))
    !                                         end if
    !                                     end if
    !                                 end do
                    
    !                                 if(ph2ei(j,k,l,norb).ne.0) then
    !                                     if(norb.eq.n)then
    !                                         if((j.eq.norb).or.(k.eq.norb))then
    !                                             if(l.ne.norb)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
    !                                                 tot =tot+&
    !                                                 (gg_1(norb-1)*(REAL(pcos1(norb)*psin2(norb))*&
    !                                                 occupancy(2,norb))*ph2ei(j,k,l,norb))
    !                                             end if
    !                                         else
    !                                             if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
    !                                                 tot = tot+&
    !                                                 (gg_1(norb-1)*(REAL(-psin1(norb)*psin2(norb))*&
    !                                                 occupancy(2,norb))*ph2ei(j,k,l,norb))
    !                                             end if 
    !                                         end if
    !                                     else
    !                                         tot = tot+&
    !                                         (gg_1(norb-1)*REAL(conjg(z1jk(j,k,2,norb))*z2l(l,1,norb))*ph2ei(j,k,l,norb))
    !                                     end if
    !                                 end if
    !                                 temp(n)=tot
    !                             end do
    !                             total=total+temp
    !                             cycle
    !                         else if(equal.eq.3)then
    !                             temp=0
    !                             vmultr=real(z1jk(j,k,:,:))*real(z2l(l,:,:))
    !                             do n=1, norb
    !                                 vmult_dd=vmultr
    !                                 if(j.eq.k)then
    !                                     temp(n)=0
    !                                     cycle
    !                                 else if((j.eq.n).or.(k.eq.n))then
    !                                     if(l.eq.n)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
    !                                         vmult_dd(1,n)=0
    !                                         vmult_dd(2,n)=real(psin1(n)*pcos2(n))*occupancy(2,n)
    !                                     else !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(0j)
    !                                         vmult_dd(1,n)=0
    !                                         vmult_dd(2,n)=-real(psin1(n)*psin2(n))*occupancy(2,n)
    !                                     end if
    !                                 else if((j.ne.n).or.(k.ne.n))then
    !                                     if(l.eq.n)then !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
    !                                         vmult_dd(1,n)=0
    !                                         vmult_dd(2,n)=real(pcos1(n)*pcos2(n))*occupancy(2,n)
    !                                     else    !before diff alive:a^(a)_(1j)*a^(b)_(1j) dead:a^(a)_(0j)*a^(b)_(0j)
    !                                         vmult_dd(1,n)=real(psin1(n)*pcos2(n))*occupancy(1,n)
    !                                         vmult_dd(2,n)=-real(pcos1(n)*psin2(n))*occupancy(2,n)
    !                                     end if
    !                                 end if
                    
    !                                 gg_2(1:norb)=(0.0,0.0)
    !                                 hh_2(1:norb)=(0.0,0.0)
    !                                 gmax2=norb
    !                                 gg_2(1)=(vmult_dd(2,1))-(vmult_dd(1,1))
                    
    !                                 do p=2, norb
    !                                     gg_2(p)=gg_2(p-1)*((vmult_dd(2,p)-vmult_dd(1,p)))
    !                                     if(gg_2(p)==(0.0,0.0))then
    !                                         gmax2=p
    !                                         EXIT 
    !                                     end if
    !                                 end do
                                    
                    
    !                                 hmin2=0
    !                                 hh_2(norb) = (vmult_dd(2,norb)+vmult_dd(1,norb))
                            
    !                                 do p=(norb-1),1,-(1)
    !                                     hh_2(p)=hh_2(p+1)*(vmult_dd(2,p)+vmult_dd(1,p))
    !                                     if(hh_2(p)==(0.0,0.0))then
    !                                         hmin2=p
    !                                         EXIT 
    !                                     end if
    !                                 end do
    !                                 tot=0.0
    !                                 if((gmax2 < hmin2))then
    !                                     temp(n)=tot
    !                                     cycle
    !                                 end if
                    
    !                                 if(ph2ei(j,k,l,1).ne.0) then
    !                                     if(n.eq.1)then
    !                                         if((j.eq.1).or.(k.eq.1))then
    !                                             if(l.ne.1)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
    !                                                 tot = tot+&
    !                                                 ((REAL(psin2(1)*pcos1(1))*occupancy(2,1))*hh_2(2)*ph2ei(j,k,l,1))
    !                                             end if
    !                                         else
    !                                             if(l.ne.1)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
    !                                                 tot =tot+&
    !                                                 ((REAL(pcos2(1)*pcos1(1))*occupancy(2,1))*hh_2(2)*ph2ei(j,k,l,1))
    !                                             end if 
    !                                         end if
    !                                     else
    !                                         tot = tot+&
    !                                         (REAL(conjg(z1jk(j,k,2,1))*z2l(l,1,1))*hh_2(2)*ph2ei(j,k,l,1))
    !                                     end if
    !                                 end if
    !                                 do p=2,norb-1
    !                                     if(ph2ei(j,k,l,p).ne.0.0) then
    !                                         if(p.eq.n)then
    !                                             if((j.eq.p).or.(k.eq.p))then
    !                                                 if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
    !                                                     tot=tot+&
    !                                                     (gg_2(p-1)*(REAL(psin2(p)*pcos1(p))*&
    !                                                     occupancy(2,p))*hh_2(p+1)*ph2ei(j,k,l,p))
    !                                                 end if
    !                                             else
    !                                                 if(l.ne.p)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
    !                                                     tot = tot+&
    !                                                     (gg_2(p-1)*(REAL(pcos2(p)*pcos1(p))*&
    !                                                     occupancy(2,p))*hh_2(p+1)*ph2ei(j,k,l,p))
    !                                                 end if 
    !                                             end if
    !                                         else
    !                                             tot = tot+&
    !                                             (gg_2(p-1)*REAL(conjg(z1jk(j,k,2,p))*z2l(l,1,p))*hh_2(p+1)*ph2ei(j,k,l,p))
    !                                         end if
    !                                     end if
    !                                 end do
                    
    !                                 if(ph2ei(j,k,l,norb).ne.0) then
    !                                     if(norb.eq.n)then
    !                                         if((j.eq.norb).or.(k.eq.norb))then
    !                                             if(l.ne.norb)then !before diff alive:0 dead:a^(a)_(1j)*a^(b)_(1j)
    !                                                 tot= tot+&
    !                                                 (gg_2(norb-1)*(REAL(psin2(norb)*pcos1(norb))*&
    !                                                 occupancy(2,norb))*ph2ei(j,k,l,norb))
    !                                             end if
    !                                         else
    !                                             if(l.ne.k)then    !before diff alive:0 dead:a^(a)_(0j)*a^(b)_(1j)
    !                                                 tot = tot+&
    !                                                 (gg_2(norb-1)*(REAL(pcos2(norb)*pcos1(norb))*&
    !                                                 occupancy(2,norb))*ph2ei(j,k,l,norb))
    !                                             end if 
    !                                         end if
    !                                     else
    !                                         tot = tot+&
    !                                         (gg_2(norb-1)*REAL(conjg(z1jk(j,k,2,norb))*z2l(l,1,norb))*ph2ei(j,k,l,norb))
    !                                     end if
    !                                 end if

    !                                 temp(n)=tot

    !                             end do
    !                             total=total+temp
    !                             cycle
    !                         end if
                           
    !                     end do
    !                 end do
    !             end do
    !             !$omp end target teams distribute parallel do

               
    !         end if
    !     end if

    !     if(equal.eq.1)then
    !         h2etot_diff_bra=total
    !         h2etot_diff_ket=total
    !     else if(equal.eq.2)then
    !         h2etot_diff_bra=total
    !     else if(equal.eq.3)then
    !         h2etot_diff_ket=total
    !     end if

    !     h2etot=h2etot*0.5
    !     h2etot_diff_bra = h2etot_diff_bra*0.5
    !     h2etot_diff_ket = h2etot_diff_ket*0.5
    !     return

    ! end subroutine two_elec_part_gpu
    
    ! ! Top level routine to allocate hamiltonian and overlap matrices 
    ! subroutine hamgen_gpu(ham,zstore,elecs,size,verb)

    !     implicit none

        
    !     type(hamiltonian), intent(inout)::ham 
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     type(elecintrgl),intent(in)::elecs
    !     integer,allocatable,dimension(:,:,:)::occupancy_an
    !     integer,allocatable,dimension(:,:,:,:)::occupancy_an_cr,occupancy_2an
    !     complex(kind=8),allocatable, dimension(:,:,:)::passback
    !     integer, allocatable,dimension(:)::IPIV1
    !     integer,intent(in)::size,verb
    !     complex(kind=8),allocatable,dimension(:)::WORK1
    !     real(kind=8),dimension(ndet,ndet,norb)::temp2
        
    !     real(kind=8),dimension(norb,norb)::ph1ei
    !     real(kind=8),dimension(norb,norb,norb,norb)::ph2ei
    !     real(kind=8) :: phnuc

    !     complex(kind=8), dimension(ndet,norb)::psin
    !     complex(kind=8), dimension(ndet,norb)::pcos
    !     real(kind=8),dimension(ndet,norb)::pphi

    !     complex(kind=8), dimension(ndet,ndet)::phjk
    !     complex(kind=8), dimension(ndet,ndet)::povrlp
    !     complex(kind=8), dimension(ndet,ndet)::pinv
    !     complex(kind=8), dimension(ndet,ndet)::pkinvh

    !     real(kind=8), dimension(ndet,ndet,norb)::pdiff_hjk
    !     real(kind=8), dimension(ndet,ndet,norb)::pdiff_ovrlp
    !     real(kind=8), dimension(ndet,ndet,ndet,norb)::pdiff_invh
    !     integer:: j,k,l,ierr

       
    !     if (errorflag .ne. 0) return

    !     allocate(occupancy_an(norb,2,norb),stat=ierr)
    !     if(ierr==0) allocate(occupancy_2an(norb,norb,2,norb),stat=ierr)
    !     if(ierr==0) allocate(occupancy_an_cr(norb,norb,2,norb),stat=ierr)
    !     if(ierr==0) allocate(passback(norb,2,norb),stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in occupancy vector allocation . ierr had value ", ierr
    !         errorflag=1
    !         return
    !     end if 
    !     allocate(IPIV1(size),stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
    !         errorflag=1
    !     end if 
    !     if (ierr==0) allocate(WORK1(size),stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
    !         errorflag=1
    !     end if  
       
        
    !     occupancy_an=1
    !     occupancy_2an=1
    !     occupancy_an_cr=1
    !     !$omp parallel
    !     !$omp do
    !     do j=1,norb
    !         occupancy_an(j,1,j)=0
    !         do l=j-1, 1, -1
    !             occupancy_an(j,1,l)=-1
    !         end do
    !     end do
    !     !$omp end do
    !     !$omp  do
    !     do j=1,norb
    !         do k=1, norb
    !             occupancy_2an(j,k,:,:)=occupancy_an(j,:,:)
    !             occupancy_2an(j,k,2,k)=occupancy_2an(j,k,1,k)
    !             occupancy_2an(j,k,1,k)=0
    !             occupancy_an_cr(j,k,:,:)=occupancy_an(j,:,:)
    !             occupancy_an_cr(j,k,1,k)=occupancy_an_cr(j,k,2,k)
    !             occupancy_an_cr(j,k,2,k)=0
    !             do l=k-1,1,-1
    !                 occupancy_2an(j,k,1,l)=occupancy_2an(j,k,1,l)*(-1)
    !                 occupancy_an_cr(j,k,1,l)=occupancy_an_cr(j,k,1,l)*(-1)
    !             end do
    !         end do
    !     end do
    !     !$omp end do
    !     !$omp end parallel 
        
    
    !     ph1ei=elecs%h1ei
    !     ph2ei=elecs%h2ei
    !     phnuc=elecs%hnuc 

    !     phjk=ham%hjk
    !     povrlp=ham%ovrlp
    !     pinv=ham%inv
    !     pkinvh=ham%kinvh
    !     if(GDflg.eq.'y')then
    !         pdiff_hjk=ham%diff_hjk
    !         pdiff_ovrlp=ham%diff_ovrlp
    !         pdiff_invh=ham%diff_invh
    !     end if
        

    !     do j=1,ndet
    !         psin(j,:)=zstore(j)%sin(:)
    !         pcos(j,:)=zstore(j)%cos(:)
    !         pphi(j,:)=zstore(j)%phi(:)
    !     end do
        
       
    !     !$omp target data map(to:ph1ei,ph2ei,phnuc,psin,pcos,pphi,occupancy_2an,&
    !     !$omp occupancy_an_cr,occupancy_an,passback,errorflag,ierr,GDflg) &
    !     !$omp & map(tofrom:phjk,povrlp) if(GDflg.eq.'y') map(tofrom:pdiff_hjk,pdiff_ovrlp)
    !     do j=1, size
    !         call he_row_gpu(ph1ei,ph2ei,phnuc,psin,pcos,pphi,occupancy_2an,&
    !         & occupancy_an_cr,occupancy_an,passback,phjk,povrlp,pdiff_hjk,pdiff_ovrlp,j,size)
    !         if(verb.eq.1)then
    !             write(6,"(a,i0,a)") "Hamiltonian row ",j, " completed"
    !         end if  
    !     end do
    !     !$omp end target data
        
    !     ham%hjk=phjk
    !     ham%ovrlp=povrlp
    !     if(GDflg.eq.'y') then 
    !     ham%diff_hjk=pdiff_hjk
    !     ham%diff_ovrlp=pdiff_ovrlp
    !     end if
         

    !     pinv=povrlp
        
       
    !     if (ierr==0) call ZGETRF(size,size,pinv,size,IPIV1,ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)")"Error in ZGETRF",ierr
    !     end if
    !     if (ierr==0) call ZGETRI(size,pinv,size,IPIV1,WORK1,size,ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)")"Error in ZGETRF",ierr
    !     end if
        
    !     ham%inv=pinv
    !     !$omp parallel
    !     !$omp workshare
    !     pkinvh=matmul(pinv,phjk)
    !     !$omp end workshare
    !     !$omp end parallel
    !     ham%kinvh=pkinvh
        
    !     if(GDflg.eq.'y')then
    !         !$omp parallel shared(temp2,ham) private(k,l)
    !         !$omp do simd collapse(2)
    !         do k=1, ndet
    !             do l=1, ndet
    !                 if(l.eq.2)then
    !                     temp2(k,2,:)=matmul(REAL(ham%inv(k,:)),ham%diff_ovrlp(2,:,:))
    !                 else
    !                     temp2(k,l,:)=real(ham%inv(k,l))*ham%diff_ovrlp(2,l,:)
    !                 end if
    !             end do
    !         end do
    !         !$omp end do simd
    !         !$omp do simd collapse(2)
    !         do k=1, ndet
    !             do l=1, ndet
    !                 ham%diff_invh(2,k,l,:)=matmul(transpose(temp2(k,:,:)),real(ham%kinvh(:,l)))*(-1)
    !             end do
    !         end do
    !         !$omp end do simd
    !         !$omp end parallel
    !     end if
        

    !     deallocate(occupancy_an,stat=ierr)
    !     if(ierr==0) deallocate(passback,stat=ierr)
    !     if(ierr==0) deallocate(occupancy_2an,stat=ierr)
    !     if(ierr==0) deallocate(occupancy_an_cr,stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in passback vector deallocation . ierr had value ", ierr
    !         errorflag=1
    !     end if

    !     if (ierr==0) deallocate(IPIV1,stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in IPIV vector deallocation . ierr had value ", ierr
    !         errorflag=1
    !     end if

    !     if (ierr==0) deallocate(WORK1,stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in WORK vector deallocation . ierr had value ", ierr
    !         errorflag=1
    !     end if
        
    !     return
       

    ! end subroutine hamgen_gpu
    
    
!END MODULE ham