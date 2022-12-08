MODULE ham

    use globvars
    use alarrays
    use operators
    use omp_lib

    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators
    subroutine he_row_gpu(ham,zstore,elecs,row,size,occupancy_2an,occupancy_an_cr,occupancy_an,passback)
        
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
          
         
            call one_elec_part_gpu(zstore(row)%sin,zstore(row)%cos,z2l,h1etot,occupancy_an_cr,&
                elecs%h1ei,h1etot_diff_bra,h1etot_diff_ket,zstore(m)%sin,zstore(m)%cos,equal)
            
           
            z2l=z2l*occupancy_an
            !!$omp flush(z2l)
         
            call two_elec_part_gpu(zstore(row)%sin,zstore(row)%cos,z1jk,z2l,h2etot,occupancy_2an,occupancy_an,&
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
                        overlap_diff = diff_overlap_gpu(zstore(row)%sin,zstore(row)%cos,zstore(m)%sin,zstore(m)%cos,equal)
                        ham%diff_ovrlp(row,m,:) = overlap_diff
                    end if
                else if(m.eq.2)then
                    ham%diff_hjk(m,row,:)=h1etot_diff_ket+h2etot_diff_ket
                    overlap_diff = diff_overlap_gpu(zstore(row)%sin,zstore(row)%cos,zstore(m)%sin,zstore(m)%cos,equal)
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

    end subroutine he_row_gpu

    subroutine one_elec_part_gpu(zs1sin,zs1cos,z2l,h1etot,occupancy,h1ei,h1etot_diff_bra,h1etot_diff_ket,zs2sin,zs2cos,equal)

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
            !$omp target teams distribute parallel do simd collapse(2) &
            !$omp & map(to:h1ei(:,:),occupancy(:,:,:,:),z2l(:,:,:),zs1sin(:),zs1cos(:),zs2sin(:),zs2cos(:),equal,len) &
            !$omp & map(tofrom:temp(:,:),h1etot_diff(:,:,:)) map(alloc:prod(len),chng_prod(len),temp_prod(len),zomt(2,len)) &
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
            !$omp end target teams distribute parallel do simd         
        else
            !$omp target teams distribute parallel do simd collapse(2) &
            !$omp & map(to:h1ei(:,:),occupancy(:,:,:,:),z2l(:,:,:),zs1sin(:),zs1cos(:),len) &
            !$omp & map(tofrom:temp(:,:)) map(alloc:zomt(2,len)) &
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
            !$omp end target teams distribute parallel do simd
        end if

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

    end subroutine one_elec_part_gpu

    ! complex(kind=8) function one_elec_body_gpu(len,zs1sin,zs1cos,z2l,occupancy,h1ei,k)

    !     implicit none
        
    !     complex(kind=8),dimension(:)::zs1sin,zs1cos
    !     integer,dimension(:,:),intent(in)::occupancy
    !     complex(kind=8),dimension(:,:),intent(in)::z2l
    !     real(kind=8), intent(in)::h1ei
    !     integer,intent(in)::k,len
    !     complex(kind=8),allocatable,dimension(:,:)::zomt
    !     complex(kind=8),allocatable,dimension(:)::zomt2
    !     integer::ierr
    !     !$omp declare target
    !     !write(0,"(a)") " one elec body allocate"
    !     allocate(zomt(2,len),stat=ierr)
    !     allocate(zomt2(len),stat=ierr)
      
    !     zomt=z2l
    !     zomt(1,k)=zomt(2,k)
    !     zomt(2,k)=cmplx(0.0,0.0)
    !     zomt=zomt*occupancy
    !     one_elec_body_gpu=product((conjg(zs1sin)*zomt(1,:))+(conjg(zs1cos)*zomt(2,:)))*h1ei
        
    !     deallocate(zomt,stat=ierr)
    !     deallocate(zomt2,stat=ierr)
    !     return

    ! end function one_elec_body_gpu

    ! function one_elec_body_grad_gpu(len,zs1sin,zs1cos,zs2sin,zs2cos,z2l,occupancy,h1ei,k,j,equal)

    !     implicit none
    !     !$omp declare target
    !     complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
    !     integer,dimension(:,:),intent(in)::occupancy
    !     complex(kind=8),dimension(:,:),intent(in)::z2l
    !     real(kind=8), intent(in)::h1ei
    !     integer,intent(in)::k,j,equal,len
    !     real(kind=8),dimension(len)::one_elec_body_grad_gpu
    !     complex(kind=8),allocatable,dimension(:,:)::zomt
    !     integer::ierr

    !     !write(0,"(a)") " one elec body grDallocate"
    !     ! allocate(one_elec_body_grad_gpu(len),stat=ierr)
    !     allocate(zomt(2,len),stat=ierr)
    !     zomt=z2l
    !     zomt(1,k)=zomt(2,k)
    !     zomt(2,k)=cmplx(0.0,0.0)
    !     zomt=zomt*occupancy
    !     one_elec_body_grad_gpu = diff_overlap_cran_gpu(len,zs1sin,zs1cos,zs2sin,zs2cos,equal,zomt,j,k,occupancy)*h1ei
    !     !write(0,"(a)") " one elec body GRAD done"

    !     deallocate(zomt,stat=ierr)
    !     return

    ! end function one_elec_body_grad_gpu
        
    subroutine two_elec_part_gpu(zs1sin,zs1cos,z1jk,z2l,h2etot,occupancy_2an,occupancy_an,h2ei,&
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
            tot(j,:,:) = two_elec_part_body_gpu(norb,zs1sin,zs2sin,z2l,z1jk(j,:,:,:),h2ei(j,:,:,:),j)
            if(equal.lt.4)then
                h2etot_diff(j,:,:,:) = two_elec_part_grad_gpu(norb,zs1sin,zs1cos,zs2sin,zs2cos,z2l,z1jk(j,:,:,:),&
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

    end subroutine two_elec_part_gpu

    function two_elec_part_body_gpu(len,zs1sin,zs2sin,z2l,z1jk,h2ei,j)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs2sin
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:),intent(in)::z1jk
        real(kind=8), dimension(:,:,:), intent(in)::h2ei
        complex(kind=8),dimension(:,:),allocatable::vmult,tot
        complex(kind=8),dimension(len,len)::two_elec_part_body_gpu
        complex(kind=8),allocatable,dimension(:)::gg,hh
        integer,intent(in)::len
        integer::j,k,l,p,jspin,gmax,hmin,ierr

        ierr=0

        allocate(tot(len,len),stat=ierr)
        ! allocate(two_elec_part_body_gpu(len,len),stat=ierr)
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
            two_elec_part_body_gpu=(0.0,0.0)
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
          
            !$omp target teams distribute parallel do simd &
            !$omp & map(alloc:hh(len),gg(len),vmult(2,len)) map(tofrom:tot(k,:)) &
            !$omp & map(to:h2ei,z1jk(k,:,:),z2l,zs2sin,gmax,hmin,len,jspin) &
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
            !$omp end target teams distribute parallel do simd
        end do


        deallocate(vmult)
        deallocate(gg,stat=ierr)
        deallocate(hh,stat=ierr)

        two_elec_part_body_gpu=tot

        deallocate(tot,stat=ierr)
        
        return

    end function two_elec_part_body_gpu

    function two_elec_part_grad_gpu(len,zs1sin,zs1cos,zs2sin,zs2cos,z2l,z1jk,h2ei,occupancy_2an,occupancy_an,j,equal)

        implicit none

        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        complex(kind=8),dimension(:,:,:),intent(in)::z2l
        complex(kind=8),dimension(:,:,:),intent(in)::z1jk
        real(kind=8), dimension(:,:,:), intent(in)::h2ei
        integer,dimension(:,:,:),intent(in)::occupancy_2an
        integer,dimension(:,:,:),intent(in)::occupancy_an
        integer,intent(in)::j,equal,len
        real(kind=8),dimension(len,len,len)::two_elec_part_grad_gpu
        integer,allocatable,dimension(:,:)::occupancy
        integer::k,l,jspin

        allocate(occupancy(2,len))
       
        two_elec_part_grad_gpu=0.0

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
                two_elec_part_grad_gpu(k,l,:) = z_an_z3_diff_gpu(len,z1jk(k,:,:),z2l(l,:,:),h2ei(k,l,:),&
                equal,occupancy,j,k,l,zs1sin,zs1cos,zs2sin,zs2cos)
            end do
        end do

        deallocate(occupancy)
        
        return

    end function two_elec_part_grad_gpu

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
        !$omp do simd collapse(2)
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
            call he_row_gpu(ham,zstore,elecs,j,size,occupancy_2an,occupancy_an_cr,occupancy_an,passback)
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

    ! function diff_overlap_cran_gpu(len,zs1sin,zs1cos,zs2sin,zs2cos,dtype,zomt,annihilate2,create2,occupancy)

    !     implicit none
    !     !$omp declare target
    !     integer,intent(in)::len
    !     complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
    !     real(kind=8),dimension(len)::diff_overlap_cran_gpu
    !     complex(kind=8),dimension(:,:)::zomt
    !     integer,intent(in)::dtype,annihilate2,create2
    !     integer,dimension(:,:),intent(in)::occupancy
    !     real(kind=8),allocatable,dimension(:)::bra_prod,ket_prod,prod
    !     integer::j,k
    !     real(kind=8),allocatable,dimension(:)::temp,temp2

    !     !if (errorflag .ne. 0) return
       
        
    !     allocate(prod(len))
    !     prod=real(((conjg(zs1sin)*zomt(1,:)))+((conjg(zs1cos)*zomt(2,:))))
    !     diff_overlap_cran_gpu=0
    !     ! Differntial of overlap of the same ZS 0 unless an operator has acted on ZS
    !     if(dtype.eq.1)then
    !         allocate(bra_prod(len))
    !         if(annihilate2.eq.create2)then
    !             if((real(zs2cos(annihilate2)).eq.0).and.(real(zs2sin(annihilate2)).eq.1).or.(real(zs2sin(annihilate2)).eq.0))then
    !                 diff_overlap_cran_gpu(annihilate2)=0
    !             else
    !                 bra_prod=prod           !dead amplitude is zero
    !                 bra_prod(annihilate2)=real(2*zs1sin(annihilate2)*zs1cos(annihilate2)*occupancy(1,annihilate2))
                    
    !                 diff_overlap_cran_gpu(annihilate2)=product(bra_prod)   
    !             end if 
    !         else if(annihilate2.ne.create2)then
    !             bra_prod=prod
    !             if((real(zs2cos(annihilate2)).eq.0).and.(real(zs2sin(annihilate2)).eq.1)) then
    !                 bra_prod(annihilate2)=-1*occupancy(2,annihilate2)
    !             else if((real(zs2sin(annihilate2)).eq.0).and.(real(zs2cos(annihilate2)).eq.1)) then
    !                 bra_prod(annihilate2)=occupancy(2,annihilate2)
    !             else
    !                 bra_prod(annihilate2)=real(((zs1cos(annihilate2)**2)-(zs1sin(annihilate2)**2))*occupancy(2,annihilate2))
    !             end if
    !             diff_overlap_cran_gpu(annihilate2)=product(bra_prod)
                
    !             bra_prod=prod
    !             if((real(zs2cos(create2)).eq.0).and.(real(zs2sin(create2)).eq.1)) then
    !                 bra_prod(create2)=-1*occupancy(1,create2)
    !             else if((real(zs2sin(create2)).eq.0).and.(real(zs2cos(create2)).eq.1)) then
    !                 bra_prod(create2)=occupancy(1,create2)
    !             else
    !                 bra_prod(create2)=real(((zs1cos(create2)**2)-(zs1sin(create2)**2))*occupancy(1,create2))
    !             end if               !dead amplitude is zero
    !             diff_overlap_cran_gpu(create2)=product(bra_prod) 
    !         end if
    !         deallocate(bra_prod)
    !     else if(dtype.eq.2)then 
    !         allocate(temp(len))
    !         allocate(bra_prod(len))
    !         bra_prod=real(zs1cos*zs2sin*occupancy(1,:)-zs1sin*zs2cos*occupancy(2,:))            
    !         if(annihilate2.eq.create2)then  !dead amplitude is zero
    !             bra_prod(annihilate2)=real(zs1cos(annihilate2)*zs2sin(annihilate2)*occupancy(1,annihilate2))
    !         else                   
    !             bra_prod(annihilate2)=-real(zs1sin(annihilate2)*zs2sin(annihilate2))*occupancy(2,annihilate2)
    !             bra_prod(create2)=real(zs1cos(create2)*zs2cos(create2))*occupancy(1,create2) 
    !         end if
    !         !!$omp parallel do simd shared(prod,bra_prod) private(temp,j)
    !         do j=1,len
    !             temp=prod
    !             temp(j)=bra_prod(j)
    !             diff_overlap_cran_gpu(j)=product(temp)
                
    !         end do 
    !         deallocate(temp)
    !         deallocate(bra_prod)
    !         !!$omp end parallel do simd 
    !     else if(dtype.eq.3)then
    !         allocate(temp2(len))
    !         allocate(ket_prod(len))
    !         ket_prod=real(zs1sin*zs2cos*occupancy(1,:)-zs1cos*zs2sin*occupancy(2,:))                  
    !         if(annihilate2.eq.create2)then  !dead amplitude is zero
    !             ket_prod(annihilate2)=real(zs1sin(annihilate2)*zs2cos(annihilate2)*occupancy(1,annihilate2))
    !         else                   
    !             ket_prod(annihilate2)=real(zs1cos(annihilate2)*zs2cos(annihilate2))*occupancy(2,annihilate2)!alive amplitude is zero
    !             ket_prod(create2)=-real(zs1sin(create2)*zs2sin(create2))*occupancy(1,create2) !dead amplitude is zero
    !         end if
    !         !!$omp parallel do simd  shared(prod,ket_prod) private(temp2,j)
    !         do j=1,len
    !             temp2=prod
    !             temp2(j)=ket_prod(j)
    !             diff_overlap_cran_gpu(j)=product(temp2)
    !         end do
    !         !!$omp end parallel do simd   
    !         deallocate(temp2)
    !         deallocate(ket_prod)
    !     end if

    !     deallocate(prod)
        
    !     RETURN


    ! end function diff_overlap_cran_gpu

     ! computes the vector of values formed by the derivative of the overlap with respect to each orbital. 
    ! Does not have capability to deal with states where creation and annihilation operators have acted
    function diff_overlap_gpu(zs1sin,zs1cos,zs2sin,zs2cos,dtype)
    
        implicit none
      
        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        integer,intent(in)::dtype
        real(kind=8),dimension(norb)::diff_overlap_gpu
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
            !!$omp target teams distribute parallel do simd &
            !!$omp map(to:prod,bra_prod) map(tofrom:diff_overlap_gpu) map(alloc:temp(norb)) &
            !$omp parallel do simd shared(prod,bra_prod) private(temp)
            do j=1,norb
                temp=prod
                temp(j)=bra_prod(j)
                diff_overlap_gpu(j)=product(temp)
            end do
            !$omp end parallel do simd  
            !!$omp end target teams distribute parallel do simd  
            deallocate(bra_prod)
            deallocate(temp)
        else if(dtype.eq.3)then
            allocate(ket_prod(norb))
            allocate(temp2(norb))
            ket_prod=real(conjg(zs1sin)*zs2cos)-real(conjg(zs1cos)*zs2sin)
            !!$omp target teams distribute parallel do simd &
            !!$omp map(to:prod,ket_prod) map(tofrom:diff_overlap_gpu) map(alloc:temp2(norb)) &
            !$omp parallel do simd shared(prod,ket_prod) private(temp2)
            do j=1,norb
                temp2=prod
                temp2(j)=ket_prod(j)
                diff_overlap_gpu(j)=product(temp2)
            end do
            !$omp end parallel do simd  
            !!$omp end target teams distribute parallel do simd  
            deallocate(ket_prod)
            deallocate(temp2)
        end if
        deallocate(prod)
        return
        
    end function diff_overlap_gpu

    ! complex(kind=8) function z_an_z3_gpu(len,z1,z2,vec)

    !     implicit none
    !     !!$omp declare target
    !     integer, intent(in)::len
    !     complex(kind=8),dimension(:,:),intent(in)::z1,z2
    !     complex(kind=8),dimension(:,:),allocatable::vmult
    !     real(kind=8),dimension(:),intent(in)::vec
    !     complex(kind=8),allocatable,dimension(:)::gg,hh
    !     complex(kind=8)::tot
    !     integer::j,gmax,hmin,ierr

    !     if (errorflag .ne. 0) return
    !     ierr=0

    !     allocate(gg(len),stat=ierr)
    !     allocate(hh(len),stat=ierr)
    !     allocate(vmult(2,len),stat=ierr)
    !     if (ierr/=0) then
    !         write(0,"(a,i0)") "Error in vmult allocation . ierr had value ", ierr
    !         errorflag=1
    !         return
    !     end if 
        
    !     vmult=conjg(z1)*(z2)
        
    
    !     gg(1:len)=(0.0,0.0)
    !     hh(1:len)=(0.0,0.0)
    !     gmax=len
    !     gg(1)=vmult(2,1)-vmult(1,1)

    !     do j=2, len
    !         gg(j)=gg(j-1)*(vmult(2,j)-vmult(1,j))
    !         if(gg(j)==(0.0,0.0))then
    !             gmax=j
    !             EXIT 
    !         end if
    !     end do
        
    !     hmin=0
    !     hh(len) = vmult(2,len)+vmult(1,len)
    !     do j=(len-1),1,-(1)
    !         hh(j)=hh(j+1)*(vmult(2,j)+vmult(1,j))
    !         if(hh(j)==(0.0,0.0))then
    !             hmin=j
    !             EXIT 
    !         end if
    !     end do


    !     tot=(0.0,0.0)
    !     if (gmax < hmin) then
    !         z_an_z3_gpu=tot
    !         return
    !     end if

    !     if(vec(1).ne.0) then
    !         tot = tot+(conjg(z1(2,1))*z2(1,1)*hh(2)*vec(1))
    !     end if

    !     !$omp parallel do simd reduction(+:tot) shared(hh,vec,gg,z1,z2)  
    !     do j=2,len-1
    !         if(vec(j).ne.0.0) then
    !             tot = tot+ (gg(j-1)*conjg(z1(2,j))*z2(1,j)*hh(j+1)*vec(j))
    !         end if
    !     end do
    !     !$omp end parallel do simd 

    !     if(vec(len).ne.0) then
    !         tot = tot +(gg(len-1)*conjg(z1(2,len))*z2(1,len)*vec(len))
    !     end if

    !     deallocate(vmult)
    !     deallocate(gg,stat=ierr)
    !     deallocate(hh,stat=ierr)
    
    !     z_an_z3_gpu=tot
    
    !     return

    ! end function z_an_z3_gpu

    function z_an_z3_diff_gpu(len,z1,z2,vec,dtype,occupancy,annihilate1,annihilate1_2,annihilate2,zs1sin,zs1cos,zs2sin,zs2cos)

        implicit none
        integer, intent(in)::len
        complex(kind=8),dimension(:,:),intent(in)::z1,z2
        real(kind=8),dimension(:,:),allocatable::vmult_dd,vmult_1d, vmult_2d,vmult
        real(kind=8),dimension(len),intent(in)::vec
        integer,dimension(:,:),intent(in)::occupancy
        integer, intent(in)::dtype,annihilate1,annihilate1_2,annihilate2
        complex(kind=8),dimension(:)::zs1sin,zs1cos,zs2sin,zs2cos
        real(kind=8),dimension(len)::z_an_z3_diff_gpu
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
       
        z_an_z3_diff_gpu=0
     
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
            !$omp target teams distribute parallel do &
            !$omp map(to:vmult,vec,occupancy,zs1cos,zs1sin,z2,hh_0,z1,gg_0,annihilate1,annihilate1_2,annihilate2,&
            !$omp & breakflag,tot0,gmax0,hmin0,len) &
            !$omp & map(tofrom:temp0) map(alloc:vmult_dd(2,len)) &
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
            !$omp end target teams distribute parallel do

            z_an_z3_diff_gpu=temp0(:)
            z_an_z3_diff_gpu=temp0(:)
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
            !$omp target teams distribute parallel do &
            !$omp map(to:vmult,vec,occupancy,zs1cos,zs1sin,zs2sin,zs2cos,z2,hh_1,gg_1,z1,annihilate1,annihilate1_2,&
            !$omp & annihilate2,tot1,gmax1,hmin1,len) &
            !$omp & map(tofrom:temp) map(alloc:vmult_1d(2,len)) &
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
            !$omp end target teams distribute parallel do
            z_an_z3_diff_gpu=temp(:)
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
            !$omp target teams distribute parallel do &
            !$omp map(to:vmult,vec,occupancy,zs1cos,zs1sin,zs2sin,zs2cos,z2,hh_2,gg_2,z1,annihilate1,annihilate1_2,&
            !$omp & annihilate2,tot2,gmax2,hmin2,len) &
            !$omp & map(tofrom:temp2) map(alloc:vmult_2d(2,len)) &
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
            !$omp end target teams distribute parallel do
     
            z_an_z3_diff_gpu=temp2(:)
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

    end function z_an_z3_diff_gpu

END MODULE ham