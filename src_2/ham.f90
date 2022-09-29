MODULE ham 

    use globvars
    use alarrays
    use operators

    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators
    subroutine he_row(ham,zstore,elecs,row,size)
        
        implicit none
        type(hamiltonian), intent(inout)::ham
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::row,size
        type(zombiest),allocatable, dimension(:,:)::z1jk
        type(zombiest),allocatable, dimension(:)::z2l
        type(zombiest)::zomt
        integer::j,k,l,m,p,jspin,ierr
        complex(kind=8)::h1etot, h2etot,temp
        real(kind=8),dimension(norb)::h1etot_diff_bra,h2etot_diff_bra,diff_temp_bra, diff_temp2_bra
        real(kind=8),dimension(norb)::h1etot_diff_ket,h2etot_diff_ket,diff_temp_ket, diff_temp2_ket  

        if (errorflag .ne. 0) return
        ierr = 0

        call alloczf(zomt)
        call alloczs2d(z1jk,norb)
        call alloczs(z2l,norb)
        h1etot=(0.0,0.0)
        h2etot=(0.0,0.0)
        temp=(0.0,0.0)
        h1etot_diff_bra(1:norb)=0.0
        h1etot_diff_ket(1:norb)=0.0
        diff_temp_bra(1:norb)=0.0
        diff_temp_ket(1:norb)=0.0
        diff_temp2_bra(1:norb)=0.0
        diff_temp2_ket(1:norb)=0.0
        !$omp parallel private(j,k,l,jspin,zomt,z2l,h1etot,h2etot,temp, &
        !$omp   h1etot_diff_bra,h2etot_diff_bra, diff_temp_bra, diff_temp2_bra, &
        !$omp   h1etot_diff_ket,h2etot_diff_ket,diff_temp_ket, diff_temp2_ket) shared(z1jk,zstore,row,size)
        !$omp do
        do j=1, norb
            do k=1, norb
                zomt=zstore(row)
                call an(zomt,j)
                call an(zomt,k)
                z1jk(j,k)=zomt
            end do
        end do
        !$omp end do
        !$omp flush(z1jk)
        !$omp do
        do m=row,size
            h1etot=(0.0,0.0)
            h2etot=(0.0,0.0)
            temp=(0.0,0.0)
            h1etot_diff_bra(1:norb)=0.0
            h1etot_diff_ket(1:norb)=0.0
            diff_temp_bra(1:norb)=0.0
            diff_temp_ket(1:norb)=0.0
            diff_temp2_bra(1:norb)=0.0
            diff_temp2_ket(1:norb)=0.0
            
            do l=1, norb
                zomt=zstore(m)
                call an(zomt,l)
                z2l(l)=zomt
            end do
            
            do j=1, norb
                do k=1, norb
                    temp=(0.0,0.0)
                    zomt=z2l(j)
                    call cr(zomt,k)
                    temp = (overlap(zstore(row),zomt)*elecs%h1ei(j,k))
                    h1etot=h1etot+temp
                    if(GDflg.eq.'y')then
                        if(m.eq.row)then
                            diff_temp_bra = diff_overlap(zstore(row)%diffalive,zstore(row)%diffdead,zomt%diffalive,zomt%diffdead)
                            diff_temp_ket = diff_temp_bra
                            h1etot_diff_bra= h1etot_diff_bra + (diff_temp_bra*elecs%h1ei(j,k))
                            h1etot_diff_ket= h1etot_diff_ket + (diff_temp_ket*elecs%h1ei(j,k))
                        else 
                            diff_temp_bra = diff_overlap(zstore(row)%diffalive,zstore(row)%diffdead,&
                                                            REAL(zomt%alive),REAL(zomt%dead))
                            diff_temp_ket = diff_overlap(REAL(zstore(row)%alive),REAL(zstore(row)%dead),&
                                                                zomt%diffalive,zomt%diffdead)
                            h1etot_diff_bra= h1etot_diff_bra + (diff_temp_bra*elecs%h1ei(j,k))
                            h1etot_diff_ket= h1etot_diff_ket + (diff_temp_ket*elecs%h1ei(j,k))
                        end if
                    end if
                end do
            end do
            
            do j=1, norb
                if(zstore(row)%alive(j)==(0.0,0.0))then
                    CYCLE
                end if

                if(modulo(j,2)==0)then
                    jspin=2
                else
                    jspin=1
                end if

                do k=1, norb
                    if(iszero(z1jk(j,k)).eqv..true.)then
                        CYCLE
                    end if
                    do l=jspin, norb, 2
                        if(zstore(m)%alive(l)==(0.0,0.0))then
                            CYCLE
                        end if
                        h2etot = h2etot + z_an_z3(z1jk(j,k),z2l(l),elecs%h2ei(j,k,l,:))
                    end do
                end do
            end do
            h2etot=h2etot*0.5
            if(GDflg.eq.'y')then
                do j=1, norb
                    do k=1, norb
                        do l=1, norb 
                            zomt=z2l(l)
                            if(m.eq.row)then
                                diff_temp_bra = diff_overlap(z1jk(j,k)%diffalive,&
                                    z1jk(j,k)%diffdead,zomt%diffalive,zomt%diffdead)
                                diff_temp_ket=diff_temp_bra
                            else 
                                diff_temp_bra = diff_overlap(z1jk(j,k)%diffalive,&
                                    z1jk(j,k)%diffdead,REAL(zomt%alive),REAL(zomt%dead))
                                diff_temp_ket = diff_overlap(REAL(z1jk(j,k)%alive),&
                                    REAL(z1jk(j,k)%dead),zomt%diffalive,zomt%diffdead)
                            end if
                            do p=1, norb
                                diff_temp2_bra=diff_temp_bra
                                diff_temp2_ket=diff_temp_ket
                                if(m.eq.row)then
                                    diff_temp2_bra(p)=z1jk(j,k)%diffdead(p)*zomt%diffalive(p)
                                    diff_temp2_ket(p)=z1jk(j,k)%diffdead(p)*zomt%diffalive(p)
                                    h2etot_diff_bra= h2etot_diff_bra + (diff_temp2_bra*elecs%h2ei(j,k,l,p))
                                    h2etot_diff_ket= h2etot_diff_ket + (diff_temp2_bra*elecs%h2ei(j,k,l,p))
                                else 
                                    diff_temp2_bra(p) = z1jk(j,k)%diffdead(p)*REAL(zomt%alive(p))
                                    h2etot_diff_bra= h2etot_diff_bra + (diff_temp2_bra*elecs%h2ei(j,k,l,p))
                                    diff_temp2_ket(p) = REAL(z1jk(j,k)%dead(p))*zomt%diffalive(p)
                                    h2etot_diff_ket= h2etot_diff_ket + (diff_temp2_ket*elecs%h2ei(j,k,l,p))
                                end if
                            end do
                        end do
                    end do
                end do
                h2etot_diff_bra=h2etot_diff_bra*0.5
                h2etot_diff_ket=h2etot_diff_ket*0.5
            end if
            
            ham%ovrlp(row,m)=overlap(zstore(row),zstore(m))
            ham%ovrlp(m,row)= ham%ovrlp(row,m)
            ham%hjk(row,m)=h1etot+h2etot+(elecs%hnuc*ham%ovrlp(row,m))
            ham%hjk(m,row)=ham%hjk(row,m)
            if(GDflg.eq.'y') then
                if(m.eq.row)then
                    ham%diff_ovrlp_bra(row,m,:) = diff_overlap(zstore(row)%diffalive,zstore(row)%diffdead,&
                                            zstore(m)%diffalive,zstore(m)%diffdead)
                    ham%diff_ovrlp_ket(row,m,:) =ham%diff_ovrlp_bra(row,m,:)
                else 
                    ham%diff_ovrlp_bra(row,m,:) = diff_overlap(zstore(row)%diffalive,zstore(row)%diffdead,&
                                            REAL(zstore(m)%alive),REAL(zstore(m)%dead))
                    ham%diff_ovrlp_ket(row,m,:) = diff_overlap(REAL(zstore(row)%alive),REAL(zstore(row)%dead),&
                                            zstore(m)%diffalive,zstore(m)%diffdead)     
                end if
                ham%diff_ovrlp_bra(m,row,:)=ham%diff_ovrlp_bra(row,m,:)
                ham%diff_ovrlp_ket(m,row,:)=ham%diff_ovrlp_ket(row,m,:)
            
                ham%diff_hjk_bra(row,m,:)=h1etot_diff_bra+h2etot_diff_bra
                ham%diff_hjk_bra(m,row,:)=ham%diff_hjk_bra(row,m,:)
                ham%diff_hjk_ket(row,m,:)=h1etot_diff_ket+h2etot_diff_ket
                ham%diff_hjk_ket(m,row,:)=ham%diff_hjk_ket(row,m,:)

            end if
        end do
        !$omp end do
        !$omp end parallel

        call dealloczs2d(z1jk)
        call dealloczs(z2l)
        call dealloczf(zomt)

        return

    end subroutine he_row


        
    
    ! Function to generate 1 electron hamiltonian element part
    complex(kind=8) function h1et(z1,z2,elecs)

        implicit none

        type(zombiest),intent(in)::z1,z2
        type(elecintrgl),intent(in)::elecs
        type(zombiest)::zomt
        integer::j,k,ierr
        complex(kind=8)::tot
        
        if (errorflag .ne. 0) return
        ierr = 0

        call alloczf(zomt)
        tot=(0.0,0.0)
        !$omp parallel private(j,k,zomt) shared(tot,elecs,z1,z2)
        !$omp do reduction(+:tot) 
        do j=1, norb
            do k=1, norb
                zomt=z2
        
                call an(zomt,j)
                call cr(zomt,k)
                tot = tot + (overlap(z1,zomt)*elecs%h1ei(j,k))
            end do
        end do
        !$omp end do
        !$omp end parallel
        call dealloczf(zomt)
        
        h1et=tot
        ! print*,h1et
        return

    end function h1et

    ! Function to generate 2 electron hamiltonian element part
    complex(kind=8) function h2et(z1,z2,elecs)

        implicit none

        type(zombiest),intent(in)::z1,z2
        type(elecintrgl),intent(in)::elecs
        type(zombiest),allocatable, dimension(:,:)::z1jk
        type(zombiest),allocatable, dimension(:)::z2l
        type(zombiest)::zomt
        integer::j,k,l,jspin, ierr
        complex(kind=8)::tot

        if (errorflag .ne. 0) return
        ierr = 0

      
        call alloczs2d(z1jk,norb)
        call alloczs(z2l,norb)
        call alloczf(zomt)
        tot=(0.0,0.0)

        !$omp parallel shared(z1jk,z1,z2,tot) private(j,k,l,jspin,zomt)
        !$omp do
        do j=1, norb
            do k=1, norb
                zomt=z1
                call an(zomt,j)
                call an(zomt,k)
                z1jk(j,k)=zomt
        
            end do
        end do
        !$omp end do NOWAIT
        !$omp do
        do l=1, norb
            zomt=z2
    
            call an(zomt,l)
            z2l(l)=zomt
            
        end do
        !$omp end do
        
        !$omp do reduction(+:tot)
        do j=1, norb
            if(z1%alive(j)==(0.0,0.0))then
                CYCLE
            end if

            if(modulo(j,2)==0)then
                jspin=2
            else
                jspin=1
            end if

            do k=1, norb
                if(iszero(z1jk(j,k)).eqv..true.)then
                    CYCLE
                end if
                do l=jspin, norb, 2
                    if(z2%alive(l)==(0.0,0.0))then
                        CYCLE
                    end if
                    tot = tot + z_an_z3(z1jk(j,k),z2l(l),elecs%h2ei(j,k,l,:))
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel

        call dealloczs2d(z1jk)
        call dealloczs(z2l)
        call dealloczf(zomt)

        h2et=tot*0.5
        ! print*,h2et
        return

    end function h2et

    
    ! Function to generate indivdual hamiltonian elements
    ! Old version not currently used but useful to keep
    complex(kind=8) function hamval(z1,z2,elecs,ovrl)

        implicit none

        type(zombiest),intent(in)::z1,z2
        type(elecintrgl),intent(in)::elecs
        complex(kind=8), intent(in)::ovrl

        hamval= h1et(z1,z2,elecs)+h2et(z1,z2,elecs)+((elecs%hnuc)*ovrl)
        return
    
    end function hamval



    ! Top level routine to allocate hamiltonian and overlap matrices 
    subroutine hamgen(ham,zstore,elecs,size)

        implicit none

        type(hamiltonian), intent(inout)::ham 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer, allocatable,dimension(:)::IPIV1
        complex(kind=8),allocatable,dimension(:)::WORK1
        ! real(kind=8),dimension(norb)::temp
        real(kind=8),dimension(ndet,ndet,norb)::temp2
        ! real(kind=8),dimension(ndet)::vector
        real(kind=8),dimension(ndet,norb)::matrix
        integer:: j,k,l,size,ierr


        do j=1, size
                call he_row(ham,zstore,elecs,j,size)
            write(6,"(a,i0,a)") "Hamiltonian row ",j, " completed"
        end do

        ham%inv=ham%ovrlp
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

        call ZGETRF(size,size,ham%inv,size,IPIV1,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if
        call ZGETRI(size,ham%inv,size,IPIV1,WORK1,size,ierr)
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

        !$omp parallel
        !$omp workshare
        ham%kinvh=matmul(ham%inv,ham%hjk)
        !$omp end workshare
        !$omp end parallel

        ! To find the derrivative of the inverse of a matrix we use the identity
        ! d(omega^-1)= omega^-1d(omega)omega^-1. Since this result is only used when being multiplied with the hamiltonian 
        ! we will actually calculate omega^-1d(omega)omega^-1*H 
        ! j decides which zombie state we differentiate by
        ! omega^-1d(omega) is first calcualted. creating an ndet*ndet*norb array 
        ! this tensor is then 
        if(GDflg.eq.'y')then
            print*,ham%diff_hjk_bra(1,2,:) 
            do j=1, ndet
                ! do l=1, ndet
                    do k=1, ndet
                        if(j.eq.k)then
                            ! vector=REAL(ham%ovrlp(k,:))
                            ! matrix=ham%diff_ovrlp_ket(:,j,:)
                            temp2(k,j,:)=matmul(REAL(ham%ovrlp(k,:)),ham%diff_ovrlp_ket(:,j,:))   !matmul(vector,matrix)
                        else
                            do l=1, ndet
                                temp2(l,k,:)=REAL(ham%ovrlp(l,j))*ham%diff_ovrlp_bra(j,k,:)
                            end do   
                        end if
                    end do
                do k=1, ndet
                    ! temp(1:norb)=0
                    matrix=temp2(k,:,:)
                    do l=1, ndet
                        ham%diff_invh(j,k,l,:)=matmul(matrix,REAL(ham%kinvh(1:ndet,l)))
                    end do
                end do
            end do
        end if
        return
        
    end subroutine hamgen



END MODULE ham