MODULE ham 

    use mod_types
    use globvars
    use alarrays
    use omp_lib
    use dnad
    
    contains

    !Level 0 Hamiltonian Routine
    ! Subroutine that controls and calcualtes all of the hamiltonian variables 
    subroutine hamgen(haml,zstore,elecs,size,verb)

        implicit none 

        type(hamiltonian), intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::size,verb
        integer, allocatable,dimension(:)::IPIV1
        real(kind=8),allocatable,dimension(:)::WORK1
        integer::ierr

      
        

        if (errorflag .ne. 0) return
        ierr=0
        
   
        !$acc data copyout(haml%hjk(1:size,1:size),haml%ovrlp(1:size,1:size),haml%diff_hjk(1:2*norb,1:size,1:size),&
        !$acc & haml%diff_ovrlp(1:2*norb,1:size,1:size))
        call haml_ovrlp_comb(haml,zstore,elecs,size,verb)
        ! call haml_ovrlp_comb_gpu(haml,zstore,elecs,size,verb)
        !$acc end data
        haml%inv=haml%ovrlp
        allocate(WORK1(size),IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV or WORK1 vector allocation . ierr had value ", ierr
            errorflag=1
        end if 

        Call dgetrf(size, size, haml%inv, size, IPIV1, ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF ",ierr
        end if
        if (ierr==0) call dgetri(size,haml%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRI ",ierr
        end if

        deallocate(WORK1,IPIV1)

        call DGEMM("N","N",size,size,size,1.d0,haml%inv,size,haml%hjk,size,0.d0,haml%kinvh,size)

        return 

    end subroutine hamgen

    !##############################################################################################################################
    
    
    !Level 1 Routines to make the Hamiltonian and Overlap matrices

    subroutine haml_ovrlp_comb(haml,zstore,elecs,size,verb)
        implicit none

        type(hamiltonian), intent(inout)::haml 
        ! real(kind=8),dimension(:,:),intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(dual2),dimension(0:2*norb)::z1d
        integer,intent(in)::verb,size
        integer::j,k,l,n,p,ierr
        type(dual2)::ovlptot,hamtot
        real(wp)::tempx,tmp,temph_2,temph,tmph
        real(wp),dimension(2*norb)::temp_dx,temph_dx2,temph_dx
        if (errorflag .ne. 0) return 
        ierr=0

        !!$omp parallel do &
        !!$omp & private(j,k,ovlptot,hamtot,z1d) &
        !!$omp & shared(elecs,zstore,haml,size) 
        !$acc data create(z1d,ovlptot,hamtot,tempx,temp_dx(1:2*norb),tmp,temph_dx2(1:2*norb),temph_dx(1:2*norb),temph_2,temph,tmph)&
        !$acc & present(elecs,zstore,norb)
        !$acc parallel loop gang private(ovlptot,hamtot,z1d,k,j,l,n,tempx,temp_dx,tmp,tmph,temph_2,temph_dx,temph_dx2)
        do j=1,size
            z1d =typ2_2_typ1(zstore(j)%val)

            ! hamtot=haml_vals_2(z1d,z1d,elecs)+(elecs%hnuc)
            temph_2=0.0d0; temph_dx2=0.0d0
            !$acc loop reduction(+:temph_2) reduction(+:temph_dx2) private(temph,temph_dx,l,j,k)
            do p=1,elecs%num
                temph=elecs%integrals(p)
                !$acc loop reduction(*:temph)
                do n=1, norb
                    temph=z1d(n)%x*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%x+&
                    z1d(n+norb)%x*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%x
                end do
                temph_2=temph_2+temph
                temph_dx=elecs%integrals(p)
                !$acc loop collapse(2) reduction(*:temph_dx)
                do l=1,norb
                    do n=1, norb
                        if(l.ne.n)then
                            tmph=(z1d(n)%x*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%x+&
                            z1d(n+norb)%x*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%x)
                            temph_dx(l)=temp_dx(l)*tmph
                            temph_dx(l+norb)=temp_dx(l+norb)*tmph
                        else
                            tmph= ((z1d(n)%x*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%dx(l)+&
                            z1d(n)%dx(l)*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%x)+&
                            (z1d(n+norb)%x*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%dx(l)+&
                            z1d(n+norb)%dx(l)*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%x))
                            temph_dx(l)=temp_dx(l)*tmph
                            tmph=((z1d(n)%x*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%dx(l+norb)+&
                            z1d(n)%dx(l+norb)*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%x)+&
                            (z1d(n+norb)%x*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%dx(l+norb)+&
                            z1d(n+norb)%dx(l+norb)*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%x))
                            temph_dx(l+norb)=temp_dx(l+norb)*tmph
                        end if 
                    end do 
                end do 
                temph_dx2=temph_dx2+temph_dx
            end do
            hamtot%x=temph_2
            hamtot%dx=temph_dx2

            haml%hjk(j,j)=hamtot%x
            haml%diff_hjk(:,j,j)=hamtot%dx(1:norb)

            haml%ovrlp(j,j)=1.0d0
            haml%diff_ovrlp(:,j,j)=0.0d0 
           
            do k=j+1,size
                tempx=1.0d0
                temp_dx=1.0d0
                !$acc loop worker reduction(*:tempx) 
                do l=1,norb
                    tempx=tempx*((z1d(l)%x*zstore(k)%val(l)%x)+(z1d(l+norb)%x*zstore(k)%val(norb+l)%x))
                end do 
                !$acc loop worker collapse(2) reduction(*:temp_dx)
                do l=1,norb
                    do n=1,norb
                        if(l.ne.n)then 
                            tmp=(z1d(n)%x * zstore(k)%val(n)%x) + (z1d(n+norb)%x * zstore(k)%val(n+norb)%x)
                            temp_dx(l)=temp_dx(l)*tmp
                            temp_dx(l+norb)=temp_dx(l+norb)*tmp
                        else
                            tmp=(z1d(l)%x*zstore(k)%val(l)%dx(l))+(z1d(l+norb)%x*zstore(k)%val(norb+l)%dx(l))

                            temp_dx(l)=temp_dx(l)*tmp
                            tmp=(z1d(l)%dx(l+norb)*zstore(k)%val(l)%x)+(z1d(l+norb)%dx(l+norb)*zstore(k)%val(norb+l)%x)
                            temp_dx(l+norb)=temp_dx(l+norb)*tmp
                        end if
                    end do
                end do
                ovlptot%x=tempx
                ovlptot%dx=temp_dx
                ! ovlptot=overlap_2(z1d,zstore(k)%val)
                haml%ovrlp(j,k)=ovlptot%x; haml%ovrlp(k,j)=haml%ovrlp(j,k)
                haml%diff_ovrlp(:,k,j)=ovlptot%dx(1:norb)
                haml%diff_ovrlp(:,j,k)=ovlptot%dx(1+norb:2*norb)

                temph_2=0.0d0; temph_dx2=0.0d0
                !$acc loop reduction(+:temph_2) reduction(+:temph_dx2) private(temph,temph_dx,l,j,k)
                do p=1,elecs%num
                    temph=elecs%integrals(p)
                    !$acc loop reduction(*:temph)
                    do n=1, norb
                        temph=z1d(n)%x*elecs%neg_a(n,p)*zstore(k)%val(elecs%alive(n,p))%x+&
                        z1d(n+norb)%x*elecs%neg_d(n,p)*zstore(k)%val(elecs%dead(n,p))%x
                    end do
                    temph_2=temph_2+temph
                    temph_dx=elecs%integrals(p)
                    !$acc loop collapse(2) reduction(*:temph_dx)
                    do l=1,norb
                        do n=1, norb
                            if(l.ne.n)then
                                tmph=(z1d(n)%x*elecs%neg_a(n,p)*zstore(k)%val(elecs%alive(n,p))%x+&
                                z1d(n+norb)%x*elecs%neg_d(n,p)*zstore(k)%val(elecs%dead(n,p))%x)
                        
                                temph_dx(l)=temp_dx(l)*tmph
                                
                                temph_dx(l+norb)=temp_dx(l+norb)*tmph
                            else
                                tmph= ((z1d(n)%x*elecs%neg_a(n,p)*zstore(k)%val(elecs%alive(n,p))%dx(l)+&
                                z1d(n)%dx(l)*elecs%neg_a(n,p)*zstore(k)%val(elecs%alive(n,p))%x)+&
                                (z1d(n+norb)%x*elecs%neg_d(n,p)*zstore(k)%val(elecs%dead(n,p))%dx(l)+&
                                z1d(n+norb)%dx(l)*elecs%neg_d(n,p)*zstore(k)%val(elecs%dead(n,p))%x))
                                
                                temph_dx(l)=temp_dx(l)*tmph
                                tmph=((z1d(n)%x*elecs%neg_a(n,p)*zstore(k)%val(elecs%alive(n,p))%dx(l+norb)+&
                                z1d(n)%dx(l+norb)*elecs%neg_a(n,p)*zstore(k)%val(elecs%alive(n,p))%x)+&
                                (z1d(n+norb)%x*elecs%neg_d(n,p)*zstore(k)%val(elecs%dead(n,p))%dx(l+norb)+&
                                z1d(n+norb)%dx(l+norb)*elecs%neg_d(n,p)*zstore(k)%val(elecs%dead(n,p))%x))
                                
                                temph_dx(l+norb)=temp_dx(l+norb)*tmph
                            end if 
                        end do 
                    end do 
                    temph_dx2=temph_dx2+temph_dx(1:2*norb)
                end do
                hamtot%x=temph_2
                hamtot%dx=temph_dx2

                !$acc wait
                ! hamtot=haml_vals_2(z1d,zstore(k)%val,elecs)+(ovlptot*elecs%hnuc)
                haml%hjk(j,k)=hamtot%x; haml%hjk(k,j)=haml%hjk(j,k)
                haml%diff_hjk(:,k,j)=hamtot%dx(1:norb)
                haml%diff_hjk(:,j,k)=hamtot%dx(1+norb:2*norb)
            end do 
            ! if(verb.eq.1)then
                ! write(6,"(a,i0,a)") "hamliltonian column ",j, " completed"
            ! end if 
        end do
        !$acc end data
        !!$omp end parallel do

    end subroutine haml_ovrlp_comb

    

    !##############################################################################################################################

    !Level 2 routines to make an Overlap and Hamiltonian matrix column

    ! Calcualates a column of a hamiltonian Start specifies the row the column
    ! is started to be calcualted 

    subroutine haml_ovrlp_column(temp,zstore,size,elecs,row)

        implicit none
        type(zombiest),dimension(1:ndet),intent(in)::zstore
        type(grad_do),intent(inout)::temp
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::row,size
        type(dual2),dimension(0:2*norb)::z1d
        type(dual2)::ovlptot,hamtot
        real(wp)::tempx,tmp,temph_2,temph,tmph
        real(wp),dimension(2*norb)::temp_dx,temph_dx2,temph_dx
        integer::j,k,l,p,n

        if (errorflag .ne. 0) return
        z1d = typ2_2_typ1(temp%zom%val)
        !$acc data copyin(z1d,size) copyout(temp%hjk(1:size,row),temp%ovrlp(1:size,row),temp%diff_ovrlp_1,temp%diff_ovrlp_2) &
        !$acc & create(ovlptot,hamtot,tempx,temp_dx(1:2*norb),tmp,temph_dx2(1:2*norb),temph_dx(1:2*norb),temph_2,temph,tmph) &
        !$acc & present(elecs,zstore,norb)
        !!$omp parallel do  default(none) &
        !!$omp & private(j,ovlptot,hamtot) &
        !!$omp & shared(elecs,zstore,temp,z1d,row,norb,size)
        !$acc parallel loop async gang private(ovlptot,hamtot,temp_dx,tempx,tmp,tmph,temph_2,temph_dx,temph_dx2)
        do j=1,size
            if (j.ne.row) then
                tempx=1.0d0
                temp_dx=1.0d0
                !$acc loop worker reduction(*:tempx) 
                do l=1,norb
                    tempx=tempx*((z1d(l)%x*zstore(j)%val(l)%x)+(z1d(l+norb)%x*zstore(j)%val(norb+l)%x))
                end do 
                !$acc loop worker collapse(2) reduction(*:temp_dx)
                do l=1,norb
                    do k=1,norb
                        if(l.ne.k)then
                            tmp=((z1d(k)%x * zstore(j)%val(k)%x) + (z1d(k+norb)%x * zstore(j)%val(k+norb)%x))
                            temp_dx(l)=temp_dx(l)*tmp
                            temp_dx(l+norb)=temp_dx(l+norb)*tmp
                        else
                            tmp=(z1d(l)%x*zstore(j)%val(l)%dx(l))+(z1d(l+norb)%x*zstore(j)%val(norb+l)%dx(l))
                            temp_dx(l)=temp_dx(l)*tmp
                            tmp=(z1d(l)%dx(l+norb)*zstore(j)%val(l)%x)+(z1d(l+norb)%dx(l+norb)*zstore(j)%val(norb+l)%x)
                            temp_dx(l+norb)=temp_dx(l+norb)*tmp
                        end if
                    end do
                end do
                ovlptot%x=tempx
                ovlptot%dx=temp_dx
                ! ovlptot=overlap_2(z1d,zstore(j)%val)
                temp%ovrlp(row,j)=ovlptot%x; temp%ovrlp(row,j)%dx=ovlptot%dx(1:norb)
                ! temp%ovrlp(j,row)=temp%ovrlp(row,j)
                temp%diff_ovrlp_1(:,j)=ovlptot%dx(1:norb)
                temp%diff_ovrlp_2(:,j)=ovlptot%dx(1+norb:2*norb)

                temph_2=0.0d0; temph_dx2=0.0d0
                !$acc loop reduction(+:temph_2) reduction(+:temph_dx2) private(temph,temph_dx,l,j,k)
                do p=1,elecs%num
                    temph=elecs%integrals(p)
                    !$acc loop reduction(*:temph)
                    do n=1, norb
                        temph=z1d(n)%x*elecs%neg_a(n,p)*zstore(j)%val(elecs%alive(n,p))%x+&
                        z1d(n+norb)%x*elecs%neg_d(n,p)*zstore(j)%val(elecs%dead(n,p))%x
                    end do
                    temph_2=temph_2+temph
                    temph_dx=elecs%integrals(p)
                    !$acc loop collapse(2) reduction(*:temph_dx)
                    do l=1,norb
                        do n=1, norb
                            if(l.ne.n)then
                                tmph=(z1d(n)%x*elecs%neg_a(n,p)*zstore(j)%val(elecs%alive(n,p))%x+&
                                z1d(n+norb)%x*elecs%neg_d(n,p)*zstore(j)%val(elecs%dead(n,p))%x)
                                temph_dx(l)=temp_dx(l)*tmph
                                temph_dx(l+norb)=temp_dx(l+norb)*tmph
                            else
                                tmph= ((z1d(n)%x*elecs%neg_a(n,p)*zstore(j)%val(elecs%alive(n,p))%dx(l)+&
                                z1d(n)%dx(l)*elecs%neg_a(n,p)*zstore(j)%val(elecs%alive(n,p))%x)+&
                                (z1d(n+norb)%x*elecs%neg_d(n,p)*zstore(j)%val(elecs%dead(n,p))%dx(l)+&
                                z1d(n+norb)%dx(l)*elecs%neg_d(n,p)*zstore(j)%val(elecs%dead(n,p))%x))
                                temph_dx(l)=temp_dx(l)*tmph
                                tmph=((z1d(n)%x*elecs%neg_a(n,p)*zstore(j)%val(elecs%alive(n,p))%dx(l+norb)+&
                                z1d(n)%dx(l+norb)*elecs%neg_a(n,p)*zstore(j)%val(elecs%alive(n,p))%x)+&
                                (z1d(n+norb)%x*elecs%neg_d(n,p)*zstore(j)%val(elecs%dead(n,p))%dx(l+norb)+&
                                z1d(n+norb)%dx(l+norb)*elecs%neg_d(n,p)*zstore(j)%val(elecs%dead(n,p))%x))
                                
                                temph_dx(l+norb)=temp_dx(l+norb)*tmph
                            end if 
                        end do 
                    end do 
                    temph_dx2=temph_dx2+temph_dx
                end do
                hamtot%x=temph_2
                hamtot%dx=temph_dx2

                !$acc wait
                ! hamtot=haml_vals_2(z1d,zstore(j)%val,elecs)+(ovlptot*elecs%hnuc)
                temp%hjk(row,j)%x=hamtot%x; temp%hjk(row,j)%dx=hamtot%dx(1:norb)
                ! temp%hjk(j,row)=temp%hjk(row,j)

                temp%diff_hjk_1(:,j)=hamtot%dx(1:norb)
                temp%diff_hjk_2(:,j)=hamtot%dx(1+norb:2*norb)
            else 
                temp%ovrlp(row,row)=1.0d0
                temp%ovrlp(row,row)%dx=0.0d0
                temp%diff_ovrlp_1(:,j)=0.0d0

                temph_2=0.0d0; temph_dx2=0.0d0
                !$acc loop reduction(+:temph_2) reduction(+:temph_dx2) private(temph,temph_dx,l,j,k)
                do p=1,elecs%num
                    temph=elecs%integrals(p)
                    !$acc loop reduction(*:temph)
                    do n=1, norb
                        temph=z1d(n)%x*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%x+&
                        z1d(n+norb)%x*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%x
                    end do
                    temph_2=temph_2+temph
                    temph_dx=elecs%integrals(p)
                    !$acc loop collapse(2) reduction(*:temph_dx)
                    do l=1,norb
                        do n=1, norb
                            if(l.ne.n)then
                                tmph=(z1d(n)%x*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%x+&
                                z1d(n+norb)%x*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%x)
                                temph_dx(l)=temp_dx(l)*tmph
                                temph_dx(l+norb)=temp_dx(l+norb)*tmph
                            else
                                tmph= ((z1d(n)%x*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%dx(l)+&
                                z1d(n)%dx(l)*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%x)+&
                                (z1d(n+norb)%x*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%dx(l)+&
                                z1d(n+norb)%dx(l)*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%x))

                                temph_dx(l)=temp_dx(l)*tmph
                                tmph=((z1d(n)%x*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%dx(l+norb)+&
                                z1d(n)%dx(l+norb)*elecs%neg_a(n,p)*z1d(elecs%alive(n,p))%x)+&
                                (z1d(n+norb)%x*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%dx(l+norb)+&
                                z1d(n+norb)%dx(l+norb)*elecs%neg_d(n,p)*z1d(elecs%dead(n,p))%x))
                                
                                temph_dx(l+norb)=temp_dx(l+norb)*tmph
                            end if 
                        end do 
                    end do 
                    temph_dx2=temph_dx2+temph_dx
                end do
                hamtot%x=temph_2
                hamtot%dx=temph_dx2
    
                ! hamtot=haml_vals_2(z1d,z1d,elecs)+(elecs%hnuc)
                temp%hjk(row,row)%x=hamtot%x; temp%hjk(row,row)%dx=hamtot%dx(1:norb)

                temp%diff_hjk_1(:,j)=hamtot%dx(1:norb)
            end if 
        end do
        !$acc wait
        !!$omp end parallel do 
        !$acc end data
        temp%ovrlp(:,row)=temp%ovrlp(row,:); temp%hjk(:,row)=temp%hjk(row,:)
        return

    end subroutine haml_ovrlp_column

    !##############################################################################################################################

    !Level 3 routine to make individual value in overlap and hamiltonian matrix

    ! calculates indvidual hamiltonian elements taking in two Zombie states and a set of 
    ! creation and annihilation operations
    function haml_vals(z1d,z2d,elecs) result(ham_tot)

        implicit none 
        type(dual2),dimension(0:),intent(in)::z1d,z2d
        type(elecintrgl),intent(in)::elecs
        type(dual2)::ham_tot
        type(dual2)::ov
        integer::j,k
        
        if (errorflag .ne. 0) return

       
        ham_tot=0.0d0
        !!$omp simd
        do j=1,elecs%num
            ov=elecs%integrals(j)
            do k=1, norb
                ov=ov*((z1d(k)*z2d(elecs%alive(k,j))*elecs%neg_a(k,j))+(z1d((k+norb))*z2d((elecs%dead(k,j)))*elecs%neg_d(k,j)))
            end do
            ham_tot=ham_tot+ov
        end do
        !!$omp end simd
        return 
      
    end function haml_vals

    ! calculates individual overlaps where no creation and annihilation operations are needed
    function overlap_1(z1d,z2d) result(ovrlp_tot)

        implicit none
        type(dual2),dimension(0:)::z1d,z2d
        type(dual2)::ovrlp_tot
        integer::j

        if (errorflag .ne. 0) return

        ovrlp_tot=1.0d0
      
        !!$omp simd
        do j=1,norb
            ovrlp_tot=ovrlp_tot*((z1d(j)*z2d(j))+(z1d(j+norb)*z2d(norb+j)))
        end do
        !!$omp end simd
      
        return 
    end function overlap_1

    !##############################################################################################################################
    
    
   
    function overlap_2(z1d,z2d) result(ovrlp_tot)
       
        implicit none
        type(dual2),dimension(0:)::z1d,z2d
        type(dual2)::ovrlp_tot
        real(wp)::temp
        real(wp),dimension(2*norb)::temp_dx
        integer::j,k

        if (errorflag .ne. 0) return

        ! ovrlp_tot=1.0d0
        !!$omp simd
        !$acc data present(z1d,z2d,norb,ovrlp_tot) create(temp,temp_dx(1:2*norb)) copyin(norb)
        temp=1.0d0
        temp_dx=1.0d0
        !$acc loop worker reduction(*:temp) 
        do j=1,norb
            temp=temp*((z1d(j)%x*z2d(j)%x)+(z1d(j+norb)%x*z2d(norb+j)%x))
        end do 
        !$acc loop worker collapse(2) reduction(*:temp_dx)
        do j=1,norb
            do k=1,norb
                if(j.ne.k)then 
                    temp_dx(j)=temp_dx(j)*((z1d(k)%x * z2d(k)%x) + (z1d(k+norb)%x * z2d(k+norb)%x))
                    temp_dx(j+norb)=temp_dx(j+norb)*((z1d(k)%x * z2d(k)%x) + (z1d(k+norb)%x * z2d(k+norb)%x))
                else
                    temp_dx(j)=temp_dx(j)*(z1d(j)%x*z2d(j)%dx(j))+(z1d(j+norb)%x*z2d(norb+j)%dx(j))
                    temp_dx(j+norb)=temp_dx(j+norb)*(z1d(j)%dx(j+norb)*z2d(j)%x)+(z1d(j+norb)%dx(j+norb)*z2d(norb+j)%x)
                end if
            end do
        end do
        ovrlp_tot%x=temp
        ovrlp_tot%dx=temp_dx
        !$acc end data
       


        ! do j=1,norb
            ! temp=z1d(j)%x*z2d(j)%x+z1d(j+norb)%x*z2d(norb+j)%x
            
            ! do k=1,2*norb
                ! ovrlp_tot%dx(k)=((z1d(j)%x*z2d(j)%dx(k)+z1d(j)%dx(k)*z2d(j)%x)+&
                                ! (z1d(j+norb)%x*z2d(norb+j)%dx(k)+z1d(j+norb)%dx(k)*z2d(norb+j)%x))*ovrlp_tot%x + &
                                ! (z1d(j)%x*z2d(j)%x+z1d(j+norb)%x*z2d(norb+j)%x)*ovrlp_tot%dx(k)
                                ! (temp)*ovrlp_tot%dx(k)
                ! ovrlp_tot%dx(k+norb)=((z1d(j)%x*z2d(j)%dx(k+norb)+z1d(j)%dx(k+norb)*z2d(j)%x)+&
                !                 (z1d(j+norb)%x*z2d(norb+j)%dx(k+norb)+z1d(j+norb)%dx(k+norb)*z2d(norb+j)%x))*ovrlp_tot%x + &
                !                 (temp)*ovrlp_tot%dx(k+norb)
            ! end do
            ! ovrlp_tot%x=ovrlp_tot%x*(z1d(j)%x*z2d(j)%x+z1d(j+norb)%x*z2d(norb+j)%x) !(temp)
        ! end do
       
        !!$omp end simd
      
        return 
    end function overlap_2

    function haml_vals_2(z1d,z2d,elecs) result(ham_tot)
        
        implicit none 
        type(dual2),dimension(0:),intent(in)::z1d,z2d
        type(elecintrgl),intent(in)::elecs
        type(dual2)::ham_tot
        type(dual2)::ov
        integer::j,k,l
        real(wp)::temp,temp_2
        real(wp),dimension(2*norb)::temp_dx,temp_dx2
        if (errorflag .ne. 0) return

       
        
      
        !$acc data present(z1d,z2d,elecs,norb,ham_tot) create(norb,temp,temp_2,temp_dx(1:2*norb),temp_dx2(1:2*norb))
        temp=0.0d0; temp_2=0.0d0
        !$acc loop reduction(+:temp_2) reduction(+:temp_dx2) private(temp,temp_dx,l,j,k)
        do j=1,elecs%num
            temp=elecs%integrals(j)
            !$acc loop reduction(*:temp)
            do k=1, norb
                temp=z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x+z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x
            end do
            temp_2=temp_2+temp
            temp_dx=elecs%integrals(j)
            !$acc loop collapse(2) reduction(*:temp_dx)
            do l=1,norb
                do k=1, norb
                    if(l.ne.k)then
                        temp_dx(l)=temp_dx(l)*(z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x+&
                        z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x)

                        temp_dx(l+norb)=temp_dx(l+norb)*(z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x+&
                        z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x)  
                    else
                        temp_dx(l)=temp_dx(l)*((z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%dx(l)+&
                            z1d(k)%dx(l)*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x)+&
                            (z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%dx(l)+&
                            z1d(k+norb)%dx(l)*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x))

                        temp_dx(l+norb)=temp_dx(l+norb)*((z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%dx(l+norb)+&
                            z1d(k)%dx(l+norb)*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x)+&
                            (z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%dx(l+norb)+&
                            z1d(k+norb)%dx(l+norb)*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x))
                    end if 
                end do 
            end do 
            temp_dx2=temp_dx2+temp_dx
        end do
        ham_tot%x=temp_2
        ham_tot%dx=temp_dx2
        !$acc end data  
               

        ! do j=1,elecs%num
            ! ov=elecs%integrals(j)
           
            ! do k=1, norb
                ! temp=z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x+z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x
                
                ! do l=1,2*norb !k
                !     ov%dx(l)=((z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%dx(l)+&
                !     z1d(k)%dx(l)*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x)+&
                !     (z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%dx(l)+&
                !     z1d(k+norb)%dx(l)*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x))*ov%x + &
                !     (z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x+z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x)*ov%dx(l)
                    ! (temp)*ov%dx(l)

                    ! ov%dx(l+1)=((z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%dx(l+1)+&
                    ! z1d(k)%dx(l+1)*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x)+&
                    ! (z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%dx(l+1)+&
                    ! z1d(k+norb)%dx(l+1)*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x))*ov%x + &
                !     ! (temp)*ov%dx(l+1)
                ! end do
                ! ov%x=ov%x*z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x+z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x
                !temp
            ! end do
            ! ham_tot=ham_tot+ov
        ! end do
        
        !!$omp end simd
        return 
      
    end function haml_vals_2



    ! function overlap_gpu(z1d,z2d,zstore,elecs) result(overlap)

    !     implicit none 
    !     type(dual2)::z1d(:),z2d(:)
    !     type(dual2)::overlap
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     type(elecintrgl),intent(in)::elecs
    !     real(wp)::ovlptot_var
    !     real(wp),dimension(norb)::ovlptot_dx
    !     integer::j,k,l
    !     ! type(dim3)::threads_per_block, ovrlp_dx_block_num

    !     if (errorflag .ne. 0) return

    !     ! threads_per_block= dim3(16,16)
    !     ! ovrlp_dx_block_num=ceiling(norb/threads_per_block%x,norb/threads_per_block%y)
        
      
    !     ovlptot_var=1.0d0
    !     !$cuf kernel do <<<*,*>>> reduce(*:ovlptot_var)
    !     do j=1,norb
    !         ovlptot_var=ovlptot_var*((z1d(j)%x * z2d(j)%x) + (z1d(j+norb)%x * z2d(j+norb)%x))
    !     end do
    !     ! call overlp_x<<<1, norb>>>(z1d, z2d, norb, ovlptot_var)
    !     ovlptot_dx=1.0d0
    !     !$cuf kernel do <<<*,*>>> reduce(*:ovlptot_dx)
    !     do j=1,norb
    !         do l=1,norb
    !             if(j.ne.l)then 
    !                 ovlptot_dx(j)=ovlptot_dx(j)*((z1d(l)%x * z2d(l)%x) + (z1d(l+norb)%x * z2d(l+norb)%x))
    !             else
    !                 ovlptot_dx(j)=ovlptot_dx(j)*(z1d(j)%x*z2d(j)%dx(j)+z1d(j)%dx(j)*z2d(j)%x)+&
    !                         (z1d(j+norb)%x*z2d(norb+j)%dx(j)+z1d(j+norb)%dx(j)*z2d(norb+j)%x)
    !             end if
    !         end do 
    !     end do
    !     ! call overlp_dx<<<ovrlp_dx_block_num, threads_per_block>>>(z1d, z2d, norb, ovlptot_dx)
        
    !     overlap%x=ovlptot_var
    !     overlap%dx(1:norb)=ovlptot_dx(1:norb)   
       
    ! end function overlap_gpu

    ! function haml_gpu(z1d,z2d,elecs) result(ham_tot)

    !     type(elecintrgl),intent(in)::elecs
    !     type(dual2)::haml_tot
    !     real(wp),dimension(elecs%num)::temp_result_x
    !     real(wp),dimension(elecs%num,norb)::temp_result_dx
    !     type(dim3)::threads_per_block_x,threads_per_block_dx,haml_grid,haml_grid_dx
    !     integer::j,k

    !     if (errorflag .ne. 0) return
    !     ham_tot=0.0d0
      

    !     ! threads_per_block_x= dim3(32,32)
    !     ! haml_grid=dim3(ceiling(elecs%num/threads_per_block_x%x),ceiling(norb/threads_per_block_x%y))
    !     ! threads_per_block_dx= dim3(4,16,16)
    !     ! haml_grid_dx=dim3(ceiling(elecs%num/threads_per_block_dx%x),&
    !     !     ceiling(norb/threads_per_block_dx%y),ceiling(norb/threads_per_block_dx%z))

    !     temp_result_x=elecs%integrals
        
    !     !$cuf kernel do(2) <<<*,*>>> reduce(*:temp_result_x)
    !     do j=1,elecs%num
    !         do l=1,norb
    !             temp_result_x(j)=temp_result_x(j)*(z1d(l)%x*elecs%neg_a(l,j)*z2d(elecs%alive(l,j))%x+&
    !             z1d(l+norb)%x*elecs%neg_d(l,j)*z2d(elecs%dead(l,j))%x)
    !         end do 
    !     end do

    !     ! call haml_x<<<haml_grid, threads_per_block_x>>>(z1d, z2d, norb, temp_result_x)

    !     ham_tot%x=0.0d0
    !     !$cuf kernel do <<<*,*>>> reduce(+:ham_tot)
    !     do j=1,elecs%num
    !         ham_tot%x=ham_tot%x+temp_result_x(j)
    !     end do
      
    !     !$cuf kernel do(2) <<<*,*>>>
    !     do k=1,norb
    !         do j=1,elecs%num
    !             temp_result_dx(j,k)=elecs%integrals(j)
    !         end do 
    !     end do

    !     !$cuf kernel do(3) <<<*,*>>> reduce(*:temp_result_dx)
    !     do j=1,elecs%num
    !         do l=1,norb
    !             do k=1, norb
    !                 if(l.ne.k)then
    !                     temp_result_dx(j,l)=temp_result_dx(j,l)*(z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x+&
    !                     z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x)  
    !                 else
    !                     temp_result_dx(j,l)=temp_result_dx(j,l)*((z1d(k)%x*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%dx(l)+&
    !                         z1d(k)%dx(l)*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))%x)+&
    !                         (z1d(k+norb)%x*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%dx(l)+&
    !                         z1d(k+norb)%dx(l)*elecs%neg_d(k,j)*z2d(elecs%dead(k,j))%x))
    !                 end if 
    !             end do 
    !         end do
    !     end do 

    !     ! call haml_dx<<<haml_grid_dx, threads_per_block_dx>>>(z1d, z2d, elecs, norb, temp_result_dx)
    !     ham_tot%dx=0.d0
    !     !$cuf kernel do(2) <<<*,*>>> reduce(+:ham_tot)
    !     do k=1,norb
    !         do j=1,elecs%num
    !             ham_tot%dx(k)=ham_tot%dx(k)+temp_result_dx(j,k)
    !         end do 
    !     end do

        

    ! end function haml_gpu


    ! subroutine haml_ovrlp_comb_gpu(haml,zstore,elecs,size,verb)
    !     implicit none

    !     type(hamiltonian), intent(inout)::haml 
    !     ! real(kind=8),dimension(:,:),intent(inout)::haml 
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     type(elecintrgl),intent(in)::elecs
    !     type(dual2),dimension(0:2*norb)::z1d
    !     integer,intent(in)::verb,size
    !     integer::j,k,ierr
    !     type(dual2)::ovlptot,hamtot
    
    !     if (errorflag .ne. 0) return 
    !     ierr=0

    !     !$omp parallel do simd&
    !     !$omp & private(j,k,ovlptot,hamtot,z1d) &
    !     !$omp & shared(elecs,zstore,haml) 
    !     do j=1,size
    !         z1d =typ2_2_typ1(zstore(j)%val)

    !         hamtot=haml_gpu(z1d,z1d,elecs)+(elecs%hnuc)
    !         haml%hjk(j,j)=hamtot%x
    !         haml%diff_hjk(:,j,j)=hamtot%dx(1:norb)

    !         haml%ovrlp(j,j)=1.0d0
    !         haml%diff_ovrlp(:,j,j)=0.0d0 
           
    !         do k=j+1,size
    !             ovlptot=1.0d0
    !             ovlptot=overlap_gpu(z1d,store(k)%val,zstore,elecs) 
    !             haml%ovrlp(j,k)=ovlptot%x; haml%ovrlp(k,j)=haml%ovrlp(j,k)
    !             haml%diff_ovrlp(:,k,j)=ovlptot%dx(1:norb)
    !             haml%diff_ovrlp(:,j,k)=ovlptot%dx(1+norb:2*norb)

    !             hamtot=haml_gpu(z1d,store(k)%val,elecs)+(ovlptot*elecs%hnuc)
    !             haml%hjk(j,k)=hamtot%x; haml%hjk(k,j)=haml%hjk(j,k)
    !             haml%diff_hjk(:,k,j)=hamtot%dx(1:norb)
    !             haml%diff_hjk(:,j,k)=hamtot%dx(1+norb:2*norb)
    !         end do 
    !         if(verb.eq.1)then
    !             write(6,"(a,i0,a)") "hamliltonian column ",j, " completed"
    !         end if 
    !     end do
    !     !$omp end parallel do simd

    ! end subroutine haml_ovrlp_comb_gpu



    ! subroutine haml_ovrlp_column_gpu(temp,zstore,size,elecs,row)

    !     implicit none
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     type(grad_do),intent(inout)::temp
    !     type(elecintrgl),intent(in)::elecs
    !     integer,intent(in)::row,size
    !     type(dual2),dimension(0:2*norb)::z1d
    !     type(dual2)::ovlptot,hamtot
    !     integer::j

    !     if (errorflag .ne. 0) return
    !     z1d = typ2_2_typ1(temp%zom%val)
    !     !$omp parallel do  default(none) &
    !     !$omp & private(j,ovlptot,hamtot) &
    !     !$omp & shared(elecs,zstore,temp,z1d,row,norb) 
    !     do j=1,size
    !         if (j.ne.row) then
    !             ovlptot=overlap_gpu(z1d,store(j)%val,zstore,elecs)
    !             temp%ovrlp(row,j)=ovlptot%x; temp%ovrlp(row,j)%dx=ovlptot%dx(1:norb)
    !             temp%ovrlp(j,row)=temp%ovrlp(row,j)

    !             temp%diff_ovrlp_1(:,j)=ovlptot%dx(1:norb)
    !             temp%diff_ovrlp_2(:,j)=ovlptot%dx(1+norb:2*norb)

    !             hamtot=haml_gpu(z1d,store(j)%val,elecs)+(ovlptot*elecs%hnuc)
    !             temp%hjk(row,j)%x=hamtot%x; temp%hjk(row,j)%dx=hamtot%dx(1:norb)
    !             temp%hjk(j,row)=temp%hjk(row,j)

    !             temp%diff_hjk_1(:,j)=hamtot%dx(1:norb)
    !             temp%diff_hjk_2(:,j)=hamtot%dx(1+norb:2*norb)
    !         else 
    !             temp%ovrlp(row,row)=1.0d0
    !             temp%ovrlp(row,row)%dx=0.0d0
    !             temp%diff_ovrlp_1(:,j)=0.0d0
    !             hamtot=haml_gpu(z1d,z1d,elecs)+(elecs%hnuc)
    !             temp%hjk(row,row)%x=hamtot%x; temp%hjk(row,row)%dx=hamtot%dx(1:norb)

    !             temp%diff_hjk_1(:,j)=hamtot%dx(1:norb)
    !         end if 
    !     end do
    !     ! $omp end parallel do 

    !     return

    ! end subroutine haml_ovrlp_column_gpu

  


END MODULE ham
