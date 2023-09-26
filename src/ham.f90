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
       
        call haml_ovrlp_comb(haml,zstore,elecs,size,verb)
        
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
    
   
    !##############################################################################################################################

    
    !Level 1 Routines to make the Hamiltonian and Overlap matrices

    subroutine haml_ovrlp_comb(haml,zstore,elecs,size,verb)
        implicit none

        type(hamiltonian), intent(inout)::haml 
        ! real(kind=8),dimension(:,:),intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::verb,size
        integer::j,k,ierr
        real(wp)::ovlptot,hamtot
    
        if (errorflag .ne. 0) return 
        ierr=0

        !$omp parallel do &
        !$omp & private(j,k,ovlptot,hamtot) &
        !$omp & shared(elecs,zstore,haml) 
        do j=1,size
           
            hamtot=haml_vals_2(zstore(j)%val%x,zstore(j)%val%x,elecs)+(elecs%hnuc)
            haml%hjk(j,j)=hamtot
           
            haml%ovrlp(j,j)=1.0d0
            haml%diff_ovrlp(:,j,j)=0.0d0 
           
            do k=j+1,size
                ovlptot=overlap_2(zstore(j)%val%x,zstore(k)%val%x)
                haml%ovrlp(j,k)=ovlptot; haml%ovrlp(k,j)=haml%ovrlp(j,k)
                hamtot=haml_vals_2(zstore(j)%val%x,zstore(k)%val%x,elecs)+(ovlptot*elecs%hnuc)
                haml%hjk(j,k)=hamtot; haml%hjk(k,j)=haml%hjk(j,k)
              
            end do 
            if(verb.eq.1)then
                write(6,"(a,i0,a)") "hamliltonian column ",j, " completed"
            end if 
        end do
        !$omp end parallel do

    end subroutine haml_ovrlp_comb

    

    !##############################################################################################################################

    !Level 2 routines to make an Overlap and Hamiltonian matrix column

    ! Calcualates a column of a hamiltonian Start specifies the row the column
    ! is started to be calcualted 

    subroutine haml_ovrlp_column(temp,zstore,size,elecs,row)

        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        type(grad_do),intent(inout)::temp
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::row,size
       
        real(wp)::ovlptot,hamtot
        integer::j

        if (errorflag .ne. 0) return
      
        !$omp parallel do  default(none) &
        !$omp & private(j,ovlptot,hamtot) &
        !$omp & shared(elecs,zstore,temp,row,norb) 
        do j=1,size
            if (j.ne.row) then
                ovlptot=overlap_2(temp%zom%val%x,zstore(j)%val%x)
                temp%ovrlp(row,j)=ovlptot
                temp%ovrlp(j,row)=temp%ovrlp(row,j)

                hamtot=haml_vals_2(temp%zom%val%x,zstore(j)%val%x,elecs)+(ovlptot*elecs%hnuc)
                temp%hjk(row,j)=hamtot
                temp%hjk(j,row)=temp%hjk(row,j)

            else 
                temp%ovrlp(row,row)=1.0d0
                hamtot=haml_vals_2(temp%zom%val%x,temp%zom%val%x,elecs)+(elecs%hnuc)
                temp%hjk(row,row)=hamtot
            end if 
        end do
        !$omp end parallel do 

        return

    end subroutine haml_ovrlp_column

    subroutine haml_ovrlp_column_grad(haml,zstore,size,elecs,row,orb)

        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        type(hamiltonian),intent(inout)::haml
        type(elecintrgl),intent(in)::elecs
        integer,intent(in)::row,size,orb
        real(wp),dimension(norb)::ovlptot,hamtot
        integer::j,st,fn

        if (errorflag .ne. 0) return
        if(orb==0)then
            st=1; fn=norb
            !$omp parallel do  default(none) &
            !$omp & private(j,ovlptot,hamtot) &
            !$omp & shared(elecs,zstore,row,norb,orb,st,fn,haml) 
            do j=1,size
                if (j.ne.row) then
                    ovlptot=overlap_2_grad(zstore(row)%val,zstore(j)%val%x,st,fn)
                
                    haml%diff_ovrlp(:,j,row)=ovlptot
                
                    hamtot=haml_vals_2_grad(zstore(row)%val,zstore(j)%val,elecs,st,fn,1)+(ovlptot*elecs%hnuc)

                    haml%diff_hjk(:,j,row)=hamtot
                
                else 
                    haml%diff_ovrlp(:,j,row)=0.0d0
                
                    hamtot=haml_vals_2_grad(zstore(row)%val,zstore(j)%val,elecs,st,fn,0)
                    haml%diff_hjk(:,j,row)=hamtot
                
                end if 
            end do
            !$omp end parallel do 
        else 
            st=orb; fn=orb
            !$omp parallel do  default(none) &
            !$omp & private(j,ovlptot,hamtot) &
            !$omp & shared(elecs,zstore,row,norb,orb,st,fn,haml) 
            do j=1,size
                if (j.ne.row) then
                    ovlptot=overlap_2_grad(zstore(row)%val,zstore(j)%val%x,st,fn)
                
                    haml%diff_ovrlp(orb,j,row)=ovlptot(orb)
                
                    hamtot=haml_vals_2_grad(zstore(row)%val,zstore(j)%val,elecs,st,fn,1)+(ovlptot*elecs%hnuc)

                    haml%diff_hjk(orb,j,row)=hamtot(orb)
                
                else 
                    haml%diff_ovrlp(orb,j,row)=0.0d0
                
                    hamtot=haml_vals_2_grad(zstore(row)%val,zstore(row)%val,elecs,st,fn,0)
                    haml%diff_hjk(orb,j,row)=hamtot(orb)
                
                end if 
            end do
            !$omp end parallel do 
        end if
      
        

        return

    end subroutine haml_ovrlp_column_grad


    !##############################################################################################################################

    !Level 3 routine to make individual value in overlap and hamiltonian matrix

   

    !##############################################################################################################################
    
    
   
    function overlap_2(z1d,z2d) result(ovrlp_tot)

        implicit none
        real(wp),dimension(0:)::z1d,z2d
        real(wp)::ovrlp_tot
        integer::j

        if (errorflag .ne. 0) return

       
        ovrlp_tot=1.0d0
       !$omp simd
        do j=1,norb
            ovrlp_tot=ovrlp_tot*(z1d(j)*z2d(j)+z1d(j+norb)*z2d(norb+j))
        end do
        !$omp end simd
      
        return 
    end function overlap_2

    function haml_vals_2(z1d,z2d,elecs) result(ham_tot)

        implicit none 
        real(wp),dimension(0:),intent(in)::z1d,z2d
        type(elecintrgl),intent(in)::elecs
        real(wp)::ham_tot
        real(wp)::ov
        real(wp)::temp
        integer::j,k
        
        if (errorflag .ne. 0) return

       
        ham_tot=0.0d0
        !$omp simd
        do j=1,elecs%num
            ov=elecs%integrals(j)
            do k=1, norb
                ov=ov*(z1d(k)*elecs%neg_a(k,j)*z2d(elecs%alive(k,j))+z1d(k+norb)*elecs%neg_d(k,j)*z2d(elecs%dead(k,j)))
            end do
            ham_tot=ham_tot+ov
        end do
        !$omp end simd
        return 
      
    end function haml_vals_2

    function overlap_2_grad(z1d,z2d,st,fn) result(ovrlp_tot)

        implicit none
        type(dual),dimension(0:)::z1d
        real(wp),dimension(0:)::z2d
        integer::st,fn
        real(wp),dimension(norb)::ovrlp_tot
        real(wp)::temp
        integer::j,k

        if (errorflag .ne. 0) return

       
        ovrlp_tot=1.0d0
       !$omp simd
        do k=st,fn
            do j=1,norb
                if(j.ne.k)then
                    ovrlp_tot(k)=ovrlp_tot(k)*(z1d(j)%x*z2d(j)+z1d(j+norb)%x*z2d(norb+j))
                else 
                    ovrlp_tot(k)=ovrlp_tot(k)*(z1d(j)%dx(k)*z2d(j)+z1d(j+norb)%dx(k)*z2d(norb+j))
                end if
            end do
        end do
        !$omp end simd
      
        return 
    end function overlap_2_grad

    function haml_vals_2_grad(z1d,z2d,elecs,st,fn,typ) result(ham_tot)

        implicit none 
        type(dual),dimension(0:)::z1d,z2d
        integer::st,fn,typ
        type(elecintrgl),intent(in)::elecs
        real(wp),dimension(norb)::ham_tot
        real(wp)::ov
        integer::j,k,l
        
        if (errorflag .ne. 0) return

       
        ham_tot=0.0d0
        if(typ==1)then
            do k=st, fn
                do j=1,elecs%num
                    ov=elecs%integrals(j)
                    do l=1,norb
                        if(l.ne.k)then 
                            ov=ov*(z1d(l)%x*elecs%neg_a(l,j)*z2d(elecs%alive(l,j))%x+&
                            z1d(l+norb)%x*elecs%neg_d(l,j)*z2d(elecs%dead(l,j))%x)
                        else
                            ov=ov*(z1d(l)%dx(k)*elecs%neg_a(l,j)*z2d(elecs%alive(l,j))%x+&
                            z1d(l+norb)%dx(l)*elecs%neg_d(l,j)*z2d(elecs%dead(l,j))%x)
                        end if
                    end do
                    ham_tot(k)=ham_tot(k)+ov
                end do
            end do
        else
            do k=st, fn
                do j=1,elecs%num
                    ov=elecs%integrals(j)
                    do l=1,norb
                        if(l.ne.k)then 
                            ov=ov*(z1d(l)%x*elecs%neg_a(l,j)*z2d(elecs%alive(l,j))%x+&
                            z1d(l+norb)%x*elecs%neg_d(l,j)*z2d(elecs%dead(l,j))%x)
                        else
                            ov=ov*(z1d(l)%dx(k)*elecs%neg_a(l,j)*z2d(elecs%alive(l,j))%x+&
                            z1d(l)%x*elecs%neg_a(l,j)*z2d(elecs%alive(l,j))%dx(k)+&
                            z1d(l+norb)%dx(l)*elecs%neg_d(l,j)*z2d(elecs%dead(l,j))%x+&
                            z1d(l+norb)%x*elecs%neg_d(l,j)*z2d(elecs%dead(l,j))%dx(k))
                        end if
                    end do
                    ham_tot(k)=ham_tot(k)+ov
                end do
            end do
        end if
      
        return 
      
    end function haml_vals_2_grad

     ! calculates indvidual hamiltonian elements taking in two Zombie states and a set of 
    ! creation and annihilation operations
    ! function haml_vals(z1d,z2d,elecs) result(ham_tot)

    !     implicit none 
    !     type(dual),dimension(0:),intent(in)::z1d,z2d
    !     type(elecintrgl),intent(in)::elecs
    !     type(dual2)::ham_tot
    !     type(dual2)::ov
    !     integer::j,k
        
    !     if (errorflag .ne. 0) return

       
    !     ham_tot=0.0d0
    !     !$omp simd
    !     do j=1,elecs%num
    !         ov=elecs%integrals(j)
    !         do k=1, norb
    !             ov=ov*((z1d(k)*z2d(elecs%alive(k,j))*elecs%neg_a(k,j))+(z1d((k+norb))*z2d((elecs%dead(k,j)))*elecs%neg_d(k,j)))
    !         end do
    !         ham_tot=ham_tot+ov
    !     end do
    !     !$omp end simd
    !     return 
      
    ! end function haml_vals

    ! calculates individual overlaps where no creation and annihilation operations are needed
    ! function overlap_1(z1d,z2d) result(ovrlp_tot)

    !     implicit none
    !     type(dual2),dimension(0:)::z1d,z2d
    !     type(dual2)::ovrlp_tot
    !     integer::j

    !     if (errorflag .ne. 0) return

    !     ovrlp_tot=1.0d0
      
    !     !$omp simd
    !     do j=1,norb
    !         ovrlp_tot=ovrlp_tot*((z1d(j)*z2d(j))+(z1d(j+norb)*z2d(norb+j)))
    !     end do
    !     !$omp end simd
      
    !     return 
    ! end function overlap_1

    ! function overlap_gpu(z1d,z2d,zstore,elecs) result(overlap)

    !     implicit none 
    !     type(dual2),dimension(0:),intent(in)::z1d,z2d
    !     type(dual2)::overlap
    !     type(zombiest),dimension(:),intent(in)::zstore
    !     type(elecintrgl),intent(in)::elecs
    !     real(wp),device,::ovlptot_var,hamtot_var
    !     real(wp),device,dimension(norb)::ovlptot_dx,hamtot_dx
    !     integer::j,k
    !     type(dim3)::threads_per_block, ovrlp_dx_block_num

    !     if (errorflag .ne. 0) return

    !     threads_per_block= dim3(16,16)
    !     ovrlp_dx_block_num=ceiling(norb/threads_per_block%x,norb/threads_per_block%y)
      
    !     ovlptot_var=1.0d0 
    !     call overlp_x<<<1, norb>>>(z1d, z2d, ovlptot_var)
    !     ovlptot_dx=1.0d0
    !     call overlp_dx<<<ovrlp_dx_block_num, threads_per_block>>>(z1d, z2d, ovlptot_dx)
        
    !     overlap%x=ovlptot_var
    !     overlap%dx(1:norb)=ovlptot_dx(1:norb)   
       
    ! end function overlap_gpu

    ! function haml_gpu(z1d,z2d,elecs) result(ham_tot)

    !     type(elecintrgl),intent(in)::elecs
    !     type(dual2)::haml_tot
    !     real(wp),device,dimension(elecs%num)::temp_result_x
    !     real(wp),device::htot_x
    !     real(wp),device,dimension(:)::htot_dx
    !     real(wp),device,dimension(elecs%num,norb)::temp_result_dx
    !     type(dim3)::threads_per_block_x,threads_per_block_dx,haml_grid,haml_grid_dx
    !     integer::j,k

    !     if (errorflag .ne. 0) return
    !     ham_tot=0.0d0
      

    !     threads_per_block_x= dim3(32,32)
    !     haml_grid=dim3(ceiling(elecs%num/threads_per_block_x%x),ceiling(norb/threads_per_block_x%y))
    !     threads_per_block_dx= dim3(4,16,16)
    !     haml_grid_dx=dim3(ceiling(elecs%num/threads_per_block_dx%x),&
    !         ceiling(norb/threads_per_block_dx%y),ceiling(norb/threads_per_block_dx%z))

    !     temp_result_x=elecs%integrals

    !     call haml_x<<<haml_grid, threads_per_block_x>>>(z1d, z2d, temp_result_x)

    !     htot_x=0.0d0
    !     !$cuf kernel do <<<*,*>>> reduce(+:htot_x)
    !     do j=1,elecs%num
    !         ham_tot%x=ham_tot%x+temp_result_x(j)
    !     end do
      
    !     !$cuf kernel do(2) <<<*,*>>>
    !     do k=1,norb
    !         do j=1,elecs%num
    !             temp_result_dx(j,k)=elecs%integrals(j)
    !         end do 
    !     end do

    !     call haml_dx<<<haml_grid_dx, threads_per_block_dx>>>(z1d, z2d, elecs, temp_result_dx)
    !     htot_dx=0.d0
    !     !$cuf kernel do(2) <<<*,*>>>
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
