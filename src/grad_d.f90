Module grad_d

    use globvars
    use infnan_mod
    contains

    ! Level 1 subroutine calcualtes the derivative when the d vector is normalised
    ! The value of the derrivative for each d_(k) is stored when differentiating by each ZS_(j)
    ! Array is stored d-vector compondent(k), differentiation wih respect to zombie state j,
    ! then differentiate W.R.T each coeficient in j
    subroutine d_normalise_diff(dvec,haml,step,diff_state,orb)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        integer, intent(in)::step,diff_state,orb
        real(kind=8),dimension(norb)::diff_ovrlp_cmpnt, diff_d_cmpnt,total
        real(kind=8)::factor
        integer::k,orblim,orbsrt

        if (errorflag .ne. 0) return
        if(diff_state.eq.0) return
        
        diff_ovrlp_cmpnt(:)=0.0

        if(orb.eq.0)then
            orbsrt=1
            orblim=norb
        else 
            orbsrt=orb
            orblim=orb
        end if

       
        call diff_of_norm_ovrlp_cmpndt(dvec,haml,diff_ovrlp_cmpnt,diff_state,orbsrt,orblim)

        if(step.eq.0)then
            do k=1, ndet
                factor=(dvec%d(k))/(2*((dvec%norm)**3))
                dvec%d_diff(k,diff_state,:)=factor*diff_ovrlp_cmpnt
            end do
        else
            call diff_of_norm_d_cmpndt(dvec,haml,diff_d_cmpnt,diff_state,orbsrt,orblim)
            total=(diff_d_cmpnt+diff_ovrlp_cmpnt)/(2*(dvec%norm))
            do k=1, ndet !for each d component
               !for differentiation w.r.t to each ZS
                dvec%d_diff(k,diff_state,:)=(((dvec%norm)*dvec%d_diff(k,diff_state,:))-((dvec%d(k))*total))/((dvec%norm)**2)    
                ! end do
            end do
        end if
      

        return

    end subroutine d_normalise_diff

    ! Level 2 subroutine calcualting compondent of the derivaive found at normalisation. 
    ! Finds sum_{lk} (d^{l*}*d^{k}*overlap_{lk}*d(overlap_{lk}))/(overlap_{lk}). 
    ! This value is the same for each d_(i).
    subroutine  diff_of_norm_ovrlp_cmpndt(dvec,haml,diff_norm_cmpndt,diff_state,orbsrt,orblim)

        implicit none
        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state,orbsrt,orblim
        real(kind=8),dimension(norb),intent(inout)::diff_norm_cmpndt
        real(kind=8),dimension(ndet)::temp
        ! real(kind=8),dimension(norb)::temp
        integer::l,j,p

        if (errorflag .ne. 0) return


      
        temp=0.0
        do j=orbsrt, orblim
            temp=0.0
            do l=1, ndet
                do p=1,ndet
                    temp(l)=temp(l)+(dvec%d(p)*haml%diff_ov_dov(diff_state,p,j,l))
                end do
            end do
            do p=1,ndet 
                diff_norm_cmpndt(j)=diff_norm_cmpndt(j)+dvec%d(p)*temp(p)
            end do
        end do
                    
           
        return
        
    end subroutine diff_of_norm_ovrlp_cmpndt

    ! Level 2 subroutine calcualting compondent of the derivaive found at normalisation. 
    ! Finds sum_{lk} (d(d^{l*})*d^{l*}*|d^{k}*overlap_{lk}|/|(d_{l}|). Since we are dealing with
    ! We exploit the fact that diff(overlap) W.R.T to each Zombie state k is a matrix of zeros 
    ! except along column and row k.
    ! This value is the same for each d_(i). 
    subroutine diff_of_norm_d_cmpndt(dvec,haml,diff_norm_cmpndt,diff_state,orbsrt,orblim)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        real(kind=8),dimension(norb),intent(inout)::diff_norm_cmpndt
        integer,intent(in)::diff_state,orbsrt,orblim
        real(kind=8),dimension(norb)::temp
        real(kind=8),dimension(ndet)::ovrd
        integer::j,l,k,p

        if (errorflag .ne. 0) return

  
        temp=0.0
        ovrd=abs(matmul(dvec%d,(haml%ovrlp)))
        do j=orbsrt, orblim
            do l=1,ndet 
                temp(j) = temp(j) + ((dvec%d_diff(l,diff_state,j)*dvec%d(l))*ovrd(l)/abs(dvec%d(l)))
            end do
        end do
        ovrd=0.0
        do j=orbsrt, orblim
            do k=1,ndet 
                do p=1,ndet
                    ovrd(k)=ovrd(k)+((dvec%d(p)*dvec%d_diff(p,diff_state,j)*abs(haml%ovrlp(p,k)))/abs(dvec%d(p)))
                end do
            end do 
            do p=1,ndet
                temp(j)=temp(j) + (abs(dvec%d(p))*ovrd(p))
            end do
        end do 

        diff_norm_cmpndt=temp
        

        return

    end subroutine diff_of_norm_d_cmpndt

    ! Level 1 subroutine calcualting compondent of the derivaive found when taking an imaginary timestep. 
    ! The time step takes d_({k}p-1) to d_({k}p).
    ! The component caused by normalisation of d_({j}p-1) has already been calculated.
    ! d(d_({k}_p)= d_({k}_p-1)- sum_{l}[ omega^-1_{kl}*d(omega_{kl})*H_{kl}*d_({l}p-1)/norm +
    ! omega^-1_{kl}*d(H_{kl})*d_({l}p-1)/norm + omega^-1_{kl}*H_{kl}*d(d_({l}p-1))]* dbeta
    ! k over all values in d, j differentiation with respect to ZS j
    subroutine timestep_diff(dvec,haml,db,diff_state,orb)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        real, intent(in):: db
        integer,intent(in)::diff_state,orb
        real(kind=8),dimension(ndet,norb)::diff_ts_invo, diff_ts_ham, diff_ts_d
        integer::k,orblim,orbsrt

        if (errorflag .ne. 0) return
        if(diff_state.eq.0) return

        diff_ts_invo=0
        diff_ts_ham=0
        diff_ts_d=0
        if(orb.eq.0)then
            orbsrt=1
            orblim=norb
        else 
            orbsrt=orb
            orblim=orb
        end if
   
       
        
        call timestep_diff_invovrlp_cmpnt(dvec,haml,diff_ts_invo,diff_state)
        call timestep_diff_ham_cmpnt(dvec,haml,diff_ts_ham,diff_state,orbsrt,orblim)
        call timestep_diff_d_cmpnt(dvec,haml,diff_ts_d,diff_state,orbsrt,orblim)
        do k=1,ndet !Over each d_{k}
            !Differentiate with respect to each ZS_{j}
            dvec%d_diff(k,diff_state,:)=dvec%d_diff(k,diff_state,:)-((diff_ts_invo(k,:)+diff_ts_ham(k,:)+diff_ts_d(k,:))*db)
        end do
      

        return

    end subroutine timestep_diff

    ! Level 2 subroutine calcualting compondent of the time step derivaive differentiating the inverse overlap matrix. 
    ! The time step takes d_({k}p-1) to d_({k}p). 
    ! The component caused by normalisation of d_({k}p-1) has already been calculated.
    ! sum_{l}[ omega^-1_{kl}*d(omega_{kl})*H_{kl}*d_({l}p-1)/norm ]

    subroutine timestep_diff_invovrlp_cmpnt(dvec,haml,ts_diff_cmpnt,diff_state)
        
        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state
        real(kind=8),dimension(:,:),intent(inout)::ts_diff_cmpnt
        integer::k,l

        if (errorflag .ne. 0) return

        do k=1, ndet !d_{k}
           !Find dependence on jth ZS
                do l=1, ndet
                    ts_diff_cmpnt(k,:)= ts_diff_cmpnt(k,:)+ (haml%diff_invh(diff_state,k,:,l)*(dvec%d(l)))
                end do
        end do

        return


    end subroutine timestep_diff_invovrlp_cmpnt

    ! Level 2 subroutine calcualting compondent of the time step derivaive differentiating the hamiltonian matrix. 
    ! The time step takes d_({k}p-1) to d_({k}p). 
    ! The component caused by normalisation of d_({k}p-1) has already been calculated.
    !  sum_{l}[ omega^-1_{kl}*d(H_{kl})*d_({l}p-1)/norm ]
    ! k specifies d_{k}, d(ham)/dj is zero except in jth column and row l is used to sum over all d_{l} compondents, 
    ! j specifies the dependence jth ZS
    ! So when j=k can move along the row k of d(ham)/dj.
    subroutine timestep_diff_ham_cmpnt(dvec,haml,ts_diff_cmpnt,diff_state,orbsrt,orblim)
        
        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state,orbsrt,orblim
        real(kind=8),dimension(:,:),intent(inout)::ts_diff_cmpnt
        integer::j,k,p
        ! real(kind=8),dimension(ndet)::temp

        if (errorflag .ne. 0) return

        ! temp=0
        do j=orbsrt, orblim
            do k=1,ndet 
                do p=1,ndet 
                    ts_diff_cmpnt(k,j)=ts_diff_cmpnt(k,j)+(haml%diff_in_dhjk(diff_state,p,j,k)*dvec%d(p))
                end do 
            end do 
        end do

    
        return

    end subroutine timestep_diff_ham_cmpnt

    ! Level 2 subroutine calcualting compondent of the time step derivaive differentiating d vectors. 
    ! The time step takes d_({j}p-1) to d_({j}p). 
    ! The component caused by normalisation of d_({j}p-1) has already been calculated.
    ! sum_{k}[ omega^-1_{jk}*H_{jk}*d(d_({k}p-1))]
    ! k specifies d_{k}, l is used to sum over all d_{l} compondents, j specifies the dependence jth ZS
    subroutine timestep_diff_d_cmpnt(dvec,haml,ts_diff_cmpnt,diff_state,orbsrt,orblim)
        
        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state,orbsrt,orblim
        real(kind=8),dimension(:,:),intent(inout)::ts_diff_cmpnt
        integer::j,k,l

        if (errorflag .ne. 0) return
        do j=orbsrt, orblim
            do k=1,ndet !Extracting compondent for d_{k}
                do l=1,ndet !Sum over all compondents of d
                    ts_diff_cmpnt(k,j)= ts_diff_cmpnt(k,j) + ((haml%kinvh(k,l))*dvec%d_diff(l,diff_state,j))
                end do
            end do
        end do

        return
            
    end subroutine timestep_diff_d_cmpnt

    subroutine final_grad(dvec,haml,grad_fin,diff_state,orb)!,strt)

        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::haml
        type(grad),intent(inout)::grad_fin
        integer,intent(in)::diff_state,orb!,strt
        integer::j,p,ierr,bgn,end
        real(kind=8),dimension(:),allocatable::dh_temp!,dh_temp_hess
        real(kind=8),dimension(:),allocatable::dham
        real(kind=8)::ov
        ! real(kind=8),allocatable,dimension(:,:)::temp
        ! integer, allocatable,dimension(:)::IPIV1
        ! real(kind=8),allocatable,dimension(:)::WORK1
        ! if(strt.eq.0)then
        !     bgn=1
        !     end=norb
        ! else if(strt.eq.1)then 
        !     bgn=1
        !     end=norb/2
        ! elseif(strt.eq.2)then  
        !     bgn=(norb/2)+1
        !     end=norb
        ! end if

        if (errorflag .ne. 0) return

        ierr=0
        allocate(dh_temp(ndet),dham(ndet),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"In dh_temp, dham allocation",ierr
            errorflag=1
        end if
       
        ! grad_fin%vars(diff_state,:)=0
        dham=2*matmul(dvec%d,haml%hjk)
        
        if(orb.eq.0)then
            grad_fin%vars(diff_state,:)=0
            do j=1,norb!bgn,end
            ! do j=1, norb
                dh_temp=dvec%d*haml%diff_hjk(diff_state,j,:)
                ov=0   
             
                do p=1,ndet
                    ov=ov+(dvec%d(p)*haml%diff_hjk(diff_state,j,p))
                end do
              
                dh_temp(diff_state)=ov
                ov=0
                
                do p=1,ndet
                    ov=ov+(dvec%d(p)*dh_temp(p)+(dham(p)*dvec%d_diff(p,diff_state,j)))
                    ! grad_fin%vars(diff_state,j)=grad_fin%vars(diff_state,j)+dvec%d(p)*dh_temp(p)+&
                        ! (dham(p)*dvec%d_diff(p,diff_state,j))
                end do
                grad_fin%vars(diff_state,j)=ov
            end do
           
        else 
            dh_temp=dvec%d*haml%diff_hjk(diff_state,orb,:)   
            ov=0
            do p=1,ndet
                ov=ov+(dvec%d(p)*haml%diff_hjk(diff_state,orb,p))
            end do

            dh_temp(diff_state)=ov
            ov=0
            do p=1,ndet
                ov=ov+dvec%d(p)*dh_temp(p)+(dham(p)*dvec%d_diff(p,diff_state,orb))
            end do
            grad_fin%vars(diff_state,orb)=ov
        end if
        
        deallocate(dh_temp,dham,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"In dh_temp, dham deallocation",ierr
            errorflag=1
        end if

        return 
       

        ! if(orb.eq.0)then
        !     if(typ.eq.0) then
        !         ierr=0
        !         allocate(temp(norb,norb))
        !         temp=0
            
        !         do j=1,norb
        !             do k=j,norb
        !                 dh_temp_hess=dvec%d*haml%hess_hjk(diff_state,j,k,:)   
        !                 dh_temp_hess(diff_state)=0

        !                 do p=1,ndet
        !                     dh_temp_hess(diff_state)=dh_temp_hess(diff_state)+(dvec%d(p)*haml%hess_hjk(diff_state,j,k,p))
        !                 end do

        !                 do p=1,ndet
        !                     temp(j,k)=temp(j,k)+dvec%d(p)*dh_temp_hess(p)
        !                 end do
        !                 temp(k,j)=temp(j,k)
        !             end do 
                    
        !         end do

        !         ierr=0
            
        !         allocate(WORK1(norb),IPIV1(norb),stat=ierr)
        !         if (ierr/=0) then
        !             write(0,"(a,i0)") "Error in IPIV or WORK1 vector allocation . ierr had value ", ierr
        !             errorflag=1
        !         end if 

        !         Call dgetrf(norb, norb, temp, norb, IPIV1, ierr)
        !         if (ierr/=0) then
        !             write(0,"(a,i0)")"Error in DGETRF",ierr
        !         end if
        !         if (ierr==0) call dgetri(norb,temp,norb,IPIV1,WORK1,norb,ierr)
        !         if (ierr/=0) then
        !             write(0,"(a,i0)")"Error in DGETRF ",ierr
        !         end if

        !         grad_fin%vars_hess(diff_state,:)=0
        !         grad_fin%hessian(diff_state,:,:)=temp
        !         ! grad_fin%vars_hess(diff_state,:)=grad_fin%vars(diff_state,:)
            
        !         do j=1,norb
        !             do k=1,norb
        !                 grad_fin%vars_hess(diff_state,j)=grad_fin%vars_hess(diff_state,j)+temp(j,k)*grad_fin%vars(diff_state,k)
        !             end do
        !             if(is_nan(grad_fin%vars_hess(diff_state,j)).eqv..true.)then
        !                 ! print*,'ingradd'
        !                 grad_fin%vars_hess(diff_state,:)=grad_fin%vars(diff_state,:)
        !                 Exit
        !             end if
        !         end do
        !         ! grad_fin%hess_sum(diff_state)=1
                
        !         ! do j=1,norb 
        !         !     grad_fin%hess_sum(diff_state)= grad_fin%hess_sum(diff_state)*temp(j,j)
        !         !     if(j.ne.IPIV1(j))then
        !         !         grad_fin%hess_sum(diff_state)= grad_fin%hess_sum(diff_state)*(-1)
        !         !     end if 
        !         ! end do 



        !         ! grad_fin%hess_sum(diff_state)= sum(temp)
        !         deallocate(WORK1,IPIV1,temp)
        !     else
        !         grad_fin%vars_hess(diff_state,:)=0
        !         do j=1,norb
        !             do k=1,norb
        !                 grad_fin%vars_hess(diff_state,j)=grad_fin%vars_hess(diff_state,j)+&
        !                 grad_fin%hessian(diff_state,j,k)*grad_fin%vars(diff_state,k)
        !             end do
        !             if(is_nan(grad_fin%vars_hess(diff_state,j)).eqv..true.)then
        !                 ! print*,'ingradd'
        !                 grad_fin%vars_hess(diff_state,:)=grad_fin%vars(diff_state,:)
        !                 Exit
        !             end if
        !         end do

        !     end if 
        ! else
            ! grad_fin%vars_hess(diff_state,:)=grad_fin%vars(diff_state,:)
        ! end if 
          
        ! return
    end subroutine final_grad

    

    END MODULE grad_d