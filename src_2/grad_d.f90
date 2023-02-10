Module grad_d

    use globvars
 
    contains

    ! Level 1 subroutine calcualtes the derivative when the d vector is normalised
    ! The value of the derrivative for each d_(k) is stored when differentiating by each ZS_(j)
    ! Array is stored d-vector compondent(k), differentiation wih respect to zombie state j,
    ! then differentiate W.R.T each coeficient in j
    subroutine d_normalise_diff(dvec,haml,step,diff_state)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        integer, intent(in)::step,diff_state
        real(kind=8),dimension(norb)::diff_ovrlp_cmpnt, diff_d_cmpnt,total
        real(kind=8)::factor
        integer::k

        if (errorflag .ne. 0) return
        if(diff_state.eq.0) return
        diff_ovrlp_cmpnt(:)=0.0

       
        call diff_of_norm_ovrlp_cmpndt(dvec,haml,diff_ovrlp_cmpnt,diff_state)

        if(step.eq.0)then
            do k=1, ndet
                factor=(dvec%d(k))/(2*((dvec%norm)**3))
                dvec%d_diff(k,diff_state,:)=factor*diff_ovrlp_cmpnt
            end do
        else
            call diff_of_norm_d_cmpndt(dvec,haml,diff_d_cmpnt,diff_state)
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
    subroutine  diff_of_norm_ovrlp_cmpndt(dvec,haml,diff_norm_cmpndt,diff_state)

        implicit none
        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state
        real(kind=8),dimension(norb),intent(inout)::diff_norm_cmpndt
        real(kind=8),dimension(ndet)::temp
        integer::l,j,p

        if (errorflag .ne. 0) return


      
        temp=0.0
        do j=1,norb
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
    subroutine diff_of_norm_d_cmpndt(dvec,haml,diff_norm_cmpndt,diff_state)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        real(kind=8),dimension(norb),intent(inout)::diff_norm_cmpndt
        integer,intent(in)::diff_state
        real(kind=8),dimension(norb)::temp
        real(kind=8),dimension(ndet)::ovrd
        integer::j,l,k,p

        if (errorflag .ne. 0) return

  
        temp=0.0
        ovrd=abs(matmul(dvec%d,(haml%ovrlp)))
        do j=1,norb
            do l=1,ndet 
                temp(j) = temp(j) + ((dvec%d_diff(l,diff_state,j)*dvec%d(l))*ovrd(l)/abs(dvec%d(l)))
            end do
        end do
        ovrd=0.0
        do j=1,norb
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
    subroutine timestep_diff(dvec,haml,db,diff_state)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::haml
        real, intent(in):: db
        integer,intent(in)::diff_state
        real(kind=8),dimension(ndet,norb)::diff_ts_invo, diff_ts_ham, diff_ts_d
        integer::k

        if (errorflag .ne. 0) return
        if(diff_state.eq.0) return

        diff_ts_invo=0
        diff_ts_ham=0
        diff_ts_d=0
   

        
        call timestep_diff_invovrlp_cmpnt(dvec,haml,diff_ts_invo,diff_state)
        call timestep_diff_ham_cmpnt(dvec,haml,diff_ts_ham,diff_state)
        call timestep_diff_d_cmpnt(dvec,haml,diff_ts_d,diff_state)
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
    subroutine timestep_diff_ham_cmpnt(dvec,haml,ts_diff_cmpnt,diff_state)
        
        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state
        real(kind=8),dimension(:,:),intent(inout)::ts_diff_cmpnt
        integer::j,k,p
        real(kind=8),dimension(ndet)::temp

        if (errorflag .ne. 0) return

        temp=0
        do j=1,norb 
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
    subroutine timestep_diff_d_cmpnt(dvec,haml,ts_diff_cmpnt,diff_state)
        
        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::haml
        integer,intent(in)::diff_state
        real(kind=8),dimension(:,:),intent(inout)::ts_diff_cmpnt
        integer::j,k,l

        if (errorflag .ne. 0) return
        do j=1,norb
            do k=1,ndet !Extracting compondent for d_{k}
                do l=1,ndet !Sum over all compondents of d
                    ts_diff_cmpnt(k,j)= ts_diff_cmpnt(k,j) + ((haml%kinvh(k,l))*dvec%d_diff(l,diff_state,j))
                end do
            end do
        end do

        return
            
    end subroutine timestep_diff_d_cmpnt

    subroutine final_grad(dvec,haml,grad_fin,diff_state,d_diff_flg)

        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::haml
        type(grad),intent(inout)::grad_fin
        integer,intent(in)::diff_state,d_diff_flg
        integer::j,k,p
        real(kind=8),dimension(ndet)::temp2,dh_temp
        real(kind=8),dimension(ndet)::dham

        if (errorflag .ne. 0) return

       
        grad_fin%vars(diff_state,:)=0
        
        if(d_diff_flg.eq.0)then
            do j=1, norb
                dh_temp=dvec%d*haml%diff_hjk(diff_state,j,:)   
                dh_temp(diff_state)=0
                do p=1,ndet
                    dh_temp(diff_state)=dh_temp(diff_state)+(dvec%d(p)*haml%diff_hjk(diff_state,j,p))
                end do
                
                do p=1,ndet
                    grad_fin%vars(diff_state,j)=grad_fin%vars(diff_state,j)+dvec%d(p)*dh_temp(p)
                end do
            end do
        else
            dham=2*matmul(dvec%d,haml%hjk)
            do j=1, norb

                dh_temp=dvec%d*haml%diff_hjk(diff_state,j,:)   
                dh_temp(diff_state)=0
                do p=1,ndet
                    dh_temp(diff_state)=dh_temp(diff_state)+(dvec%d(p)*haml%diff_hjk(diff_state,j,p))
                end do
                  
                do p=1,ndet
                    grad_fin%vars(diff_state,j)=grad_fin%vars(diff_state,j)+dvec%d(p)*dh_temp(p)
                end do

                do p=1,ndet
                    grad_fin%vars(diff_state,j)=grad_fin%vars(diff_state,j)+(dvec%d(p)*dh_temp(p))+&
                    (dham(p)*dvec%d_diff(p,diff_state,j))
                end do
            end do

        end if 
          
        ! print*,grad_fin%vars(diff_state,:)
        return
    end subroutine final_grad

    ! subroutine final_grad_gpu(pvars,phjk,pdiff_hjk,d_diff,d,diff_state,d_diff_flg)

    !     implicit none

    !     real(kind=8),dimension(:,:),intent(inout)::pvars
    !     complex(kind=8), dimension(:,:),intent(inout)::phjk
    !     real(kind=8), dimension(:,:,:),intent(inout)::pdiff_hjk
    !     complex(kind=8), dimension(:)::d
    !     real(kind=8), dimension(:,:,:)::d_diff
    !     integer,intent(in)::diff_state,d_diff_flg
    !     integer::j,l
    !     real(kind=8),dimension(norb)::temp1,temp2
    !     real(kind=8),dimension(ndet)::dham

    !     if (errorflag .ne. 0) return
    !     dham=matmul(REAL(d),REAL(phjk))
    !     j=diff_state !do j=1, ndet !Each ZS{j} dependence

    !     !!$omp target map(to:d,dham,j,d_diff) map(alloc:temp1(norb),temp2(norb))
    !     temp1(:)=0
    !     temp2(:)=0
    !     !$omp parallel do
    !     do l=1, ndet
    !         temp1 = temp1 + dham(l)*d_diff(l,diff_state,:)
    !         if(l.eq.diff_state)then
    !             temp2=temp2+(real(d(l))*matmul(real(d),pdiff_hjk(diff_state,:,:)))
    !         else 
    !             temp2=temp2+(real(d(l)*d(diff_state))*pdiff_hjk(diff_state,l,:))
    !         end if
    !     end do
    !     !$omp end parallel do
    !     if(d_diff_flg.eq.0)then 
    !         temp1=0 
    !     end if
    !     pvars(diff_state,:)=(2*temp1)+temp2
    !     !$omp target update to(pvars)
    !     !!$omp end target

    !     return
    ! end subroutine final_grad_gpu

    END MODULE grad_d