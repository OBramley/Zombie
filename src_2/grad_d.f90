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
        integer::j,k

        if (errorflag .ne. 0) return
        if(diff_state.eq.0) return
        diff_ovrlp_cmpnt(:)=0.0

       
        call diff_of_norm_ovrlp_cmpndt(dvec,haml,diff_ovrlp_cmpnt,diff_state)

        if(step.eq.0)then
            do k=1, ndet
                factor=(Real(dvec%d(k)))/(2*((dvec%norm)**3))
                j=diff_state!do j=1,ndet
                    dvec%d_diff(k,j,:)=factor*diff_ovrlp_cmpnt
            end do
        else
            call diff_of_norm_d_cmpndt(dvec,haml,diff_d_cmpnt,diff_state)
            total=(diff_d_cmpnt+diff_ovrlp_cmpnt)/(2*(dvec%norm))
            do k=1, ndet !for each d component
                j=diff_state!do j=1,ndet !fir differentiation w.r.t to each ZS
                    dvec%d_diff(k,j,:)=(((dvec%norm)*dvec%d_diff(k,j,:))-(REAL(dvec%d(k))*total))/((dvec%norm)**2)    
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
        real(kind=8),dimension(norb)::temp
        integer::j,k

        if (errorflag .ne. 0) return

        
        j=diff_state!do j=1, ndet !Each derivative
        temp(1:norb)=0.0
        do k=1, ndet !Over d_{k}
            ! do k=1,ndet
            if(k.eq.j)then
                temp=temp+real(dvec%d(k))*matmul(real(dvec%d*(haml%ovrlp(j,:)/abs(haml%ovrlp(j,:)))),haml%diff_ovrlp(j,:,:))
            else 
                temp=temp+real((dvec%d(k)*dvec%d(j)*haml%ovrlp(j,k))/abs(haml%ovrlp(j,k)))*haml%diff_ovrlp(j,k,:)
            end if                
        end do
             
            diff_norm_cmpndt=temp
            
        !end do

     
       
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
        integer::j,l,k

        if (errorflag .ne. 0) return

        j=diff_state !do j=1,ndet !over each derivative j
            temp(1:norb)=0.0
            do l=1,ndet !over each d_{l}
                do k=1,ndet !over each d_{k}
                    temp = temp + ((dvec%d_diff(l,j,:)*REAL(dvec%d(l)*abs(dvec%d(k)*haml%ovrlp(l,k))))/abs(REAL(dvec%d(l))))
                end do
            end do
            diff_norm_cmpndt=2*temp
        !end do
            

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
        real(kind=8),dimension(ndet,ndet,norb)::diff_ts_invo, diff_ts_ham, diff_ts_d
        integer::j,k

        if (errorflag .ne. 0) return
        if(diff_state.eq.0) return

        diff_ts_invo(:,:,:)=0
        diff_ts_ham(:,:,:)=0
        diff_ts_d(:,:,:)=0
   

        
        call timestep_diff_invovrlp_cmpnt(dvec,haml,diff_ts_invo,diff_state)
        call timestep_diff_ham_cmpnt(dvec,haml,diff_ts_ham,diff_state)
        call timestep_diff_d_cmpnt(dvec,haml,diff_ts_d,diff_state)
        do k=1,ndet !Over each d_{k}
            j=diff_state!do j=1,ndet !Differentiate with respect to each ZS_{j}
                dvec%d_diff(k,j,:)=dvec%d_diff(k,j,:)-(diff_ts_invo(k,j,:)+diff_ts_ham(k,j,:)+diff_ts_d(k,j,:))*db
         
            !end do
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
        real(kind=8),dimension(ndet,ndet,norb),intent(inout)::ts_diff_cmpnt
        integer::j,k,l

        if (errorflag .ne. 0) return

        do k=1, ndet !d_{k}
            j=diff_state!do j=1, ndet !Find dependence on jth ZS
                do l=1, ndet
                    ts_diff_cmpnt(k,j,:)= ts_diff_cmpnt(k,j,:)+ (haml%diff_invh(j,k,l,:)*REAL(dvec%d(l)))
                end do
            !end do
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
        real(kind=8),dimension(ndet,ndet,norb),intent(inout)::ts_diff_cmpnt
        integer::j,k

        if (errorflag .ne. 0) return

        do k=1, ndet !d_{k}
            j=diff_state!do j=1,ndet ! Find dependence on jth zs
                if(k.eq.j)then !if j==k there's a complete row in diff_ham
                    ts_diff_cmpnt(k,j,:)=ts_diff_cmpnt(k,j,:)+matmul(real(haml%inv(j,:)*dvec%d),haml%diff_hjk(j,:,:))
                else 
                    ts_diff_cmpnt(k,j,:)=ts_diff_cmpnt(k,j,:)+(real(haml%inv(j,k)*dvec%d(j))*haml%diff_hjk(j,k,:))
                end if
            !end do
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
        real(kind=8),dimension(ndet,ndet,norb),intent(inout)::ts_diff_cmpnt
        integer::j,k,l

        if (errorflag .ne. 0) return

        do k=1,ndet !Extracting compondent for d_{k}
            j=diff_state!do j=1, ndet !Sort dependence 
                do l=1,ndet !Sum over all compondents of d
                    ts_diff_cmpnt(k,j,:)= ts_diff_cmpnt(k,j,:) + (REAL(haml%kinvh(k,l))*dvec%d_diff(l,j,:))
                end do
            !end do
        end do

        return
            
    end subroutine timestep_diff_d_cmpnt

    subroutine final_grad(dvec,haml,grad_fin,diff_state,d_diff_flg)

        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::haml
        type(grad),intent(inout)::grad_fin
        integer,intent(in)::diff_state,d_diff_flg
        integer::j,l
        real(kind=8),dimension(norb)::temp1,temp2
        real(kind=8),dimension(ndet)::dham

        if (errorflag .ne. 0) return

        dham=matmul(REAL(dvec%d),REAL(haml%hjk))
        j=diff_state !do j=1, ndet !Each ZS{j} dependence
        temp1(:)=0
        temp2(:)=0
        do l=1, ndet
            temp1 = temp1 + dham(l)*dvec%d_diff(l,diff_state,:)
            ! print*,temp1
            if(l.eq.diff_state)then
                temp2=temp2+(real(dvec%d(l))*matmul(real(dvec%d),haml%diff_hjk(diff_state,:,:)))
            else 
                temp2=temp2+(real(dvec%d(l)*dvec%d(diff_state))*haml%diff_hjk(diff_state,l,:))
            end if
            ! print*,temp2
        end do
        
        if(d_diff_flg.eq.0)then 
            temp1=0 
        end if
        grad_fin%vars(diff_state,:)=(2*temp1)+temp2
        

        return
    end subroutine final_grad

    subroutine final_grad_gpu(pvars,phjk,pdiff_hjk,d_diff,d,diff_state,d_diff_flg)

        implicit none

        real(kind=8),dimension(:,:), pointer,intent(inout)::pvars
        complex(kind=8), dimension(:,:), pointer,intent(inout)::phjk
        real(kind=8), dimension(:,:,:), pointer,intent(inout)::pdiff_hjk
        complex(kind=8), dimension(:)::d
        real(kind=8), dimension(:,:,:)::d_diff
        integer,intent(in)::diff_state,d_diff_flg
        integer::j,l
        real(kind=8),dimension(norb)::temp1,temp2
        real(kind=8),dimension(ndet)::dham

        if (errorflag .ne. 0) return
        dham=matmul(REAL(d),REAL(phjk))
        j=diff_state !do j=1, ndet !Each ZS{j} dependence

        !!$omp target map(to:d,dham,j,d_diff) map(alloc:temp1(norb),temp2(norb))
        temp1(:)=0
        temp2(:)=0
        !$omp parallel do
        do l=1, ndet
            temp1 = temp1 + dham(l)*d_diff(l,diff_state,:)
            if(l.eq.diff_state)then
                temp2=temp2+(real(d(l))*matmul(real(d),pdiff_hjk(diff_state,:,:)))
            else 
                temp2=temp2+(real(d(l)*d(diff_state))*pdiff_hjk(diff_state,l,:))
            end if
        end do
        !$omp end parallel do
        if(d_diff_flg.eq.0)then 
            temp1=0 
        end if
        pvars(diff_state,:)=(2*temp1)+temp2
        !$omp target update to(pvars)
        !!$omp end target

        return
    end subroutine final_grad_gpu

    END MODULE grad_d