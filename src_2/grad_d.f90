Module grad_d

    use globvars

    contains

    ! Level 1 subroutine calcualtes the derivative when the d vector is normalised
    ! The value of the derrivative for each d_(k) is stored when differentiating by each ZS_(j)
    ! Array is stored d-vector compondent(k), differentiation wih respect to zombie state j,
    ! then differentiate W.R.T each coeficient in j
    subroutine d_normalise_diff(dvec,ham,step)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        integer, intent(in)::step
        real(kind=8),dimension(ndet,norb)::diff_ovrlp_cmpnt, diff_d_cmpnt,total
        real(kind=8)::factor
        integer::j,k

        if (errorflag .ne. 0) return

        diff_ovrlp_cmpnt(:,:)=0.0
        call diff_of_norm_ovrlp_cmpndt(dvec,ham,diff_ovrlp_cmpnt)

        if(step.eq.0)then
            do k=1, ndet
                factor=(Real(dvec%d(k)))/(2*((dvec%norm)**3))
                do j=1,ndet
                    dvec%d_diff(k,j,:)=factor*diff_ovrlp_cmpnt(j,:)
                end do
            end do
        else
            call diff_of_norm_d_cmpndt(dvec,ham,diff_d_cmpnt)
            total=(diff_d_cmpnt+diff_ovrlp_cmpnt)/(2*(dvec%norm))
            do k=1, ndet
                do j=1,ndet
                    dvec%d_diff(k,j,:)=(((dvec%norm)*dvec%d_diff(k,j,:))-(REAL(dvec%d(k))*total(j,:)))/((dvec%norm)**2)
                end do
            end do
        end if

        return

    end subroutine d_normalise_diff

    ! Level 2 subroutine calcualting compondent of the derivaive found at normalisation. 
    ! Finds sum_{lk} (d^{l*}*d^{k}*overlap_{lk}*d(overlap_{lk}))/(overlap_{lk}). 
    ! This value is the same for each d_(i).
    subroutine  diff_of_norm_ovrlp_cmpndt(dvec,ham,diff_norm_cmpndt)

        implicit none
        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        real(kind=8),dimension(ndet,norb),intent(inout)::diff_norm_cmpndt
        real(kind=8),dimension(norb)::temp
        integer::j,k,l

        if (errorflag .ne. 0) return

        
        do j=1, ndet !Each derivative
            temp(1:norb)=0.0
            do l=1, ndet !Over d_{l}
                if(l.eq.j)then !The row with values 
                    do k=1,ndet
                    temp=temp+(REAL(dvec%d(l)*dvec%d(k)*ham%ovrlp(l,k))*ham%diff_ovrlp_bra(l,k,:))/abs(REAL(ham%ovrlp(l,k)))
                
                    end do
                else ! when differentiating by j only values in jth column
                    temp=temp + (REAL(dvec%d(l)*dvec%d(j)*ham%ovrlp(l,j))*ham%diff_ovrlp_ket(l,j,:))/abs(REAL(ham%ovrlp(l,j)))
                end if
                diff_norm_cmpndt(j,:)=temp
            end do

        end do
       
        return
        
    end subroutine diff_of_norm_ovrlp_cmpndt

    ! Level 2 subroutine calcualting compondent of the derivaive found at normalisation. 
    ! Finds sum_{lk} (d(d^{l*})*d^{l*}*|d^{k}*overlap_{lk}|/|(d_{l}|). Since we are dealing with
    ! We exploit the fact that diff(overlap) W.R.T to each Zombie state k is a matrix of zeros 
    ! except along column and row k.
    ! This value is the same for each d_(i). 
    subroutine diff_of_norm_d_cmpndt(dvec,ham,diff_norm_cmpndt)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        real(kind=8),dimension(ndet,norb),intent(inout)::diff_norm_cmpndt
        real(kind=8),dimension(norb)::temp
        integer::j,l,k

        if (errorflag .ne. 0) return

        do j=1,ndet !over each derivative j
            temp(1:norb)=0.0
            do l=1,ndet !over each d_{l}
                do k=1,ndet !over each d_{k}
                    temp = temp + ((dvec%d_diff(l,j,:)*REAL(dvec%d(l)*abs(dvec%d(k)*ham%ovrlp(l,k))))/abs(REAL(dvec%d(l))))
                end do
            end do
            diff_norm_cmpndt(j,:)=2*temp
        end do
                

        ! do k=1, ndet
        !     temp(1:norb)=0.0
        !     ovrlpd=REAL(dvec%d)*abs(matmul(REAL(ham%ovrlp),REAL(dvec%d)))
        !     do l=1, ndet
        !             temp=temp+((ovrlpd(l)*dvec%d_diff(l,j,:))/abs(REAL(dvec%d(l))))    
        !             ! ((*REAL(dvec%d(l))*abs(REAL(dvec%d(k))*ham%ovrlp(l,k)))/(abs(REAL(dvec%d(l)))))
        !         end do
        !     diff_norm_cmpndt(j,:)=2*temp
        ! end do

        return

    end subroutine diff_of_norm_d_cmpndt

    ! Level 1 subroutine calcualting compondent of the derivaive found when taking an imaginary timestep. 
    ! The time step takes d_({k}p-1) to d_({k}p).
    ! The component caused by normalisation of d_({j}p-1) has already been calculated.
    ! d(d_({k}_p)= d_({k}_p-1)- sum_{l}[ omega^-1_{kl}*d(omega_{kl})*H_{kl}*d_({l}p-1)/norm +
    ! omega^-1_{kl}*d(H_{kl})*d_({l}p-1)/norm + omega^-1_{kl}*H_{kl}*d(d_({l}p-1))]* dbeta
    ! k over all values in d, j differentiation with respect to ZS j
    subroutine timestep_diff(dvec,ham,db)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        real, intent(in):: db
        real(kind=8),dimension(ndet,ndet,norb)::diff_ts_invo, diff_ts_ham, diff_ts_d
        integer::j,k

        if (errorflag .ne. 0) return

        diff_ts_invo(:,:,:)=0
        diff_ts_ham(:,:,:)=0
        diff_ts_d(:,:,:)=0
        ! invd=REAL(matmul(ham%inv,dvec%d))
        ! hamd=REAL(matmul(ham%hjk,dvec%d))

        
        call timestep_diff_invovrlp_cmpnt(dvec,ham,diff_ts_invo)
        call timestep_diff_ham_cmpnt(dvec,ham,diff_ts_ham)
        call timestep_diff_d_cmpnt(dvec,ham,diff_ts_d)
        do k=1,ndet !Over each d_{k}
            do j=1,ndet !Differentiate with respect to each ZS_{j}
                dvec%d_diff(k,j,:)=dvec%d_diff(k,j,:)-(diff_ts_invo(k,j,:)+diff_ts_ham(k,j,:)+diff_ts_d(k,j,:))*db
            end do
        end do


        return

    end subroutine timestep_diff

    ! Level 2 subroutine calcualting compondent of the time step derivaive differentiating the inverse overlap matrix. 
    ! The time step takes d_({k}p-1) to d_({k}p). 
    ! The component caused by normalisation of d_({k}p-1) has already been calculated.
    ! sum_{l}[ omega^-1_{kl}*d(omega_{kl})*H_{kl}*d_({l}p-1)/norm ]

    subroutine timestep_diff_invovrlp_cmpnt(dvec,ham,ts_diff_cmpnt)
        
        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::ham
        real(kind=8),dimension(ndet,ndet,norb),intent(inout)::ts_diff_cmpnt
        integer::j,k,l

        if (errorflag .ne. 0) return

        do k=1, ndet !d_{k}
            do j=1, ndet !Find dependence on jth ZS
                do l=1, ndet
                    ts_diff_cmpnt(k,j,:)= ts_diff_cmpnt(k,j,:)+ (ham%diff_invh(j,k,l,:)*REAL(dvec%d(l)))
                end do
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
    subroutine timestep_diff_ham_cmpnt(dvec,ham,ts_diff_cmpnt)
        
        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::ham
        real(kind=8),dimension(ndet,ndet,norb),intent(inout)::ts_diff_cmpnt
        integer::j,k,l

        if (errorflag .ne. 0) return

        do k=1, ndet !d_{k}
            do j=1,ndet ! Find dependence on jth zs
                if(k.eq.j)then !if j==k there's a complete row in diff_ham
                    do l=1, ndet !sum over all values in d
                        ts_diff_cmpnt(k,j,:)=ts_diff_cmpnt(k,j,:)+(REAL(ham%inv(j,l))*REAL(dvec%d(l))*ham%diff_hjk_bra(j,l,:))
                    end do
                else !else there is a single value in column j of row k
                    ts_diff_cmpnt(k,j,:)=ts_diff_cmpnt(k,j,:)+(REAL(ham%inv(k,j))*REAL(dvec%d(j))**ham%diff_hjk_ket(k,j,:))
                end if
            end do
        end do 
    
        return

    end subroutine timestep_diff_ham_cmpnt

    ! Level 2 subroutine calcualting compondent of the time step derivaive differentiating d vectors. 
    ! The time step takes d_({j}p-1) to d_({j}p). 
    ! The component caused by normalisation of d_({j}p-1) has already been calculated.
    ! sum_{k}[ omega^-1_{jk}*H_{jk}*d(d_({k}p-1))]
    ! k specifies d_{k}, l is used to sum over all d_{l} compondents, j specifies the dependence jth ZS
    subroutine timestep_diff_d_cmpnt(dvec,ham,ts_diff_cmpnt)
        
        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::ham
        real(kind=8),dimension(ndet,ndet,norb),intent(inout)::ts_diff_cmpnt
        integer::j,k,l

        
        do k=1,ndet !Extracting compondent for d_{k}
            do j=1, ndet !Sort dependence 
                do l=1,ndet !Sum over all compondents of d
                    ts_diff_cmpnt(k,j,:)= ts_diff_cmpnt(k,j,:) + (REAL(ham%kinvh(k,l))*dvec%d_diff(l,j,:))
                end do
            end do
        end do

        return
            
    end subroutine timestep_diff_d_cmpnt

    subroutine final_grad(dvec,ham,grad_fin)

        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::ham
        type(grad),intent(inout)::grad_fin
        integer::j,k,l
        real(kind=8),dimension(norb)::temp1,temp2
        real(kind=8),dimension(ndet)::dham

        dham=matmul(REAL(dvec%d),REAL(ham%hjk))
        do j=1, ndet !Each ZS{j} dependence
            temp1(:)=0
            do k=1, ndet
                temp1 = temp1 + dham(k)*dvec%d_diff(k,j,:)
            end do
            grad_fin%vars(j,:)=2*temp1
        end do

        do j=1,ndet
            temp2(:)=0
            do l=1, ndet
                if(l.eq.j)then
                    do k=1,ndet
                        temp2=temp2 + (REAL(dvec%d(l)*dvec%d(k))*ham%diff_hjk_bra(j,l,:))
                    end do
                else 
                    temp2 = temp2 + (REAL(dvec%d(j)*dvec%d(l))*ham%diff_hjk_ket(l,j,:))
                end if
            end do
            grad_fin%vars(j,:)=grad_fin%vars(j,:)+temp2
        end do

        print*,"here"
        print*,grad_fin%vars(1,1:norb)

    end subroutine final_grad

    subroutine zombie_alter(zstore,grad_fin,gamma)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        type(grad),intent(in)::grad_fin
        real(kind=8),intent(in)::gamma
        integer::j

        do j=1, ndet
            zstore(j)%phi=zstore(j)%phi-(gamma*(grad_fin%vars(j,:)))
            zstore(j)%sin=sin(cmplx(zstore(j)%phi,0.0d0,kind=8))
            zstore(j)%cos=cos(cmplx(zstore(j)%phi,0.0d0,kind=8))
        end do


    end subroutine zombie_alter


END MODUle grad_d