Module grad_d

    use globvars

    contains

    ! Level 1 subroutine calcualtes the derivative when the d vector is normalised
    ! The value of the derrivative for each d_(j) is stored when differentiating by each ZS_(k)
    ! Array is stored d-vector compondent(j), differentiation wih respect to zombie state k,
    ! then differentiate W.R.T each coeficient in k
    subroutine d_normalise_diff(dvec,ham,step)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        integer, intent(in)::step
        real(kind=8),dimension(ndet,norb)::diff_ovrlp_cmpnt, diff_d_cmpnt,total
        real(kind=8)::factor
        integer::j,k,l

        if (errorflag .ne. 0) return

        diff_ovrlp_cmpnt(:,:)=0.0
        call diff_of_norm_ovrlp_cmpndt(dvec,ham,diff_ovrlp_cmpnt)

        if(step.eq.0)then
            do j=1, ndet
                factor=(Real(dvec%d(j)))/(2*((dvec%norm)**3))
                do k=1,ndet
                    dvec%d_diff(j,k,:)=factor*diff_ovrlp_cmpnt(j,:)
                end do
            end do
        else
            call diff_of_norm_d_cmpndt(dvec,ham,diff_d_cmpnt)
            total=(diff_d_cmpnt+diff_ovrlp_cmpnt)/(2*(dvec%norm))
            do j=1, ndet
                do k=1,ndet
                    dvec%d_diff(j,k,:)=(((dvec%norm)*dvec%d_diff(j,k,:))-(REAL(dvec%d(j))*total(j,:)))/((dvec%norm)**2)
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
        real(kind=8),dimension(ndet,norb)::ldifovrlp
        real(kind=8),dimension(norb)::temp
        real(kind=8),dimension(ndet)::ddovrlp,kd,kovrlp
        integer::j,k,l

        if (errorflag .ne. 0) return

        
        do j=1, ndet
            temp(1:norb)=0.0
            kd=REAL(dvec%d)
            kd(j)=1.0
            kovrlp=abs(REAL(ham%ovrlp(j,:)))
            kovrlp(j)=1
            ddovrlp=(REAL(dvec%d)*kd*REAL(ham%ovrlp(j,:))*kovrlp)/(abs(REAl(ham%ovrlp(j,:)))*kovrlp)
            ldifovrlp=ham%diff_ovrlp(j,:,:)
            ldifovrlp(j,:)=1.0
            do k=1, ndet
                temp=temp+ ddovrlp(k)*ham%diff_ovrlp(j,k,:)*ldifovrlp     
            end do
            diff_norm_cmpndt(j,:)=temp
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
        real(kind=8),dimension(ndet)::ovrlpd
        real(kind=8),dimension(norb)::temp
        integer::j,l

        if (errorflag .ne. 0) return

        do j=1, ndet
            temp(1:norb)=0.0
            ovrlpd=REAL(dvec%d)*abs(matmul(REAL(ham%ovrlp),REAL(dvec%d)))
            do l=1, ndet
                    temp=temp+((ovrlpd(l)*dvec%d_diff(l,j,:))/abs(REAL(dvec%d(l))))    
                    ! ((*REAL(dvec%d(l))*abs(REAL(dvec%d(k))*ham%ovrlp(l,k)))/(abs(REAL(dvec%d(l)))))
                end do
            diff_norm_cmpndt(j,:)=2*temp
        end do

        return

    end subroutine diff_of_norm_d_cmpndt

    ! Level 1 subroutine calcualting compondent of the derivaive found when taking an imaginary timestep. 
    ! The time step takes d_({j}p-1) to d_({j}p).
    ! The component caused by normalisation of d_({j}p-1) has already been calculated.
    ! d(d_({j}_p)= d_({j}_p-1)- sum_{k}[ omega^-1_{jk}*d(omega_{jk})*H_{jk}*d_({k}p-1)/norm +
    ! omega^-1_{jk}*d(H_{jk})*d_({k}p-1)/norm + omega^-1_{jk}*H_{jk}*d(d_({k}p-1))]* dbeta
    ! j over all values in d, k differentiation with respect to ZS k
    subroutine timestep_diff(dvec,ham,db)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        real, intent(in):: db
        real(kind=8),dimension(norb)::temp1,temp2,temp3
        real(kind=8),dimension(ndet)::invd, hamd
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
        do j=1,ndet
            do k=1,ndet
                dvec%d_diff(j,k,:)=dvec%d_diff(j,k,:)-(diff_ts_invo(j,k,:)+diff_ts_ham(j,k,:)+diff_ts_d(j,k,:))*db
            end do
        end do


        return

    end subroutine timestep_diff

    ! Level 2 subroutine calcualting compondent of the time step derivaive differentiating the inverse overlap matrix. 
    ! The time step takes d_({j}p-1) to d_({j}p). 
    ! The component caused by normalisation of d_({j}p-1) has already been calculated.
    ! sum_{k}[ omega^-1_{jk}*d(omega_{jk})*H_{jk}*d_({k}p-1)/norm ]

    subroutine timestep_diff_invovrlp_cmpnt(dvec,ham,ts_diff_cmpnt)
        
        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::ham
        real(kind=8),dimension(ndet,ndet,norb),intent(inout)::ts_diff_cmpnt
        integer::j,k

    end subroutine timestep_diff_invovrlp_cmpnt

    ! Level 2 subroutine calcualting compondent of the time step derivaive differentiating the hamiltonian matrix. 
    ! The time step takes d_({j}p-1) to d_({j}p). 
    ! The component caused by normalisation of d_({j}p-1) has already been calculated.
    !  sum_{k}[ omega^-1_{jk}*d(H_{jk})*d_({k}p-1)/norm ]
    ! j specifies d_{j}, d(ham)/dk is zero except in column k and row k  l is used to sum over all d_{l} compondents, k specifies the dependence kth ZS
    ! So when j=k can move along the row j of d(ham)/dk. Else only one none zero value per row in the kth column. 
    subroutine timestep_diff_ham_cmpnt(dvec,ham,ts_diff_cmpnt)
        
        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::ham
        real(kind=8),dimension(ndet,ndet,norb),intent(inout)::ts_diff_cmpnt
        integer::j,k,l

        do j=1, ndet
            do k=1, ndet
                if(k.eq.j)then 
                    do l=1, ndet
                        ts_diff_cmpnt(j,k,:)=ts_diff_cmpnt(j,k,:)+(REAL(ham%inv(j,l))*REAL(dvec%d(l))*ham%diff_hjk(k,l,:))
                    end do
                else 
                    ts_diff_cmpnt(j,j,:)=ts_diff_cmpnt(j,k,:)+(REAL(ham%inv(k,j))*REAL(dvec%d(j))**ham%diff_hjk(k,j,:))
                end if
            end do
        end do 
    
        return

    end subroutine timestep_diff_ham_cmpnt

    ! Level 2 subroutine calcualting compondent of the time step derivaive differentiating d vectors. 
    ! The time step takes d_({j}p-1) to d_({j}p). 
    ! The component caused by normalisation of d_({j}p-1) has already been calculated.
    ! sum_{k}[ omega^-1_{jk}*H_{jk}*d(d_({k}p-1))]
    ! j specifies d_{j}, l is used to sum over all d_{l} compondents, k specifies the dependence kth ZS
    subroutine timestep_diff_d_cmpnt(dvec,ham,ts_diff_cmpnt)
        
        implicit none

        type(dvector),intent(in)::dvec
        type(hamiltonian),intent(in)::ham
        real(kind=8),dimension(ndet,ndet,norb),intent(inout)::ts_diff_cmpnt
        real(kind=8),dimension(norb)::temp
        integer::j,k,l

        do j=1,ndet !Extracting compondent for d_{j}
            do l=1,ndet !Sum over all compondents of d
                do k=1, ndet !Sort dependence 
                    ts_diff_cmpnt(j,k,:)= ts_diff_cmpnt(j,k,:) + (REAL(ham%kinvh(j,l))*dvec%d_diff(l,k,:))
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
        integer::j,k
        real(kind=8),dimension(norb)::temp1,temp2,temp3
        

        do j=1, ndet
            temp1(:)=0
            temp2(:)=0
            temp3(:)=0
            do k=1, ndet
                temp1=temp1+(REAL(dvec%d(j))*REAL(dvec%d(k))*ham%diff_hjk(j,k,:))
                temp2=temp2+(REAL(dvec%d(j))*REAL(ham%hjk(j,k))*dvec%d_diff(k,:))
                temp3=temp3+(REAL(dvec%d(k))*REAL(ham%hjk(j,k))*dvec%d_diff(j,:))
            end do
            grad_fin%vars(j,:)=temp1+temp2+temp3
        end do
        print*,"here"
        print*,grad_fin%vars(1,1:norb)

    end subroutine final_grad


END MODUle grad_d