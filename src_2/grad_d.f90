Module grad_d

    use globvars
    ! use 

    contains

    ! Calcualtes the derivative when the d vector is normalised
    subroutine d_normalise_diff(dvec,ham,step)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        integer, intent(in)::step
        real(kind=8),dimension(ndet,norb)::diff_norm
        integer::j,k

        if (errorflag .ne. 0) return

        

        if(step.eq.0)then
            do j=1, ndet
                diff_norm(j,:)=0.0
                do k=1, ndet
                    diff_norm(j,:)=diff_norm(j,:)+(ham%diff_ovrlp(j,k,:)/abs(REAl(ham%ovrlp(j,k))))
                end do
                dvec%d_diff(j,:)=(REAL(dvec%d(j))*(dvec%norm)*diff_norm(j,:))
            end do
        else
            diff_norm(:,:)=0.0
            call diff_of_norm(dvec,ham,diff_norm)

            do j=1, ndet
                dvec%d_diff(j,:)=((dvec%norm*dvec%d_diff(j,:))-(Real(dvec%d(j))*diff_norm(j,:)))/((dvec%norm)**2)
            end do
        end if

        return

    end subroutine d_normalise_diff

    ! Level 2 subroutine calcualting the compondent of the derivaive found at normalisation 
    ! caused by differentiating the overlap matrix for each j norm^2 *sum_{m} d(overlap_{jm})/(overlap_jm)
    subroutine  diff_of_norm(dvec,ham,diff_norm)

        implicit none
        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        real(kind=8),dimension(ndet,norb),intent(inout)::diff_norm
        real(kind=8),dimension(norb)::temp
        integer::j,k

        if (errorflag .ne. 0) return

        
        do j=1, ndet
            temp(1:norb)=0.0
            do k=1, ndet
                temp=temp+(ham%diff_ovrlp(j,k,:)/abs(REAl(ham%ovrlp(j,k))))
            end do

            diff_norm(j,:)=(1/(dvec%norm))*(&
                (dvec%d_diff(j,:)*(REAL(dvec%d(j))*abs(dot_product(dvec%d,ham%ovrlp(j,:)))))/(abs(dvec%d(j)))+&
                (dvec%norm**2)*temp)
        end do
       
        return
        
    end subroutine

    subroutine timestep_diff(dvec,ham,db)

        implicit none

        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        real, intent(in):: db
        real(kind=8),dimension(norb)::temp1,temp2,temp3
        real(kind=8),dimension(ndet)::invd, hamd
        integer::j,k

        if (errorflag .ne. 0) return

        invd=REAL(matmul(ham%inv,dvec%d))
        hamd=REAL(matmul(ham%hjk,dvec%d))

        do j=1, ndet
            temp1=0.0
            temp2=0.0
            temp3=0.0
            do k=1, ndet
                temp1=temp1 + (real(dvec%d(k))*ham%diff_inv(j,k,:))
                temp2=temp2 + (invd(k)*ham%diff_hjk(j,k,:))
                temp3=temp3 + REAL(ham%kinvh(k,j))*dvec%d_diff(j,:)
            end do
            dvec%d_diff(j,:)=dvec%d_diff(j,:)-(temp1+temp2+temp3)*db
        end do


        return


    end subroutine timestep_diff

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