Module grad_d

    use globvars
    ! use 

    contains

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

    subroutine  diff_of_norm(dvec,ham,diff_norm)

        implicit none
        type(dvector),intent(inout)::dvec
        type(hamiltonian),intent(in)::ham
        ! integer, intent(in)::step
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


END MODUle grad_d