MODULE clean

    use globvars
    use alarrays
    use zom
    use ham
    use outputs

    contains

    subroutine clean_setup(cstore,nume,cleanham,elecs)

        implicit none

        type(zombiest),dimension(:),allocatable,intent(inout)::cstore
        type(hamiltonian), allocatable, intent(inout)::cleanham
        type(elecintrgl),intent(inout)::elecs 
        integer, intent(in)::nume
        integer, allocatable, dimension(:,:)::combs,combs2
        integer::j,k,l,ierr,total,total2,checker


        total=int(gamma(dble(norb+1))/(gamma(dble(nume+1))*gamma(dble(norb-nume)+1)))
        allocate(combs(total,nel),stat=ierr)
        if(ierr==0)  allocate (combs2(total,nel),stat=ierr)
            if(ierr/=0) then
                write(0,"(a,i0)") "Error in combination matrix allocation. ierr had value ", ierr
                errorflag=1
                return
            end if

        call combinations(norb,nume,combs)
        
        l=1
        total2=0
        do j=1, total
            checker=0
            do k=1,nume
                checker=checker+modulo(combs(j,k),2)
            end do
            if(checker==((nume/2)-spin))then
                combs2(l,:)=combs(j,:)
                total2=total2+1
                l=l+1
            end if
        end do

        call alloczs(cstore,total2)
        do k=1, total2
            call zomhf(cstore(k),combs2(k,:))
        end do

        deallocate(combs,stat=ierr)
        if(ierr==0)  deallocate(combs2,stat=ierr)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in combination matrix deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        call allocham(cleanham,total2)
        call hamgen(cleanham,cstore,elecs,total2)
        call matrixwriter(ham%hjk,ndet,"cleanham.csv")

        return


    end subroutine clean_setup


END MODULE clean

