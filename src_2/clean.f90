MODULE clean

    use globvars
    use alarrays
    use zom
    use ham
    use outputs
    use operators

    contains

    subroutine clean_setup(cstore,nume,cleanham,elecs,clean_ndet)

        implicit none

        type(zombiest),dimension(:),allocatable,intent(out)::cstore
        type(hamiltonian), allocatable, intent(out)::cleanham
        type(elecintrgl),intent(in)::elecs 
        integer, intent(in)::nume
        integer, intent(out)::clean_ndet
        integer, allocatable, dimension(:,:)::combs,combs2
        integer::j,k,l,ierr,total,total2,checker

        if (errorflag .ne. 0) return

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
        call matrixwriter(cleanham%hjk,ndet,"cleanham.csv")
        clean_ndet=total2

        return


    end subroutine clean_setup

    subroutine cleaner(zstore,cleanzom,dvec, dvec_clean, cleantot,norm)

        implicit none
        type(zombiest),dimension(:),intent(in)::zstore,cleanzom
        type(dvector),intent(in)::dvec
        type(dvector),intent(inout):: dvec_clean
        integer,intent(in)::cleantot
        complex(kind=8),intent(out)::norm
        complex(kind=8)::ovrlp1, ovrlp2
        integer::j,k,l
        
        if (errorflag .ne. 0) return
        
        norm=(0.0d0,0.0d0)
        do j=1,cleantot
            do k=1, ndet
                ovrlp1=overlap(cleanzom(j),zstore(k))
                dvec_clean%d(j)=dvec_clean%d(j)+(dvec%d(k)*ovrlp1)
                do l=1, ndet
                    ovrlp2=overlap(zstore(l),cleanzom(j))
                    norm = norm + (conjg(dvec%d(l))*dvec%d(k)*ovrlp2*ovrlp1)
                end do
            end do
        end do

        return

    end subroutine cleaner

END MODULE clean

