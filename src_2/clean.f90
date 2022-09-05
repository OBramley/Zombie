MODULE clean

    use globvars
    use alarrays
    use zom
    use ham
    use outputs
    use operators

    contains

    subroutine clean_setup(cstore,nume,cleanham,elecs,clean_ndet,zstore)

        implicit none

        type(zombiest),dimension(:),allocatable,intent(out)::cstore
        type(hamiltonian), intent(out)::cleanham
        type(elecintrgl),intent(in)::elecs 
        type(zombiest),dimension(:),intent(in)::zstore
        integer, intent(in)::nume
        integer, intent(out)::clean_ndet
        type(zombiest),dimension(:),allocatable::cstoretemp
        integer, allocatable, dimension(:,:)::combs,combs2
        integer, allocatable, dimension(:)::magovrlp
        integer::j,k,ierr,total,total2,total3,checker
        complex(kind=8)::magnitude


        if (errorflag .ne. 0) return

        total=choose(norb,nume)
        allocate(combs(total,nume),stat=ierr)
        if(ierr==0)  allocate (combs2(total,nume),stat=ierr)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in combination matrix allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        write(6,"(a,i0)") 'Total combinations ',total
        ! combs2(1:total,1:norb)=0

        ! The occupational combiantions for the correct number of electrons are found 
        call combinations(norb,nume,combs,total)
       
        ! These compinations are then reduced down to only those with the correct spin
        total2=0
        !$omp parallel private(j,k,checker) shared(combs,combs2,total2)
        !$omp do
        do j=1, total
            checker=0
            do k=1,nume
                checker=checker+modulo(combs(j,k),2)
            end do
            if(checker==((nume/2)-spin))then
                !$omp critical
                total2=total2+1
                combs2(total2,:)=combs(j,:)
                !$omp end critical
            end if
        end do
        !$omp end do
        !$omp end parallel
       
        write(6,"(a,i0)") 'Combinations with corect spin ',total2
        ! A temporary set of Zombie states is created with the HF states with the correct number of electrons and spin
        call alloczs(cstoretemp,total2)
    
        
       
        !$omp parallel private(k) shared(cstore,combs2)
        !$omp do
        do j=1, total2
            call zomhfc(cstoretemp(j),combs2(j,:))
        end do
        !$omp end do
        !$omp end parallel

        
        deallocate(combs,stat=ierr)
        if(ierr==0)  deallocate(combs2,stat=ierr)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in combination matrix deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if


        ! This section reduces the size of cleaning matrix further by finding the overlap between each cleaning 
        ! state and the biased basis. Only cleanign states with overlaps higher than the set threshold are accepted
        ! and placed into the final cleaning set to produce the cleaning Hamiltonian.
        if(total2>500)then
            total3=0
            allocate(magovrlp(total2),stat=ierr)
            if(ierr/=0) then
                write(0,"(a,i0)") "Error in magovrlp allocation. ierr had value ", ierr
                errorflag=1
                return
            end if

            !$omp parallel private(j,k,magnitude) shared(total3,magovrlp,cstoretemp,zstore)
            !$omp do
            do j=1, total2
                magnitude=(0.0,0.0)
                do k=1, ndet
                    magnitude=magnitude+overlap(cstoretemp(j),zstore(k))
                end do
                if(REAL(magnitude)>0.00001) then
                    !$omp critical
                    total3=total3+1
                    magovrlp(total3)=j
                    !$omp end critical
                end if
            end do
            !$omp end do
            !$omp end parallel

            call alloczs(cstore,total3)

            !$omp parallel private(j) shared(cstoretemp, cstore, magovrlp)
            !$omp do
            do j=1, total3
                cstore(j)=cstoretemp(magovrlp(j))
            end do
            !$omp end do
            !$omp end parallel

            deallocate(magovrlp,stat=ierr)
            if(ierr/=0) then
                write(0,"(a,i0)") "Error in magovrlp deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if
            write(6,"(a,i0)") 'Combinations with corect spin and enough contribution ',total3
            clean_ndet=total3
            
        else
            call alloczs(cstore,total2)
            cstore=cstoretemp
            clean_ndet=total2 
        end if

        call dealloczs(cstoretemp)
        call allocham(cleanham,clean_ndet)
        call hamgen(cleanham,cstore,elecs,clean_ndet)
        call matrixwriter(cleanham%hjk,clean_ndet,"clean_ham.csv")
        
        do j=1,clean_ndet
            call zombiewriter_c(cstore(j),j)
        end do
        return


    end subroutine clean_setup

    subroutine zomhfc(zom,occ)

        implicit none
        type(zombiest),intent(inout)::zom 
        integer, dimension(:), intent(in)::occ
        integer::j

        if (errorflag .ne. 0) return

        zom%dead(1:norb)=(1.0d0,0.0d0)

        do j=1, size(occ)
            zom%alive(occ(j))=(1.0d0,0.0d0)
            zom%dead(occ(j))=(0.0d0,0.0d0)
        end do

        return

    end subroutine zomhfc

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

