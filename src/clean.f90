MODULE clean

    use globvars
    use alarrays
    use zom
    use outputs
    ! use operators
    use readpars
    use ham

    contains

    subroutine sd_anal(zstore,nume,dvec,pass)

        implicit none 

        type(zombiest),dimension(:),allocatable::cstore
        type(zombiest),dimension(:),intent(in)::zstore
        type(dvector),intent(in)::dvec
        integer, intent(in)::nume,pass
        integer, allocatable, dimension(:,:)::combs,combs2
        integer, allocatable, dimension(:)::position
        real(kind=8), allocatable, dimension(:)::magovrlp
        integer::j,k,l,ierr,total,total2,checker
        ! complex(kind=8)::norm
        real(kind=8)::norm
        logical,allocatable,dimension(:)::excld
        ! complex(kind=8)::ovrlp1, ovrlp2
        real(kind=8)::ovrlp1, ovrlp2


        total=choose(norb,nume)
        allocate(combs(total,nume),stat=ierr)
        if(ierr==0)  allocate (combs2(total,nume),stat=ierr)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in combination matrix allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        write(6,"(a,i0)") 'Total combinations ',total


        ! The occupational combiantions for the correct number of electrons are found 
        call combinations(norb,nume,combs,total)
    
        
        total2=0
        do j=1, total
            checker=0
            do k=1,nume
                checker=checker+modulo(combs(j,k),2)
            end do
            if(checker==((nume/2)-spin))then
                total2=total2+1
                combs2(total2,:)=combs(j,:)
            end if
        end do
        write(6,"(a,i0)") 'Total combinations with correct spin ',total2
        call alloczs(cstore,total2)

        allocate(magovrlp(total2),stat=ierr)
        if(ierr==0) allocate(position(total2),stat=ierr)
        if(ierr==0) allocate(excld(total2),stat=ierr)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in magovrlp allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
       
        deallocate(combs,stat=ierr)
   
        do j=1, total2
            call zomhfc(cstore(j),combs2(j,:))
        end do

        ! norm=(0.0d0,0.0d0)
        ! magovrlp=(0.0,0.0)
        norm=0.0
        magovrlp=0.0
        
        excld=.TRUE.
        !!$omp parallel do collapse(2) reduction(+:norm)
        do j=1,total2
            do k=1,ndet 
                ! ovrlp1=overlap(cstore(j),zstore(k))*dvec%d(k)
                ovrlp1=product((cstore(j)%sin*zstore(k)%sin)+(cstore(j)%cos*zstore(k)%cos))*dvec%d(k)
                ! magovrlp(j)=magovrlp(j)+real(ovrlp1)
                magovrlp(j)=magovrlp(j)+ovrlp1
                do l=1,ndet
                    ovrlp2=product((cstore(j)%sin*zstore(l)%sin)+(cstore(j)%cos*zstore(l)%cos))
                    ! ovrlp2=overlap(zstore(l),cstore(j))
                    ! norm=norm + (conjg(dvec%d(l))*ovrlp1*ovrlp2)
                    norm=norm + (dvec%d(l)*ovrlp1*ovrlp2)

                end do
            end do
        end do
        !!$omp end parallel do 
        write(6,"(a,e25.17e3)") 'The norm for states with correct spin and electrons is ',real(norm)
        
 
        if(pass.eq.1)then 
            open(unit=9,file='slt_ovrlp.csv',status="new", iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening slt_ovrlp.csv. ierr had value ", ierr
                errorflag=1
                return
            end if

            write(9,'(a,e25.17e3)') 'Norm before gradient descent is : ', real(norm)
            close(9)
        else
            
            do j=1, total2
                position(j)=minloc(magovrlp,1,excld)
                excld(position(j))=.FALSE.
            end do
            open(unit=9,file='slt_ovrlp.csv',status="old",access='append',iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening slt_ovrlp.csv. ierr had value ", ierr
                errorflag=1
                return
            end if

            write(9,'(a,e25.17e3)') 'Norm after gradient descent is : ', real(norm)
            do j=total2, 1, -1 
                write(9,'(e25.17e3,a,*(i0, :", "))') magovrlp(position(j)),' ', (combs2(position(j),k),k=nume,1,-1)
            end do

            close(9)
        end if 
        call dealloczs(cstore)
        deallocate(position,stat=ierr)
        deallocate(excld,stat=ierr)
        deallocate(magovrlp,stat=ierr)

        return 


    end subroutine sd_anal

    subroutine clean_setup(cstore,nume,cleanham,elecs,clean_ndet,zstore,an_cr,an2_cr2)

        implicit none

        type(zombiest),dimension(:),allocatable,intent(inout)::cstore
        type(hamiltonian), intent(inout)::cleanham
        integer, intent(inout)::clean_ndet
        type(elecintrgl),intent(in)::elecs 
        type(oprts),intent(in)::an_cr,an2_cr2
        type(zombiest),dimension(:),intent(in)::zstore
        integer, intent(in)::nume
        type(zombiest),dimension(:),allocatable::cstoretemp
        integer, allocatable, dimension(:,:)::combs,combs2,combsfix
        integer, allocatable, dimension(:)::magovrlp
        integer::j,k,ierr,total,total2,total3,totalf,checker
        ! complex(kind=8)::magnitude
        real(kind=8)::magnitude


        if (errorflag .ne. 0) return

      
        

        total=choose(norb,nume)
        allocate(combs(total,nume),stat=ierr)
        if(ierr==0)  allocate (combs2(total,nume),stat=ierr)
        if(ierr==0)  allocate (combsfix(total,nume),stat=ierr)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in combination matrix allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        write(6,"(a,i0)") 'Total combinations ',total


        ! The occupational combiantions for the correct number of electrons are found 
        call combinations(norb,nume,combs,total)
        
        
        totalf=0
        total2=0
        total3=0
        !$omp parallel private(j,k,checker) &
        !$omp shared(combs,combsfix,combs2,total2,totalf)
        
        ! This additional bit of code allows specific orbitals to be fixed as occupied
        !$omp do
        do j=1, total
            if((ANY(combs(j,:).eq.1)).and.(ANY(combs(j,:).eq.2)))then
                !$omp critical
                totalf=totalf+1
                combsfix(totalf,:)=combs(j,:)
                !$omp end critical
            end if
        end do
        !$omp end do

        ! This can be turn off by making combsfix=combs and totalf=total
        ! totalf=total
        ! combsfix=combs

        !$OMP barrier
        !$omp master
        write(6,"(a,i0)") 'Combinations with the first orbital fixed as occupied ',totalf
        ! A temporary set of Zombie states is created with the HF states with the correct number of electrons and spin
        !$omp end master

        ! These combinations are then reduced down to only those with the correct spin

        !$omp do
        do j=1, totalf
            checker=0
            do k=1,nume
                checker=checker+modulo(combs(j,k),2)
            end do
            if(checker==((nume/2)-spin))then
                !$omp critical
                total2=total2+1
                combs2(total2,:)=combsfix(j,:)
                !$omp end critical
            end if
        end do
        !$omp end do
        !$omp end parallel

    
        write(6,"(a,i0)") 'Combinations with corect spin ',total2
        ! A temporary set of Zombie states is created with the HF states with the correct number of electrons and spin
        call alloczs(cstoretemp,total2)
        
        allocate(magovrlp(total2),stat=ierr)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in magovrlp allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
    
        deallocate(combs,stat=ierr)
        if(ierr==0)deallocate(combsfix,stat=ierr)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in combination matrix deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if
    
        !$omp parallel private(j,k,magnitude) &
        !$omp shared(combs2, total2,total3,magovrlp,cstoretemp,zstore,clean_ndet) 
        
        !$omp do
        do j=1, total2
            call zomhfc(cstoretemp(j),combs2(j,:))
        end do
        !$omp end do
    
        
        ! This section reduces the size of cleaning matrix further by finding the overlap between each cleaning 
        ! state and the biased basis. Only cleanign states with overlaps higher than the set threshold are accepted
        ! and placed into the final cleaning set to produce the cleaning Hamiltonian.
        if(total2>99)then
            !$omp do
            do j=1, total2
                magnitude=0.0
                ! magnitude=(0.0,0.0)
                do k=1, ndet
                    magnitude=magnitude+product((cstore(j)%sin*zstore(k)%sin)+(cstore(j)%cos*zstore(k)%cos))    
                    ! overlap(cstoretemp(j),zstore(k))
                end do
                ! if(abs(REAL(magnitude))>0.0005) then
                if(abs(magnitude)>0.0005) then
                    !$omp critical
                    total3=total3+1
                    magovrlp(total3)=j
                    !$omp end critical
                end if
            end do
            !$omp end do

            !$omp master
            write(6,"(a,i0)") 'Combinations with corect spin and enough contribution ',total3
            clean_ndet=total3
            !$omp end master
            
        else
            !$omp master
            clean_ndet=total2
            !$omp end master
        end if
        !$omp end parallel
        
        call alloczs(cstore,clean_ndet)
        if(total2>99)then
            !$omp parallel do
            do j=1, total3
                cstore(j)=cstoretemp(magovrlp(j))
            end do
            !$omp end parallel do
        else 
            cstore=cstoretemp
        end if

        call dealloczs(cstoretemp)
    
        deallocate(magovrlp,stat=ierr)
        if(ierr==0)deallocate(combs2,stat=ierr)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in magovrlp deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        

        do j=1,clean_ndet
            call zombiewriter_c(cstore(j),j)
        end do
       
        
        call allocham(cleanham,clean_ndet,1)
        call hamgen(cleanham,cstore,elecs,ndet,an_cr,an2_cr2,1)
        ! call hamgen(cleanham,cstore,elecs,clean_ndet,1)
        call matrixwriter(cleanham%hjk,clean_ndet,"data/clean_ham.csv")
        
        
        return


    end subroutine clean_setup

    subroutine zomhfc(zoms,occ)

        implicit none
        type(zombiest),intent(inout)::zoms 
        integer, dimension(:), intent(in)::occ
        integer::j

        if (errorflag .ne. 0) return

        ! zoms%dead(1:norb)=(1.0d0,0.0d0)
        ! zoms%cos(1:norb)=(1.0d0,0.0d0)
        ! zoms%sin(1:norb)=(0.0d0,0.0d0)
        ! zoms%alive(1:norb)=(0.0d0,0.0d0)
        ! zoms%phi(1:norb)=0.0

        zoms%dead(1:norb)=1.0d0
        zoms%cos(1:norb)=1.0d0
        zoms%sin(1:norb)=0.0d0
        zoms%alive(1:norb)=0.0d0
        zoms%phi(1:norb)=0.0

        do j=1, size(occ)
            ! zoms%alive(occ(j))=(1.0d0,0.0d0)
            ! zoms%dead(occ(j))=(0.0d0,0.0d0)
            ! zoms%sin(occ(j))=(1.0d0,0.0d0)
            ! zoms%cos(occ(j))=(0.0d0,0.0d0)
            ! zoms%phi(occ(j))=0.5*pirl
            zoms%alive(occ(j))=1.0d0
            zoms%dead(occ(j))=0.0d0
            zoms%sin(occ(j))=1.0d0
            zoms%cos(occ(j))=0.0d0
            zoms%phi(occ(j))=0.5*pirl
        end do

        return

    end subroutine zomhfc

    subroutine cleaner(zstore,cleanzom,dvec, dvec_clean, cleantot,norm)

        implicit none
        type(zombiest),dimension(:),intent(in)::zstore,cleanzom
        type(dvector),intent(in)::dvec
        type(dvector),intent(inout):: dvec_clean
        integer,intent(in)::cleantot
        ! complex(kind=8),intent(out)::norm
        ! complex(kind=8)::ovrlp1, ovrlp2
        real(kind=8),intent(out)::norm
        real(kind=8)::ovrlp1, ovrlp2
        integer::j,k,l
        
        if (errorflag .ne. 0) return
        
        ! norm=(0.0d0,0.0d0)
        norm=0.0
      
        do j=1,cleantot
            do k=1, ndet
                ! ovrlp1=overlap(cleanzom(j),zstore(k))
                ovrlp1=product((cleanzom(j)%sin*zstore(k)%sin)+(cleanzom(j)%cos*zstore(k)%cos))
                dvec_clean%d(j)=dvec_clean%d(j)+(dvec%d(k)*ovrlp1)
                
                do l=1, ndet
                    ! ovrlp2=overlap(zstore(l),cleanzom(j))
                    ovrlp2=product((cleanzom(j)%sin*zstore(l)%sin)+(cleanzom(j)%cos*zstore(l)%cos))
                    ! norm = norm + (conjg(dvec%d(l))*dvec%d(k)*ovrlp2*ovrlp1)
                    norm = norm + (dvec%d(l))*dvec%d(k)*ovrlp2*ovrlp1
                end do
            end do
        end do

     

        return

    end subroutine cleaner

    subroutine clean_read(cstore,cleanham,clean_ndet,elecs,an_cr,an2_cr2)
        
        implicit none

        type(zombiest),dimension(:),allocatable,intent(inout)::cstore
        type(hamiltonian), intent(inout)::cleanham
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer, intent(inout):: clean_ndet
        integer::pyscfc
     
        if (errorflag .ne. 0) return
    
        pyscfc=1
   
        if(pyscfc.eq.1)then
            call pyscf_clean(cstore,clean_ndet,nel)
            call allocham(cleanham,clean_ndet,1)
            call hamgen(cleanham,cstore,elecs,ndet,an_cr,an2_cr2,1)
            ! call hamgen(cleanham,cstore,elecs,clean_ndet,1)
            call matrixwriter(cleanham%hjk,clean_ndet,"data/clean_ham.csv")
        else
            clean_ndet = lines_clean(clean_ndet)        
            call alloczs(cstore,clean_ndet)
            call read_zombie_c(cstore,clean_ndet)
            call allocham(cleanham,clean_ndet,1)
            call read_ham_c(cleanham,clean_ndet)
        end if

    
        return 

    end subroutine clean_read

    integer function lines_clean(nlines)
        implicit none

        integer, intent(INOUT):: nlines
        integer:: ierr
        
        if (errorflag .ne. 0) return

        ierr=0
        nlines=0
        open(unit=204, file='data/clean_ham.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a,i0)") 'Error in opening clean_ham.csv file',ierr
            errorflag = 1
            return
        end if

        do 
            read(204,*, iostat=ierr)
            if(ierr<0)then
                ! write(0,"(a,i0)") "nlines has value ", nlines
                lines_clean=nlines
                close(204)
                return
            else if (ierr/=0) then
                write(0,"(a,i0)") "Error in counting h1ea rows. ierr had value ", ierr
                errorflag=1
                return
            end if
            nlines=nlines+1
        end do

        return 

    end function lines_clean


    subroutine pyscf_clean(cstore, clean_ndet,nume)

        implicit none

        type(zombiest),dimension(:),allocatable,intent(inout)::cstore
        integer, intent(inout)::clean_ndet
        integer, intent(in)::nume
        integer, allocatable, dimension(:,:)::combs
        integer::Line2, Line3, Line4, Line5, Line6,Line7
        real::Line1
        integer:: ierr, nlines,j

        if (errorflag .ne. 0) return
        ierr=0

        open(unit=204, file='data/FCIconfigs_equilibrium.txt',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a,i0)") 'Error in opening pyscf cleaning file ',ierr
            errorflag = 1
            return
        end if

        nlines=0
        do 
            read(204,*, iostat=ierr)
            if(ierr<0)then
                ! write(0,"(a,i0)") "nlines has value ", nlines
                clean_ndet=nlines
                close(204)
                exit
            else if (ierr/=0) then
                write(0,"(a,i0)") "Error in counting pyscf cleaning file rows. ierr had value ", ierr
                errorflag=1
                return
            end if
            nlines=nlines+1
        end do
        ierr=0
        print*, clean_ndet
        call alloczs(cstore,clean_ndet)
        allocate(combs(clean_ndet,nume),stat=ierr)
        if (ierr.ne.0) then
            write(0,"(a,i0)") 'Error in combination matrix allocation',ierr
            errorflag = 1
            return
        end if
      
        open(unit=204, file='data/FCIconfigs_equilibrium.txt',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a,i0)") 'Error in opening pyscf cleaning file',ierr
            errorflag = 1
            return
        end if
    
        do j=1,clean_ndet
            read(204,*, iostat=ierr) Line1, Line2, Line3, Line4, Line5, Line6,Line7
            combs(j,1)=Line2
            combs(j,2)=Line3
            combs(j,3)=Line4
            combs(j,4)=Line5
            combs(j,5)=Line6
            combs(j,6)=Line7
        end do

        close(204)
     
        !$omp parallel shared(combs,cstore) private(j)
        !$omp do
        do j=1,clean_ndet
            combs(j,1:3)=(combs(j,1:3)*2)-1
            combs(j,4:6)=(combs(j,4:6)*2)
            ! print*, combs(j,:)
        end do
        !$omp end do

        !$omp do
        do j=1, clean_ndet
            call zomhfc(cstore(j),combs(j,:)) 
        end do
        !$omp end do
        !$omp end parallel
   
        return


    end subroutine

END MODULE clean

