MODULE zom 

    use globvars

    contains



    subroutine gen_ran_zs(zstore,num)

        implicit none

        type(zombiest),dimension(:),intent(inout)::zstore
        integer, intent(in)::num
        integer::j,k
        DOUBLE PRECISION, external::ZBQLU01
        real(kind=8)::dummy
        ! real(kind=8),dimension(norb)::rands
        ! complex(kind=8)::ctemp
        if (errorflag .ne. 0) return

        
        if(imagflg=='n') then
            do j=1,num
                ! call random_number(rands)
                do k=1,norb
                    dummy=-1
                    !$omp critical
                    do while((dummy.lt.0))
                       dummy=2*pirl*ZBQLU01(1)
                    !    if((dummy.gt.0.5*pirl))then
                    !         dummy=-1
                    !    end if
                    end do
                    !$omp end critical
                    zstore(j)%phi(k)=dummy
                end do 
                zstore(j)%sin=sin(cmplx(zstore(j)%phi,0.0d0,kind=8))
                zstore(j)%cos=cos(cmplx(zstore(j)%phi,0.0d0,kind=8))
            end do
            if(rhf_1=='y') then
                zstore(1)%phi(1:nel)=0.5*pirl
                zstore(1)%phi(nel+1:)=0
                zstore(1)%sin=0
                zstore(1)%cos=1
                zstore(1)%sin(1:nel)=cmplx(1,0.0d0,kind=8)
                zstore(1)%cos(1:nel)=cmplx(0,0.0d0,kind=8)
            end if 
        else if(imagflg=='y')then
            do j=1,num
                do k=1,norb
                    zstore(j)%phi(k)=2*pirl*ZBQLU01(1)
                    zstore(j)%img=2*pirl*ZBQLU01(1)
                end do
                zstore(j)%sin=sin(zstore(j)%phi*exp(i*zstore(j)%img))
                zstore(j)%cos=cos(cmplx(zstore(j)%phi,0.0d0,kind=8))
            end do
            if(rhf_1=='y') then
                zstore(1)%alive(1:nel)=(1.0d0,0.0d0)
                zstore(1)%dead(1:nel)=(0.0d0,0.0d0)
                zstore(1)%alive((nel+1):norb)=(0.0d0,0.0d0)
                zstore(1)%dead((nel+1):norb)=(1.0d0,0.0d0)
            end if 
        end if

        return

    end subroutine gen_ran_zs

    subroutine gen_hf_zs(zstore)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        integer::j,k,total,ierr,count
        integer, allocatable, dimension(:,:)::combs
        integer, dimension(ndet,norb)::combs2

        if (errorflag .ne. 0) return
        ierr=0

        combs2(1:ndet,1:norb)=0
        
        zstore(1)%dead(1:norb)=(1.0d0,0.0d0)
        count=2
        !$omp parallel shared(count,zstore,combs2,errorflag) private(combs,j,k,total,ierr)
        !$omp do
        do j=1,norb
            if(errorflag==1) then
                cycle
            end if
            total=choose(norb,j)
            allocate(combs(total,j),stat=ierr)
            if(ierr/=0) then
                write(0,"(a,i0)") "Error in combination matrix allocation. ierr had value ", ierr
                errorflag=1
                cycle
            end if
            call combinations(norb,j,combs,total)
            
            !$omp critical 
            do k=1, total
                combs2(count,1:j)=combs(k,:)
                count=count+1
            end do
            !$omp end critical
          
            deallocate(combs,stat=ierr)
            if(ierr/=0) then
                write(0,"(a,i0)") "Error in combination matrix deallocation. ierr had value ", ierr
                errorflag=1
                cycle
            end if
        end do
        !$omp end do
        !$omp do
        do j=1,ndet
            if(errorflag==1) then
                cycle
            end if
            call zomhf(zstore(j),combs2(j,:))
        end do
        !$omp end do
        !$omp end parallel


        return  
    
    end subroutine gen_hf_zs

    function choose(n, k)

    ! code adapted from https://rosettacode.org/wiki/Combinations#Fortran accessed 11:45 18/07/2022
       
        implicit none

        integer :: choose
        integer, intent(in) :: n, k
        integer :: jmax, j, jmin
     
    
        if ( (n < 0 ) .or. (k < 0 ) ) then
           write(0, *) "negative in choose"
           choose = 0
        else
           if ( n < k ) then
              choose = 0
           else if ( n == k ) then
              choose = 1
           else
              jmax = max(k, n-k)
              jmin = min(k, n-k)
              choose = 1
              do j = jmax+1, n
                 choose = choose * j
              end do
              do j = 2, jmin
                 choose = choose / j
              end do
           end if
        end if
        
    end function choose

    subroutine combinations(n_max,m_max,final,tot)

     ! code adapted from https://rosettacode.org/wiki/Combinations#Fortran accessed 11:45 18/07/2022
        
        implicit none

        type comb_result
            integer, dimension(:), allocatable :: combs
        end type comb_result

        integer, intent(in)::n_max,m_max,tot
        integer,dimension(:,:),intent(inout)::final
        type(comb_result), dimension(:), pointer :: co
        integer :: j, k, s, jx, kx, t,ierr
    
        allocate(co(0:tot-1),stat=ierr)
        do j = 0, tot-1
            allocate(co(j)%combs(0:m_max-1))
         end do
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in co matrix allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        do j = 0, tot-1
            jx = j; kx = m_max
            do s = 0, n_max-1
               if ( kx == 0 ) exit
               t = choose(n_max-(s+1), kx-1)
               if ( jx < t ) then
                  co(j)%combs(kx-1) = s
                  kx = kx - 1
               else
                  jx = jx - t
               end if
            end do
        end do

        do j=0, tot-1
            do k = m_max-1,0,-1
                final(j+1,k+1)=(co(j)%combs(k)+1)
            end do
            deallocate(co(j)%combs)
            if(ierr/=0) then
                write(0,"(a,i0)") "Error in combs matrix deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if 
        end do
        deallocate(co)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in co matrix deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if 

        return
        
    end subroutine combinations


    subroutine zomhf(zom,occ)

        implicit none
        type(zombiest),intent(inout)::zom 
        integer, dimension(:), intent(in)::occ
        integer::j

        if (errorflag .ne. 0) return

        zom%dead(1:norb)=(1.0d0,0.0d0)
        zom%alive(1:norb)=(0.0d0,0.0d0)
        zom%cos(1:norb)=(1.0d0,0.0d0)
        zom%sin(1:norb)=(0.0d0,0.0d0)
        zom%phi(1:norb)=0.0

        do j=1, norb
            if(occ(j)==0)then
                return
            end if
            zom%sin(occ(j))=(1.0d0,0.0d0)
            zom%cos(occ(j))=(0.0d0,0.0d0)
            zom%phi(occ(j))=0.5*pirl
            zom%alive(occ(j))=(1.0d0,0.0d0)
            zom%dead(occ(j))=(0.0d0,0.0d0)
        end do

        return

    end subroutine zomhf

    subroutine gen_biased_zs(zstore)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        DOUBLE PRECISION, external::ZBQLNOR,ZBQLUAB,ZBQLU01
        real(kind=8)::mu((norb/2)),sig(norb/2)
        real(kind=8)::val
        integer::j,k

        if (errorflag .ne. 0) return
        ! mu and sigma values for Li2 10 spin orbitals used in OG paper
        mu=(/0.25,0.25,0.25,0.0,0.0/)
        sig=(/0.0001,0.0001,0.175,0.351,0.120/)

        ! mu and sigma values for BH with 38 spin orbtials
        ! mu=(/0.25,0.213632469,0.193380738,0.001262455,0.000505343,0.00062495,0.000530594,9.57371E-06,0.000169358,3.27753E-05, &
        ! 0.004644281,0.000396432,0.000387224,5.15685E-05,0.004644276,0.000396434,0.000387213,5.16551E-05,9.58165E-06/)
        ! sig(1:norb/2)=0.15
        ! mu(1:(nel/2))=0.25
        ! call musig(mu,sig)
        
        if(imagflg=='n') then
            !$omp parallel shared(zstore) private(j,k)
            !$omp do
            do j=1, ndet
                !$omp critical
                do k=1,norb/2
                    val=-1
                    do while(val.lt.0)
                        ! val=2*pirl*sig(k)*ZBQLU01(1)    ! 
                         val=2*pirl*ZBQLNOR(mu(k),sig(k))
                        ! if((val.gt.0.5*pirl))then
                        !     val=-1
                        ! end if
                    end do
                    zstore(j)%phi(2*k-1)=val
                    val=-1
                    do while(val.lt.0)
                        ! val=2*pirl*sig(k)*ZBQLU01(1)
                        val=2*pirl*ZBQLNOR(mu(k),sig(k))
                    
                        ! if((val.gt.0.5*pirl))then
                        !     val=-1
                        ! end if
                    end do
                    zstore(j)%phi(2*k)=val
                    ! print*,val(2*k-1),val(2*k)
                    ! val(2*k-1)=2*pirl*mu(k)*exp(-ZBQLUAB(0,0.1))
                    ! val(2*k)=2*pirl*mu(k)*exp(-ZBQLUAB(0,0.1))
                end do
                !$omp end critical
                zstore(j)%sin=sin(cmplx(zstore(j)%phi,0.0d0,kind=8))
                zstore(j)%cos=cos(cmplx(zstore(j)%phi,0.0d0,kind=8))
            end do
            !$omp end do
            !$omp end parallel
            if(rhf_1=='y') then
                zstore(1)%phi(1:nel)=0.5*pirl
                zstore(1)%phi(nel+1:)=0
                zstore(1)%sin=0
                zstore(1)%cos=1
                zstore(1)%sin(1:nel)=cmplx(1,0.0d0,kind=8)
                zstore(1)%cos(1:nel)=cmplx(0,0.0d0,kind=8)
            end if 
            ! do j=1, ndet
            !     zstore(j)%sin(:)=zstore(j)%sin(:)*100
            !     ! zstore(j)%cos(:)=zstore(j)%cos(:)*100
            ! end do
        else if(imagflg=='y')then
            print*,"not yet written"
        end if
        return

    end subroutine gen_biased_zs

    subroutine musig(mu,sig)

        implicit none
        real(kind=8),dimension(:),intent(inout)::mu,sig
        integer::alive,j
        real(kind=8)::asrt,aend,dsrt


        alive=int(nel/2)
        mu(1:alive)=0.25
        mu(alive+1:)=0

        asrt=0.0   !0.0001
        ! asrt=0
        aend=0.90     !  0.200
        dsrt=0.90/((norb/2)-1)     !      0.08
        !dend= 0.0006

        do j=1,(norb/2)
            sig(j)=exp(-0.7*asrt)-0.5
            asrt=asrt+dsrt
        end do
        
        ! do j=0, (alive-1)
        !     sig(j+1)=((aend-asrt)/(alive-1)*j)+asrt
        ! end do

        ! do j=0, (((norb/2)-alive)-1)
        !     sig(j+1+alive)=((dend-dsrt)/(((norb/2)-alive)-1)*j)+dsrt
        ! end do

        return

    end subroutine

    subroutine genzf(zstore,num)
        
        implicit none
        type(zombiest), dimension(:), intent(inout)::zstore
        integer, intent(in)::num

        if (errorflag .ne. 0) return
        
        select case(zst)
            case('HF')
                call gen_hf_zs(zstore)
            case('RN')
                call gen_ran_zs(zstore,num)
            case('BB')
                call gen_biased_zs(zstore)
            case default
                write(0,"(a)") "Error! Initial zombie type method not recognised!"
                write(0,"(a)") "This should have been caught at the read stage!"
                errorflag = 1
                return
        end select 
        return
    
    end subroutine genzf

END MODULE zom