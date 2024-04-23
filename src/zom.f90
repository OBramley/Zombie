MODULE zom 

    use mod_types
    use randgen
    use globvars
    use infnan_mod
    contains



    subroutine gen_ran_zs(zstore,num)

        implicit none

        type(zombiest),dimension(:),intent(inout)::zstore
        integer, intent(in)::num
        integer::j,k

        if (errorflag .ne. 0) return

        if(imagflg=='n') then
            do j=1,num
                do k=1,norb
                    zstore(j)%phi(k)=2*pirl*(ZBQLU01(1)) 
                end do
                call val_set(zstore(j))
            end do
            if(rhf_1=='y') then

                zstore(1)%phi(1:nel)=0.5*pirl
                zstore(1)%phi(nel+1:)=0
                zstore(1)%val(1:norb)=0 
                zstore(1)%val(norb+1:)=1 
                zstore(1)%val(1:nel)=1
                zstore(1)%val(norb+1:norb+nel)=0
                ! if((gramflg.eq.'y').and.(GDflg.eq.'y'))then
                    ! zstore(1)%phi(nel)=0
                    ! zstore(1)%val(nel)=0
                    ! zstore(1)%val(nel+norb)=1
                !     zstore(1)%phi(nel+1)=0.5*pirl
                !     zstore(1)%val(nel+1)=1
                !     zstore(1)%val(nel+1+norb)=0
                ! end if 
            end if 
            
        else if(imagflg=='y')then
            ! do j=1,num
            !     do k=1,norb
            !         call random_number(r)
            !         zstore(j)%phi(k)=2*pirl*r   !ZBQLU01(1)
            !         call random_number(r)
            !         zstore(j)%img=2*pirl*r  !ZBQLU01(1)
            !     end do
            !     zstore(j)%sin=sin(zstore(j)%phi*exp(i*zstore(j)%img))
            !     zstore(j)%cos=cos(cmplx(zstore(j)%phi,0.0d0,kind=8))
            ! end do
            ! if(rhf_1=='y') then
            !     zstore(1)%alive(1:nel)=(1.0d0,0.0d0)
            !     zstore(1)%dead(1:nel)=(0.0d0,0.0d0)
            !     zstore(1)%alive((nel+1):norb)=(0.0d0,0.0d0)
            !     zstore(1)%dead((nel+1):norb)=(1.0d0,0.0d0)
            ! end if 
        end if

        return

    end subroutine gen_ran_zs

    subroutine gen_hf_zs(zstore)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        integer::count,j
        integer::total,k
        integer, allocatable, dimension(:,:)::combs
        integer, dimension(ndet,norb)::combs2
        integer::ierr=0

        if (errorflag .ne. 0) return
      

        combs2(1:ndet,1:norb)=0
    
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
                write(stderr,"(a,i0)") "Error in combination matrix allocation. ierr had value ", ierr
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
                write(stderr,"(a,i0)") "Error in combination matrix deallocation. ierr had value ", ierr
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
        integer:: jmax, j, jmin
     
    
        if ( (n < 0 ) .or. (k < 0 ) ) then
           write(stderr, *) "negative in choose"
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

        integer, intent(in)::tot
        integer, intent(in)::n_max,m_max
        integer,dimension(:,:),intent(inout)::final
        type(comb_result), dimension(:), pointer :: co
        integer::j, k, jx, t
        integer :: ierr,s,kx

     
    
        allocate(co(0:tot-1),stat=ierr)
        do j = 0, tot-1
            allocate(co(j)%combs(0:m_max-1))
         end do
        if(ierr/=0) then
            write(stderr,"(a,i0)") "Error in co matrix allocation. ierr had value ", ierr
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
                write(stderr,"(a,i0)") "Error in combs matrix deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if 
        end do
        
        deallocate(co)
        if(ierr/=0) then
            write(stderr,"(a,i0)") "Error in co matrix deallocation. ierr had value ", ierr
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


        zom%val(1+norb:2*norb)=1.0d0
        zom%val(1:norb)=0.0d0
        zom%phi(1:norb)=0

        do j=1, norb
            if(occ(j)==0)then
                return
            end if
            
            zom%val(occ(j))=1.0d0
            zom%val(norb+occ(j))=0.0d0
            zom%phi(occ(j))=0.5*pirl
  
        end do

        return

    end subroutine zomhf

    subroutine biased_func(z1)
        implicit none
        type(zombiest),intent(inout)::z1
        integer::k
        z1%phi=0.0001
        if(nel.gt.4)then
            z1%phi(1:4)=0.5*pirl
            if(modulo(nel,2)==0) then
                do k=5,nel+4
                    z1%phi(k)=0.5*pirl*ZBQLU01(1)
                end do
                ! do k=nel+5,norb
                !     z1%phi(k)=z1%phi(k)*ZBQLU01(1)
                ! end do
            else
                do k=5,nel+5
                    z1%phi(k)=0.5*pirl*ZBQLU01(1)
                end do
                ! do k=nel+6,norb
                !     z1%phi(k)=z1%phi(k)*ZBQLU01(1)
                ! end do
            end if
        else if(nel.eq.4)then
            z1%phi(1:2)=0.5*pirl
            if(modulo(nel,2)==0) then
                do k=3,nel+4
                    z1%phi(k)=0.5*pirl*ZBQLU01(1)
                end do
                ! do k=nel+5,norb
                !     z1%phi(k)=z1%phi(k)*ZBQLU01(1)
                ! end do
            else
                do k=3,nel+5
                    z1%phi(k)=0.5*pirl*ZBQLU01(1)
                end do
                ! do k=nel+6,norb
                !     z1%phi(k)=z1%phi(k)*ZBQLU01(1)
                ! end do
            end if
        else 
            if(modulo(nel,2)==0) then
                do k=1,nel+4
                    z1%phi(k)=0.5*pirl*ZBQLU01(1)
                end do
                ! do k=nel+5,norb
                !     z1%phi(k)=z1%phi(k)*ZBQLU01(1)
                ! end do
            else
                do k=1,nel+5
                    z1%phi(k)=0.5*pirl*ZBQLU01(1)
                end do
                ! do k=nel+6,norb
                !     z1%phi(k)=z1%phi(k)*ZBQLU01(1)
                ! end do
            end if
        end if 
        return 

    end subroutine biased_func

    subroutine gen_biased_zs(zstore)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        ! real(wp)::mu((norb/2)),sig(norb/2)
        ! real(wp)::val
        integer::j,k
       

        if (errorflag .ne. 0) return
 
        ! call musig(mu,sig)
        if(imagflg=='n') then
            do j=1, ndet
                call biased_func(zstore(j))

                ! !$omp critical
                ! do k=1,norb/2
                !     val=-1
                !     do while(val.lt.0)
                !         val=2*pirl*random_normal(mu(k),sig(k)) 
                !     end do
                !     if((is_nan(val).eqv..true.))then
                !         val=2*pirl*(ZBQLU01(1))
                !     end if
                    
                !     zstore(j)%phi(2*k-1)=val
                !     val=-1
                !     do while(val.lt.0)
                !         val=2*pirl*random_normal(mu(k),sig(k))
                !     end do
                !     if((is_nan(val).eqv..true.))then
                !         val=2*pirl*(ZBQLU01(1)) 
                !     end if 
                !     zstore(j)%phi(2*k)=val
                ! end do
                ! !$omp end critical
               
                call val_set(zstore(j))
                ! zstore(j)%val(1:norb) = sin(zstore(j)%phi)
                ! zstore(j)%val(norb+1:2*norb)=cos(zstore(j)%phi)
             
            end do
            if(rhf_1=='y') then
                zstore(1)%phi(1:nel)=0.5*pirl
                zstore(1)%phi(nel+1:)=0
                zstore(1)%val(1:norb)=0 
                zstore(1)%val(norb+1:)=1 
                zstore(1)%val(1:nel)=1
                zstore(1)%val(norb+1:norb+nel)=0
            end if 
        else if(imagflg=='y')then
            print*,"not yet written"
        end if

       

        return

    end subroutine gen_biased_zs

    subroutine musig(mu,sig)

        implicit none
        real(wp),dimension(:),intent(inout)::mu,sig
        integer::alive,j
        real(wp)::asrt,aend,dsrt,dend,val


        alive=int(nel/2)
        mu(1:alive)=0.25
        mu(alive+1:)=0

        val=(norb/10)

        asrt=0.0001/ceiling(val)
        aend=0.17/ceiling(val)   
        dsrt=0.35/ceiling(val)  
        dend= 0.15/ceiling(val)

       
        
        do j=0, (alive-1)
            sig(j+1)=((aend-asrt)/(alive-1)*j)+asrt
        end do

        do j=0, (((norb/2)-alive)-1)
            ! mu(j+1+alive)=((dend-dsrt)/(((norb/2)-alive)-1)*j)+dsrt
            sig(j+1+alive)=((dend-dsrt)/(((norb/2)-alive)-1)*j)+dsrt
        end do
 
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
                write(stderr,"(a)") "Error! Initial zombie type method not recognised!"
                write(stderr,"(a)") "This should have been caught at the read stage!"
                errorflag = 1
                return
        end select 
        return
    
    end subroutine genzf


    ! FUNCTION random_normal(MU,SIGMA)
    !     !
    !     !       Returns a random number Normally distributed with mean
    !     !       MU and standard deviation |SIGMA|, using the Box-Muller
    !     !       algorithm
    !     !
    !     DOUBLE PRECISION THETA,R,ZBQLNOR,random_normal,PI,MU,SIGMA
    !     DOUBLE PRECISION SPARE
    !     INTEGER STATUS
    !     SAVE STATUS,SPARE,PI
    !     DOUBLE PRECISION p
    !     DATA STATUS /-1/
        
    !     IF (STATUS.EQ.-1) PI = 4.0D0*DATAN(1.0D0)

    !     IF (STATUS.LE.0) THEN
    !     call random_number(p)
    !     THETA = 2.0D0*PI*p 
    !     call random_number(p)

    !     R = DSQRT( -2.0D0*DLOG(p) )

    !     ZBQLNOR = (R*DCOS(THETA))
    !     SPARE = (R*DSIN(THETA))
    !     STATUS = 1
    !     ELSE
    !     ZBQLNOR = SPARE
    !     STATUS = 0
    !     ENDIF

    !     ZBQLNOR = MU + (SIGMA*ZBQLNOR)
    !     random_normal=ZBQLNOR

    ! END

END MODULE zom