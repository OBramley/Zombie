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
                    zstore(j)%phi(k)=0.5*pirl*(ZBQLU01()) 
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
            !         zstore(j)%phi(k)=2*pirl*r   !ZBQLU01()
            !         call random_number(r)
            !         zstore(j)%img=2*pirl*r  !ZBQLU01()
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
        ! do j=1,norb
        do j=nel, nel 
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
            combs2(1:ndet,1:nel)=combs(1:ndet,:)
            ! !$omp critical 
            ! do k=1, total
            !     combs2(count,1:j)=combs(k,:)
            !     count=count+1
            ! end do
            ! !$omp end critical
          
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
        zom%val(1:norb)=0.0d0 !1.0d-16
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
        integer::k,mult,a1,a2,a3,a4,a5,a6,a7,a8
        real(wp)::step

       
        z1%phi=0
  
        if(nel.gt.10)then
            a1=5 !1s2s
            a2=11 !2p
            a3=13 !3s
            a4=19 !3p
            a5=29 !3d
            a6=31 !4s
            a7=37 !4p  
            a8=47 !4d
        else if(nel.gt.4)then
            a1=3 !1s
            a2=5 !2s
            a3=11 !2p
            a4=13 !3s
            a5=19 !3p
            a6=29 !3d  
            a7=31 !4s
            a8=37 !4p
        else if(nel.gt.2)then
            a1=2 !1s1
            a2=5 !1s22s
            a3=11 !2p
            a4=13 !3s
            a5=19 !3p
            a6=29 !3d  
            a7=31 !4s
            a8=37 !4p
        else
            a1=0
            a2=3 !1s
            a3=5 !2s
            a4=11 !2p
            a5=13 !3s
            a6=19 !3p
            a7=29 !3d
            a8=37 !4s
    
        end if
       
        do k=1,norb
            if(k .lt.a1)then 
                z1%phi(k)=0.5*pirl
            else if( k .lt. a2 )then 
                z1%phi(k)= value_maker(1) 
            else if( k .lt. a3 )then
                z1%phi(k)=value_maker(2)
            else if( k .lt. a4 )then 
                z1%phi(k)=value_maker(3) 
            else if( k .lt. a5 )then 
                ! z1%phi(k)=value_maker(4) 
                z1%phi(k)=1.0d-10*ZBQLU01()
            else if( k .lt. a6 )then 
                ! z1%phi(k)=value_maker(5) 
                z1%phi(k)=1.0d-12*ZBQLU01()
            else if( k .lt. a7 )then 
                ! z1%phi(k)=value_maker(6)
                z1%phi(k)=0.0d0 !1.0d-13*ZBQLU01()
            else if( k .lt. a8 )then 
                ! z1%phi(k)=value_maker(7)
                z1%phi(k)=1.0d-11*ZBQLU01()
            end if 
        end do 

        return 

    end subroutine biased_func

    function value_maker(num) result(val)
        implicit none
        integer::num
        real(wp)::val,val2
    
        val=0
        do while(val == 0) 
            val2=abs((num-ZBQLU01())/(nel*ZBQLU01())*nel)!*ZBQLU01()*2))
            val=0.5*pirl*exp(-val2)
            if((is_nan(val).eqv..true.).or.(is_inf(val).eqv..true.))then
                val=0
            end if
        end do
    end function value_maker

    subroutine gen_biased_zs(zstore)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        real(wp)::mu((norb/2)),sig((norb/2))
        ! real(wp)::val
        integer::j
       

        if (errorflag .ne. 0) return
 
        if(imagflg=='n') then
            do j=1, ndet
                call biased_func(zstore(j))
                call val_set(zstore(j))
                
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


END MODULE zom