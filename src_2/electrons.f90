MODULE electrons

use globvars
use alarrays

contains

! Program to generate one and two electron integrals from PyScf
! Some of this code has been adapted from George Booth
subroutine electronintegrals(elecs,an_cr,an2_cr2)

    implicit none

    type(elecintrgl), intent(inout)::elecs
    type(oprts),intent(inout)::an_cr,an2_cr2
    real::e1in,e2in
    integer:: ierr,e1,e2
    

    if (errorflag .ne. 0) return


    ierr = 0
    open(unit=129, file='integrals/h1e.csv',status='old',iostat=ierr)
    if (ierr.ne.0) then
        write(0,"(a)") 'Error in opening h1ec.csv file'
        errorflag = 1
        return
    end if
    read(129,*, iostat=ierr)e1in
    if (ierr .ne. 0) then
        write(0,"(a)") 'Error reading h1e'
        errorflag = 1
        return
      end if
    close(129)

    e1=int(e1in)

    open(unit=131, file='integrals/h2e.csv',status='old',iostat=ierr)
    if (ierr.ne.0) then
        write(0,"(a)") 'Error in opening h2e.csv file'
        errorflag = 1
        return
    end if
    read(131,*, iostat=ierr)e2in
    if (ierr .ne. 0) then
        write(0,"(a)") 'Error reading h1e'
        errorflag = 1
        return
      end if
    close(131)

    e2=int(e2in)

    call allocintgrl(elecs,e1,e2)
    
    call  one_electrons(elecs,an_cr,e1)
    call  two_electrons(elecs,an2_cr2,e2)

    

    ! nlines = lines(nlines)
    ! call spattospin1(elecs,nlines)
    ! call spattospin2(elecs,nlines)

    open(unit=128, file='integrals/hnuc.csv',status='old',iostat=ierr)
    if (ierr.ne.0) then
        write(0,"(a)") 'Error in opening hnuc.csv file'
        errorflag = 1
        return
    end if
    read(128,*, iostat=ierr)elecs%hnuc
    if (ierr .ne. 0) then
        write(0,"(a)") 'Error reading Hnuc'
        errorflag = 1
        return
      end if
    close(128)
    
    write(6,"(a)") "1 & 2 electron integrals successfully generated"

    return

end subroutine electronintegrals

! 
subroutine two_electrons(elecs,an2_cr2,e2)

    implicit none 
    type(elecintrgl), intent(inout)::elecs
    type(oprts),intent(inout)::an2_cr2
    integer,intent(in)::e2
    real(kind=8),dimension(e2+1,5)::read_in
    integer::ierr,l,k,an1,an2,cr1,cr2,j

    
    open(unit=131, file='integrals/h2e.csv',status='old',iostat=ierr)
    if (ierr.ne.0) then
        write(0,"(a)") 'Error in opening hnuc.csv file'
        errorflag = 1
        return
    end if
   
    do j=1,e2+1
        read(131,*) (read_in(j,k),k=1,5)
    end do 
    

    close(131)
   
    elecs%h2ei=read_in(2:,1)

    call alloc_oprts(an2_cr2,e2)

   
    do l=1,e2
        ! cr1=int(read_in(l+1,2)) !3    1:cr1,2:cr2,3:an1,4:an2
        ! cr2=int(read_in(l+1,3)) !4
        ! an1=int(read_in(l+1,4)) !2
        ! an2=int(read_in(l+1,5)) !1
        cr1=int(read_in(l+1,2)) !3
        cr2=int(read_in(l+1,3)) !4
        an1=int(read_in(l+1,4)) 
        an2=int(read_in(l+1,5)) 
        !annihilation
        an2_cr2%dead(an1,l)=an2_cr2%alive(an1,l)
        an2_cr2%neg_dead(an1,l)=an2_cr2%neg_alive(an1,l)
        an2_cr2%alive(an1,l)=0
        an2_cr2%neg_alive(:an1-1,l)=int(an2_cr2%neg_alive(:an1-1,l)*(-1),kind=1)
        !annihilation
        an2_cr2%dead(an2,l)=an2_cr2%alive(an2,l)
        an2_cr2%neg_dead(an2,l)=an2_cr2%neg_alive(an2,l)
        an2_cr2%alive(an2,l)=0
        an2_cr2%neg_alive(:an2-1,l)=int(an2_cr2%neg_alive(:an2-1,l)*(-1),kind=1) 
        !creation 
        an2_cr2%alive(cr1,l)=an2_cr2%dead(cr1,l)
        an2_cr2%neg_alive(cr1,l)=an2_cr2%neg_dead(cr1,l)
        an2_cr2%dead(cr1,l)=0
        an2_cr2%neg_alive(:cr1-1,l)=int(an2_cr2%neg_alive(:cr1-1,l)*(-1),kind=1)
        !creation 
        an2_cr2%alive(cr2,l)=an2_cr2%dead(cr2,l)
        an2_cr2%neg_alive(cr2,l)=an2_cr2%neg_dead(cr2,l)
        an2_cr2%dead(cr2,l)=0
        an2_cr2%neg_alive(:cr2-1,l)=int(an2_cr2%neg_alive(:cr2-1,l)*(-1),kind=1)
        
        
        
    end do 

    if(GDflg.eq.'y')then
        do l=1,e2
            cr2=int(read_in(l+1,2))
            cr1=int(read_in(l+1,3))
            an2=int(read_in(l+1,4))
            an1=int(read_in(l+1,5))
            ! print*,'*************************'
            ! print*,an2,an1,cr1,cr2
            do j=1,norb 
                if((an2.eq.j).or.(cr2.eq.j).or.(an1.eq.j).or.(cr1.eq.j))then
                    an2_cr2%dcnt(0,j)=int(an2_cr2%dcnt(0,j)+1,kind=2)
                    an2_cr2%dcnt(an2_cr2%dcnt(0,j),j)=int(l,kind=2)
                    if((an2.eq.j).or.(an1.eq.j))then
                        if((cr1.eq.j).or.(cr2.eq.j))then
                            ! print*,j,' equals ', cr2,' or ',cr1, ' and one equals ', an2, an1
                            an2_cr2%alive_diff(j,:,l)=an2_cr2%alive(:,l)
                            an2_cr2%dead_diff(j,:,l)=an2_cr2%dead(:,l)
                            an2_cr2%neg_alive_diff(j,:,l)=an2_cr2%neg_alive(:,l)
                            an2_cr2%neg_dead_diff(j,:,l)=an2_cr2%neg_dead(:,l)
                            an2_cr2%alive_diff(j,j,l)=int(j+norb,kind=2)
                            an2_cr2%neg_alive_diff(j,j,l)=int(2*an2_cr2%neg_alive_diff(j,j,l),kind=1)

                        else 
                            ! print*,j,' equals ', an2,' or ',an1, ' but not ', cr2, cr1
                            an2_cr2%alive_diff(j,:,l)=an2_cr2%alive(:,l)
                            an2_cr2%dead_diff(j,:,l)=an2_cr2%dead(:,l)
                            an2_cr2%neg_alive_diff(j,:,l)=an2_cr2%neg_alive(:,l)
                            an2_cr2%neg_dead_diff(j,:,l)=an2_cr2%neg_dead(:,l)
                            an2_cr2%alive_diff(j,j,l)=int(j,kind=2)
                            an2_cr2%dead_diff(j,j,l)=int(j+norb,kind=2)
                            an2_cr2%neg_alive_diff(j,j,l)=int(an2_cr2%neg_dead_diff(j,j,l)*(-1),kind=1)
                        end if
                    else
                        an2_cr2%dcnt(0,j)=int(an2_cr2%dcnt(0,j)+1,kind=2)
                        ! print*,j,' equals ', cr2,' or ',cr1, ' but not ', an2, an1
                        an2_cr2%alive_diff(j,:,l)=an2_cr2%alive(:,l)
                        an2_cr2%dead_diff(j,:,l)=an2_cr2%dead(:,l)
                        an2_cr2%neg_alive_diff(j,:,l)=an2_cr2%neg_alive(:,l)
                        an2_cr2%neg_dead_diff(j,:,l)=an2_cr2%neg_dead(:,l)
                        an2_cr2%alive_diff(j,j,l)=int(j,kind=2)
                        an2_cr2%dead_diff(j,j,l)=int(j+norb,kind=2)
                        an2_cr2%neg_dead_diff(j,j,l)=an2_cr2%neg_alive_diff(j,j,l)
                        an2_cr2%neg_alive_diff(j,j,l)=int(an2_cr2%neg_alive_diff(j,j,l)*(-1),kind=1)
                    end if
                end if 
                
            end do
        end do
      
    end if

    return

end subroutine two_electrons


subroutine one_electrons(elecs,an_cr,e1)

    implicit none 
    type(elecintrgl), intent(inout)::elecs
    type(oprts),intent(inout)::an_cr
    integer,intent(in)::e1
    real(kind=8),dimension(e1+1,3)::read_in
    integer::ierr,l,k,an,cr,j

    
    open(unit=129, file='integrals/h1e.csv',status='old',iostat=ierr)
    if (ierr.ne.0) then
        write(0,"(a)") 'Error in opening hnuc.csv file'
        errorflag = 1
        return
    end if

    do j=1,e1+1
        read(129,*) (read_in(j,k),k=1,3)
    end do 
    close(129)

    
    elecs%h1ei=read_in(2:,1)
    
    
    call alloc_oprts(an_cr,e1)

    do l=1,e1
        an=int(read_in(l+1,2))
        cr=int(read_in(l+1,3))
        !annihilation
        an_cr%dead(an,l)=an_cr%alive(an,l)
        an_cr%neg_dead(an,l)=an_cr%neg_alive(an,l)
        an_cr%alive(an,l)=0
        an_cr%neg_alive(:an-1,l)=int(an_cr%neg_alive(:an-1,l)*(-1),kind=1)
        !creation 
        an_cr%alive(cr,l)=an_cr%dead(cr,l)
        an_cr%neg_alive(cr,l)=an_cr%neg_dead(cr,l)
        an_cr%dead(cr,l)=0
        an_cr%neg_alive(:cr-1,l)=int(an_cr%neg_alive(:cr-1,l)*(-1),kind=1)
    end do 
    if(GDflg.eq.'y')then
        do l=1,e1
            an=int(read_in(l+1,2))
            cr=int(read_in(l+1,3))
            do j=1,norb 
                if((an.eq.j).or.(cr.eq.j))then
                    an_cr%dcnt(0,j)=int(an_cr%dcnt(0,j)+1,kind=2)
                    an_cr%dcnt(an_cr%dcnt(0,j),j)=int(l,kind=2)
                    if(an.eq.cr)then 
                        an_cr%alive_diff(j,:,l)=an_cr%alive(:,l)
                        an_cr%dead_diff(j,:,l)=an_cr%dead(:,l)
                        an_cr%neg_alive_diff(j,:,l)=an_cr%neg_alive(:,l)
                        an_cr%neg_dead_diff(j,:,l)=an_cr%neg_dead(:,l)
                        an_cr%alive_diff(j,j,l)=int(j+norb,kind=2)
                        an_cr%neg_alive_diff(j,j,l)=int(2*an_cr%neg_alive_diff(j,j,l),kind=1)
                    else if(j.eq.an)then
                        an_cr%alive_diff(j,:,l)=an_cr%alive(:,l)
                        an_cr%dead_diff(j,:,l)=an_cr%dead(:,l)
                        an_cr%neg_alive_diff(j,:,l)=an_cr%neg_alive(:,l)
                        an_cr%neg_dead_diff(j,:,l)=an_cr%neg_dead(:,l)
                        an_cr%alive_diff(j,j,l)=int(j,kind=2)
                        an_cr%dead_diff(j,j,l)=int(j+norb,kind=2)
                        an_cr%neg_alive_diff(j,j,l)=int(an_cr%neg_dead_diff(j,j,l)*(-1),kind=1)
                    else if(j.eq.cr)then
                        an_cr%alive_diff(j,:,l)=an_cr%alive(:,l)
                        an_cr%dead_diff(j,:,l)=an_cr%dead(:,l)
                        an_cr%neg_alive_diff(j,:,l)=an_cr%neg_alive(:,l)
                        an_cr%neg_dead_diff(j,:,l)=an_cr%neg_dead(:,l)
                        an_cr%alive_diff(j,j,l)=int(j,kind=2)
                        an_cr%dead_diff(j,j,l)=int(j+norb,kind=2)
                        an_cr%neg_dead_diff(j,j,l)=an_cr%neg_alive_diff(j,j,l)
                        an_cr%neg_alive_diff(j,j,l)=int(an_cr%neg_alive_diff(j,j,l)*(-1),kind=1)
                    end if
                end if 
            end do
        end do 
    
    end if
    
    
    return


end subroutine one_electrons


!integer function lines(nlines)
!     implicit none

!     integer, intent(INOUT):: nlines
!     integer:: ierr
    
!     ierr=0
!     nlines=0
!     open(unit=129, file='integrals/h1ea.csv',status='old',iostat=ierr)
!     if (ierr.ne.0) then
!         write(0,"(a,i0)") 'Error in opening h1ea.csv file',ierr
!         errorflag = 1
!         return
!     end if

!     do 
!         read(129,*, iostat=ierr)
!         if(ierr<0)then
!             ! write(0,"(a,i0)") "nlines has value ", nlines
!             lines=nlines
!             close(129)
!             return
!         else if (ierr/=0) then
!             write(0,"(a,i0)") "Error in counting h1ea rows. ierr had value ", ierr
!             errorflag=1
!             return
!         end if
!         nlines=nlines+1
!     end do

    
!     return 


! end function lines


! subroutine spattospin1(elecs,nlines)

!     implicit none

!     type(elecintrgl), intent(inout)::elecs
!     integer, intent(in)::nlines
!     real(kind=8), dimension(:,:),allocatable ::h1ea
!     integer:: ierr, j, k,jj,kk, nspao

!     ierr=0 
!     allocate (h1ea(nlines,nlines), stat=ierr)
!     if (ierr/=0) then
!         write(0,"(a,i0)") "Error in h1ea allocation. ierr had value ", ierr
!         errorflag=1
!         return
!     end if


!     open(unit=130, file='integrals/h1ea.csv',status='old',iostat=ierr)
!     if (ierr.ne.0) then
!         write(0,"(a)") 'Error in opening h1ea.csv file'
!         errorflag = 1
!         return
!     end if
  
!     do j=1, nlines
!         read(130,*,iostat=ierr) (h1ea(j,k),k=1,nlines)
!     end do

!     close(130)


!     nspao = int(norb/2)

!     !$omp parallel shared(elecs,h1ea) private(j,k,jj,kk)
!     !$omp do 
!     do j=1,nspao
!         do k=1, nspao
!             jj=(j*2)-1
!             kk=(k*2)-1
!             ! if(h1ea(j,k).ne.0) h1ea(j,k)=1
!             elecs%h1ei(jj,kk)=h1ea(j,k)
!             elecs%h1ei(jj+1,kk+1)=elecs%h1ei(jj,kk)
!         end do
!     end do
!     !$omp end do
!     !$omp end parallel
    
!     deallocate(h1ea, stat=ierr)
!     if (ierr/=0) then
!         write(0,"(a,i0)") "Error in h1ea deallocation. ierr had value ", ierr
!         errorflag=1
!         return
!     end if

!     return

! end subroutine spattospin1

! subroutine spattospin2(elecs,nlines)

!     implicit none

!     type(elecintrgl), intent(inout)::elecs
!     integer, intent(in)::nlines
!     real(kind=8), dimension(:,:,:,:),allocatable ::h2ea
!     integer:: ierr, j, k, l, m, jj, kk, ll, mm, nspao
!     character(len=4)::val1,val2

!     ierr=0 
!     allocate (h2ea(nlines,nlines,nlines,nlines), stat=ierr)
!     if (ierr/=0) then
!         write(0,"(a,i0)") "Error in h2ea allocation. ierr had value ", ierr
!         errorflag=1
!         return
!     end if

!     nspao = int(norb/2)
!     do j=1,nlines
!         do k=1, nlines
!             write(val1,'(i0)')j
!             write(val2,'(i0)')k
!             open(unit=(131+j+k), file='integrals/h2ea_'//trim(val1)//'_'//trim(val2)//'.csv', status='old',iostat=ierr)
!             if (ierr.ne.0) then
!                 write(0,"(a)") 'Error in opening h2ea_'//trim(val1)//'_'//trim(val2)//'.csv file'
!                 errorflag = 1
!                 return
!             end if
!             do l=1, nlines
!                 read((131+j+k),*) (h2ea(j,k,l,m),m=1,nlines)
!             end do
!             close(131+j+k)
!         end do
!     end do

!     !$omp parallel shared(elecs,h2ea) private(j,k,l,m,jj,kk,ll,mm) 
!     !$omp do
!     do j=1,nspao
!         jj=(2*j)-1
!         do k=1, nspao
!             kk=(2*k)-1
!             do l=1, nspao
!                 ll=(2*l)-1
!                 do m=1, nspao
!                     mm=(2*m)-1
!                     ! if(h2ea(j,k,l,m).ne.0) h2ea(j,k,l,m)=1
!                     elecs%h2ei(jj,ll,kk,mm)=h2ea(j,k,l,m)
!                     elecs%h2ei(jj+1,ll,kk+1,mm)=h2ea(j,k,l,m)
!                     elecs%h2ei(jj,ll+1,kk,mm+1)=h2ea(j,k,l,m)
!                     elecs%h2ei(jj+1,ll+1,kk+1,mm+1)=h2ea(j,k,l,m)
!                 end do
!             end do
!         end do
!     end do
!     !$omp end do
!     !$omp end parallel

    
!     deallocate(h2ea, stat=ierr)
!     if (ierr/=0) then
!         write(0,"(a,i0)") "Error in h2ea deallocation. ierr had value ", ierr
!         errorflag=1
!         return
!     end if

!     return
                    
! end subroutine spattospin2
    





END MODULE electrons