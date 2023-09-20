MODULE electrons

    use globvars
    use alarrays
    
    contains
    
    ! Program to generate one and two electron integrals from PyScf
    ! Some of this code has been adapted from George Booth
    subroutine electronintegrals(elecs)
    
        implicit none
    
        type(elecintrgl), intent(inout)::elecs
        type(oprts)::an_cr,an2_cr2
        real(wp), dimension(:), allocatable::h1ei
        real(wp), dimension(:), allocatable::h2ei
        integer::e1,e2,cnt
        integer:: ierr,j
        
    
        if (errorflag .ne. 0) return
       
    
        ierr = 0
        
        write(6,"(a)") "Starting one electron integral allocation"
        call  one_electrons(h1ei,e1,an_cr)
        write(6,"(a)") "completed one electron integral allocation"
       
        write(6,"(a)") "Starting two electron integral allocation"
        call  two_electrons(h2ei,e2,an2_cr2)
        write(6,"(a)") "completed two electron integral allocation"

        elecs%num=e1+e2

        allocate(elecs%integrals(elecs%num),stat=ierr)
        allocate(elecs%alive(norb,elecs%num),stat=ierr)
        allocate(elecs%dead(norb,elecs%num),stat=ierr)
        allocate(elecs%neg_a(norb,elecs%num),stat=ierr)
        allocate(elecs%neg_d(norb,elecs%num),stat=ierr)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in electron integral  allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        cnt=1
        do j=1,e1
            elecs%integrals(cnt)=h1ei(j)
            elecs%alive(:,cnt)=an_cr%alive(:,j)
            elecs%dead(:,cnt)=an_cr%dead(:,j)
            elecs%neg_a(:,cnt)=an_cr%neg_alive(:,j)
            elecs%neg_d(:,cnt)=an_cr%neg_dead(:,j)
            cnt=cnt+1
        end do
        
        
        do j=1,e2
            elecs%integrals(cnt)=h2ei(j)
            elecs%alive(:,cnt)=an2_cr2%alive(:,j)
            elecs%dead(:,cnt)=an2_cr2%dead(:,j)
            elecs%neg_a(:,cnt)=an2_cr2%neg_alive(:,j)
            elecs%neg_d(:,cnt)=an2_cr2%neg_dead(:,j)
            cnt=cnt+1
        end do
        

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

        call dealloc_oprts(an_cr)
        call dealloc_oprts(an2_cr2)
        deallocate(h1ei,h2ei,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in h1ei  deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        
        write(6,"(a)") "1 & 2 electron integrals successfully generated"
    
        return
    
    end subroutine electronintegrals
    
    
    subroutine two_electrons(h2ei,e2,an2_cr2)
    
        implicit none 
        real(wp), dimension(:), intent(inout),allocatable::h2ei
        type(oprts),intent(inout)::an2_cr2
        integer,intent(inout)::e2
        real(kind=8),dimension(:,:),allocatable::read_in
        integer::ierr,l,k,an1,an2,cr1,cr2,j,max,len
       
    
        ierr=0
        len=norb*norb*norb*norb
        allocate(read_in(len,5),stat=ierr)
       
        open(unit=131, file='integrals/h2e.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a)") 'Error in opening hnuc.csv file'
            errorflag = 1
            return
        end if

        read(131,*) (read_in(1,k),k=1,5)
        e2=int(read_in(1,1))
        do j=2,e2+1
            read(131,*) (read_in(j,k),k=1,5)
        end do 
       
        close(131)
      

        allocate (h2ei(e2),stat=ierr)
        if (ierr/=0) then
        write(0,"(a,i0)") "Error in 2 electron integral  allocation. ierr had value ", ierr
        errorflag=1
        return
        end if
       
        h2ei=read_in(2:,1)
        h2ei=h2ei*0.5_wp
        call alloc_oprts(an2_cr2,e2)
       
    
        do l=1,e2
            cr1=int(read_in(l+1,2)) 
            cr2=int(read_in(l+1,3)) 
            an2=int(read_in(l+1,4)) 
            an1=int(read_in(l+1,5))
           
            !annihilation
            an2_cr2%dead(an1,l)=an2_cr2%alive(an1,l)
            an2_cr2%neg_dead(an1,l)=an2_cr2%neg_alive(an1,l)
            an2_cr2%alive(an1,l)=0
            an2_cr2%neg_alive(:an1-1,l)=int(an2_cr2%neg_alive(:an1-1,l)*(-1),kind=int8)
            !annihilation
            an2_cr2%dead(an2,l)=an2_cr2%alive(an2,l)
            an2_cr2%neg_dead(an2,l)=an2_cr2%neg_alive(an2,l)
            an2_cr2%alive(an2,l)=0
            an2_cr2%neg_alive(:an2-1,l)=int(an2_cr2%neg_alive(:an2-1,l)*(-1),kind=int8) 
            !creation 
            an2_cr2%alive(cr1,l)=an2_cr2%dead(cr1,l)
            an2_cr2%neg_alive(cr1,l)=an2_cr2%neg_dead(cr1,l)
            an2_cr2%dead(cr1,l)=0
            an2_cr2%neg_alive(:cr1-1,l)=int(an2_cr2%neg_alive(:cr1-1,l)*(-1),kind=int8)
            !creation 
            an2_cr2%alive(cr2,l)=an2_cr2%dead(cr2,l)
            an2_cr2%neg_alive(cr2,l)=an2_cr2%neg_dead(cr2,l)
            an2_cr2%dead(cr2,l)=0
            an2_cr2%neg_alive(:cr2-1,l)=int(an2_cr2%neg_alive(:cr2-1,l)*(-1),kind=int8)
    
        end do

        deallocate(read_in, stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in read_in  deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return

    
    end subroutine two_electrons
    
    
    subroutine one_electrons(h1ei,e1,an_cr)
    
        implicit none 
        type(oprts),intent(inout)::an_cr
        real(wp), dimension(:), intent(inout),allocatable::h1ei
        integer,intent(inout)::e1
        real(kind=8),dimension(:,:),allocatable::read_in
        integer::ierr,l,k,an,cr,j,max,len

        
        ierr=0

        len=norb*norb
        allocate(read_in(len,3),stat=ierr)
    
        open(unit=129, file='integrals/h1e.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a)") 'Error in opening hnuc.csv file'
            errorflag = 1
            return
        end if
        read(129,*) (read_in(1,k),k=1,3)
        e1=int(read_in(1,1))
       
        do j=2,e1+1
            read(129,*) (read_in(j,k),k=1,3)
        end do 
        close(129)
        
       
        allocate (h1ei(e1), stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in electron integral  allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        h1ei=read_in(2:,1)
       
        call alloc_oprts(an_cr,e1)
          
        do l=1,e1
            an=int(read_in(l+1,2))
            cr=int(read_in(l+1,3))
            !annihilation
            an_cr%dead(an,l)=an_cr%alive(an,l)
            an_cr%neg_dead(an,l)=an_cr%neg_alive(an,l)
            an_cr%alive(an,l)=0
            an_cr%neg_alive(:an-1,l)=int(an_cr%neg_alive(:an-1,l)*(-1),kind=int8)
            !creation 
            an_cr%alive(cr,l)=an_cr%dead(cr,l)
            an_cr%neg_alive(cr,l)=an_cr%neg_dead(cr,l)
            an_cr%dead(cr,l)=0
            an_cr%neg_alive(:cr-1,l)=int(an_cr%neg_alive(:cr-1,l)*(-1),kind=int8)
        end do 
         
        deallocate(read_in, stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in read_in  deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return
         
        return
    
    
    end subroutine one_electrons
    
    END MODULE electrons