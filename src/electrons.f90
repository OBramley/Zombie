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
        integer,dimension(:,:),allocatable::orbital_choice2,orbital_choice2_temp,orbital_choice
        integer(int16),dimension(:,:), allocatable::alive
        integer(int16),dimension(:,:), allocatable::dead
        integer(int8),dimension(:,:), allocatable::neg_a
        integer(int8),dimension(:,:), allocatable::neg_d
        integer::e1,e2,cnt,j,k,cnt2,neg
        integer::ierr=0
    
        if (errorflag .ne. 0) return
       
        write(stdout,"(a)") "Starting one electron integral allocation"
        call  one_electrons(h1ei,e1,an_cr)
        write(stdout,"(a)") "completed one electron integral allocation"
       
        write(stdout,"(a)") "Starting two electron integral allocation"
        call  two_electrons(h2ei,e2,an2_cr2)
        write(stdout,"(a)") "completed two electron integral allocation"

        elecs%num=e1+e2

        allocate(elecs%integrals(elecs%num),stat=ierr)
        allocate(alive(norb,elecs%num),stat=ierr)
        allocate(dead(norb,elecs%num),stat=ierr)
        allocate(neg_a(norb,elecs%num),stat=ierr)
        allocate(neg_d(norb,elecs%num),stat=ierr)
        if(ierr/=0) then
            write(stderr,"(a,i0)") "Error in electron integral  allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        cnt=1
        do j=1,e1
            elecs%integrals(cnt)=h1ei(j)
            alive(:,cnt)=an_cr%alive(:,j)
            dead(:,cnt)=an_cr%dead(:,j)
            neg_a(:,cnt)=an_cr%neg_alive(:,j)
            neg_d(:,cnt)=an_cr%neg_dead(:,j)
            cnt=cnt+1
        end do
        
        do j=1,e2
            elecs%integrals(cnt)=h2ei(j)
            alive(:,cnt)=an2_cr2%alive(:,j)
            dead(:,cnt)=an2_cr2%dead(:,j)
            neg_a(:,cnt)=an2_cr2%neg_alive(:,j)
            neg_d(:,cnt)=an2_cr2%neg_dead(:,j)
            cnt=cnt+1
        end do
        

        open(unit=128, file='integrals/hnuc.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(stderr,"(a)") 'Error in opening hnuc.csv file'
            errorflag = 1
            return
        end if
        read(128,*, iostat=ierr)elecs%hnuc
        if (ierr .ne. 0) then
            write(stderr,"(a)") 'Error reading Hnuc'
            errorflag = 1
            return
          end if
        close(128)

        call dealloc_oprts(an_cr)
        call dealloc_oprts(an2_cr2)
        deallocate(h1ei,h2ei,stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in h1ei  deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        !0 - no change; 1 - annihilation (cos*sin in dead); 2 - creation (sin*cos in alive); 
        !3 - both creation and annihilation (sin*sin); 4 - negative alive; -5 - negative dead; -6 - both alive and dead negative 
        allocate(orbital_choice(norb,elecs%num),stat=ierr)
        allocate(elecs%orbital_choice(norb,elecs%num),stat=ierr)
        orbital_choice=0
        do j=1,elecs%num
            neg=0
            do k=1,norb
                if(alive(k,j).eq.k)then
                    if(dead(k,j).eq.(k+norb))then
                        if(neg_a(k,j).eq.1)then
                            if(neg_d(k,j).eq.1)then
                                orbital_choice(k,j)=0
                            else 
                                write(stdout,"(a)") "Error in electron integral allocations cannot have negaitve &
                                & dead and occupied alive"
                               errorflag=1 
                               return
                            end if
                        else
                            if(neg_d(k,j).eq.1)then
                                orbital_choice(k,j)=4
                            else 
                               write(stdout,"(a)") "Error in electron integral allocations cannot have negaitve dead and alive"
                               errorflag=1 
                               return
                            end if
                        end if 
                    else if(dead(k,j).eq.0)then
                        if(neg_a(k,j).eq.1)then
                            orbital_choice(k,j)=3
                        else
                            orbital_choice(k,j)=3
                            neg=neg+1
                        end if 
                    end if 
                else if(alive(k,j).eq.(k+norb))then
                    if(dead(k,j).eq.0)then
                        if(neg_a(k,j).eq.1)then
                            orbital_choice(k,j)=2
                        else
                            orbital_choice(k,j)=2
                            neg=neg+1
                        end if
                    end if
                else if(alive(k,j).eq.0)then
                    if(dead(k,j).eq.k)then
                        if(neg_d(k,j).eq.1)then
                            orbital_choice(k,j)=1
                        else
                            orbital_choice(k,j)=1
                            neg=neg+1
                        end if 
                    end if 
                end if
            end do 
            if(modulo(neg,2).ne.0)then
                elecs%integrals(j)=elecs%integrals(j)*(-1)
            end if
        end do 
        allocate(elecs%orbital_choice3(norb))
        elecs%orbital_choice3=counter2(orbital_choice)
        do k=1,norb 
            elecs%orbital_choice(k,:)=orbital_choice(elecs%orbital_choice3(k),:)
        end do 
        deallocate(orbital_choice,stat=ierr)
        allocate(orbital_choice2(0:norb,elecs%num*2),stat=ierr)
        if(ierr/=0) then
            write(stderr,"(a,i0)") "Error in orbital_choice2  allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
   
        cnt2=0
        orbital_choice2=0
        do k=1,norb
            call sorts(elecs%orbital_choice,elecs%integrals,norb,elecs%num,k)
        end do 
      
        cnt=elecs%orbital_choice(1,1)
        do j=2,elecs%num
            if(elecs%orbital_choice(1,j).ne.cnt)then
                orbital_choice2(0,1)=orbital_choice2(0,1)+1
                orbital_choice2(1,orbital_choice2(0,1))=j-1
                cnt=elecs%orbital_choice(1,j)
            end if
        end do
        orbital_choice2(0,1)=orbital_choice2(0,1)+1
        orbital_choice2(1,orbital_choice2(0,1))=elecs%num
        do k=2, norb
            cnt=elecs%orbital_choice(k,1)
            cnt2=1
            do j=2,elecs%num
                if((elecs%orbital_choice(k,j).ne.cnt).or.(j.gt.orbital_choice2(k-1,cnt2)))then
                    orbital_choice2(0,k)=orbital_choice2(0,k)+1
                    orbital_choice2(k,orbital_choice2(0,k))=j-1
                    cnt=elecs%orbital_choice(k,j)
                    if((j.gt.orbital_choice2(k-1,cnt2)))then
                        cnt2=cnt2+1
                    end if
                end if
                
            end do
            orbital_choice2(0,k)=orbital_choice2(0,k)+1
            orbital_choice2(k,orbital_choice2(0,k))=elecs%num
        end do
        cnt=maxval(orbital_choice2(0,1:norb))
      
        allocate(orbital_choice2_temp(0:norb,cnt),stat=ierr)
        if(ierr/=0) then
            write(stderr,"(a,i0)") "Error in orbital_choice2_temp  allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        orbital_choice2_temp=0
        do k=1,norb
            cnt=1
            do j=1,orbital_choice2(0,k)
                if(elecs%orbital_choice(k,orbital_choice2(k,j)).ne.0)then
                    orbital_choice2_temp(0,k)=orbital_choice2_temp(0,k)+1
                    orbital_choice2_temp(k,(orbital_choice2_temp(0,k)*2)-1)=cnt
                    orbital_choice2_temp(k,(orbital_choice2_temp(0,k)*2))=orbital_choice2(k,j)
                    cnt=orbital_choice2(k,j)+1
                else 
                    cnt=orbital_choice2(k,j)+1
                end if
            end do
        end do
        cnt=2*maxval(orbital_choice2_temp(0,1:norb))
        allocate(elecs%orbital_choice2(0:norb,cnt),stat=ierr)
        if(ierr/=0) then
            write(stderr,"(a,i0)") "Error in eelcs%orbital_choice2  allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        elecs%orbital_choice2=orbital_choice2_temp(0:norb,1:cnt)
       
        deallocate(alive,stat=ierr)
        deallocate(dead,stat=ierr)
        deallocate(neg_a,stat=ierr)
        deallocate(neg_d,stat=ierr)
        deallocate(orbital_choice2,stat=ierr)
        deallocate(orbital_choice2_temp,stat=ierr)
        if(ierr/=0) then
            write(stderr,"(a,i0)") "Error in orbital_choice2  deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        
        write(stdout,"(a)") "1 & 2 electron integrals successfully generated"
        
        
        return
    
    end subroutine electronintegrals

    subroutine sorts(arr,arr2,N,M,col)
        implicit none
        integer, intent(in) :: M, N,col
        integer, dimension(:,:),intent(inout) :: arr
        real(wp), dimension(:),intent(inout) :: arr2
        integer, dimension(N,M)::sort
        real(wp), dimension(M):: sort2
        integer :: l, j, k, cnt,cnt2
        integer,dimension(5)::check_val
        integer::chck
        cnt=0
        cnt2=1

        if(col==1)then 
            check_val=counter(arr(1,:))
            do l=1,5
                chck=check_val(l)
                do j=1,M
                    if(arr(col,j).eq.chck)then
                        cnt=cnt+1
                        sort(:,cnt)=arr(:,j)
                        sort2(cnt)=arr2(j)
                    end if
                end do
            end do
        else 
            do while(cnt2.lt.M)
                do k=cnt2,M
                    if(arr(col-1,k).ne.arr(col-1,cnt2))then   
                        exit 
                    end if
                end do 
                check_val=counter(arr(col,cnt2:k-1))
                do l=1,5
                    chck=check_val(l)
                    do j=cnt2,k-1
                        if(arr(col,j).eq.chck)then
                            cnt=cnt+1
                            sort(:,cnt)=arr(:,j)
                            sort2(cnt)=arr2(j)
                        end if
                    end do
                end do
                cnt2=cnt+1 
            end do
        end if
        if(cnt.ne.M)then
            if(cnt.eq.M-1)then
                sort(:,M)=arr(:,M)
                sort2(M)=arr2(M)
            else
                write(stderr,"(a)") "Error in sorting"
                errorflag=1
                return
            end if
        end if 
        arr=sort
        arr2=sort2
       
    end subroutine sorts

    function counter(arr) result(indices)
        implicit none
        integer, dimension(:),intent(in)::arr
        integer,dimension(5)::vals
        integer,dimension(5)::indices
        integer::j
        j=1

        do j=0,4
            vals(j+1)=count(arr(:)==j)
        end do
     
        do j=5,1,-1
           indices(j)=maxloc(vals,1)
            vals(indices(j))=-10
        end do 
        indices=indices-1
    end function counter

    function counter2(arr) result(indices)
        implicit none
        integer, dimension(:,:),intent(in)::arr
        integer,dimension(5)::vals
        integer,dimension(norb)::vals_2
        integer,dimension(norb)::indices
        integer::j,k

        do k=1,norb
            do j=0,4
                vals(j+1)=count(arr(k,:)==j)
            end do
            vals_2(k)=maxval(vals)
        end do
    
        do j=norb,1,-1
            indices(j)=maxloc(vals_2,1)
            vals_2(indices(j))=-10
        end do 
     
    end function counter2

    subroutine two_electrons(h2ei,e2,an2_cr2)
    
        implicit none 
        real(wp), dimension(:), intent(inout),allocatable::h2ei
        type(oprts),intent(inout)::an2_cr2
        integer,intent(inout)::e2
        real(wp),dimension(:,:),allocatable::read_in
        integer::l,k,an1,an2,cr1,cr2,j,len
        integer::ierr=0
    
        len=norb*norb*norb*norb
        allocate(read_in(len,5),stat=ierr)
       
        open(unit=131, file='integrals/h2e.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(stderr,"(a)") 'Error in opening hnuc.csv file'
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
        write(stderr,"(a,i0)") "Error in 2 electron integral  allocation. ierr had value ", ierr
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
            write(stderr,"(a,i0)") "Error in read_in  deallocation. ierr had value ", ierr
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
        real(wp),dimension(:,:),allocatable::read_in
        integer::l,k,an,cr,j,len
        integer::ierr=0
        

        len=norb*norb
        allocate(read_in(len,3),stat=ierr)
    
        open(unit=129, file='integrals/h1e.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(stderr,"(a)") 'Error in opening hnuc.csv file'
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
            write(stderr,"(a,i0)") "Error in electron integral  allocation. ierr had value ", ierr
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
            write(stderr,"(a,i0)") "Error in read_in  deallocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        return
         
        return
    
    
    end subroutine one_electrons
    
    END MODULE electrons    