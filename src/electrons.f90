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
        integer:: ierr,j,k
        
    
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
        allocate(elecs%ali_dead(2*norb,elecs%num),stat=ierr)
        allocate(elecs%negs(2*norb,elecs%num),stat=ierr)
        if(ierr/=0) then
            write(0,"(a,i0)") "Error in electron integral  allocation. ierr had value ", ierr
            errorflag=1
            return
        end if

        cnt=1
        do j=1,e1
            elecs%integrals(cnt)=h1ei(j)
            elecs%ali_dead(1:norb,cnt)=an_cr%alive(:,j)
            elecs%ali_dead(norb+1:,cnt)=an_cr%dead(:,j)
            elecs%negs(1:norb,cnt)=an_cr%neg_alive(:,j)
            elecs%negs(norb+1:,cnt)=an_cr%neg_dead(:,j)
            ! do k=1,norb
                ! if(an_cr%alive(k,j)==0)then
                !     elecs%ali_dead((2*k)-1,cnt)=0
                ! else if(an_cr%alive(k,j)>norb)then 
                !     elecs%ali_dead((2*k)-1,cnt)=int((an_cr%alive(k,j)-norb)*2,kind=int16)
                ! else
                !     elecs%ali_dead((2*k)-1,cnt)=int((an_cr%alive(k,j)*2)-1,kind=int16)
                ! end if

                ! if(an_cr%dead(k,j)==0)then
                !     elecs%ali_dead((2*k),cnt)=0
                ! else if(an_cr%dead(k,j)<norb)then 
                !     elecs%ali_dead((2*k),cnt)=int((an_cr%dead(k,j)*2)-1,kind=int16)
                ! else
                !     elecs%ali_dead((2*k),cnt)=int((an_cr%dead(k,j)-norb)*2,kind=int16)
                ! end if     

                ! elecs%negs((2*k)-1,cnt)=an_cr%neg_alive(k,j)
                ! elecs%negs((2*k),cnt)=an_cr%neg_dead(k,j)
            ! end do

            cnt=cnt+1
        end do
        
        
        do j=1,e2
            elecs%integrals(cnt)=h2ei(j)
            elecs%ali_dead(1:norb,cnt)=an2_cr2%alive(:,j)
            elecs%ali_dead(norb+1:,cnt)=an2_cr2%dead(:,j)
            elecs%negs(1:norb,cnt)=an2_cr2%neg_alive(:,j)
            elecs%negs(norb+1:,cnt)=an2_cr2%neg_dead(:,j)
            ! do k=1,norb
            !     if(an2_cr2%alive(k,j)==0)then
            !         elecs%ali_dead((2*k)-1,cnt)=0
            !     else if(an2_cr2%alive(k,j)>norb)then 
            !         elecs%ali_dead((2*k)-1,cnt)=int((an2_cr2%alive(k,j)-norb)*2,kind=int16)
            !     else
            !         elecs%ali_dead((2*k)-1,cnt)=int((an2_cr2%alive(k,j)*2)-1,kind=int16)
            !     end if
            !     if(an2_cr2%dead(k,j)==0)then
            !         elecs%ali_dead((2*k),cnt)=0
            !     else if(an2_cr2%dead(k,j)<norb)then 
            !         elecs%ali_dead((2*k),cnt)=int((an2_cr2%dead(k,j)*2)-1,kind=int16)
            !     else
            !         elecs%ali_dead((2*k),cnt)=int((an2_cr2%dead(k,j)-norb)*2,kind=int16)
            !     end if   
            !     elecs%negs((2*k)-1,cnt)=an2_cr2%neg_alive(k,j)
            !     elecs%negs((2*k),cnt)=an2_cr2%neg_dead(k,j)  
            ! end do
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

        ! if(GDflg.eq.'y')then
        !     allocate(temp_diff(0:e2,norb),stat=ierr)
        !     if (ierr/=0) then
        !         write(0,"(a,i0)") "Error in temp operators allocation. ierr had value ", ierr
        !         errorflag=1
        !         return
        !     end if
        !     call alloc_oprts_2(temp_an2_cr2,e2)
        !     max=0
        !     temp_diff=0
            
        !     do j=1,norb
        !         do l=1,e2
        !             cr1=int(read_in(l+1,2))
        !             cr2=int(read_in(l+1,3))
        !             an2=int(read_in(l+1,4))
        !             an1=int(read_in(l+1,5))
        !             if((an2.eq.j).or.(cr2.eq.j).or.(an1.eq.j).or.(cr1.eq.j))then
        !                 temp_diff(0,j)=temp_diff(0,j)+1
        !                 ! print*,temp_diff(0,j)
        !                 temp_diff(temp_diff(0,j),j)=l

        !                 temp_an2_cr2%alive(:,l)=an2_cr2%ham%alive(:,l)
        !                 temp_an2_cr2%dead(:,l)=an2_cr2%ham%dead(:,l)
        !                 temp_an2_cr2%neg_alive(:,l)=an2_cr2%ham%neg_alive(:,l)
        !                 temp_an2_cr2%neg_dead(:,l)=an2_cr2%ham%neg_dead(:,l)
    
        !                 if((an2.eq.j).or.(an1.eq.j))then
        !                     if((cr1.eq.j).or.(cr2.eq.j))then
        !                         temp_an2_cr2%alive(j,l)=int(j+norb,kind=2)
        !                         temp_an2_cr2%neg_alive(j,l)=int(2*temp_an2_cr2%neg_alive(j,l),kind=1)
        !                         temp_an2_cr2%neg_dead(j,l)=0
        !                     else
        !                         temp_an2_cr2%alive(j,l)=int(j,kind=2)
        !                         temp_an2_cr2%dead(j,l)=int(j+norb,kind=2)
        !                         temp_an2_cr2%neg_alive(j,l)=int(temp_an2_cr2%neg_dead(j,l)*(-1),kind=1)
        !                     end if
        !                 else
        !                     temp_an2_cr2%alive(j,l)=int(j,kind=2)
        !                     temp_an2_cr2%dead(j,l)=int(j+norb,kind=2)
        !                     temp_an2_cr2%neg_dead(j,l)=temp_an2_cr2%neg_alive(j,l)
        !                     temp_an2_cr2%neg_alive(j,l)=int(temp_an2_cr2%neg_alive(j,l)*(-1),kind=1)
        !                 end if
        !             end if
        !         end do
        !         if(max.lt.temp_diff(0,j))then 
        !             max=temp_diff(0,j)
        !         end if
        !         call alloc_oprts_2(an2_cr2%diff(j),temp_diff(0,j))
        !         do l=1,temp_diff(0,j)
        !             an2_cr2%diff(j)%alive(:,l)=temp_an2_cr2%alive(:,temp_diff(l,j))
        !             an2_cr2%diff(j)%dead(:,l)=temp_an2_cr2%dead(:,temp_diff(l,j))
        !             an2_cr2%diff(j)%neg_alive(:,l)=temp_an2_cr2%neg_alive(:,temp_diff(l,j))
        !             an2_cr2%diff(j)%neg_dead(:,l)=temp_an2_cr2%neg_dead(:,temp_diff(l,j))
        !         end do
        !     end do

        !     allocate(an2_cr2%dcnt(0:max,norb))
        !     an2_cr2%dcnt=0
        !     do j=1, norb 
        !         an2_cr2%dcnt(0:max,j)=temp_diff(0:max,j)
        !         ! print*,an2_cr2%dcnt(0:,j)
        !     end do
        !     ! an2_cr2%dcnt(0,:)=temp_diff(0,:)
        !     ! do j=1, norb 
        !     !     an2_cr2%dcnt(1:,j)=temp_diff(1:max,j)
        !     ! end do 
        !     deallocate(temp_diff,stat=ierr)
        !     deallocate(read_in,stat=ierr)
        !     if (ierr/=0) then
        !         write(0,"(a,i0)") "Error in temp operators deallocation. ierr had value ", ierr
        !         errorflag=1
        !         return
        !     end if
        !     call dealloc_oprts_2(temp_an2_cr2)
 
        ! end if
        
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
         
        ! if(GDflg.eq.'y')then
          
        !     allocate(temp_diff(0:e1,norb),stat=ierr)
        !     if (ierr/=0) then
        !         write(0,"(a,i0)") "Error in temp operators allocation. ierr had value ", ierr
        !         errorflag=1
        !         return
        !     end if
        !     temp_diff=0
        !     call alloc_oprts_2(temp_an_cr,e1) 
            
        !     max=0
        !     do j=1,norb
        !         do l=1,e1
        !             an=int(read_in(l+1,2))
        !             cr=int(read_in(l+1,3))
        !             if((an.eq.j).or.(cr.eq.j))then
        !                 temp_diff(0,j)=temp_diff(0,j)+1
        !                 temp_diff(temp_diff(0,j),j)=l

        !                 temp_an_cr%alive(:,l)=an_cr%ham%alive(:,l)
        !                 temp_an_cr%dead(:,l)=an_cr%ham%dead(:,l)
        !                 temp_an_cr%neg_alive(:,l)=an_cr%ham%neg_alive(:,l)
        !                 temp_an_cr%neg_dead(:,l)=an_cr%ham%neg_dead(:,l)
        !                 if(an.eq.cr)then 
        !                     temp_an_cr%alive(j,l)=int(j+norb,kind=2)
        !                     temp_an_cr%neg_alive(j,l)=int(2*temp_an_cr%neg_alive(j,l),kind=1)
        !                     temp_an_cr%neg_dead(j,l)=0
        !                 else if(j.eq.an)then
        !                     temp_an_cr%alive(j,l)=int(j,kind=2)
        !                     temp_an_cr%dead(j,l)=int(j+norb,kind=2)
        !                     temp_an_cr%neg_alive(j,l)=int(temp_an_cr%neg_dead(j,l)*(-1),kind=1)
        !                 else if(j.eq.cr)then
        !                     temp_an_cr%alive(j,l)=int(j,kind=2)
        !                     temp_an_cr%dead(j,l)=int(j+norb,kind=2)
        !                     temp_an_cr%neg_dead(j,l)=temp_an_cr%neg_alive(j,l)
        !                     temp_an_cr%neg_alive(j,l)=int(temp_an_cr%neg_alive(j,l)*(-1),kind=1)
        !                 end if
        !             end if
        !         end do
        !         if(max.lt.temp_diff(0,j))then 
        !             max=temp_diff(0,j)
        !         end if
        !         call alloc_oprts_2(an_cr%diff(j),temp_diff(0,j))
        !         do l=1,temp_diff(0,j)
        !             an_cr%diff(j)%alive(:,l)=temp_an_cr%alive(:,temp_diff(l,j))
        !             an_cr%diff(j)%dead(:,l)=temp_an_cr%dead(:,temp_diff(l,j))
        !             an_cr%diff(j)%neg_alive(:,l)=temp_an_cr%neg_alive(:,temp_diff(l,j))
        !             an_cr%diff(j)%neg_dead(:,l)=temp_an_cr%neg_dead(:,temp_diff(l,j))
        !         end do
        !     end do
          
        !     allocate(an_cr%dcnt(0:max,norb))
         
        !     do j=1, norb 
        !         an_cr%dcnt(0:max,j)=temp_diff(0:max,j)
        !     end do 
            
        !     deallocate(temp_diff,stat=ierr)
        !     deallocate(read_in,stat=ierr)
        !     if (ierr/=0) then
        !         write(0,"(a,i0)") "Error in temp operators deallocation. ierr had value ", ierr
        !         errorflag=1
        !         return
        !     end if
          
        !     call dealloc_oprts_2(temp_an_cr)

        
        
        ! end if
        
        
        return
    
    
    end subroutine one_electrons
    
    END MODULE electrons