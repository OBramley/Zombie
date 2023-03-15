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
        ! real::e1in,e2in
        integer:: ierr,e1,e2
        
    
        if (errorflag .ne. 0) return
       
    
        ierr = 0
        ! open(unit=129, file='integrals/h1e.csv',status='old',iostat=ierr)
        ! if (ierr.ne.0) then
        !     write(0,"(a)") 'Error in opening h1ec.csv file'
        !     errorflag = 1
        !     return
        ! end if
        ! read(129,*, iostat=ierr)e1in
        ! if (ierr .ne. 0) then
        !     write(0,"(a)") 'Error reading h1e'
        !     errorflag = 1
        !     return
        !   end if
        ! close(129)
        
      
       
        
    
        ! open(unit=131, file='integrals/h2e.csv',status='old',iostat=ierr)
        ! if (ierr.ne.0) then
        !     write(0,"(a)") 'Error in opening h2e.csv file'
        !     errorflag = 1
        !     return
        ! end if
        ! read(131,*, iostat=ierr)e2in
        ! if (ierr .ne. 0) then
        !     write(0,"(a)") 'Error reading h1e'
        !     errorflag = 1
        !     return
        !   end if
        ! close(131)

        ! e1=int(e1in)
        ! e2=int(e2in)
        
        ! call allocintgrl(elecs,e1,e2)
        write(6,"(a)") "Starting one electron integral allocation"
        call  one_electrons(elecs,an_cr)
        write(6,"(a)") "completed one electron integral allocation"
       
        write(6,"(a)") "Starting two electron integral allocation"
        call  two_electrons(elecs,an2_cr2)
        write(6,"(a)") "completed two electron integral allocation"
       
        
    
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
    
    
    subroutine two_electrons(elecs,an2_cr2)
    
        implicit none 
        type(elecintrgl), intent(inout)::elecs
        type(oprts),intent(inout)::an2_cr2
        integer::e2
        type(oprts_2)::temp_an2_cr2
        real(kind=8),dimension(norb*norb*norb*norb,5)::read_in
        integer::ierr,l,k,an1,an2,cr1,cr2,j,max
        integer,allocatable,dimension(:,:)::temp_diff
        ! integer,allocatable,dimension(:,:,:)::temp_hess

        
        ierr=0
        
        open(unit=131, file='integrals/h2e.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a)") 'Error in opening hnuc.csv file'
            errorflag = 1
            return
        end if

        read(131,*) (read_in(1,k),k=1,5)
        elecs%h2_num=int(read_in(1,1))
        e2=elecs%h2_num
        do j=2,e2+1
            read(131,*) (read_in(j,k),k=1,5)
        end do 
       
        close(131)
      

        allocate (elecs%h2ei(elecs%h2_num),stat=ierr)
        if (ierr/=0) then
        write(0,"(a,i0)") "Error in electron integral  allocation. ierr had value ", ierr
        errorflag=1
        return
        end if
       
        elecs%h2ei=read_in(2:,1)
       
        call alloc_oprts(an2_cr2,e2)
        
    
        do l=1,e2
            cr1=int(read_in(l+1,2)) 
            cr2=int(read_in(l+1,3)) 
            an2=int(read_in(l+1,4)) 
            an1=int(read_in(l+1,5))
           
            !annihilation
            an2_cr2%ham%dead(an1,l)=an2_cr2%ham%alive(an1,l)
            an2_cr2%ham%neg_dead(an1,l)=an2_cr2%ham%neg_alive(an1,l)
            an2_cr2%ham%alive(an1,l)=0
            an2_cr2%ham%neg_alive(:an1-1,l)=int(an2_cr2%ham%neg_alive(:an1-1,l)*(-1),kind=1)
            !annihilation
            an2_cr2%ham%dead(an2,l)=an2_cr2%ham%alive(an2,l)
            an2_cr2%ham%neg_dead(an2,l)=an2_cr2%ham%neg_alive(an2,l)
            an2_cr2%ham%alive(an2,l)=0
            an2_cr2%ham%neg_alive(:an2-1,l)=int(an2_cr2%ham%neg_alive(:an2-1,l)*(-1),kind=1) 
            !creation 
            an2_cr2%ham%alive(cr1,l)=an2_cr2%ham%dead(cr1,l)
            an2_cr2%ham%neg_alive(cr1,l)=an2_cr2%ham%neg_dead(cr1,l)
            an2_cr2%ham%dead(cr1,l)=0
            an2_cr2%ham%neg_alive(:cr1-1,l)=int(an2_cr2%ham%neg_alive(:cr1-1,l)*(-1),kind=1)
            !creation 
            an2_cr2%ham%alive(cr2,l)=an2_cr2%ham%dead(cr2,l)
            an2_cr2%ham%neg_alive(cr2,l)=an2_cr2%ham%neg_dead(cr2,l)
            an2_cr2%ham%dead(cr2,l)=0
            an2_cr2%ham%neg_alive(:cr2-1,l)=int(an2_cr2%ham%neg_alive(:cr2-1,l)*(-1),kind=1)
    
        end do

        
       
        if(GDflg.eq.'y')then
            allocate(temp_diff(0:e2,norb),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in temp operators allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
            call alloc_oprts_2(temp_an2_cr2,e2)
            max=0

            do j=1,norb
                do l=1,e2
                    cr1=int(read_in(l+1,2))
                    cr2=int(read_in(l+1,3))
                    an2=int(read_in(l+1,4))
                    an1=int(read_in(l+1,5))
                    if((an2.eq.j).or.(cr2.eq.j).or.(an1.eq.j).or.(cr1.eq.j))then
                        temp_diff(0,j)=int(temp_diff(0,j)+1,kind=2)
                        temp_diff(temp_diff(0,j),j)=l

                        temp_an2_cr2%alive(:,l)=an2_cr2%ham%alive(:,l)
                        temp_an2_cr2%dead(:,l)=an2_cr2%ham%dead(:,l)
                        temp_an2_cr2%neg_alive(:,l)=an2_cr2%ham%neg_alive(:,l)
                        temp_an2_cr2%neg_dead(:,l)=an2_cr2%ham%neg_dead(:,l)
    
                        if((an2.eq.j).or.(an1.eq.j))then
                            if((cr1.eq.j).or.(cr2.eq.j))then
                                temp_an2_cr2%alive(j,l)=int(j+norb,kind=2)
                                temp_an2_cr2%neg_alive(j,l)=int(2*temp_an2_cr2%neg_alive(j,l),kind=1)
                                temp_an2_cr2%neg_dead(j,l)=0
                            else
                                temp_an2_cr2%alive(j,l)=int(j,kind=2)
                                temp_an2_cr2%dead(j,l)=int(j+norb,kind=2)
                                temp_an2_cr2%neg_alive(j,l)=int(temp_an2_cr2%neg_dead(j,l)*(-1),kind=1)
                            end if
                        else
                            temp_an2_cr2%alive(j,l)=int(j,kind=2)
                            temp_an2_cr2%dead(j,l)=int(j+norb,kind=2)
                            temp_an2_cr2%neg_dead(j,l)=temp_an2_cr2%neg_alive(j,l)
                            temp_an2_cr2%neg_alive(j,l)=int(temp_an2_cr2%neg_alive(j,l)*(-1),kind=1)
                        end if
                    end if
                end do
                if(max.lt.temp_diff(0,j))then 
                    max=temp_diff(0,j)
                end if
                call alloc_oprts_2(an2_cr2%diff(j),temp_diff(0,j))
                do l=1,temp_diff(0,j)
                    an2_cr2%diff(j)%alive(:,l)=temp_an2_cr2%alive(:,temp_diff(l,j))
                    an2_cr2%diff(j)%dead(:,l)=temp_an2_cr2%dead(:,temp_diff(l,j))
                    an2_cr2%diff(j)%neg_alive(:,l)=temp_an2_cr2%neg_alive(:,temp_diff(l,j))
                    an2_cr2%diff(j)%neg_dead(:,l)=temp_an2_cr2%neg_dead(:,temp_diff(l,j))
                end do
            end do

            allocate(an2_cr2%dcnt(0:max,norb))
            an2_cr2%dcnt(0,:)=temp_diff(0,:)
            do j=1, norb 
                an2_cr2%dcnt(1:,j)=temp_diff(1:max,j)
            end do 
            deallocate(temp_diff,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in temp operators deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if

            ! allocate(temp_hess(0:e2,norb,norb),stat=ierr)
            ! if (ierr/=0) then
            !     write(0,"(a,i0)") "Error in temp operators allocation. ierr had value ", ierr
            !     errorflag=1
            !     return
            ! end if
            ! max=0

            ! do j=1,norb
            !     do k=j,norb 
            !         do l=1,e2
            !             cr1=int(read_in(l+1,2))
            !             cr2=int(read_in(l+1,3))
            !             an2=int(read_in(l+1,4))
            !             an1=int(read_in(l+1,5))
            !             if(((an2.eq.j).or.(cr2.eq.j).or.(an1.eq.j).or.(cr1.eq.j)))then
            !                 if((an2.eq.k).or.(cr2.eq.k).or.(an1.eq.k).or.(cr1.eq.k))then
            !                     temp_hess(0,j,k)=int(temp_hess(0,j,k)+1,kind=2)
            !                     temp_hess(temp_hess(0,j,k),j,k)=l

            !                     temp_an2_cr2%alive(:,l)=an2_cr2%ham%alive(:,l)
            !                     temp_an2_cr2%dead(:,l)=an2_cr2%ham%dead(:,l)
            !                     temp_an2_cr2%neg_alive(:,l)=an2_cr2%ham%neg_alive(:,l)
            !                     temp_an2_cr2%neg_dead(:,l)=an2_cr2%ham%neg_dead(:,l)
            !                     if(j.ne.k)then
            !                         if((j.eq.an1).or.(j.eq.an2))then !j=an1 or j=an2
            !                             if((j.eq.cr1).or.(j.eq.cr2))then !j=an1 or j=an2 and j=cr1 or j=cr2
            !                                 temp_an2_cr2%alive(j,l)=int(j+norb,kind=2)
            !                                 temp_an2_cr2%neg_alive(j,l)=int(2*temp_an2_cr2%neg_alive(j,l),kind=1)
            !                                 temp_an2_cr2%neg_dead(j,l)=0
            !                             else 
            !                                 temp_an2_cr2%alive(j,l)=int(j,kind=2)
            !                                 temp_an2_cr2%dead(j,l)=int(j+norb,kind=2)
            !                                 temp_an2_cr2%neg_alive(j,l)=int(temp_an2_cr2%neg_dead(j,l)*(-1),kind=1)
            !                             end if
            !                         else 
            !                             temp_an2_cr2%alive(j,l)=int(j,kind=2)
            !                             temp_an2_cr2%dead(j,l)=int(j+norb,kind=2)
            !                             temp_an2_cr2%neg_dead(j,l)=temp_an2_cr2%neg_alive(j,l)
            !                             temp_an2_cr2%neg_alive(j,l)=int(temp_an2_cr2%neg_alive(j,l)*(-1),kind=1)
            !                         end if
                                    
            !                         if((k.eq.an1).or.(k.eq.an2))then !k=an1 or k=an2
            !                             if((k.eq.cr1).or.(k.eq.cr2))then !k=an1 or k=an2 and k=cr1 or k=cr2
            !                                 temp_an2_cr2%alive(k,l)=int(k+norb,kind=2)
            !                                 temp_an2_cr2%neg_alive(k,l)=int(2*temp_an2_cr2%neg_alive(k,l),kind=1)
            !                                 temp_an2_cr2%neg_dead(k,l)=0
            !                             else 
            !                                 temp_an2_cr2%alive(k,l)=int(k,kind=2)
            !                                 temp_an2_cr2%dead(k,l)=int(k+norb,kind=2)
            !                                 temp_an2_cr2%neg_alive(k,l)=int(temp_an2_cr2%neg_dead(k,l)*(-1),kind=1)
            !                             end if
            !                         else 
            !                             temp_an2_cr2%alive(k,l)=int(k,kind=2)
            !                             temp_an2_cr2%dead(k,l)=int(k+norb,kind=2)
            !                             temp_an2_cr2%neg_dead(k,l)=temp_an2_cr2%neg_alive(k,l)
            !                             temp_an2_cr2%neg_alive(k,l)=int(temp_an2_cr2%neg_alive(k,l)*(-1),kind=1)
            !                         end if 
            !                     else
            !                         if((j.eq.an1).or.(j.eq.an2))then !j=k=an1 or j=k=an2
            !                             if((j.eq.cr1).or.(j.eq.cr2))then !j=k=an1 or j=k=an2 and j=k=cr1 or j=k=cr2 double differentiate 
            !                                 temp_an2_cr2%alive(j,l)=int(j,kind=2)
            !                                 temp_an2_cr2%dead(j,l)=int(j+norb,kind=2)
            !                                 temp_an2_cr2%neg_dead(j,l)=int(2*temp_an2_cr2%neg_alive(j,l),kind=1)  
            !                                 temp_an2_cr2%neg_alive(j,l)=int(-2*temp_an2_cr2%neg_alive(j,l),kind=1)
            !                             else !double differentiate at j annihilation
            !                                 temp_an2_cr2%alive(j,l)=int(j+norb,kind=2)
            !                                 temp_an2_cr2%neg_alive(j,l)=int(-4*temp_an2_cr2%neg_dead(j,l),kind=1)
            !                                 temp_an2_cr2%neg_dead(j,l)=0
            !                             end if
            !                         else !j=k=cr1 or j=k=cr2 so double differentiate at j creation
            !                             temp_an2_cr2%alive(j,l)=int(j+norb,kind=2)
            !                             temp_an2_cr2%neg_alive(k,l)=int(-4*temp_an2_cr2%neg_alive(k,l),kind=1)
            !                             temp_an2_cr2%neg_dead(j,l)=0
            !                         end if 



                                    
            !                     end if
            !                 end if
            !             end if
            !         end do
            !         if(max.lt.temp_hess(0,j,k))then 
            !             max=temp_hess(0,j,k)
            !         end if
            !         call alloc_oprts_2(an2_cr2%hess(j,k),temp_hess(0,j,k))
            !         do l=1,temp_hess(0,j,k)
            !             an2_cr2%hess(j,k)%alive(:,l)=temp_an2_cr2%alive(:,temp_hess(l,j,k))
            !             an2_cr2%hess(j,k)%dead(:,l)=temp_an2_cr2%dead(:,temp_hess(l,j,k))
            !             an2_cr2%hess(j,k)%neg_alive(:,l)=temp_an2_cr2%neg_alive(:,temp_hess(l,j,k))
            !             an2_cr2%hess(j,k)%neg_dead(:,l)=temp_an2_cr2%neg_dead(:,temp_hess(l,j,k))
            !         end do 
            !     end do
            ! end do 

            ! allocate(an2_cr2%hcnt(0:max,norb,norb))
            ! an2_cr2%hcnt(0,:,:)=temp_hess(0,:,:)
            ! do j=1, norb
            !     do k=1,norb 
            !         an2_cr2%hcnt(1:,j,k)=temp_hess(1:max,j,k)
            !     end do
            ! end do 
            ! deallocate(temp_hess,stat=ierr)
            ! if (ierr/=0) then
            !     write(0,"(a,i0)") "Error in temp operators deallocation. ierr had value ", ierr
            !     errorflag=1
            !     return
            ! end if

            call dealloc_oprts_2(temp_an2_cr2)
        end if
        
        return
    
    end subroutine two_electrons
    
    
    subroutine one_electrons(elecs,an_cr)
    
        implicit none 
        type(elecintrgl), intent(inout)::elecs
        type(oprts),intent(inout)::an_cr
        type(oprts_2)::temp_an_cr
        integer::e1
        real(kind=8),dimension(norb*norb,3)::read_in
        integer::ierr,l,k,an,cr,j,max
        integer,allocatable,dimension(:,:)::temp_diff
        ! integer,allocatable,dimension(:,:,:)::temp_hess
        
        ierr=0
    
        open(unit=129, file='integrals/h1e.csv',status='old',iostat=ierr)
        if (ierr.ne.0) then
            write(0,"(a)") 'Error in opening hnuc.csv file'
            errorflag = 1
            return
        end if
        read(129,*) (read_in(1,k),k=1,3)
        elecs%h1_num=int(read_in(1,1))
        e1=elecs%h1_num
        do j=2,e1+1
            read(129,*) (read_in(j,k),k=1,3)
        end do 
        close(129)
        
        elecs%h1_num=int(read_in(1,1))
        allocate (elecs%h1ei( elecs%h1_num), stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in electron integral  allocation. ierr had value ", ierr
            errorflag=1
            return
        end if
        elecs%h1ei=read_in(2:,1)
        e1=elecs%h1_num
        call alloc_oprts(an_cr,e1)
          
        do l=1,e1
            an=int(read_in(l+1,2))
            cr=int(read_in(l+1,3))
            !annihilation
            an_cr%ham%dead(an,l)=an_cr%ham%alive(an,l)
            an_cr%ham%neg_dead(an,l)=an_cr%ham%neg_alive(an,l)
            an_cr%ham%alive(an,l)=0
            an_cr%ham%neg_alive(:an-1,l)=int(an_cr%ham%neg_alive(:an-1,l)*(-1),kind=1)
            !creation 
            an_cr%ham%alive(cr,l)=an_cr%ham%dead(cr,l)
            an_cr%ham%neg_alive(cr,l)=an_cr%ham%neg_dead(cr,l)
            an_cr%ham%dead(cr,l)=0
            an_cr%ham%neg_alive(:cr-1,l)=int(an_cr%ham%neg_alive(:cr-1,l)*(-1),kind=1)
        end do 
          
        if(GDflg.eq.'y')then
            ! print*,'here'   
            allocate(temp_diff(0:e1,norb),stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in temp operators allocation. ierr had value ", ierr
                errorflag=1
                return
            end if
            ! print*,'here'   
            call alloc_oprts_2(temp_an_cr,e1) 
            ! print*,'here'   
            max=0
            do j=1,norb
                do l=1,e1
                    an=int(read_in(l+1,2))
                    cr=int(read_in(l+1,3))
                    if((an.eq.j).or.(cr.eq.j))then
                        temp_diff(0,j)=int(temp_diff(0,j)+1,kind=2)
                        temp_diff(temp_diff(0,j),j)=l

                        temp_an_cr%alive(:,l)=an_cr%ham%alive(:,l)
                        temp_an_cr%dead(:,l)=an_cr%ham%dead(:,l)
                        temp_an_cr%neg_alive(:,l)=an_cr%ham%neg_alive(:,l)
                        temp_an_cr%neg_dead(:,l)=an_cr%ham%neg_dead(:,l)
                        if(an.eq.cr)then 
                            temp_an_cr%alive(j,l)=int(j+norb,kind=2)
                            temp_an_cr%neg_alive(j,l)=int(2*temp_an_cr%neg_alive(j,l),kind=1)
                            temp_an_cr%neg_dead(j,l)=0
                        else if(j.eq.an)then
                            temp_an_cr%alive(j,l)=int(j,kind=2)
                            temp_an_cr%dead(j,l)=int(j+norb,kind=2)
                            temp_an_cr%neg_alive(j,l)=int(temp_an_cr%neg_dead(j,l)*(-1),kind=1)
                        else if(j.eq.cr)then
                            temp_an_cr%alive(j,l)=int(j,kind=2)
                            temp_an_cr%dead(j,l)=int(j+norb,kind=2)
                            temp_an_cr%neg_dead(j,l)=temp_an_cr%neg_alive(j,l)
                            temp_an_cr%neg_alive(j,l)=int(temp_an_cr%neg_alive(j,l)*(-1),kind=1)
                        end if
                    end if
                end do
                if(max.lt.temp_diff(0,j))then 
                    max=temp_diff(0,j)
                end if
                call alloc_oprts_2(an_cr%diff(j),temp_diff(0,j))
                do l=1,temp_diff(0,j)
                    an_cr%diff(j)%alive(:,l)=temp_an_cr%alive(:,temp_diff(l,j))
                    an_cr%diff(j)%dead(:,l)=temp_an_cr%dead(:,temp_diff(l,j))
                    an_cr%diff(j)%neg_alive(:,l)=temp_an_cr%neg_alive(:,temp_diff(l,j))
                    an_cr%diff(j)%neg_dead(:,l)=temp_an_cr%neg_dead(:,temp_diff(l,j))
                end do
            end do
            ! print*,'here'   
            allocate(an_cr%dcnt(0:max,norb))
         
            do j=1, norb 
                an_cr%dcnt(0:,j)=temp_diff(0:max,j)
            end do 
            deallocate(temp_diff,stat=ierr)
            if (ierr/=0) then
                write(0,"(a,i0)") "Error in temp operators deallocation. ierr had value ", ierr
                errorflag=1
                return
            end if

            ! allocate(temp_hess(0:e1,norb,norb),stat=ierr)
            ! if (ierr/=0) then
            !     write(0,"(a,i0)") "Error in temp operators allocation. ierr had value ", ierr
            !     errorflag=1
            !     return
            ! end if
            ! max=0
            ! do j=1,norb
            !     do k=j,norb
            !         do l=1,e1
            !             an=int(read_in(l+1,2))
            !             cr=int(read_in(l+1,3)) 
            !             if(((an.eq.j).or.(cr.eq.j)).and.((an.eq.k).or.(cr.eq.k)))then
            !                 temp_hess(0,j,k)=int(temp_hess(0,j,k)+1,kind=2)
            !                 temp_hess(temp_hess(0,j,k),j,k)=l

            !                 temp_an_cr%alive(:,l)=an_cr%ham%alive(:,l)
            !                 temp_an_cr%dead(:,l)=an_cr%ham%dead(:,l)
            !                 temp_an_cr%neg_alive(:,l)=an_cr%ham%neg_alive(:,l)
            !                 temp_an_cr%neg_dead(:,l)=an_cr%ham%neg_dead(:,l)


            !                 if(j.ne.k)then
            !                     temp_an_cr%alive(an,l)=int(an,kind=2)
            !                     temp_an_cr%dead(an,l)=int(an+norb,kind=2)
            !                     temp_an_cr%neg_alive(an,l)=int(temp_an_cr%neg_dead(an,l)*(-1),kind=1)

            !                     temp_an_cr%alive(cr,l)=int(cr,kind=2)
            !                     temp_an_cr%dead(cr,l)=int(cr+norb,kind=2)
            !                     temp_an_cr%neg_dead(cr,l)=temp_an_cr%neg_alive(cr,l)
            !                     temp_an_cr%neg_alive(cr,l)=int(temp_an_cr%neg_alive(cr,l)*(-1),kind=1)
            !                 else
            !                     if(an.eq.cr)then
            !                         temp_an_cr%alive(j,l)=int(j,kind=2)
            !                         temp_an_cr%dead(j,l)=int(j+norb,kind=2)
            !                         temp_an_cr%neg_alive(j,l)=int(-2*temp_an_cr%neg_alive(j,l),kind=1)
            !                         temp_an_cr%neg_dead(j,l)=int(2*temp_an_cr%neg_alive(j,l),kind=1)
            !                     else if(j.eq.an)then 
            !                         temp_an_cr%alive(an,l)=int(an+norb,kind=2)
            !                         temp_an_cr%neg_alive(an,l)=int(-4*temp_an_cr%neg_dead(an,l),kind=1)
            !                         temp_an_cr%neg_dead(an,l)=0
            !                     else if(j.eq.cr)then
            !                         temp_an_cr%alive(cr,l)=int(cr+norb,kind=2)
            !                         temp_an_cr%neg_alive(cr,l)=int(-4*temp_an_cr%neg_alive(cr,l),kind=1)
            !                         temp_an_cr%neg_dead(cr,l)=0
            !                     end if
            !                 end if 
            !             end if
            !         end do
            !         if(max.lt.temp_hess(0,j,k))then 
            !             max=temp_hess(0,j,k)
            !         end if
            !         call alloc_oprts_2(an_cr%hess(j,k),temp_hess(0,j,k))
            !         do l=1,temp_hess(0,j,k)
            !             an_cr%hess(j,k)%alive(:,l)=temp_an_cr%alive(:,temp_hess(l,j,k))
            !             an_cr%hess(j,k)%dead(:,l)=temp_an_cr%dead(:,temp_hess(l,j,k))
            !             an_cr%hess(j,k)%neg_alive(:,l)=temp_an_cr%neg_alive(:,temp_hess(l,j,k))
            !             an_cr%hess(j,k)%neg_dead(:,l)=temp_an_cr%neg_dead(:,temp_hess(l,j,k))
            !         end do
            !     end do
            ! end do 

            ! allocate(an_cr%hcnt(0:max,norb,norb))
            ! do j=1, norb 
            !     do k=1,norb 
            !         an_cr%hcnt(0:,j,k)=temp_hess(0:max,j,k)
            !     end do
            ! end do 
            ! deallocate(temp_hess,stat=ierr)
            ! if (ierr/=0) then
            !     write(0,"(a,i0)") "Error in temp operators deallocation. ierr had value ", ierr
            !     errorflag=1
            !     return
            ! end if

            call dealloc_oprts_2(temp_an_cr)
        
        end if
        
        
        return
    
    
    end subroutine one_electrons
    
    END MODULE electrons