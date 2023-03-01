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
        
        
        ierr=0
        
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
            cr1=int(read_in(l+1,2)) 
            cr2=int(read_in(l+1,3)) 
            an2=int(read_in(l+1,4)) 
            an1=int(read_in(l+1,5))
           
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
            do j=1,norb
                do l=1,e2
                    cr1=int(read_in(l+1,2))
                    cr2=int(read_in(l+1,3))
                    an2=int(read_in(l+1,4))
                    an1=int(read_in(l+1,5))
                    if((an2.eq.j).or.(cr2.eq.j).or.(an1.eq.j).or.(cr1.eq.j))then
                        an2_cr2%dcnt(0,j)=int(an2_cr2%dcnt(0,j)+1,kind=2)
                        an2_cr2%dcnt(an2_cr2%dcnt(0,j),j)=l
    
                        if((an2.eq.j).or.(an1.eq.j))then
                            if((cr1.eq.j).or.(cr2.eq.j))then
                                an2_cr2%alive_diff(j,:,l)=an2_cr2%alive(:,l)
                                an2_cr2%dead_diff(j,:,l)=an2_cr2%dead(:,l)
                                an2_cr2%neg_alive_diff(j,:,l)=an2_cr2%neg_alive(:,l)
                                an2_cr2%neg_dead_diff(j,:,l)=an2_cr2%neg_dead(:,l)
                                an2_cr2%alive_diff(j,j,l)=int(j+norb,kind=2)
                                an2_cr2%neg_alive_diff(j,j,l)=int(2*an2_cr2%neg_alive_diff(j,j,l),kind=1)
                            else
                                an2_cr2%alive_diff(j,:,l)=an2_cr2%alive(:,l)
                                an2_cr2%dead_diff(j,:,l)=an2_cr2%dead(:,l)
                                an2_cr2%neg_alive_diff(j,:,l)=an2_cr2%neg_alive(:,l)
                                an2_cr2%neg_dead_diff(j,:,l)=an2_cr2%neg_dead(:,l)
                                an2_cr2%alive_diff(j,j,l)=int(j,kind=2)
                                an2_cr2%dead_diff(j,j,l)=int(j+norb,kind=2)
                                an2_cr2%neg_alive_diff(j,j,l)=int(an2_cr2%neg_dead_diff(j,j,l)*(-1),kind=1)
                            end if
                        else
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

            do j=1,norb
                do k=j,norb 
                    do l=1,e2
                        cr1=int(read_in(l+1,2))
                        cr2=int(read_in(l+1,3))
                        an2=int(read_in(l+1,4))
                        an1=int(read_in(l+1,5))
                        if(((an2.eq.j).or.(cr2.eq.j).or.(an1.eq.j).or.(cr1.eq.j)))then
                            if((an2.eq.k).or.(cr2.eq.k).or.(an1.eq.k).or.(cr1.eq.k))then
                                ! print*,an1,an2,cr1,cr2,j,k
                                an2_cr2%hcnt(0,j,k)=int(an2_cr2%hcnt(0,j,k)+1,kind=2)
                                an2_cr2%hcnt(an2_cr2%hcnt(0,j,k),j,k)=l
                                an2_cr2%alive_hess(j,k,1:norb,l)=an2_cr2%alive(1:norb,l)
                                an2_cr2%dead_hess(j,k,1:norb,l)=an2_cr2%dead(1:norb,l)
                                an2_cr2%neg_alive_hess(j,k,1:norb,l)=an2_cr2%neg_alive(1:norb,l)
                                an2_cr2%neg_dead_hess(j,k,1:norb,l)=an2_cr2%neg_dead(1:norb,l)
                                if(j.ne.k)then
                                    if(an1.eq.j)then
                                        if(an2.eq.k)then
                                            an2_cr2%neg_dead_hess(j,k,j,l)=int(-4*an2_cr2%neg_dead_hess(j,k,j,l),kind=1)
                                            an2_cr2%neg_dead_hess(j,k,k,l)=int(-4*an2_cr2%neg_dead_hess(j,k,k,l),kind=1)
                                        else 
                                            an2_cr2%neg_dead_hess(j,k,j,l)=int(-4*an2_cr2%neg_dead_hess(j,k,j,l),kind=1)
                                            an2_cr2%neg_alive_hess(j,k,k,l)=int(-4*an2_cr2%neg_alive_hess(j,k,k,l),kind=1)
                                        end if
                                    else if(an2.eq.j)then
                                        an2_cr2%neg_dead_hess(j,k,j,l)=int(-4*an2_cr2%neg_dead_hess(j,k,j,l),kind=1)
                                        an2_cr2%neg_alive_hess(j,k,k,l)=int(-4*an2_cr2%neg_alive_hess(j,k,k,l),kind=1)
                                    else 
                                        an2_cr2%neg_alive_hess(j,k,j,l)=int(-4*an2_cr2%neg_alive_hess(j,k,j,l),kind=1)
                                        an2_cr2%neg_dead_hess(j,k,k,l)=int(-4*an2_cr2%neg_dead_hess(j,k,k,l),kind=1)
                                    end if  
                                else
                                    an2_cr2%alive_hess(j,k,j,l)=int(j,kind=2)
                                    an2_cr2%dead_hess(j,k,j,l)=int(j+norb,kind=2)
                                    an2_cr2%neg_alive_hess(j,k,j,l)=int(-2*an2_cr2%neg_alive_hess(j,k,j,l),kind=1)
                                    an2_cr2%neg_dead_hess(j,k,j,l)=int(2*an2_cr2%neg_alive_hess(j,k,j,l),kind=1)  
                                end if
                            end if
                        end if
                    end do 
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
    
        
        ierr=0
    
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
                        an_cr%dcnt(an_cr%dcnt(0,j),j)=l
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
                    do k=j,norb 
                        if(((an.eq.j).or.(cr.eq.j)).and.((an.eq.k).or.(cr.eq.k)))then
                            an_cr%hcnt(0,j,k)=int(an_cr%hcnt(0,j,k)+1,kind=2)
                            an_cr%hcnt(an_cr%hcnt(0,j,k),j,k)=l
                            if(an.eq.cr)then
                                an_cr%alive_hess(j,k,:,l)=an_cr%alive(:,l)
                                an_cr%dead_hess(j,k,:,l)=an_cr%dead(:,l)
                                an_cr%neg_alive_hess(j,k,:,l)=an_cr%neg_alive(:,l)
                                an_cr%neg_dead_hess(j,k,:,l)=an_cr%neg_dead(:,l)
                                an_cr%alive_hess(j,k,j,l)=int(j,kind=2)
                                an_cr%dead_hess(j,k,j,l)=int(j+norb,kind=2)
                                an_cr%neg_alive_hess(j,k,j,l)=int(-2*an_cr%neg_alive_hess(j,k,j,l),kind=1)
                                an_cr%neg_dead_hess(j,k,j,l)=int(2*an_cr%neg_alive_hess(j,k,j,l),kind=1) 
                            else
                                an_cr%alive_hess(j,k,:,l)=an_cr%alive(:,l)
                                an_cr%dead_hess(j,k,:,l)=an_cr%dead(:,l)
                                an_cr%neg_alive_hess(j,k,:,l)=an_cr%neg_alive(:,l)
                                an_cr%neg_dead_hess(j,k,:,l)=an_cr%neg_dead(:,l)
                                an_cr%neg_dead_hess(j,k,an,l)=int(-4*an_cr%neg_dead_hess(j,k,an,l),kind=1)
                                an_cr%neg_alive_hess(j,k,cr,l)=int(-4*an_cr%neg_alive_hess(j,k,cr,l),kind=1)
                            end if
                        end if
    
                    end do
    
    
    
    
                end do
            end do 
        
        end if
        
        
        return
    
    
    end subroutine one_electrons
    
    END MODULE electrons