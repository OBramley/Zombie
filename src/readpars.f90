MODULE readpars

    use mod_types
    use globvars
    use dnad

    contains

    subroutine readrunconds   !   Level 1 Subroutine
        implicit none
        character(LEN=100)::LINE1, LINE2, LINE3, LINE4, LINE5, LINE6, LINE7, LINE8,LINE9, LINE10
        character(LEN=100):: LINE11,LINE12, LINE13, LINE14, LINE15, LINE16, LINE17
        integer::ierr, n
        
        if (errorflag .ne. 0) return

        open(unit=140,file='rundata.csv',status='old',iostat=ierr)

        if (ierr.ne.0) then
          write(stderr,"(a)") 'Error in opening rundata.csv file'
          errorflag = 1
          return
        end if

        read(140,*,iostat=ierr)LINE1, LINE2, LINE3, LINE4, LINE5, LINE6, LINE7, LINE8,LINE9, LINE10
        read(140,*,iostat=ierr)LINE11, LINE12, LINE13,LINE14, LINE15, LINE16, LINE17
        if (ierr.ne.0) then
            write(stderr,"(a,i0)") "Error reading rundata.csv of input file",ierr
            errorflag = 1
            return
        end if
        close(140)
      
        ! print*,LINE1, LINE2, LINE3, LINE4, LINE5, LINE6, LINE7, LINE8, LINE17, LINE18
        ! print*,LINE9, LINE10, LINE11,LINE12, LINE13, LINE14, LINE15, LINE16
        n=0        
        if((LINE1(1:1).eq.'y').or.(LINE1(1:1).eq.'Y')) then
            zomgflg="y"
        else if((LINE1(1:1).eq.'n').or.(LINE1(1:1).eq.'N')) then
            zomgflg="n"
        else
            write(stderr,"(a,a)") "Error. zomflg value must be YES/NO. Read ", trim(LINE1)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
            hamgflg="y"
        else if((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
            hamgflg="n"
        else
            write(stderr,"(a,a)") "Error. hamflg value must be YES/NO. Read ", trim(LINE2)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE3(1:1)=='y').or.(LINE3(1:1).eq.'Y')) then
            propflg="y"
        else if((LINE3(1:1)=='n').or.(LINE3(1:1).eq.'N')) then
            propflg="n"
        else
            write(stderr,"(a,a)") "Error. imaginary time flag must be YES/NO. Read ", trim(LINE3)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE4,*,iostat=ierr)beta
        if(ierr/=0) then
            write(stderr,"(a,a)") "Error reading beta Read ", trim(LINE4)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE5,*,iostat=ierr)timesteps
        if(ierr/=0) then
            write(stderr,"(a,a)") "Error reading timesteps. Read ", trim(LINE5)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE6(1:1)=='y').or.(LINE6(1:1).eq.'Y')) then
            cleanflg="y"
        else if((LINE6(1:1)=='n').or.(LINE6(1:1).eq.'N')) then
            cleanflg="n"
        else if((LINE6(1:1)=='f').or.(LINE6(1:1).eq.'F')) then
            cleanflg="f"
        else
            write(stderr,"(a,a)") "Error. cleaning flag must be YES/NO/f. Read ", trim(LINE6)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE7(1:1)=='y').or.(LINE7(1:1).eq.'Y')) then
            gramflg="y"
        else if((LINE7(1:1)=='n').or.(LINE7(1:1).eq.'N')) then
            gramflg="n"
        else
            write(stderr,"(a,a)") "Error. Gram Schmidt flag must be YES/NO. Read ", trim(LINE7)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE8,*,iostat=ierr)gramnum
        if(ierr/=0) then
            write(stderr,"(a,a)") "Error reading number of GS states. Read ", trim(LINE8)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE9(1:1)=='y').or.(LINE9(1:1).eq.'Y')) then
            GDflg="y"
        else if((LINE9(1:1)=='n').or.(LINE9(1:1).eq.'N')) then
            GDflg="n"
        else
            write(stderr,"(a,a)") "Error. GDflg flag must be YES/NO. Read ", trim(LINE9)
            errorflag=1
            return
        end if
        n=n+1 
        if((LINE10(1:1)=='y').or.(LINE10(1:1).eq.'Y')) then
            rstrtflg="y"
        else if((LINE10(1:1)=='n').or.(LINE10(1:1).eq.'N')) then
            rstrtflg="n"
        else
            write(stderr,"(a,a)") "Error. Restart flag must be YES/NO. Read ", trim(LINE10)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE11,*,iostat=ierr)norb
        if(ierr/=0) then
            write(stderr,"(a,a)") "Error reading number of orbitals. Read ", trim(LINE11)
            errorflag=1
            return
        end if
        norb=norb*2
        n=n+1
        read(LINE12,*,iostat=ierr)nel
        if(ierr/=0) then
            write(stderr,"(a,a)") "Error reading number of electrons. Read ", trim(LINE12)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE13,*,iostat=ierr)spin
        if(ierr/=0) then
            write(stderr,"(a,a)") "Error reading spin. Read ", trim(LINE13)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE14,*,iostat=ierr)ndet
        if(ierr/=0) then
            write(stderr,"(a,a)") "Error reading number of zombie states. Read ", trim(LINE14)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE15(1:1)=='r').or.(LINE15(1:1).eq.'R')) then
            zst="RN"
        else if((LINE15(1:1)=='h').or.(LINE15(1:1).eq.'H')) then
            zst="HF"
        else if((LINE15(1:1)=='b').or.(LINE15(1:1).eq.'B')) then
            zst="BB"
        else
            write(stderr,"(a,a)") "Error. Zombie state type must be hf/ran/bb. Read ", trim(LINE15)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE16(1:1)=='y').or.(LINE16(1:1).eq.'Y')) then
            rhf_1="y"
        else if((LINE16(1:1)=='n').or.(LINE16(1:1).eq.'N')) then
            rhf_1="n"
        else
            write(stderr,"(a,a)") "Error. rhf_1 flag must be YES/NO. Read ", trim(LINE16)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE17(1:1)=='y').or.(LINE17(1:1).eq.'Y')) then
            imagflg="y"
        else if((LINE17(1:1)=='n').or.(LINE17(1:1).eq.'N')) then
            imagflg="n"
        else
            write(stderr,"(a,a)") "Error. imagflg flag must be YES/NO. Read ", trim(LINE17)
            errorflag=1
            return
        end if
        n=n+1

        if (n.ne.17) then
            write(stderr,"(a)") "Not all required variables read in readrunconds subroutine"
            write(stderr,"(a,i0,a)") "Read a total of ", n, "of an expected 17 parameters"
            errorflag = 1
            return
          end if

          return
        

    end subroutine readrunconds

    subroutine read_zombie(zstore,gst)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        real(wp),dimension(norb)::phi,img
        real(wp),dimension(norb*2)::ccos,csin
        integer,intent(in)::gst
        character(len=4)::num
        character(len=2)::gst_num
        integer::ierr,j,k,l,zomnum,lines,gram_st
        character(LEN=20)::filenm
        ierr=0
        lines=0
        
     
        if(imagflg=='n') then
            do j=1, ndet
                if(GDflg.eq.'y')then
                    write(num,"(i4.4)")(j+1000)
                    ! write(num,"(i4.4)")j
                else
                    write(num,"(i4.4)")j
                end if
                if(gramflg=='n')then
                    filenm="data/zombie_"//trim(num)//".csv"
                else
                    write(gst_num,"(i2.1)")gst
                    filenm="data/zom_"//trim(gst_num)//"_"//trim(num)//".csv"
                end if
                ! filenm="data/zombie_"//trim(num)//".csv"
                zomnum=500+j
                
                ! if(rstrtflg.eq.'y')then
                open(unit=zomnum,file=trim(filenm),status="old",iostat=ierr)
                if(ierr/=0)then
                    write(stderr,"(a,i0)") "Error in opening zombie state file to read in. ierr had value ", ierr
                    errorflag=1
                    return
                end if
                lines=0
                do 
                    read(zomnum,*,iostat=ierr)
                    
                    if(ierr<0)then
                        close(zomnum)
                        exit
                    else if (ierr/=0) then
                        write(stderr,"(a,i0)") "Error in counting zombie rows. ierr had value ", ierr
                        errorflag=1
                        return
                    end if
                    lines=lines+1
                end do
                close(zomnum) 
    
                ierr=0
                open(unit=zomnum,file=trim(filenm),status="old",iostat=ierr)
                if(ierr/=0)then
                    write(stderr,"(a,i0)") "Error in opening zombie state file to read in. ierr had value ", ierr
                    errorflag=1
                    return
                end if
              
                if(lines.gt.3)then
                    do k=1,lines-3
                        read(zomnum,*)
                    end do
                end if 
       

                read(zomnum,*) phi
                zstore(j)%phi=phi
                zstore(j)%val(1:norb)=sin(zstore(j)%phi)
                zstore(j)%val(norb+1:2*norb)=cos(zstore(j)%phi)
                
               
                close(zomnum)
               
            end do
        else if(imagflg=='y')then
            ! do j=1, ndet
            !     write(num,"(i4.4)")j
            !     filenm="data/zombie_"//trim(num)//".csv"
            !     zomnum=500+j
            !     open(unit=zomnum,file=trim(filenm),status="old",iostat=ierr)
            !     if(ierr/=0)then
            !         write(stderr,"(a,i0)") "Error in opening zombie state file to read in. ierr had value ", ierr
            !         errorflag=1
            !         return
            !     end if

            !     read(zomnum,*) phi
            !     read(zomnum,*) img
            !     read(zomnum,*) ccos
            !     read(zomnum,*) csin
            !     do k=1,norb
            !         zstore(j)%phi(k)=phi(k)
            !         ! zstore(j)%img(k)=img(k)
            !     end do
            !     do k=1,(norb*2),2
            !         ! zstore(j)%cos((k+1)/2)=cmplx(ccos(k),ccos(k+1),wp)
            !         ! zstore(j)%sin((k+1)/2)=cmplx(csin(k),csin(k+1),wp)
            !     end do
            !     close(zomnum)
            ! end do
        end if
        
        
        write(6,"(a)") "Zombie states succeffuly read in"
        return
        
    end subroutine read_zombie

    subroutine read_zombie_c(cstore,clean_ndet)

        implicit none
        type(zombiest),dimension(:),intent(inout)::cstore
        integer,intent(in)::clean_ndet
        integer,dimension(norb)::dead,alive
        ! real(wp),dimension(norb*2)::cdead,calive
        character(len=6)::num
        integer::ierr,j,k,zomnum
        character(LEN=28)::filenm

        ierr=0
        ! if(imagflg=='n') then
            do j=1, clean_ndet
                write(num,"(i6.6)")j
                filenm="data/clean_zombie_"//trim(num)//".csv"
                zomnum=500+j
                open(unit=zomnum,file=trim(filenm),status="old",iostat=ierr)
                if(ierr/=0)then
                    write(stderr,"(a,i0)") "Error in opening zombie state file to read in. ierr had value ", ierr
                    errorflag=1
                    return
                end if
                read(zomnum,*) dead
                read(zomnum,*) alive
                close(zomnum)
                cstore(j)%val(norb+1:2*norb)=dead
                cstore(j)%val(1:norb)=alive
                do k=1, norb
                    if(cstore(j)%val(k).eq.1)then
                        cstore(j)%phi(k)=0.5*pirl
                    else
                        cstore(j)%phi(k)=0
                    end if 
                end do
            end do
        ! else if(imagflg=='y')then
            ! do j=1, clean_ndet
            !     write(num,"(i6.6)")j
            !     filenm="data/clean_zombie_"//trim(num)//".csv"
            !     zomnum=500+j
            !     open(unit=zomnum,file=trim(filenm),status="old",iostat=ierr)
            !     if(ierr/=0)then
            !         write(stderr,"(a,i0)") "Error in opening zombie state file to read in. ierr had value ", ierr
            !         errorflag=1
            !         return
            !     end if

            !     read(zomnum,*) cdead
            !     read(zomnum,*) calive
            !     do k=1,(norb*2),2

            !         cstore(j)%dead((k+1)/2)=cmplx(cdead(k),cdead(k+1),wp)
            !         cstore(j)%alive((k+1)/2)=cmplx(calive(k),calive(k+1),wp)
            !     end do
            !     close(zomnum)
            ! end do
        ! end if

        write(6,"(a)") "Cleaning Zombie states succeffuly read in"
        return
        
    end subroutine read_zombie_c

    subroutine read_ham(ham,size)

        implicit none
        type(hamiltonian),intent(inout)::ham
        integer,intent(in)::size
        integer::ierr,j,k
        REAL(wp),dimension(size)::line
        REAL(wp),dimension(size*2)::cline
        character(LEN=100)::hamnm,ovrlpnm
        integer, allocatable,dimension(:)::IPIV1
        real(wp),allocatable,dimension(:)::WORK1
        ! complex(wp),allocatable,dimension(:)::WORK1

        if (errorflag .ne. 0) return
        ierr=0

        open(unit=140,file='rundata.csv',status='old',iostat=ierr)

        if (ierr.ne.0) then
          write(stderr,"(a)") 'Error in opening rundata.csv file'
          errorflag = 1
          return
        end if

        read(140,*)
        read(140,*)
        read(140,*,iostat=ierr)hamnm, ovrlpnm
        if (ierr.ne.0) then
            write(stderr,"(a)") "Error reading rundata.csv of input file"
            errorflag = 1
            return
        end if
        close(140)
      

        if(imagflg=='n') then
            open(unit=200,file='data/'//trim(hamnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening hamiltonian file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(200,*) line
                do k=1, size
                    ham%hjk(j,k)=line(k)
                end do
            end do
            close(200)

            open(unit=201,file='data/'//trim(ovrlpnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening overlap file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(201,*) line
                do k=1, size
                    ham%ovrlp(j,k)=line(k)
                end do
            end do
            close(201)
        else if (imagflg=='y')then
            open(unit=200,file='data/'//trim(hamnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening hamiltonian file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(200,*) cline
                do k=1, (size*2),2
                end do
            end do
            close(200)

            open(unit=201,file='data/'//trim(ovrlpnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening overlap file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(201,*) line
                do k=1, (size*2),2
                    ! ham%ovrlp(j,(k+1)/2)=cmplx(cline(k),cline(k+1),wp)
                end do
            end do
            close(201)
        end if

        ham%inv=ham%ovrlp

        allocate(IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        
        allocate(WORK1(size),stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if   

        call ZGETRF(size,size,ham%inv,size,IPIV1,ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)")"Error in ZGETRF",ierr
        end if
        call ZGETRI(size,ham%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)")"Error in ZGETRF",ierr
        end if

        deallocate(IPIV1,stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if

        deallocate(WORK1,stat=ierr)
        if (ierr/=0) then
            write(stderr,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if
        
        ham%kinvh=matmul(ham%inv,ham%hjk)


    end subroutine read_ham

    subroutine dvec_read(d,size,p,filenm)

        implicit none
        real(wp),dimension(:),intent(inout)::d
        character(LEN=13),intent(in)::filenm
        integer,intent(in)::size,p
        REAL(wp),dimension(size)::line
        REAL(wp),dimension(size*2)::cline
        integer::ierr,j,vec
        if (errorflag .ne. 0) return

        ierr=0

        vec=900+p
        open(unit=vec,file='data/'//trim(filenm),status="old",iostat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in opening dvector file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(imagflg=='n')then
            read(vec,*) line
            do j=1, size
                ! d(j)=cmplx(line(j),0.0,wp)
                d(j)=line(j)
            end do
        else if(imagflg=='y')then
            read(vec,*) cline
            do j=1, (size*2),2
                ! d(j)=cmplx(cline(j),cline(j+1),wp)
            end do
        end if
        close(vec)
        return

    end subroutine dvec_read

    subroutine read_ham_c(ham,size)

        implicit none
        type(hamiltonian),intent(inout)::ham
        integer,intent(in)::size
        integer::ierr,j,k
        REAL(wp),dimension(size)::line
        REAL(wp),dimension(size*2)::cline
        character(LEN=100)::a,b,hamnm
       

        if (errorflag .ne. 0) return
        ierr=0

        open(unit=140,file='rundata.csv',status='old',iostat=ierr)

        if (ierr.ne.0) then
          write(stderr,"(a)") 'Error in opening rundata.csv file'
          errorflag = 1
          return
        end if

        read(140,*)
        read(140,*)
        read(140,*,iostat=ierr)a,b,hamnm
        if (ierr.ne.0) then
            write(stderr,"(a)") "Error reading rundata.csv of input file"
            errorflag = 1
            return
        end if
        close(140)
      
        
        if(imagflg=='n') then
            open(unit=204,file='data/'//trim(hamnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening hamiltonian file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(204,*) line
                do k=1, size
                    ! ham%hjk(j,k)=cmplx(line(k),0.0,wp)
                    ham%hjk(j,k)=line(k)
                end do
            end do
            close(204)

        else if (imagflg=='y')then
            open(unit=204,file='data/'//trim(hamnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening hamiltonian file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(200,*) cline
                do k=1, (size*2),2
                    ! ham%hjk(j,(k+1)/2)=cmplx(cline(k),cline(k+1),wp)
                end do
            end do
            close(204)
        end if

    end subroutine read_ham_c


    subroutine restart_chk()
        implicit none

        logical::file_exists


        inquire(file="data/zombie_1001.csv",exist=file_exists)
        if(file_exists.eqv..false.)then 
            return
        else 
            zomgflg ='n'
        end if 
        
        inquire(file="epoc.csv",exist=file_exists)
        if(file_exists.eqv..True.)then
            rstrtflg='y'
        end if 

        return

    end subroutine


END MODULE readpars
