MODULE readpars

    use globvars

    contains

    subroutine readrunconds   !   Level 1 Subroutine
        implicit none
        character(LEN=100)::LINE1, LINE2, LINE3, LINE4, LINE5, LINE6, LINE7, LINE8
        character(LEN=100)::LINE9, LINE10, LINE11,LINE12, LINE13, LINE14, LINE15, LINE16
        integer::ierr, n
        
        if (errorflag .ne. 0) return

        open(unit=140,file='rundata.csv',status='old',iostat=ierr)

        if (ierr.ne.0) then
          write(0,"(a)") 'Error in opening rundata.csv file'
          errorflag = 1
          return
        end if

        read(140,*,iostat=ierr)LINE1, LINE2, LINE3, LINE4, LINE5, LINE6, LINE7, LINE8
        read(140,*,iostat=ierr)LINE9, LINE10, LINE11,LINE12, LINE13, LINE14, LINE15, LINE16
        if (ierr.ne.0) then
            write(0,"(a)") "Error reading rundata.csv of input file"
            errorflag = 1
            return
        end if
        close(140)
      

        n=0        
        if((LINE1(1:1).eq.'y').or.(LINE1(1:1).eq.'Y')) then
            zomgflg="y"
        else if((LINE1(1:1).eq.'n').or.(LINE1(1:1).eq.'N')) then
            zomgflg="n"
        else
            write(0,"(a,a)") "Error. zomflg value must be YES/NO. Read ", trim(LINE1)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
            hamgflg="y"
        else if((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
            hamgflg="n"
        else
            write(0,"(a,a)") "Error. hamflg value must be YES/NO. Read ", trim(LINE2)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE3(1:1)=='y').or.(LINE3(1:1).eq.'Y')) then
            propflg="y"
        else if((LINE3(1:1)=='n').or.(LINE3(1:1).eq.'N')) then
            propflg="n"
        else
            write(0,"(a,a)") "Error. imaginary time flag must be YES/NO. Read ", trim(LINE3)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE4,*,iostat=ierr)beta
        if(ierr/=0) then
            write(0,"(a,a)") "Error reading beta Read ", trim(LINE4)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE5,*,iostat=ierr)timesteps
        if(ierr/=0) then
            write(0,"(a,a)") "Error reading timesteps. Read ", trim(LINE5)
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
            write(0,"(a,a)") "Error. cleaning flag must be YES/NO/f. Read ", trim(LINE6)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE7(1:1)=='y').or.(LINE7(1:1).eq.'Y')) then
            gramflg="y"
        else if((LINE7(1:1)=='n').or.(LINE7(1:1).eq.'N')) then
            gramflg="n"
        else
            write(0,"(a,a)") "Error. Gram Schmidt flag must be YES/NO. Read ", trim(LINE7)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE8,*,iostat=ierr)gramnum
        if(ierr/=0) then
            write(0,"(a,a)") "Error reading number of GS states. Read ", trim(LINE8)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE9,*,iostat=ierr)norb
        if(ierr/=0) then
            write(0,"(a,a)") "Error reading number of orbitals. Read ", trim(LINE9)
            errorflag=1
            return
        end if
        norb=norb*2
        n=n+1
        read(LINE10,*,iostat=ierr)nel
        if(ierr/=0) then
            write(0,"(a,a)") "Error reading number of electrons. Read ", trim(LINE10)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE11,*,iostat=ierr)spin
        if(ierr/=0) then
            write(0,"(a,a)") "Error reading spin. Read ", trim(LINE11)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE12,*,iostat=ierr)ndet
        if(ierr/=0) then
            write(0,"(a,a)") "Error reading number of zombie states. Read ", trim(LINE12)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE13(1:1)=='r').or.(LINE13(1:1).eq.'R')) then
            zst="RN"
        else if((LINE13(1:1)=='h').or.(LINE13(1:1).eq.'H')) then
            zst="HF"
        else if((LINE13(1:1)=='b').or.(LINE13(1:1).eq.'B')) then
            zst="BB"
        else
            write(0,"(a,a)") "Error. Zombie state type must be hf/ran/bb. Read ", trim(LINE13)
            errorflag=1
            return
        end if
        n=n+1
        read(LINE14,*,iostat=ierr)bb_improv
        if(ierr/=0) then
            write(0,"(a,a)") "Error reading number ofbb_improv. Read ", trim(LINE14)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE15(1:1)=='y').or.(LINE15(1:1).eq.'Y')) then
            rhf_1="y"
        else if((LINE15(1:1)=='n').or.(LINE15(1:1).eq.'N')) then
                rhf_1="n"
        else
            write(0,"(a,a)") "Error. rhf_1 flag must be YES/NO. Read ", trim(LINE15)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE16(1:1)=='y').or.(LINE16(1:1).eq.'Y')) then
            imagflg="y"
        else if((LINE16(1:1)=='n').or.(LINE16(1:1).eq.'N')) then
            imagflg="n"
        else
            write(0,"(a,a)") "Error. imagflg flag must be YES/NO. Read ", trim(LINE16)
            errorflag=1
            return
        end if
        n=n+1

        if (n.ne.16) then
            write(0,"(a)") "Not all required variables read in readrunconds subroutine"
            write(0,"(a,i0,a)") "Read a total of ", n, "of an expected 16 parameters"
            errorflag = 1
            return
          end if

          return
        

    end subroutine readrunconds

    subroutine read_zombie(zstore)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        real(kind=8),dimension(norb)::dead,alive
        real(kind=8),dimension(norb*2)::cdead,calive
        character(len=4)::num
        integer::ierr,j,k,zomnum
        character(LEN=20)::filenm

        ierr=0

        if(imagflg=='n') then
            do j=1, ndet
                write(num,"(i4.4)")j
                filenm="data/zombie_"//trim(num)//".csv"
                zomnum=500+j
                open(unit=zomnum,file=trim(filenm),status="old",iostat=ierr)
                if(ierr/=0)then
                    write(0,"(a,i0)") "Error in opening zombie state file to read in. ierr had value ", ierr
                    errorflag=1
                    return
                end if

                read(zomnum,*) dead
                read(zomnum,*) alive
                do k=1,norb
                    zstore(j)%dead(k)=cmplx(dead(k),0.0,kind=8)
                    zstore(j)%alive(k)=cmplx(alive(k),0.0,kind=8)
                end do 
                close(zomnum)
            end do
        else if(imagflg=='y')then
            do j=1, ndet
                write(num,"(i4.4)")j
                filenm="data/zombie_"//trim(num)//".csv"
                zomnum=500+j
                open(unit=zomnum,file=trim(filenm),status="old",iostat=ierr)
                if(ierr/=0)then
                    write(0,"(a,i0)") "Error in opening zombie state file to read in. ierr had value ", ierr
                    errorflag=1
                    return
                end if

                read(zomnum,*) cdead
                read(zomnum,*) calive
                do k=1,(norb*2),2
                    zstore(j)%dead((k+1)/2)=cmplx(cdead(k),cdead(k+1),kind=8)
                    zstore(j)%alive((k+1)/2)=cmplx(calive(k),calive(k+1),kind=8)
                end do
                close(zomnum)
            end do
        end if
        write(6,"(a)") "Zombie states succeffuly read in"
        return
        
    end subroutine read_zombie

    subroutine read_zombie_c(zstore,clean_ndet)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore
        integer,intent(in)::clean_ndet
        real(kind=8),dimension(norb)::dead,alive
        real(kind=8),dimension(norb*2)::cdead,calive
        character(len=6)::num
        integer::ierr,j,k,zomnum
        character(LEN=28)::filenm

        ierr=0

        if(imagflg=='n') then
            do j=1, clean_ndet
                write(num,"(i6.6)")j
                filenm="data/clean_zombie_"//trim(num)//".csv"
                zomnum=500+j
                open(unit=zomnum,file=trim(filenm),status="old",iostat=ierr)
                if(ierr/=0)then
                    write(0,"(a,i0)") "Error in opening zombie state file to read in. ierr had value ", ierr
                    errorflag=1
                    return
                end if

                read(zomnum,*) dead
                read(zomnum,*) alive
                do k=1,norb
                    zstore(j)%dead(k)=cmplx(dead(k),0.0,kind=8)
                    zstore(j)%alive(k)=cmplx(alive(k),0.0,kind=8)
                end do 
                close(zomnum)
            end do
        else if(imagflg=='y')then
            do j=1, clean_ndet
                write(num,"(i6.6)")j
                filenm="data/clean_zombie_"//trim(num)//".csv"
                zomnum=500+j
                open(unit=zomnum,file=trim(filenm),status="old",iostat=ierr)
                if(ierr/=0)then
                    write(0,"(a,i0)") "Error in opening zombie state file to read in. ierr had value ", ierr
                    errorflag=1
                    return
                end if

                read(zomnum,*) cdead
                read(zomnum,*) calive
                do k=1,(norb*2),2
                    zstore(j)%dead((k+1)/2)=cmplx(cdead(k),cdead(k+1),kind=8)
                    zstore(j)%alive((k+1)/2)=cmplx(calive(k),calive(k+1),kind=8)
                end do
                close(zomnum)
            end do
        end if
        write(6,"(a)") "Zombie states succeffuly read in"
        return
        
    end subroutine read_zombie_c

    subroutine read_ham(ham,size)

        implicit none
        type(hamiltonian),intent(inout)::ham
        integer,intent(in)::size
        integer::ierr,j,k
        REAL(kind=8),dimension(size)::line
        REAL(kind=8),dimension(size*2)::cline
        character(LEN=100)::hamnm,ovrlpnm
        integer, allocatable,dimension(:)::IPIV1
        complex(kind=8),allocatable,dimension(:)::WORK1

        if (errorflag .ne. 0) return
        ierr=0

        open(unit=140,file='rundata.csv',status='old',iostat=ierr)

        if (ierr.ne.0) then
          write(0,"(a)") 'Error in opening rundata.csv file'
          errorflag = 1
          return
        end if

        read(140,*)
        read(140,*)
        read(140,*,iostat=ierr)hamnm, ovrlpnm
        if (ierr.ne.0) then
            write(0,"(a)") "Error reading rundata.csv of input file"
            errorflag = 1
            return
        end if
        close(140)
      

        if(imagflg=='n') then
            open(unit=200,file='data/'//trim(hamnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening hamiltonian file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(200,*) line
                do k=1, size
                    ham%hjk(j,k)=cmplx(line(k),0.0,kind=8)
                end do
            end do
            close(200)

            open(unit=201,file='data/'//trim(ovrlpnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening overlap file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(201,*) line
                do k=1, size
                    ham%ovrlp(j,k)=cmplx(line(k),0.0,kind=8)
                end do
            end do
            close(201)
        else if (imagflg=='y')then
            open(unit=200,file='data/'//trim(hamnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening hamiltonian file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(200,*) cline
                do k=1, (size*2),2
                    ham%hjk(j,(k+1)/2)=cmplx(cline(k),cline(k+1),kind=8)
                end do
            end do
            close(200)

            open(unit=201,file='data/'//trim(ovrlpnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening overlap file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(201,*) line
                do k=1, (size*2),2
                    ham%ovrlp(j,(k+1)/2)=cmplx(cline(k),cline(k+1),kind=8)
                end do
            end do
            close(201)
        end if

        ham%inv=ham%ovrlp

        allocate(IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if 
        
        allocate(WORK1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if   

        call ZGETRF(size,size,ham%inv,size,IPIV1,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if
        call ZGETRI(size,ham%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in ZGETRF",ierr
        end if

        deallocate(IPIV1,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if

        deallocate(WORK1,stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in WORK vector allocation . ierr had value ", ierr
            errorflag=1
            return
        end if
        
        ham%kinvh=matmul(ham%inv,ham%hjk)


    end subroutine read_ham

    subroutine dvec_read(d,size,p,filenm)

        implicit none
        complex(kind=8),dimension(:),intent(inout)::d
        character(LEN=13),intent(in)::filenm
        integer,intent(in)::size,p
        REAL(kind=8),dimension(size)::line
        REAL(kind=8),dimension(size*2)::cline
        integer::ierr,j,vec
        if (errorflag .ne. 0) return

        ierr=0

        vec=900+p
        open(unit=vec,file='data/'//trim(filenm),status="old",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening dvector file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(imagflg=='n')then
            read(vec,*) line
            do j=1, size
                d(j)=cmplx(line(j),0.0,kind=8)
            end do
        else if(imagflg=='y')then
            read(vec,*) cline
            do j=1, (size*2),2
                d(j)=cmplx(cline(j),cline(j+1),kind=8)
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
        REAL(kind=8),dimension(size)::line
        REAL(kind=8),dimension(size*2)::cline
        character(LEN=100)::a,b,hamnm
       

        if (errorflag .ne. 0) return
        ierr=0

        open(unit=140,file='rundata.csv',status='old',iostat=ierr)

        if (ierr.ne.0) then
          write(0,"(a)") 'Error in opening rundata.csv file'
          errorflag = 1
          return
        end if

        read(140,*)
        read(140,*)
        read(140,*,iostat=ierr)a,b,hamnm
        if (ierr.ne.0) then
            write(0,"(a)") "Error reading rundata.csv of input file"
            errorflag = 1
            return
        end if
        close(140)
      

        if(imagflg=='n') then
            open(unit=204,file='data/'//trim(hamnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening hamiltonian file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(200,*) line
                do k=1, size
                    ham%hjk(j,k)=cmplx(line(k),0.0,kind=8)
                end do
            end do
            close(204)

        else if (imagflg=='y')then
            open(unit=204,file='data/'//trim(hamnm),status="old",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening hamiltonian file. ierr had value ", ierr
                errorflag=1
                return
            end if

            do j=1, size
                read(200,*) cline
                do k=1, (size*2),2
                    ham%hjk(j,(k+1)/2)=cmplx(cline(k),cline(k+1),kind=8)
                end do
            end do
            close(204)
        end if

    end subroutine read_ham_c


END MODULE readpars
