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
        close(140)
        if (ierr.ne.0) then
            write(0,"(a)") "Error reading first line of input file"
            errorflag = 1
            return
        end if

        
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
            propflg="y"
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
            cleanflg="y"
        else
            write(0,"(a,a)") "Error. cleaning flag must be YES/NO. Read ", trim(LINE6)
            errorflag=1
            return
        end if
        n=n+1
        if((LINE7(1:1)=='y').or.(LINE7(1:1).eq.'Y')) then
            gramflg="y"
        else if((LINE7(1:1)=='n').or.(LINE7(1:1).eq.'N')) then
            gramflg="y"
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


END MODULE readpars
