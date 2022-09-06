MODULE outputs
    use globvars

    contains

    subroutine matrixwriter_real(out,size,filenm)

        implicit none

        real(kind=8),dimension(:,:),intent(in)::out
        integer,intent(in)::size
        character(LEN=*),intent(in)::filenm
        integer::ierr,j,k

        if (errorflag .ne. 0) return
        
        ierr=0

        open(unit=200,file=filenm,status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening matrix file. ierr had value ", ierr
            errorflag=1
            return
        end if

        do j=1, size
            write(200,'(*(e23.15e3 :", "))') (out(j,k),k=1,size)
        end do

        close(200)

        return

    end subroutine matrixwriter_real



    subroutine matrixwriter(out,size,filenm)

        implicit none

        complex(kind=8),dimension(:,:),intent(in)::out
        integer,intent(in)::size
        character(LEN=*),intent(in)::filenm
        integer::ierr,j,k

        if (errorflag .ne. 0) return
        
        ierr=0

        open(unit=200,file=filenm,status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening matrix file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(imagflg=='n')then
            do j=1, size
                write(200,'(*(e25.17e3 :", "))') (REAL(out(j,k)),k=1,size)
            end do
        else if (imagflg=='y')then
            do j=1, size
                write(200,'(*(e25.17e3 :", "))') ((out(j,k)),k=1,size)
            end do
        end if

        close(200)

        return

    end subroutine matrixwriter

    subroutine zombiewriter(zom,num)

        implicit none

        type(zombiest),intent(in)::zom 
        integer,intent(in)::num
        character(LEN=20)::filenm
        integer::ierr,zomnum,j
        character(LEN=4)::nums

        if (errorflag .ne. 0) return

        ierr=0

        write(nums,"(i4.4)")num

        
        filenm = "data/zombie_"//trim(nums)//".csv"

        zomnum=300+num
        
        open(unit=zomnum,file=trim(filenm),status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening zombie state file. ierr had value ", ierr
            errorflag=1
            return
        end if

        if(imagflg=='n') then
            write(zomnum,'(*(e25.17e3 :", "))') (REAL(zom%dead(j)),j=1,norb)
            write(zomnum,'(*(e25.17e3 :", "))') (REAL(zom%alive(j)),j=1,norb)
        else if(imagflg=='y') then
            write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%dead(j)),j=1,norb)
            write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%alive(j)),j=1,norb)
        end if

        close(zomnum)

        return

    end subroutine zombiewriter

    subroutine zombiewriter_c(zom,num)

        implicit none

        type(zombiest),intent(in)::zom 
        integer,intent(in)::num
        character(LEN=28)::filenm
        integer::ierr,zomnum,j
        character(LEN=6)::nums

        if (errorflag .ne. 0) return

        ierr=0

        write(nums,"(i6.6)")num

        
        filenm = "data/clean_zombie_"//trim(nums)//".csv"

        zomnum=300+num
        
        open(unit=zomnum,file=trim(filenm),status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening zombie state file. ierr had value ", ierr
            errorflag=1
            return
        end if

        if(imagflg=='n') then
            write(zomnum,'(*(e25.17e3 :", "))') (REAL(zom%dead(j)),j=1,norb)
            write(zomnum,'(*(e25.17e3 :", "))') (REAL(zom%alive(j)),j=1,norb)
        else if(imagflg=='y') then
            write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%dead(j)),j=1,norb)
            write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%alive(j)),j=1,norb)
        end if

        close(zomnum)

        return

    end subroutine zombiewriter_c

    subroutine energywriter(time,erg,filenm,j)

        implicit none

        real(kind=8), dimension(:),intent(in)::time
        complex(kind=8), dimension(:),intent(in)::erg 
        character(LEN=*),intent(in)::filenm
        integer,intent(in)::j
        integer::ergnum,ierr,k

        if (errorflag .ne. 0) return

        ergnum=400+j
        ierr=0
        open(unit=ergnum,file=trim(filenm),status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening energy output file. ierr had value ", ierr
            errorflag=1
            return
        end if

     
        write(ergnum,'(*(e25.17e3 :", "))') (time(k),k=1,timesteps+1)
        write(ergnum,'(*(e25.17e3 :", "))') (REAL(erg(k)),k=1,timesteps+1)
        if(imagflg=='y')then
            write(ergnum,'(*(e25.17e3 :", "))') (CMPLX(erg(k)),k=1,timesteps+1)
        end if

        close(ergnum)

        return

    end subroutine energywriter

    subroutine dvec_writer(d,size,p)

        implicit none
        complex(kind=8),dimension(:),intent(in)::d
        character(LEN=4)::stateno
        integer,intent(in)::size,p
        integer::ierr,j,vec
        if (errorflag .ne. 0) return
        
        ierr=0

        write(stateno,"(i4.4)")p

        vec=900+p
        open(unit=vec,file="data/dvec_"//trim(stateno)//".csv",status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening dvector file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(imagflg=='n')then
            write(vec,'(*(e25.17e3 :", "))') (REAL(d(j)),j=1,size)
        else if(imagflg=='y')then
            write(vec,'(*(1x,es25.17e3 :", "))') ((d(j)),j=1,size*2)
        end if
        close(vec)
        return

    end subroutine dvec_writer

    subroutine dvec_writer_c(d,size,p)

        implicit none
        complex(kind=8),dimension(:),intent(in)::d
        character(LEN=4)::stateno
        integer,intent(in)::size,p
        integer::ierr,j,vec
        if (errorflag .ne. 0) return
        
        ierr=0

        write(stateno,"(i4.4)")p

        vec=900+p
        open(unit=vec,file="data/clean_dvec_"//trim(stateno)//".csv",status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening dvector file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(imagflg=='n')then
            write(vec,'(*(e25.17e3 :", "))') (REAL(d(j)),j=1,size)
        else if(imagflg=='y')then
            write(vec,'(*(1x,es25.17e3 :", "))') ((d(j)),j=1,size*2)
        end if
        close(vec)
        return

    end subroutine dvec_writer_c

END MODULE outputs