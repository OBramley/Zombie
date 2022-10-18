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

    subroutine zombiewriter(zom,num,loop)

        implicit none

        type(zombiest),intent(in)::zom 
        integer,intent(in)::num,loop
        character(LEN=20)::filenm
        integer::ierr,zomnum,j,l
        character(LEN=4)::nums

        if (errorflag .ne. 0) return

        ierr=0
        
        if(GDflg.eq.'y')then
            write(nums,"(i4.4)")(num+1000)
        else
            write(nums,"(i4.4)")num
        end if

        
        filenm = "data/zombie_"//trim(nums)//".csv"

        zomnum=300+num
        
        open(unit=zomnum,file=trim(filenm),status="unknown",iostat=ierr)
        if(ierr/=0)then
            close(zomnum)
            write(0,"(a,i0)") "Error in opening zombie state file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(loop.ge.1)then
            if(imagflg=='y')then
                l=4
            else 
                l=3
            end if
            do j=1, l*(loop-1)
                read(zomnum,*)
            end do
        end if
        if(imagflg=='n') then
            write(zomnum,'(*(e25.17e3 :", "))') ((zom%phi(j)),j=1,norb)
            write(zomnum,'(*(e25.17e3 :", "))') (REAL(zom%cos(j)),j=1,norb)
            write(zomnum,'(*(e25.17e3 :", "))') (REAL(zom%sin(j)),j=1,norb)
        else if(imagflg=='y') then
            write(zomnum,'(*(e25.17e3 :", ": ))') ((zom%phi(j)),j=1,norb)
            write(zomnum,'(*(e25.17e3 :", ": ))') ((zom%img(j)),j=1,norb)
            write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%cos(j)),j=1,norb)
            write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%sin(j)),j=1,norb)
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
        
        open(unit=zomnum,file=trim(filenm),status="unknown",iostat=ierr)
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

    subroutine energywriter(time,erg,filenm,j,loop)

        implicit none

        real(kind=8), dimension(:),intent(in)::time
        complex(kind=8), dimension(:),intent(in)::erg 
        character(LEN=*),intent(in)::filenm
        integer,intent(in)::j,loop
        integer::ergnum,ierr,k,l

        if (errorflag .ne. 0) return

        ergnum=400+j
        ierr=0
        open(unit=ergnum,file=trim(filenm),status="unknown",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening energy output file. ierr had value ", ierr
            close(ergnum)
            errorflag=1
            return
        end if
        if(loop.gt.1)then
            if(imagflg=='y')then
                l=2
            else 
                l=1
            end if
          
            do k=1, l*(loop-1)
                read(ergnum,*)
            end do
        end if
        if(loop.eq.1)then
            write(ergnum,'(*(e25.17e3 :", "))') (time(k),k=1,timesteps+1)
        end if
        write(ergnum,'(*(e25.17e3 :", "))') (REAL(erg(k)),k=1,timesteps+1)
        if(imagflg=='y')then
            write(ergnum,'(*(e25.17e3 :", "))') (CMPLX(erg(k)),k=1,timesteps+1)
        end if

        close(ergnum)

        return

    end subroutine energywriter

    subroutine epoc_writer(erg,step,chng_trk)

        implicit none
        real(kind=8),intent(in)::erg 
        integer,intent(in)::step
        integer,dimension(:),intent(in)::chng_trk
        integer::epoc,ierr,k
    
        if (errorflag .ne. 0) return
        epoc=450
        ierr=0
        open(unit=epoc,file='epoc.csv',status="unknown",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(step.gt.0)then
            do k=1, step
                read(epoc,*)
            end do
        end if
        write(epoc,'(i0,",",e25.17e3,",",*(i0:", "))') step,erg,(chng_trk(k),k=2,ndet)
        close(epoc)
        return
        
    end subroutine epoc_writer


    subroutine dvec_writer(d,size,p,loop)

        implicit none
        complex(kind=8),dimension(:),intent(in)::d
        character(LEN=4)::stateno
        integer,intent(in)::size,p,loop
        integer::ierr,j,vec
        if (errorflag .ne. 0) return
        
        ierr=0

        write(stateno,"(i4.4)")p
        
        vec=900+p
        open(unit=vec,file="data/dvec_"//trim(stateno)//".csv",status="unknown",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening dvector file. ierr had value ", ierr
            errorflag=1
            close(vec)
            return
        end if
        if(loop.gt.1)then
            do j=1, (loop-1)
                read(vec,*)
            end do
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

    subroutine clean_erg_write(clean_ndet, clean_erg,clean_norm,j)

        implicit none
        integer, intent(in)::clean_ndet,j
        complex(kind=8),intent(in)::clean_erg,clean_norm
        integer::cleane,ierr

        cleane=400+j
        ierr=0

        open(unit=cleane,file='clean_energy.csv',status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening clean energy output file. ierr had value ", ierr
            errorflag=1
            return
        end if

        ! write(cleane,"(a)")"|    No. Cleaning states     |      clean energy     |     clean norm     |       energy/norm      |"
        if(imagflg=='n')then
            write(cleane,'(e25.17e3,",",e25.17e3,",",e25.17e3,",",i6)') &
                REAL(clean_erg), REAL(clean_norm), (REAL(clean_erg)/REAL(clean_norm)),clean_ndet
        end if

        close(cleane)

        return

    end subroutine clean_erg_write
        
    



END MODULE outputs