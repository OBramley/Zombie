MODULE outputs
    use globvars

    interface epoc_writer

        module procedure epoc_writer_int, epoc_writer_array,epoc_writer_array_orbital

    end interface epoc_writer

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

        ! complex(kind=8),dimension(:,:),intent(in)::out
        real(kind=8),dimension(:,:),intent(in)::out
        integer,intent(in)::size
        character(LEN=*),intent(in)::filenm
        integer::ierr,j,k
        logical :: file_exists
        if (errorflag .ne. 0) return
        
        ierr=0
        inquire(file=filenm,exist=file_exists)
        if(file_exists.eqv..false.) then
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
        else 
            write(6,"(a,a)") trim(filenm), " alread exists so not rewriting"
        end if 

        return

    end subroutine matrixwriter

    subroutine zombiewriter(zom,num,pass)

        implicit none

        type(zombiest),intent(in)::zom 
        integer,intent(in)::num,pass
        character(LEN=20)::filenm
        integer::ierr,zomnum,j
        character(LEN=4)::nums
        logical :: file_exists

        if (errorflag .ne. 0) return

        ierr=0
        
        if(GDflg.eq.'y')then
            write(nums,"(i4.4)")(num+1000)
        else
            write(nums,"(i4.4)")num
        end if

        
        filenm = "data/zombie_"//trim(nums)//".csv"

        inquire(file=filenm,exist=file_exists)
        zomnum=300+num
        if(file_exists.eqv..false.) then
            open(unit=zomnum,file=trim(filenm),status="new",iostat=ierr)
            if(ierr/=0)then
                close(zomnum)
                write(0,"(a,i0)") "Error in opening zombie state file. ierr had value ", ierr
                errorflag=1
                return
            end if
        
            if(imagflg=='n') then
                write(zomnum,'(*(e25.17e3 :", "))') ((zom%phi(j)),j=1,norb)
                ! write(zomnum,'(*(e25.17e3 :", "))') (REAL(zom%cos(j)),j=1,norb)
                ! write(zomnum,'(*(e25.17e3 :", "))') (REAL(zom%sin(j)),j=1,norb)
                write(zomnum,'(*(e25.17e3 :", "))') ((zom%cos(j)),j=1,norb)
                write(zomnum,'(*(e25.17e3 :", "))') ((zom%sin(j)),j=1,norb)
            else if(imagflg=='y') then
                write(zomnum,'(*(e25.17e3 :", ": ))') ((zom%phi(j)),j=1,norb)
                write(zomnum,'(*(e25.17e3 :", ": ))') ((zom%img(j)),j=1,norb)
                write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%cos(j)),j=1,norb)
                write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%sin(j)),j=1,norb)
            end if
            close(zomnum)
        else if(file_exists.eqv..true.) then
            open(unit=zomnum,file=trim(filenm),status="old",access='append',iostat=ierr)
            if(ierr/=0)then
                close(zomnum)
                write(0,"(a,i0)") "Error in opening zombie state file. ierr had value ", ierr
                errorflag=1
                return
            end if
           
            write(zomnum,*)' '
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
        end if
       

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
            write(zomnum,'(*(e25.17e3 :", "))') ((zom%dead(j)),j=1,norb)
            write(zomnum,'(*(e25.17e3 :", "))') ((zom%alive(j)),j=1,norb)
            ! write(zomnum,'(*(e25.17e3 :", "))') (REAL(zom%dead(j)),j=1,norb)
            ! write(zomnum,'(*(e25.17e3 :", "))') (REAL(zom%alive(j)),j=1,norb)
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
        real(kind=8), dimension(:),intent(in)::erg 
        character(LEN=*),intent(in)::filenm
        real::db
        integer,intent(in)::j
        integer::ergnum,ierr,k
        logical :: file_exists

        if (errorflag .ne. 0) return

        ergnum=400+j
        ierr=0
        db=beta/timesteps
      
        inquire(file=trim(filenm),exist=file_exists)
        if(file_exists.eqv..false.) then
            open(unit=ergnum,file=trim(filenm),status="new",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening energy file. ierr had value ", ierr
                errorflag=1
                close(ergnum)
                return
            end if
            write(ergnum,'(*(e25.17e3 :", "))') (db*(k-1),k=1,timesteps+1)
        else 
            open(unit=ergnum,file=trim(filenm),status="old",access='append',iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening energy file. ierr had value ", ierr
                errorflag=1
                close(ergnum)
                return
            end if
            write(ergnum,*)' '
        end if
        
        ! write(ergnum,'(*(e25.17e3 :", "))') (REAL(erg(k)),k=1,timesteps+1)
        write(ergnum,'(*(e25.17e3 :", "))') ((erg(k)),k=1,timesteps+1)
        if(imagflg=='y')then
            ! write(ergnum,'(*(e25.17e3 :", "))') (CMPLX(erg(k)),k=1,timesteps+1)
        end if

        close(ergnum)

        return

    end subroutine energywriter

    subroutine epoc_writer_int(erg,step,chng_trk,lr,pass)

        implicit none
        real(kind=8),intent(in)::erg,lr 
        integer,intent(in)::step,pass
        integer,intent(in)::chng_trk
        integer::epoc,ierr
        logical :: file_exists

        if (errorflag .ne. 0) return

        
        inquire(file='epoc.csv',exist=file_exists)
        epoc=450
        ierr=0
        if(file_exists.eqv..false.) then
            open(unit=epoc,file='epoc.csv',status="new",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
                errorflag=1
                return
            end if
            write(epoc,'(a,",",a,",",a,","a)') "EPOC", "Energy", "Learning rate", "Zombie state altered"
            write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",i0)') 0,erg,0.0,0
            close(epoc)
        else if(file_exists.eqv..true.) then
            open(unit=epoc,file='epoc.csv',status="old",access='append',iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
                errorflag=1
                return
            end if
            if(pass.eq.1)then
                write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",i0)') step,erg,lr,chng_trk
            else
                write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",i0)') step,erg,lr,chng_trk
            end if
            close(epoc)
        end if
        return
        
    end subroutine epoc_writer_int

    subroutine epoc_writer_array_orbital(erg,step,chng_trk,pass)

        implicit none
        real(kind=8),intent(in)::erg
        integer,intent(in)::step,pass
        integer,dimension(:),intent(in)::chng_trk
        integer::epoc,ierr,k

        if (errorflag .ne. 0) return

    
        epoc=450
        ierr=0
       
        open(unit=epoc,file='epoc.csv',status="old",access='append',iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(pass.eq.1)then
            write(epoc,'(i0,",",e25.17e3,",",a,",",*(i0:", "))') step,erg,"   ",(chng_trk(k),k=1,ndet-1)
            ! do k=1,ndet-1
            !     write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",*(i0:", "))') step,erg,(chng_trk(k),k=1,ndet-1)
            !     if(chng_trk(k).eq.0)then
            !         EXIT 
            !     end if 
            !     write(epoc,'(a,",",e25.17e3,",",e25.17e3,",",i0)') "   ",erg_dim(k),lr(k),chng_trk(k)
            ! end do
        else
            write(epoc,'(i0,",",e25.17e3,",",a,",",*(i0:", "))') step,erg,"   ",(chng_trk(k),k=1,ndet-1)
            ! do k=1,ndet-1
            !     write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",*(i0:", "))') step,erg,(chng_trk(k),k=1,ndet-1)
            !     if(chng_trk(k).eq.0)then
            !         EXIT 
            !     end if 
            !     write(epoc,'(a,",",e25.17e3,",",e25.17e3,",",i0)') "   ",erg_dim(k),lr(k),chng_trk(k)
            ! end do
            ! write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",*(i0:", "))') step,erg,lr,(chng_trk(k),k=1,ndet-1)
        end if
        close(epoc)
       
        return
        
    end subroutine epoc_writer_array_orbital

    subroutine epoc_writer_array(erg,step,chng_trk,erg_dim,lr,pass)

        implicit none
        real(kind=8),intent(in)::erg
        real(kind=8),dimension(:),intent(in)::lr,erg_dim 
        integer,intent(in)::step,pass
        integer,dimension(:),intent(in)::chng_trk
        integer::epoc,ierr,k

        if (errorflag .ne. 0) return

    
        epoc=450
        ierr=0
       
        open(unit=epoc,file='epoc.csv',status="old",access='append',iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(pass.eq.1)then
            do k=1,ndet-1
                if(chng_trk(k).eq.0)then
                    EXIT 
                end if 
                write(epoc,'(a,",",e25.17e3,",",e25.17e3,",",i0)') "   ",erg_dim(k),lr(k),chng_trk(k)
            end do
            write(epoc,'(i0,",",e25.17e3,",",a,",",*(i0:", "))') step,erg,"   ",(chng_trk(k),k=1,ndet-1)
        else
           
            do k=1,ndet-1
                if(chng_trk(k).eq.0)then
                    EXIT 
                end if 
                write(epoc,'(a,",",e25.17e3,",",e25.17e3,",",i0)') "   ",erg_dim(k),lr(k),chng_trk(k)
            end do
            write(epoc,'(i0,",",e25.17e3,",",a,",",*(i0:", "))') step,erg,"   ",(chng_trk(k),k=1,ndet-1)
            ! write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",*(i0:", "))') step,erg,lr,(chng_trk(k),k=1,ndet-1)
        end if
        close(epoc)
       
        return
        
    end subroutine epoc_writer_array


    subroutine dvec_writer(d,size,p)

        implicit none
        ! complex(kind=8),dimension(:),intent(in)::d
        real(kind=8),dimension(:),intent(in)::d
        character(LEN=4)::stateno
        integer,intent(in)::size,p
        integer::ierr,j,vec
        logical :: file_exists
        if (errorflag .ne. 0) return
        
        ierr=0

        write(stateno,"(i4.4)")p
        
        vec=900+p
        inquire(file="data/dvec_"//trim(stateno)//".csv",exist=file_exists)
        if(file_exists.eqv..false.) then
            open(unit=vec,file="data/dvec_"//trim(stateno)//".csv",status="new",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening dvector file. ierr had value ", ierr
                errorflag=1
                close(vec)
                return
            end if
        else 
            open(unit=vec,file="data/dvec_"//trim(stateno)//".csv",status="old",access='append',iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening dvector file. ierr had value ", ierr
                errorflag=1
                close(vec)
                return
            end if
            write(vec,*)' '
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
        ! complex(kind=8),dimension(:),intent(in)::d
        real(kind=8),dimension(:),intent(in)::d
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
        ! complex(kind=8),intent(in)::clean_erg,clean_norm
        real(kind=8),intent(in)::clean_erg,clean_norm
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