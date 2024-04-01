MODULE outputs
    use mod_types
    use globvars

    interface epoc_writer

        module procedure epoc_writer_int,epoc_writer_array,epoc_writer_array_orbital,&
                epoc_writer_int_gram,epoc_writer_array_gram,epoc_writer_array_orbital_gram

    end interface epoc_writer

    contains

    subroutine matrixwriter_real(out,size,filenm)

        implicit none

        real(wp),dimension(:,:),intent(in)::out
        integer,intent(in)::size
        character(LEN=*),intent(in)::filenm
        integer::j,k
        integer::ierr=0

        if (errorflag .ne. 0) return
        
        open(unit=200,file=filenm,status="new",iostat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in opening matrix file. ierr had value ", ierr
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

        real(wp),dimension(:,:),intent(in)::out
        integer,intent(in)::size
        character(LEN=*),intent(in)::filenm
        integer::j,k
        logical :: file_exists
        integer::ierr=0
        if (errorflag .ne. 0) return
        
        ierr=0
        inquire(file=filenm,exist=file_exists)
        if(file_exists.eqv..false.) then
            open(unit=200,file=filenm,status="new",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening matrix file. ierr had value ", ierr
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
            write(stdout,"(a,a)") trim(filenm), " alread exists so not rewriting"
        end if 

        return

    end subroutine matrixwriter

    subroutine zombiewriter(zom,num,gst)

        implicit none

        type(zombiest),intent(in)::zom 
        integer,intent(in)::num,gst
        character(LEN=20)::filenm
        integer::zomnum,j
        character(LEN=4)::nums
        character(len=2)::gst_num
        logical :: file_exists
        integer::ierr=0

        if (errorflag .ne. 0) return

        
        if(GDflg.eq.'y')then
            write(nums,"(i4.4)")(num+1000)
        else
            write(nums,"(i4.4)")num
        end if

        if(gst==0)then
            filenm="data/zombie_"//trim(nums)//".csv"
        else
            write(gst_num,"(i2.1)")gst
            filenm="data/zom_"//trim(gst_num)//"_"//trim(nums)//".csv"
        end if

        inquire(file=filenm,exist=file_exists)
        zomnum=300+num
        if(file_exists.eqv..false.) then
            open(unit=zomnum,file=trim(filenm),status="new",iostat=ierr)
            if(ierr/=0)then
                close(zomnum)
                write(stderr,"(a,i0)") "Error in opening zombie state file. ierr had value ", ierr
                errorflag=1
                return
            end if
        
            if(imagflg=='n') then
                write(zomnum,'(*(e25.17e3 :", "))') ((zom%phi(j)),j=1,norb)
                write(zomnum,'(*(e25.17e3 :", "))') ((zom%val(j)),j=1+norb,2*norb)
                write(zomnum,'(*(e25.17e3 :", "))') ((zom%val(j)),j=1,norb)
            else if(imagflg=='y') then
                write(zomnum,'(*(e25.17e3 :", ": ))') ((zom%phi(j)),j=1,norb)
                write(zomnum,'(*(e25.17e3 :", "))') ((zom%val(j)),j=1+norb,2*norb)
                write(zomnum,'(*(e25.17e3 :", "))') ((zom%val(j)),j=1,norb)
            end if
            close(zomnum)
        else if(file_exists.eqv..true.) then
            open(unit=zomnum,file=trim(filenm),status="old",access='append',iostat=ierr)
            if(ierr/=0)then
                close(zomnum)
                write(stderr,"(a,i0)") "Error in opening zombie state file. ierr had value ", ierr
                errorflag=1
                return
            end if
           
            write(zomnum,*)' '
            if(imagflg=='n') then
                write(zomnum,'(*(e25.17e3 :", "))') ((zom%phi(j)),j=1,norb)
                write(zomnum,'(*(e25.17e3 :", "))') ((zom%val(j)),j=1+norb,2*norb)
                write(zomnum,'(*(e25.17e3 :", "))') ((zom%val(j)),j=1,norb)
            else if(imagflg=='y') then
                ! write(zomnum,'(*(e25.17e3 :", ": ))') ((zom%phi(j)),j=1,norb)
                ! write(zomnum,'(*(e25.17e3 :", ": ))') ((zom%img(j)),j=1,norb)
                ! write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%cos(j)),j=1,norb)
                ! write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%sin(j)),j=1,norb)
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
        integer::zomnum,j
        character(LEN=6)::nums
        integer::ierr=0

        if (errorflag .ne. 0) return

        write(nums,"(i6.6)")num

        
        filenm = "data/clean_zombie_"//trim(nums)//".csv"

        zomnum=300+num
        
        open(unit=zomnum,file=trim(filenm),status="unknown",iostat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in opening zombie state file. ierr had value ", ierr
            errorflag=1
            return
        end if

        if(imagflg=='n') then
            write(zomnum,'(*(e25.17e3 :", "))') ((zom%val(j)),j=1+norb,2*norb)
            write(zomnum,'(*(e25.17e3 :", "))') ((zom%val(j)),j=1,norb)
        
        else if(imagflg=='y') then
            ! write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%dead(j)),j=1,norb)
            ! write(zomnum,'(*(1x,es25.17e3 :", "))') ((zom%alive(j)),j=1,norb)
        end if

        close(zomnum)

        return

    end subroutine zombiewriter_c

    subroutine energywriter(erg,filenm,j)

        implicit none

        real(wp), dimension(:,:),intent(in)::erg 
        character(LEN=*),intent(in)::filenm
        real::db
        integer,intent(in)::j
        integer::ergnum,k
        logical :: file_exists
        integer::ierr=0

        if (errorflag .ne. 0) return

        ergnum=400+j
     
        db=beta/timesteps
      
        inquire(file=trim(filenm),exist=file_exists)
        if(file_exists.eqv..false.) then
            open(unit=ergnum,file=trim(filenm),status="new",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening energy file. ierr had value ", ierr
                errorflag=1
                close(ergnum)
                return
            end if
            do k=1,timesteps+1
                write(ergnum,'(*(e25.17e3 :", "))') db*(k-1),erg(:,k)
            end do
            close(ergnum)
        else 
            write(stdout,"(a)") "File alread exists so not rewriting"
     
        end if
        
        return

    end subroutine energywriter

    subroutine epoc_writer_int(erg,step,chng_trk,learningr,pass)

        implicit none
        real(wp),intent(in)::erg,learningr 
        integer,intent(in)::step,pass
        integer,intent(in)::chng_trk
        logical :: file_exists
        integer::ierr=0,epoc=450
        if (errorflag .ne. 0) return

        
        inquire(file='epoc.csv',exist=file_exists)
       
        if(file_exists.eqv..false.) then
            open(unit=epoc,file='epoc.csv',status="new",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
                errorflag=1
                return
            end if
            write(epoc,'(a,",",a,",",a,","a)') "EPOC", "Energy", "Learning rate", "Zombie state altered"
            write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",i0)') 0,erg,0.0,0
            close(epoc)
        else if(file_exists.eqv..true.) then
            open(unit=epoc,file='epoc.csv',status="old",access='append',iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
                errorflag=1
                return
            end if
            if(pass.eq.1)then
                write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",i0)') step,erg,learningr,chng_trk
            else
                write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",i0)') step,erg,learningr,chng_trk
            end if
            close(epoc)
        end if
        return
        
    end subroutine epoc_writer_int

    subroutine epoc_writer_array_orbital(erg,step,learningr,chng_trk,pass)

        implicit none
        real(wp),intent(in)::erg
        integer,intent(in)::step,pass
        real(wp),intent(in)::learningr
        integer,dimension(:),intent(in)::chng_trk
        integer::k
        integer::ierr=0, epoc=450
        if (errorflag .ne. 0) return

       
        open(unit=epoc,file='epoc.csv',status="old",access='append',iostat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(pass.eq.1)then
            write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",*(i0:", "))') step,erg,learningr,(chng_trk(k),k=1,ndet-1)
         
        else
            write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",*(i0:", "))') step,erg,learningr,(chng_trk(k),k=1,ndet-1)
        
        end if
        close(epoc)
       
        return
        
    end subroutine epoc_writer_array_orbital

    subroutine epoc_writer_array(erg,step,chng_trk,erg_dim,learningr,pass)

        implicit none
        real(wp),intent(in)::erg
        real(wp),dimension(:),intent(in)::learningr,erg_dim 
        integer,intent(in)::step,pass
        integer,dimension(:),intent(in)::chng_trk
        integer::k
        integer::ierr=0, epoc=450

        if (errorflag .ne. 0) return

       
        open(unit=epoc,file='epoc.csv',status="old",access='append',iostat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(pass.eq.1)then
            do k=1,ndet-1
                if(chng_trk(k).eq.0)then
                    EXIT 
                end if 
                write(epoc,'(a,",",e25.17e3,",",e25.17e3,",",i0)') "   ",erg_dim(k),learningr(k),chng_trk(k)
            end do
            write(epoc,'(i0,",",e25.17e3,",",a,",",*(i0:", "))') step,erg,"   ",(chng_trk(k),k=1,ndet-1)
        else
           
            do k=1,ndet-1
                if(chng_trk(k).eq.0)then
                    EXIT 
                end if 
                write(epoc,'(a,",",e25.17e3,",",e25.17e3,",",i0)') "   ",erg_dim(k),learningr(k),chng_trk(k)
            end do
            write(epoc,'(i0,",",e25.17e3,",",a,",",*(i0:", "))') step,erg,"   ",(chng_trk(k),k=1,ndet-1)
        end if
        close(epoc)
       
        return
        
    end subroutine epoc_writer_array

    subroutine epoc_writer_int_gram(erg,step,chng_trk,learningr,pass,gst)

        implicit none
        real(wp),intent(in)::erg,learningr 
        integer,intent(in)::step,pass
        integer,intent(in)::chng_trk
        integer,intent(in)::gst
        character(len=2)::gst_num
        logical :: file_exists
        integer::ierr=0,epoc=450

        if (errorflag .ne. 0) return

        write(gst_num,"(i2.1)")gst
        inquire(file="epoc_"//trim(gst_num)//".csv",exist=file_exists)
       
        if(file_exists.eqv..false.) then
            open(unit=epoc,file="epoc_"//trim(gst_num)//".csv",status="new",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
                errorflag=1
                return
            end if
            write(epoc,'(a,",",a,",",a,","a)') "EPOC", "Energy", "Learning rate", "Zombie state altered"
            write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",i0)') 0,erg,0.0,0
            close(epoc)
        else if(file_exists.eqv..true.) then
            open(unit=epoc,file="epoc_"//trim(gst_num)//".csv",status="old",access='append',iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
                errorflag=1
                return
            end if
            if(pass.eq.1)then
                write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",i0)') step,erg,learningr,chng_trk
            else
                write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",i0)') step,erg,learningr,chng_trk
            end if
            close(epoc)
        end if
        return
        
    end subroutine epoc_writer_int_gram

    subroutine epoc_writer_array_orbital_gram(erg,step,learningr,chng_trk,pass,gst)

        implicit none
        real(wp),intent(in)::erg
        integer,intent(in)::step,pass
        real(wp),intent(in)::learningr
        integer,intent(in)::gst
        character(len=2)::gst_num
        integer,dimension(:),intent(in)::chng_trk
        integer::k
        integer::ierr=0, epoc=450
        if (errorflag .ne. 0) return

        write(gst_num,"(i2.1)")gst
       
        open(unit=epoc,file="epoc_"//trim(gst_num)//".csv",status="old",access='append',iostat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(pass.eq.1)then
            write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",*(i0:", "))') step,erg,learningr,(chng_trk(k),k=1,ndet-1)
         
        else
            write(epoc,'(i0,",",e25.17e3,",",e25.17e3,",",*(i0:", "))') step,erg,learningr,(chng_trk(k),k=1,ndet-1)
        
        end if
        close(epoc)
       
        return
        
    end subroutine epoc_writer_array_orbital_gram

    subroutine epoc_writer_array_gram(erg,step,chng_trk,erg_dim,learningr,pass,gst)

        implicit none
        real(wp),intent(in)::erg
        real(wp),dimension(:),intent(in)::learningr,erg_dim 
        integer,intent(in)::step,pass
        integer,dimension(:),intent(in)::chng_trk
        integer,intent(in)::gst
        character(len=2)::gst_num
        integer::k
        integer::ierr=0, epoc=450

        if (errorflag .ne. 0) return
        write(gst_num,"(i2.1)")gst
       
        open(unit=epoc,file="epoc_"//trim(gst_num)//".csv",status="old",access='append',iostat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in opening epoc output file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(pass.eq.1)then
            do k=1,ndet-1
                if(chng_trk(k).eq.0)then
                    EXIT 
                end if 
                write(epoc,'(a,",",e25.17e3,",",e25.17e3,",",i0)') "   ",erg_dim(k),learningr(k),chng_trk(k)
            end do
            write(epoc,'(i0,",",e25.17e3,",",a,",",*(i0:", "))') step,erg,"   ",(chng_trk(k),k=1,ndet-1)
        else
           
            do k=1,ndet-1
                if(chng_trk(k).eq.0)then
                    EXIT 
                end if 
                write(epoc,'(a,",",e25.17e3,",",e25.17e3,",",i0)') "   ",erg_dim(k),learningr(k),chng_trk(k)
            end do
            write(epoc,'(i0,",",e25.17e3,",",a,",",*(i0:", "))') step,erg,"   ",(chng_trk(k),k=1,ndet-1)
        end if
        close(epoc)
       
        return
        
    end subroutine epoc_writer_array_gram


    subroutine dvec_writer(d,size,p)

        implicit none
        real(wp),dimension(:),intent(in)::d
        character(LEN=4)::stateno
        integer,intent(in)::size,p
        integer::j,vec
        logical :: file_exists
        integer::ierr=0

        if (errorflag .ne. 0) return

        write(stateno,"(i4.4)")p
        
        vec=900+p
        inquire(file="data/dvec_"//trim(stateno)//".csv",exist=file_exists)
        if(file_exists.eqv..false.) then
            open(unit=vec,file="data/dvec_"//trim(stateno)//".csv",status="new",iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening dvector file. ierr had value ", ierr
                errorflag=1
                close(vec)
                return
            end if
        else 
            open(unit=vec,file="data/dvec_"//trim(stateno)//".csv",status="old",access='append',iostat=ierr)
            if(ierr/=0)then
                write(stderr,"(a,i0)") "Error in opening dvector file. ierr had value ", ierr
                errorflag=1
                close(vec)
                return
            end if
            write(vec,*)' '
        end if
       
        if(imagflg=='n')then
            write(vec,'(*(e25.17e3 :", "))') ((d(j)),j=1,size)
        else if(imagflg=='y')then
            write(vec,'(*(1x,es25.17e3 :", "))') ((d(j)),j=1,size*2)
        end if
        close(vec)
        return

    end subroutine dvec_writer

    subroutine dvec_writer_c(d,size,p)

        implicit none
        ! complex(wp),dimension(:),intent(in)::d
        real(wp),dimension(:),intent(in)::d
        character(LEN=4)::stateno
        integer,intent(in)::size,p
        integer::j,vec
        integer::ierr=0

        if (errorflag .ne. 0) return
    
        write(stateno,"(i4.4)")p

        vec=900+p
        open(unit=vec,file="data/clean_dvec_"//trim(stateno)//".csv",status="new",iostat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in opening dvector file. ierr had value ", ierr
            errorflag=1
            return
        end if
        if(imagflg=='n')then
            write(vec,'(*(e25.17e3 :", "))') ((d(j)),j=1,size)
        else if(imagflg=='y')then
            write(vec,'(*(1x,es25.17e3 :", "))') ((d(j)),j=1,size*2)
        end if
        close(vec)
        return

    end subroutine dvec_writer_c

    subroutine clean_erg_write(clean_ndet, clean_erg,clean_norm,j)

        implicit none
        integer, intent(in)::clean_ndet,j
        ! complex(wp),intent(in)::clean_erg,clean_norm
        real(wp),intent(in)::clean_erg,clean_norm
        integer::cleane
        integer::ierr=0

        cleane=400+j
    
        open(unit=cleane,file='clean_energy.csv',status="new",iostat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in opening clean energy output file. ierr had value ", ierr
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

    subroutine elec_inegrals_write(elecs)

        implicit none
        type(elecintrgl), intent(in)::elecs
        integer::j,k,choice2_dim,ierr=0

        if (errorflag .ne. 0) return

        choice2_dim=size(elecs%orbital_choice2,dim=2)
        open(unit=500,file='integrals/elec_integrals.csv',status="new",iostat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error in opening elec_integrals.csv file. ierr had value ", ierr
            errorflag=1
            return
        end if

        write(500,'(i0)') elecs%num
        do k=1,elecs%num
            write(500, '(e25.17e3,",",*(i0 : ", "))') elecs%integrals(k),(elecs%orbital_choice(j,k),j=1,norb)
        end do
        write(500, '(i0)') choice2_dim
        write(500, '(*(i0 : ", "))') (elecs%orbital_choice2(0,j),j=1,norb)
        do k=1,norb
            write(500, '(*(i0: ", "))') (elecs%orbital_choice2(k,j),j=1,2*elecs%orbital_choice2(0,k))
        end do
        write(500, '(*(i0 : ", "))') (elecs%orbital_choice3(j),j=1,norb)
        write(500, '(e25.17e3)') elecs%hnuc
        close(500)

    end subroutine elec_inegrals_write

END MODULE outputs