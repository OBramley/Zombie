MODULE gram_schmidt

    use mod_types
    use globvars
    use alarrays
    use readpars
    use ham 
    use outputs
    use zom 
    use imgtp
    use gradient_descent
    contains


    subroutine gram_schmidt_control(elect)

        implicit none

        type(elecintrgl),intent(in)::elect
        type(gram),dimension(:),allocatable::gramstore
      
        real(kind=8), dimension(:,:),allocatable::erg
        integer::ierr,j,k
        character(LEN=20)::filenm
        logical :: file_exists
        character(len=2)::gst_num


        if(errorflag.ne.0) return

        ierr=0
        if(gramwave.eq.1)then
            write(stdout,"(a)") "The Gram Schmidt in this routine requires multiple wavefunctions to be generated."
            errorflag=1
            return
        end if

        allocate(gramstore(gramnum+1),stat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error allocating zstore array ierr had value ",ierr
            errorflag=1
            return
        end if

        do j=1,gramnum+1
            call allocgram(gramstore(j),j,ndet,norb)
        end do
        write(stdout,"(a)") "Gram Schmidt orthogonalisation components allocated"
        if(zomgflg=='y')then
            call genzf(gramstore(1)%zstore,ndet) 
            do k=1,ndet
                call zombiewriter(gramstore(1)%zstore(k),k,gramstore(1)%zstore(k)%gram_num)
            end do
            rhf_1='n'
            do j=2,gramnum+1
                gramstore(j)%zstore=gramstore(1)%zstore
                ! call gen_biased_zs(gramstore(j)%zstore) 
                do k=1,ndet
                    call zombiewriter(gramstore(j)%zstore(k),k,gramstore(j)%zstore(k)%gram_num)
                end do
            end do
            write(stdout,"(a)") "Zombie states generated"
        else if (zomgflg=='n') then
            do j=1,gramnum+1
                call read_zombie(gramstore(j)%zstore,gramstore(j)%zstore(k)%gram_num)
            end do
            write(stdout,"(a)") "Zombie read in"
        end if
    
        if(hamgflg=='y')then
            write(stdout,"(a)") "To hamiltonian gen"
            do j=1,gramnum+1
                call hamgen(gramstore(j)%haml,gramstore(j)%zstore,elect,ndet,1)
                write(gst_num,"(i2.1)")j
                call matrixwriter(gramstore(j)%haml%hjk,ndet,"data/ham_"//trim(gst_num)//".csv")
                call matrixwriter(gramstore(j)%haml%ovrlp,ndet,"data/ovlp_"//trim(gst_num)//".csv")
            end do
            write(stdout,"(a)") "Hamiltonian successfully generated"
        else if (hamgflg=='n') then
            do j=1,gramnum+1
                write(gst_num,"(i2.1)")j
                call read_ham(gramstore(j)%haml,ndet,"ham_"//trim(gst_num)//".csv","ovlp_"//trim(gst_num)//".csv")
            end do
            write(stdout,"(a)") "Hamiltonian read in"
        end if
        write(stdout,"(a)") "Setting up  wavefunction overlap matrices "
        do j=1,gramnum+1
            call gram_ovrlp_fill(gramstore,j)
        end do 

        allocate(erg(gramnum+1,timesteps+1),stat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error allocating erg array ierr had value ",ierr
            errorflag=1
            return
        end if
        write(stdout,"(a)") "Imaginary time propagation for wave function 1"
        call imaginary_time(gramstore(1)%dvecs,erg(1,:),gramstore(1)%haml, ndet)
        gramstore(1)%d_ovrlp_d=dot_product(gramstore(1)%dvecs%d,matmul(gramstore(1)%haml%ovrlp,gramstore(1)%dvecs%d))
        do j=2,gramnum+1
            write(stdout,"(a,i0)") "Imaginary time propagation for wave function ",j
            call imaginary_time(gramstore,erg(j,:),ndet,j)
            gramstore(j)%d_ovrlp_d=dot_product(gramstore(j)%dvecs%d,matmul(gramstore(j)%haml%ovrlp,gramstore(j)%dvecs%d))
        end do
        
        do j=1,gramnum+1
            write(stdout,"(a,i1,a,f21.16)") "State, ",j," energy: ", erg(j,timesteps+1)
            call dvec_writer(gramstore(j)%dvecs%d,ndet,j)
        end do
        call energywriter(erg,"energy.csv",0)
        write(stdout,"(a)") "Imaginary time propagation complete"
        
        if(GDflg.eq."y")then
            do j=1,gramnum+1
                gramstore(j)%grads%prev_erg=erg(j,timesteps+1)
            end do 
            call zombie_alter_gram(gramstore,elect) 
            GDflg='n'
            do k=1,ndet
                call zombiewriter(gramstore(1)%zstore(k),k,gramstore(1)%zstore(k)%gram_num)
            end do
            call imaginary_time(gramstore(1)%dvecs,erg(1,:),gramstore(1)%haml, ndet)
            gramstore(1)%d_ovrlp_d=dot_product(gramstore(1)%dvecs%d,matmul(gramstore(1)%haml%ovrlp,gramstore(1)%dvecs%d))
            write(stdout,"(a,f21.16)") "State, 1 energy: ", erg(1,timesteps+1)
            call matrixwriter(gramstore(j)%haml%hjk,ndet,"data/ham_final_1.csv")
            call matrixwriter(gramstore(j)%haml%ovrlp,ndet,"data/ovlp_final_1.csv")
            do j=2,gramnum+1
                do k=1,ndet
                    call zombiewriter(gramstore(j)%zstore(k),k,gramstore(j)%zstore(k)%gram_num)
                end do
                call imaginary_time(gramstore(j)%dvecs,erg(j,:),gramstore(j)%haml,ndet)
                gramstore(j)%d_ovrlp_d=dot_product(gramstore(j)%dvecs%d,matmul(gramstore(j)%haml%ovrlp,gramstore(j)%dvecs%d))
                write(stdout,"(a,i1,a,f21.16)") "State, ",j," energy: ", erg(j,timesteps+1)
                write(gst_num,"(i2.1)")j
                call matrixwriter(gramstore(j)%haml%hjk,ndet,"data/ham_final_"//trim(gst_num)//".csv")
                call matrixwriter(gramstore(j)%haml%ovrlp,ndet,"data/ovlp_final_"//trim(gst_num)//".csv")
            end do
            call energywriter(erg,"energy_final.csv",0)
        end if

        deallocate(erg,stat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error deallocating erg array ierr had value ",ierr
            errorflag=1
            return
        end if

    end subroutine gram_schmidt_control

    subroutine zombie_alter_gram(gramstore,elect)

        implicit none

        type(gram),dimension(:)::gramstore
        type(elecintrgl),intent(in)::elect
        integer,dimension(:),allocatable::epoc_cnts
        character(len=2)::gst_num
        integer::cnt,j,k,epocfile
        integer::ierr=0
       
        if (errorflag .ne. 0) return
        if(ierr==0) allocate(epoc_cnts(gramnum+1),stat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error allocating picker array ierr had value ", ierr
            errorflag=1
            return
        end if
        epoc_cnt=1 !epoc counter
        if(rstrtflg.eq.'y')then 
            do j=1,gramnum+1
                epocfile=450+j
                write(gst_num,"(i2.1)")j
                open(unit=epocfile,file="epoc_"//trim(gst_num)//".csv",status="old",iostat=ierr)
                if(ierr/=0)then
                    write(stderr,"(a,i0)") "Error in opening epoc file to read in. ierr had value ", ierr
                    errorflag=1
                    return
                end if
                cnt=0
                do 
                    read(epocfile,*,iostat=ierr)
                    if(ierr<0)then
                        exit
                    else if(ierr>0)then
                        write(stderr,"(a,i0)") "Error in counting epocs. ierr had value ", ierr
                        errorflag=1
                        return
                    else 
                        cnt=cnt+1
                    end if
                    
                end do
                close(epocfile) 
                open(unit=epocfile,file="epoc_"//trim(gst_num)//".csv",status="old",iostat=ierr)
                if(ierr/=0)then
                    write(stderr,"(a,i0)") "Error in opening epoc file to read in. ierr had value ", ierr
                    errorflag=1
                    return
                end if
                do k=1,cnt-1
                    read(epocfile,*,iostat=ierr)
                end do 
                read(epocfile,*,iostat=ierr) epoc_cnts(j)
        
                close(epocfile) 

                write(stderr,"(a,i0)") "Epoc read in as ",  epoc_cnts(j)
        
                call epoc_writer(gramstore(j)%grads%prev_erg,epoc_cnt,0,0.0d0,1,j)
                epoc_cnts(j)= epoc_cnts(j)+1
            end do 
        else
            do j=1,gramnum+1
                call epoc_writer(gramstore(j)%grads%prev_erg,0,0,0.0d0,0,j)
            end do 
        end if

     
        call orbital_gd_gram_control(gramstore,elect,epoc_cnts)
        


        return

    end subroutine zombie_alter_gram

    

    subroutine gram_ovrlp_var(gramstore,state1,state2,var)
    
        implicit none
        type(gram),dimension(:)::gramstore
        integer::state1,state2,var
        integer::k
    
        if(errorflag.ne.0) return
    
        do k=1,ndet
            gramstore(state1)%wf_ovrlp(state2,k,var)=product(gramstore(state1)%zstore(k)%val(1:norb)*&
                                                    gramstore(state2)%zstore(var)%val(1:norb)+&
                                                    gramstore(state1)%zstore(k)%val(1+norb:2*norb)*&
                                                    gramstore(state2)%zstore(var)%val(1+norb:2*norb))
            gramstore(state1)%wf_ovrlp(state2,var,k)=gramstore(state1)%wf_ovrlp(state2,k,var)
        end do
            
        return
    
    end subroutine gram_ovrlp_var

    

End Module gram_schmidt