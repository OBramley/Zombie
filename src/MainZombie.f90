program MainZombie
    
    
    use globvars
    use readpars
    use alarrays
    use electrons
    use ham
    use outputs
    use imgtp
    use clean
    use zom
    use grad_d
    use gradient_descent
    use omp_lib

    implicit none
    
  
    type(zombiest), dimension(:), allocatable:: zstore, cstore
    type(dvector), dimension(:), allocatable:: dvecs, dvec_clean
    type(energy):: en
    type(elecintrgl)::elect
    type(hamiltonian)::haml, clean_haml
    type(grad)::gradients
    type(oprts)::an_cr,an2_cr2
    integer:: j, istat,clean_ndet,ierr,diff_state
    real(kind=8)::clean_norm, clean_erg
    ! complex(kind=8)::clean_norm, clean_erg
    character(LEN=4)::stateno
    character(LEN=100) :: CWD
    real(kind=8):: starttime, stoptime, runtime
    integer(kind=8):: randseed
   
    ! DOUBLE PRECISION, external::ZBQLU01

    call CPU_TIME(starttime) !used to calculate the runtime, which is output at the end
    write(6,"(a)") " ________________________________________________________________ "
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|               Zombie State Simulation Program v3.10            |"
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|________________________________________________________________|"
    write(6,"(a)") ""
    write(6,"(a)") ""
    write(6,"(a)") ""

    ierr=0
    istat=0
    call initialise
    call readrunconds
    call restart_chk

    open(unit=570, file="/dev/urandom", access="stream", &
    form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
        read(570) randseed    ! This takes the random seed from the true-random bin. If
        close(570)           ! the urandom bin does not exist the random seed is set
    else                   ! to zero which forces the date to be used
        randseed=0
    end if
    randseed =  abs(randseed)    ! Negative seed values seem to cause instability
  
    call ZBQLINI(randseed,0)   ! Generates the seed value using the UCL random library
   
    write(6,"(a)") "Random seed set"

   
    diff_state=0
    if(GDflg=="y")then
        call allocgrad(gradients,ndet,norb)
        write(6,"(a)") "Gradients allocated"
        diff_state=2
    end if

    ! generate 1 and 2 electron integrals
    if((cleanflg=="y").or.(cleanflg=="f").or.((hamgflg=='y')).or.(GDflg=='y'))then
       

        write(6,"(a)") "Setting electron"
        call electronintegrals(elect,an_cr,an2_cr2)
        write(6,"(a)") "Electrons allocated"
      
        ! generate zombie states
        call alloczs(zstore,int(ndet,kind=16))
        write(6,"(a)") "Zombie states allocated"
        if(zomgflg=='y')then
            call genzf(zstore,ndet)
            do j=1,ndet
                call zombiewriter(zstore(j),j,0)
            end do
            write(6,"(a)") "Zombie states generated"
        else if (zomgflg=='n') then
            call read_zombie(zstore)
            write(6,"(a)") "Zombie read in"
        end if
        call flush(6)
        call flush(0)
    end if

   
    if(propflg=="y")then
        ! generate Hamiltonian and overlap
        call allocham(haml,ndet,norb)
        write(6,"(a)") "Hamiltonian allocated"
        if(gramflg.eq."n")then
            call allocdv(dvecs,1,ndet,norb)
            call allocerg(en,1)
            write(6,"(a)") "d-vector and energy array allocated"
        else if(gramflg.eq."y")then
            ! write(0,"(a,i0)") "Gram  ", ierr
            call allocdv(dvecs,1+gramnum,ndet,norb)
            call allocerg(en,1+gramnum)
        else
            write(0,"(a,i0)") "Error in gramflg setting. This should have been caught ", ierr
                errorflag=1
        end if
    
        if(hamgflg=='y')then
            write(6,"(a)") "To hamiltonian gen"
            call hamgen(haml,zstore,elect,ndet,an_cr,an2_cr2,1)
            call matrixwriter(haml%hjk,ndet,"data/ham.csv")
            call matrixwriter(haml%ovrlp,ndet,"data/ovlp.csv")
            write(6,"(a)") "Hamiltonian successfully generated"
        else if (hamgflg=='n')then
            call read_ham(haml,ndet)
            write(6,"(a)") "Hamiltonian successfully read in"
        end if
         
        ! Imaginary time propagation
        write(6,"(a)") "Imaginary time propagation started"
        ! call imgtime_prop(dvecs,en,haml,diff_state,0)
        call imgtime_prop(dvecs,en,haml,0,0)
        write(6,"(a)") "Imaginary time propagation finished"
       
        if(gramflg.eq."n")then
            call dvec_writer(dvecs(1)%d,ndet,0)
            
            call energywriter(en%t,en%erg(1,:),"energy.csv",0)
        else if(gramflg.eq."y")then
            do j=1, 1+gramnum
                write(stateno,"(i4.4)")j
                call dvec_writer(dvecs(j)%d,ndet,j)
                call energywriter(en%t,en%erg(j,:),"energy_state_"//trim(stateno)//".csv",j)
            end do
        end if
        
        
        if(GDflg.eq."y")then
            
            gradients%prev_erg=en%erg(1,timesteps+1)
            write(6,"(a,f20.16)") "Initial energy: ", gradients%prev_erg
           
            if(rstrtflg.eq.'n')then 
                ! call epoc_writer(gradients%prev_erg,0,chng_trk,0.0d0,0)
                call epoc_writer(gradients%prev_erg,0,0,0.0d0,0)
            end if
          
            call zombie_alter(zstore,gradients,haml,elect,en,dvecs,an_cr,an2_cr2)
            
            GDflg='n'
            do j=1,ndet
                call zombiewriter(zstore(j),j,0)
            end do
            dvecs(1)%d=(0.0,0.0)
            call imgtime_prop(dvecs,en,haml,0,0)
            
            write(6,"(a,f21.16)") "Final energy: ", en%erg(1,timesteps+1)
            ! write(6,"(a,f21.16)") "Final energy: ", real(en%erg(1,timesteps+1))
            call energywriter(en%t,en%erg(1,:),"energy_final.csv",0)
            call matrixwriter(haml%hjk,ndet,"data/ham_final.csv")
            call matrixwriter(haml%ovrlp,ndet,"data/ovlp_final.csv")
           
            ! call sd_anal(zstore,nel,dvecs(1),2)
        end if
     
        call deallocerg(en)
        write(6,"(a)") "Energy deallocated"
        call deallocham(haml)
        write(6,"(a)") "Hamiltonian deallocated"
        if(cleanflg=="n")then
            call deallocdv(dvecs)
            write(6,"(a)") "d-vector deallocated"
            !if(hamgflg=='y')then
            call dealloczs(zstore)
            write(6,"(a)") "Zombie states deallocated"
            call deallocintgrl(elect)
            write(6,"(a)") "Electron integrals deallocated"
            call dealloc_oprts(an_cr)
            call dealloc_oprts(an2_cr2)
            write(6,"(a)") "creation and annihilation operators deallocated"
            !end if
        end if

        call flush(6)
        call flush(0)

    else if((propflg=="n"))then
        if((cleanflg=="y").or.(cleanflg=="f"))then
            if(gramflg.eq."n")then
                call allocdv(dvecs,1,ndet,norb)
                call dvec_read(dvecs(1)%d,ndet,0,'dvec_0000.csv')
                write(6,"(a)") "d-vector read in"
            else if(gramflg.eq."y")then
                ! write(0,"(a,i0)") "Gram  ", ierr
                call allocdv(dvecs,1+gramnum,ndet,norb)
                do j=1, 1+gramnum
                    write(stateno,"(i4.4)")j
                    call dvec_read(dvecs(j)%d,ndet,j,"dvec_"//trim(stateno)//".csv")
                end do
                write(6,"(a)") "d-vectors read in"
            else
                write(0,"(a,i0)") "Error in gramflg setting. This should have been caught ", ierr
                    errorflag=1
            end if
        else if((cleanflg=="n").and.(hamgflg=='y'))then
            write(6,"(a)") "The program if here has done nothing except read in some values and then deallocate them"
            call dealloczs(zstore)
            write(6,"(a)") "Zombie states deallocated"
            call deallocintgrl(elect)
            write(6,"(a)") "Electron integrals deallocated"
            call dealloc_oprts(an_cr)
            call dealloc_oprts(an2_cr2)
            write(6,"(a)") "creation and annihilation operators deallocated"
        end if
    end if
    
    if((cleanflg=="y").or.(cleanflg=="f"))then
        if(cleanflg=="y")then
            call clean_setup(cstore,nel,clean_haml,elect,clean_ndet,zstore,an_cr,an2_cr2)
            write(6,"(a)") "Cleaning hamiltonian generated"
        else if(cleanflg=="f")then
            call clean_read(cstore,clean_haml,clean_ndet,elect,an_cr,an2_cr2)
            write(6,"(a)") "Cleaning hamiltonian read in"
        end if
        
        call allocdv(dvec_clean,1,clean_ndet,1)
        call cleaner(zstore,cstore,dvecs(1),dvec_clean(1),clean_ndet,clean_norm)
        ! clean_erg=dot_product(dvec_clean(1)%d,matmul(clean_haml%hjk,dvec_clean(1)%d))
        clean_erg=ergcalc(clean_haml%hjk,dvec_clean(1)%d)
        write(6,"(a)") "Cleaning process complete"
        call clean_erg_write(clean_ndet,clean_erg,clean_norm,99)
        call dvec_writer_c(dvec_clean(1)%d,clean_ndet,0)
        call deallocdv(dvec_clean)
        write(6,"(a)") "Cleaning d-vector dealocated"
        call deallocham(clean_haml)
        write(6,"(a)") "Cleaning hamiltonian dealocated"
        call dealloczs(cstore)
        write(6,"(a)") "Cleaning Zombie states dealocated"
        call deallocdv(dvecs)
        write(6,"(a)") "d-vector deallocated"
        call dealloczs(zstore)
        write(6,"(a)") "Zombie states deallocated"
        !if((cleanflg=="y").or.((hamgflg=='y')))then
        call deallocintgrl(elect)
        write(6,"(a)") "Electron integrals deallocated"
        !end if
        call dealloc_oprts(an_cr)
        call dealloc_oprts(an2_cr2)
        write(6,"(a)") "creation and annihilation operators deallocated"
    end if

    if(GDflg=="y")then
        call deallocgrad(gradients)
    end if

    write(6,"(a)") "All values deallocated"

    if (errorflag .ne. 0) then
        write(6,"(a)") "Program terminated early."
        write(6,"(a,i0)") "errorflag value is ", errorflag
    end if
    
    call CPU_TIME(stoptime)
    runtime = stoptime-starttime
    call getcwd(CWD)

    if (errorflag.eq.0) write(6,"(a,a)") 'Successfully Executed Zombie states Program in ', trim(CWD)
    if (errorflag.ne.0) write(6,"(a,a)") 'Unsuccessfully Executed Zombie states  Program in ', trim(CWD)
  
    if (runtime/3600.0d0 .gt. 1.0d0)then
        runtime = runtime/3600.0d0
        write(6,"(a,es12.5,a)") 'Time taken : ', runtime, ' hours'
    else if (runtime/60.0d0 .gt. 1.0d0)then
        runtime = runtime/60.0d0
        write(6,"(a,es12.5,a)") 'Time taken : ', runtime , ' mins'
    else
        write(6,"(a,es12.5,a)") 'Time taken : ', runtime, ' seconds'
    end if


    call flush(6)
    call flush(0)

    stop


end program MainZombie



