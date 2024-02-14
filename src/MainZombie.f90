program MainZombie
    
    use mod_types
    use randgen
    use globvars
    use readpars
    use alarrays
    use electrons
    use ham
    use outputs
    use imgtp
    use clean
    use zom
    use gradient_descent
    use operators
    ! use neural_network
    use omp_lib

    implicit none
    
  
    type(zombiest), dimension(:), allocatable:: zstore, cstore
    type(dvector):: dvecs
    type(elecintrgl)::elect
    type(hamiltonian)::haml
    integer:: j,clean_ndet
    real(wp), dimension(:),allocatable::erg
    type(hamiltonian)::clean_haml
    type(dvector)::dvec_clean
    real(wp)::clean_norm, clean_erg,num1,num2
    character(LEN=100) :: CWD
    real(wp):: starttime, stoptime, runtime
    integer::ierr=0, istat=0
    ! type(neural_network_layer),dimension(:),allocatable::neural_net

    call CPU_TIME(starttime) !used to calculate the runtime, which is output at the end
    write(stdout,"(a)") " ________________________________________________________________ "
    write(stdout,"(a)") "|                                                                |"
    write(stdout,"(a)") "|                                                                |"
    write(stdout,"(a)") "|               Zombie State Simulation Program v3.10            |"
    write(stdout,"(a)") "|                                                                |"
    write(stdout,"(a)") "|________________________________________________________________|"
    write(stdout,"(a)") ""
    write(stdout,"(a)") ""
    write(stdout,"(a)") ""

  
    istat=0
    call initialise
    call readrunconds
    call restart_chk
    if(rstrtflg.eq.'y')then
        ndet=zom_count()
    end if 
    if(randseed.eq.0)then
        open(unit=570, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
            read(570) randseed    ! This takes the random seed from the true-random bin. If
            close(570)           ! the urandom bin does not exist the random seed is set
        else                   ! to zero which forces the date to be used
            randseed=0
        end if
        randseed =  abs(randseed)    ! Negative seed values seem to cause instability
    end if
    call ZBQLINI(randseed,0)   ! Generates the seed value using the UCL random library
    write(stdout,"(a)") "Random seed set"
   
    ! generate 1 and 2 electron integrals
    if((cleanflg=="y").or.(cleanflg=="f").or.((hamgflg=='y')).or.(GDflg=='y'))then
       
        write(stdout,"(a)") "Setting electron"
        call electronintegrals(elect)
        write(stdout,"(a)") "Electrons allocated"
        
        if (gramflg.eq."n") then 
            ! generate zombie states
            call alloczs(zstore,ndet)
            write(stdout,"(a)") "Zombie states allocated"
            if(zomgflg=='y')then
                call genzf(zstore,ndet) 
                do j=1,ndet
                    call zombiewriter(zstore(j),j,0)
                end do
                write(stdout,"(a)") "Zombie states generated"
            else if (zomgflg=='n') then
                call read_zombie(zstore,0)
             
                write(stdout,"(a)") "Zombie read in"
            end if
            call flush(6)
            call flush(0)
        end if 
    end if

    if(propflg=="y")then
        if (gramflg.eq."n") then 
            ! generate Hamiltonian and overlap
            call allocham(haml,ndet)
            write(stdout,"(a)") "Hamiltonian allocated"
            call allocdv(dvecs,ndet)
            allocate(erg(timesteps+1),stat=ierr)
            if(ierr/=0)then 
                errorflag=1
                write(stderr,"(a,i0)") "Error in erg allocation. ierr had value ", ierr
            end if
            write(stdout,"(a)") "d-vector and energy array allocated"
        
            if(hamgflg=='y')then
                write(stdout,"(a)") "To hamiltonian gen"
                call hamgen(haml,zstore,elect,ndet,1)
                ! call neural_network_control(10000,elect,haml,zstore,neural_net)
                call matrixwriter(haml%hjk,ndet,"data/ham.csv")
                call matrixwriter(haml%ovrlp,ndet,"data/ovlp.csv")
                write(stdout,"(a)") "Hamiltonian successfully generated"
            else if (hamgflg=='n')then
                call read_ham(haml,ndet)
                write(stdout,"(a)") "Hamiltonian successfully read in"
            end if
            ! Imaginary time propagation
            write(stdout,"(a)") "Imaginary time propagation started"
            call imaginary_time_prop2(dvecs,erg,haml,ndet)
           
            write(stdout,"(a)") "Imaginary time propagation finished"
       
            write(stdout,"(a,f21.16)") "Initial energy: ", erg(timesteps+1)
            
            call dvec_writer(dvecs%d,ndet,0)
            call energywriter(erg,"energy.csv",0)
       
            if(GDflg.eq."y")then
                deallocate(erg,stat=ierr)
                call zombie_alter(zstore,haml,elect,dvecs)!,neural_net)
                GDflg='n'
                num1=0.0d0;num2=0.0d0
                do j=1,ndet
                    call zombiewriter(zstore(j),j,0)
                    num1=num1+zstore(j)%num_strt
                    zstore(j)%num=numf(zstore(j),zstore(j))
                    print*,zstore(j)%num_strt, zstore(j)%num
                    num2=num2+zstore(j)%num
                end do
                print*,num1/ndet,num2/ndet
                dvecs%d=0.0d0
                allocate(erg(timesteps+1),stat=ierr)
                call imaginary_time_prop2(dvecs,erg,haml,ndet)
                
                write(stdout,"(a,f21.16)") "Final energy: ", erg(timesteps+1)
                call energywriter(erg,"energy_final.csv",0)
                call matrixwriter(haml%hjk,ndet,"data/ham_final.csv")
                call matrixwriter(haml%ovrlp,ndet,"data/ovlp_final.csv")
                
                ! call sd_anal(zstore,nel,dvecs(1),2)
            end if
     
            deallocate(erg,stat=ierr)
            if(ierr/=0)then
                errorflag=1
                write(stderr,"(a,i0)") "Error in erg deallocation. ierr had value ", ierr
            end if
            write(stdout,"(a)") "Energy deallocated"
            call deallocham(haml)
            write(stdout,"(a)") "Hamiltonian deallocated"
            if(cleanflg=="n")then
                call deallocdv(dvecs)
                write(stdout,"(a)") "d-vector deallocated"
                !if(hamgflg=='y')then
                call dealloczs(zstore)
                write(stdout,"(a)") "Zombie states deallocated"
                call deallocintgrl(elect)
                write(stdout,"(a)") "Electron integrals deallocated"
                !end if
            end if

            call flush(6)
            call flush(0)
        else if(gramflg.eq."n") then 
            ! call gram_schmidt_control(elect,ndet)
        else
            write(stderr,"(a,i0)") "Error in gramflg setting. This should have been caught ", ierr
            errorflag=1
        end if 

    else if((propflg=="n"))then
        if((cleanflg=="y").or.(cleanflg=="f"))then
            if(gramflg.eq."n")then
                call allocdv(dvecs,ndet)
                call dvec_read(dvecs%d,ndet,0,'dvec_0000.csv')
                write(stdout,"(a)") "d-vector read in"
            else if(gramflg.eq."y")then
                ! call gram_schmidt_control(elect,ndet)
            else
                write(stderr,"(a,i0)") "Error in gramflg setting. This should have been caught ", ierr
                    errorflag=1
            end if
        else if((cleanflg=="n").and.(hamgflg=='y'))then
            write(stdout,"(a)") "The program if here has done nothing except read in some values and then deallocate them"
            call dealloczs(zstore)
            write(stdout,"(a)") "Zombie states deallocated"
            call deallocintgrl(elect)
            write(stdout,"(a)") "Electron integrals deallocated"
    
        end if
    end if
    
    if((cleanflg=="y").or.(cleanflg=="f"))then
        if(cleanflg=="y")then
            call clean_setup(cstore,nel,clean_haml,elect,clean_ndet,zstore)
            write(stdout,"(a)") "Cleaning hamiltonian generated"
        else if(cleanflg=="f")then
            call clean_read(cstore,clean_haml,clean_ndet,elect)
            write(stdout,"(a)") "Cleaning hamiltonian read in"
        end if
        
        call allocdv(dvec_clean,clean_ndet)
        call cleaner(zstore,cstore,dvecs,dvec_clean,clean_ndet,clean_norm)
        ! clean_erg=dot_product(dvec_clean(1)%d,matmul(clean_haml%hjk,dvec_clean(1)%d))
        clean_erg=ergcalc(clean_haml%hjk,dvec_clean%d)
        write(stdout,"(a)") "Cleaning process complete"
        call clean_erg_write(clean_ndet,clean_erg,clean_norm,99)
        call dvec_writer_c(dvec_clean%d,clean_ndet,0)
        call deallocdv(dvec_clean)
        write(stdout,"(a)") "Cleaning d-vector dealocated"
        call deallocham(clean_haml)
        write(stdout,"(a)") "Cleaning hamiltonian dealocated"
        call dealloczs(cstore)
        write(stdout,"(a)") "Cleaning Zombie states dealocated"
        call deallocdv(dvecs)
        write(stdout,"(a)") "d-vector deallocated"
        call dealloczs(zstore)
        write(stdout,"(a)") "Zombie states deallocated"
        !if((cleanflg=="y").or.((hamgflg=='y')))then
        call deallocintgrl(elect)
        write(stdout,"(a)") "Electron integrals deallocated"
        !end if
    end if


    write(stdout,"(a)") "All values deallocated"

    if (errorflag .ne. 0) then
        write(stdout,"(a)") "Program terminated early."
        write(stdout,"(a,i0)") "errorflag value is ", errorflag
    end if
    
    call CPU_TIME(stoptime)
    runtime = stoptime-starttime
    call getcwd(CWD)

    if (errorflag.eq.0) write(stdout,"(a,a)") 'Successfully Executed Zombie states Program in ', trim(CWD)
    if (errorflag.ne.0) write(stdout,"(a,a)") 'Unsuccessfully Executed Zombie states  Program in ', trim(CWD)
  
    if (runtime/3600.0d0 .gt. 1.0d0)then
        runtime = runtime/3600.0d0
        write(stdout,"(a,es12.5,a)") 'Time taken : ', runtime, ' hours'
    else if (runtime/60.0d0 .gt. 1.0d0)then
        runtime = runtime/60.0d0
        write(stdout,"(a,es12.5,a)") 'Time taken : ', runtime , ' mins'
    else
        write(stdout,"(a,es12.5,a)") 'Time taken : ', runtime, ' seconds'
    end if


    call flush(6)
    call flush(0)

    stop


end program MainZombie



