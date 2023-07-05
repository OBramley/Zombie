program bbi

    use globvars
    use readpars
    use alarrays
    use electrons
    use zom
    use ham
    use imgtp
    use outputs
    use comparison


    implicit none


    type(zombiest), dimension(:), allocatable:: zstore, zstore_act,zstore_temp
    type(dvector), dimension(:), allocatable:: dvecs
    type(energy):: en
    type(elecintrgl)::elect
    type(hamiltonian)::haml,haml_act
    real(kind=8):: starttime, stoptime, runtime, erg
    DOUBLE PRECISION, external::ZBQLUAB, ZBQLU01, ZBQLNOR
    character(LEN=100) :: CWD
    real(kind=8),dimension(:,:),allocatable::avalchange,ergchange
    real(kind=8),dimension(:),allocatable::eresult,aresult
    integer:: j,k,l,passes, istat,ierr,iters
    integer(kind=8):: randseed

    real(kind=8)::mu(19),sig(19)
    real(kind=8),dimension(38)::val
   

    call CPU_TIME(starttime) !used to calculate the runtime, which is output at the end
    write(6,"(a)") " ________________________________________________________________ "
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|             Zombie State bias improve Program v2.51            |"
    write(6,"(a)") "|                                                                |"
    write(6,"(a)") "|________________________________________________________________|"
    write(6,"(a)") ""
    write(6,"(a)") ""
    write(6,"(a)") ""

    
    ierr=0
    istat=0
    call initialise
    call readrunconds

    call allocintgrl(elect)
    call electronintegrals(elect)
    
    passes=1
    iters=20
    
    allocate(ergchange(((iters+1)*passes)+1,norb),stat=ierr)
    if(ierr==0)allocate(avalchange((iters*passes)+1,norb),stat=ierr)
    if (ierr/=0) then
        write(0,"(a,i0)") "Error in results allocation. ierr had value ", ierr
        errorflag=1
        return
    end if
    ergchange(:,:)=0.0
    avalchange(:,:)=0.0

    open(unit=570, file="/dev/urandom", access="stream", &
    form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
        read(570) randseed    ! This takes the random seed from the true-random bin. If
        close(570)            ! the urandom bin does not exist the random seed is set
    else                       ! to zero which forces the date to be used
        randseed=0
    end if

    randseed = abs(randseed)    ! Negative seed values seem to cause instability

    call ZBQLINI(randseed,0)   ! Generates the seed value using the UCL random library
    write(6,"(a)") "Random seed set"


    val(1:norb)=-1

    call alloczs(zstore,ndet)
    if(zomgflg=='y') then
        call musig(mu,sig)
        do j=1, ndet
            !$omp critical
            do k=1,norb/2
                do while(val(2*k-1).lt.0)
                    val(2*k-1)=2*pirl*ZBQLNOR(mu(k),sig(k))
                end do
                do while(val(2*k).lt.0)
                    val(2*k)=2*pirl*ZBQLNOR(mu(k),sig(k))
                end do
            end do
            !$omp end critical
            do k=1, norb
                zstore(j)%alive(k)=cmplx(sin(val(k)),0.0d0,kind=8)
                zstore(j)%dead(K)=cmplx(cos(val(k)),0.0d0,kind=8)
            end do
        end do
        if(rhf_1=='y') then
            zstore(1)%alive(1:nel)=(1.0d0,0.0d0)
            zstore(1)%dead(1:nel)=(0.0d0,0.0d0)
            zstore(1)%alive((nel+1):norb)=(0.0d0,0.0d0)
            zstore(1)%dead((nel+1):norb)=(1.0d0,0.0d0)
        end if 
        ! call genzf(zstore,ndet)
        do j=1,ndet
            call zombiewriter(zstore(j),j)
        end do
    else 
        call read_zombie(zstore)
    end if

    write(6,"(a)") "Initial Zombie states generated"
    do j=1,norb
        avalchange(1,j)=(dasin(REAL(zstore(2)%alive(j))))!/(2*pirl))
        ! print*,achange(1,j)
    end do

    
    call allocham(haml,ndet)
    call hamgen(haml,zstore,elect,ndet)
    call matrixwriter(haml%hjk,ndet,"data/ham.csv")
    call matrixwriter(haml%ovrlp,ndet,"data/ovrl.csv")
    ! call matrixwriter(haml%inv,ndet,"inv.csv")
    write(6,"(a)") "Initial Hamiltonian generated"
    call allocdv(dvecs,1,ndet)
    call allocerg(en,1)
    call imgtime_prop(dvecs,en,haml)

    erg=REAL(en%erg(1,timesteps+1))
    ergchange(1,1:norb)=erg
    write(6,"(a,e25.17e3)") "Starting energy: ", erg

    ! call random_ord(norb,mixed)
    ! print*,mixed(1:norb)

    call alloczs(zstore_act,ndet)
    call alloczs(zstore_temp,ndet)
    call allocham(haml_act,ndet)
    zstore_act=zstore
    zstore_temp=zstore
    haml_act=haml


    allocate(eresult(iters),stat=ierr)
    if(ierr==0) allocate(aresult(iters),stat=ierr)
    if (ierr/=0) then
        write(0,"(a,i0)") "Error in temporary results allocation. ierr had value ", ierr
        errorflag=1
        return
    end if
    

    aresult(1:iters)=0.0
    eresult(1:iters)=0.0
    ! zstore(2)%alive(1)=0.5
    ! zstore(2)%dead(1)=0.8660254038


    
   ! !$omp do ordered
    do l=1, passes
        write(6,"(a,i0)") "Starting Pass ", l
        !$omp parallel shared(zstore_temp,zstore,haml,elect,l,avalchange,ergchange,iters) &
        !$omp private(zstore_act,haml_act,eresult,aresult)
        !$omp do 
        do j=1, norb
            zstore_act=zstore
            haml_act=haml
            ! !$omp critical
            ! write(6,"(a,i0)") "Orbital ",j
            ! !$omp end critical
            call compare(zstore_act,erg,elect,haml_act,eresult,aresult,j,iters)
            !$omp critical
            zstore_temp(2)%alive(j)=zstore_act(2)%alive(j)
            zstore_temp(2)%dead(j)=zstore_act(2)%dead(j)
            ergchange(2+((iters+1)*(l-1)):(1+iters)+((iters+1)*(l-1)),j)=eresult(1:iters)
            avalchange((2+(iters*(l-1))):((1+iters)+(iters*l-1)),j)=aresult(1:iters)
            !$omp end critical
        end do
        !$omp end do
        !$omp end parallel 
      
        aresult(1:iters)=0.0
        eresult(1:iters)=0.0
        zstore=zstore_temp
        call hamgen(haml,zstore,elect,ndet)
        call imgtime_prop(dvecs,en,haml)
        zstore_act=zstore
        haml_act=haml
        erg=REAL(en%erg(1,timesteps+1))
        ergchange(((2+iters)+((l-1)*(iters+1))),1:norb)=erg
        write(6,"(a,e25.17e3)") "New energy: ", erg
    end do


    write(6,"(a)") "All passes complete"
    do j=1,ndet
        call zombiewriter_c(zstore(j),j)
    end do
    call matrixwriter(haml%hjk,ndet,"data/finalham.csv")
    call matrixwriter(haml%ovrlp,ndet,"data/finalovrl.csv")
    open(unit=200,file="ergchange.csv",status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening matrix file. ierr had value ", ierr
            errorflag=1
            return
        end if
        
        do j=1, ((iters+1)*passes)+1
            write(200,'(*(e25.17e3 :", "))') ((ergchange(j,k)),k=1,norb)
        end do
    close(200)

    open(unit=201,file="achange.csv",status="new",iostat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error in opening matrix file. ierr had value ", ierr
            errorflag=1
            return
        end if
        
        do j=1, (iters*passes)+1
            write(201,'(*(e25.17e3 :", "))') ((avalchange(j,k)),k=1,norb)
        end do
    close(201)

    call deallocerg(en)
    call deallocham(haml)
    call dealloczs(zstore)
    call deallocdv(dvecs)
    call deallocham(haml_act)
    call dealloczs(zstore_act)
    call dealloczs(zstore_temp)

    deallocate(ergchange,stat=ierr)
    if(ierr==0)deallocate(avalchange,stat=ierr)
    if(ierr==0)deallocate(aresult,stat=ierr)
    if(ierr==0)deallocate(eresult,stat=ierr)
    if (ierr/=0) then
        write(0,"(a,i0)") "Error in results deallocation. ierr had value ", ierr
        errorflag=1
        return
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

end program 
