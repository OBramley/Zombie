program testbed
    use globvars
    use readpars
    use alarrays
    use electrons
    use ham
    use outputs
    use imgtp
    use operators
    use zom 

    implicit none


    type(zombiest), dimension(:), allocatable:: zstore, cstore
    type(zs2),dimension(:),allocatable::zstore2
    type(dvector), dimension(:), allocatable:: dvecs, dvec_clean
    type(energy):: en, en_clean
    type(elecintrgl)::elect
    type(hamiltonian)::haml, clean_haml
    integer:: j, k,l,m, istat, clean_ndet,ierr
    complex(kind=8)::clean_norm, clean_erg
    character(LEN=4)::stateno
    character(LEN=100) :: CWD
    real(kind=8):: starttime, stoptime, runtime
    integer(kind=8):: randseed
    DOUBLE PRECISION, external::ZBQLU01
    integer,dimension(:,:,:),allocatable::occupy
    real(kind=8),dimension(:,:),allocatable::shift
    type(zombiest),allocatable, dimension(:)::z2l
    type(zs2)::zcran
    type(zombiest)::zomt
    real(kind=8),dimension(:),allocatable::z1_alive,z1_dead,z2_alive, z2_dead
    ! 
    ! call CPU_TIME(stoptime)
    ! runtime = stoptime-starttime
    ! call getcwd(CWD)

    ! if (errorflag.eq.0) write(6,"(a,a)") 'Successfully Executed Zombie states Program in ', trim(CWD)
    ! if (errorflag.ne.0) write(6,"(a,a)") 'Unsuccessfully Executed Zombie states  Program in ', trim(CWD)
  
    ! if (runtime/3600.0d0 .gt. 1.0d0)then
    !     runtime = runtime/3600.0d0
    !     write(6,"(a,es12.5,a)") 'Time taken : ', runtime, ' hours'
    ! else if (runtime/60.0d0 .gt. 1.0d0)then
    !     runtime = runtime/60.0d0
    !     write(6,"(a,es12.5,a)") 'Time taken : ', runtime , ' mins'
    ! else
    !     write(6,"(a,es12.5,a)") 'Time taken : ', runtime, ' seconds'
    ! end if

    ierr=0
    istat=0
    call initialise

    open(unit=570, file="/dev/urandom", access="stream", &
    form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
        read(570) randseed    ! This takes the random seed from the true-random bin. If
        close(570)           ! the urandom bin does not exist the random seed is set
    else                   ! to zero which forces the date to be used
        randseed=0
    end if

    randseed = abs(randseed)    ! Negative seed values seem to cause instability

    call ZBQLINI(randseed,0)   ! Generates the seed value using the UCL random library
    write(6,"(a)") "Random seed set"

    norb=38
    ndet=100
    nel=6
    spin=0
    zst='RN'
    imagflg='n'

    
    print*,norb
    call alloczs(zstore,ndet)
    call alloczs2(zstore2,ndet)
    call gen_ran2(zstore2,ndet)
    ! call genzf(zstore,ndet)
    do j=1, ndet
        zstore(j)%alive=sin(zstore2(j)%phi)
        zstore(j)%dead=cos(zstore2(j)%phi)
    end do
    runtime=0.0
    call CPU_TIME(starttime)
    call alloczf(zomt)
    call alloczs(z2l,norb)
    
    do l=1,ndet
        do m=l,ndet
            do j=1, norb
                zomt=zstore(m)
                call an(zomt,j)
                z2l(j)=zomt
            end do

            do j=1,norb
                do  k=1, norb
                    zomt=z2l(j)
                    call cr(zomt,k)
                ! call cr(zomt,10)
                ! print*, REAL(zstore(j)%alive)
                ! print*, REAL(zstore(j)%dead)  
                    runtime = runtime+REAL(overlap(zstore(l),zomt))
                end do
            end do
        end do
    end do
    print*,runtime
    ! do j=1, ndet
    !     do k=1 ,ndet
    !         print*, overlap(zstore(j),zstore(k))
    !     end do
    ! end do

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

    ! call dealloczs(zstore)

    runtime=0.0
    call CPU_TIME(starttime)

    print*,norb

    call alloczf2(zcran)
    
    allocate(occupy((norb*norb),2,norb))
    allocate(shift(norb*norb,norb))
    occupy=1
    shift=0
    l=0
    
    do j=1,norb
        do k=1,norb
            l=l+1
            occupy(l,2,j)=1
            occupy(l,1,j)=0
            shift(l,j)=shift(l,j)-0.5*pirl
            occupy(l,1,1:j-1)=occupy(l,1,1:j-1)*(-1)
            occupy(l,1,k)=1
            occupy(l,2,k)=0
            shift(l,k)=shift(l,k)+0.5*pirl
            occupy(l,1,1:k-1)=occupy(l,1,1:k-1)*(-1)
        end do
    end do

    ! print*, occupy(1,1,:)
    ! print*, occupy(1,2,:)
    ! print*, shift(1,:)
    ! occupy(2,4)=0
    ! occupy(2,10)=0
    ! occupy(1,1:3)=(-1)*occupy(1,1:3)
    ! occupy(1,1:9)=(-1)*occupy(1,1:9)
    ! shift(4)=0.5*pirl
    ! shift(10)=0.5*pirl
    ! do j=1, ndet
    !     do k=1, norb
    !         zstore2(j)%phi(k)=2*pirl*ZBQLU01(1)
    !     end do
    ! end do
    allocate(z1_alive(norb))
    allocate(z1_dead(norb))
    do l=1,ndet
        z1_alive=sin(zstore2(l)%phi(:))
        z1_dead=cos(zstore2(l)%phi(:))
        do m=l, ndet
            do j=1, (norb*norb)
                zcran%phi=zstore2(m)%phi+shift(j,:)
                ! print*, zcran%phi
                zcran%alive=(zstore2(m)%alive(1:norb))*occupy(j,1,1:norb)
                zcran%dead=(zstore2(m)%dead(1:norb))*occupy(j,2,1:norb)
                ! print*, zcran%alive
                ! print*, zcran%dead
                ! call cr_2(zstore2(j),2)
                ! print*,  zstore2(j)%alive*sin(zstore2(j)%phi)
                ! print*,  zstore2(j)%dead*cos(zstore2(j)%phi)
                ! cr_2(zstore2(1),zstore2(j),2)
                ! runtime = runtime+overlap_two(zstore2(l),zcran)
                runtime = runtime+overlap_with_cran(z1_alive,z1_dead,zcran)
            end do
        end do
    end do
    print*, runtime

    ! do j=1, ndet
    !     do k=1 ,ndet
    !         print*, overlap2(zstore(j),zstore(k))
    !     end do
    ! end do

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





    ! call allocintgrl(elect)
    ! call electronintegrals(elect)

    



    ! call dealloczs(zstore)

end program