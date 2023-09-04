MODULE gram_schmidt

    use globvars
    use alarrays
    use readpars
    use ham 
    use outputs
    use zom 
    use imgtp
    use gradient_descent
    contains


    subroutine gram_schmidt_control(elect,ndet,an_cr,an2_cr2,nbf)

        implicit none

        type(elecintrgl),intent(in)::elect
        type(oprts),intent(in)::an_cr,an2_cr2
        type(zombiest), dimension(:,:), allocatable:: zstore
        type(dvector),dimension(:),allocatable:: dvecs
        type(hamiltonian),dimension(:),allocatable::haml
        real(kind=8), dimension(:,:),allocatable::erg
        integer::ierr,j,k
        character(LEN=20)::filenm
        logical :: file_exists
        character(len=2)::gst_num


        if(errorflag.ne.0) return

        ierr=0


        if(GDflg.eq."y")then
            ! Allocate 2D zombie array
            allocate(zstore(gramnum+1,ndet),stat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error allocating zstore array ierr had value ",ierr
                errorflag=1
                return
            end if
            ! Allocate 2D array of Hamiltonian objects
            allocate(haml(gramnum+1),stat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error allocating hamiltonian array ierr had value ",ierr
                errorflag=1
                return
            end if
            ! Allocate individual Zombie states
            do k=1,gramnum+1
                do j=1,ndet
                    call alloczf(zstore(k,j))
                end do
            end do
            ! Generate Zombie stattes
            if(zomgflg=='y')then
                call genzf(zstore(1,:),ndet)
                do k=1,gramnum+1
                    if(k>1)then
                        zstore(k,:)=zstore(1,:)
                    end if
                    do j=1,ndet
                        call zombiewriter(zstore(k,j),j,k)
                    end do
                end do
                write(6,"(a)") "Zombie states generated"
                ! Allocate first hamiltonian 
                call allocham(haml(1),ndet,norb)
                write(6,"(a)") "To hamiltonian gen"
                ! Generate hamiltonian values
                call hamgen(haml(1),zstore(1,:),elect,ndet,an_cr,an2_cr2,1)
                 ! Write hamiltonian to file 
                call matrixwriter(haml%hjk,ndet,"data/ham.csv")
                call matrixwriter(haml%ovrlp,ndet,"data/ovlp.csv")
                write(6,"(a)") "Hamiltonian successfully generated"
                ! Allocate excited state hamiltonians 
                do k=2,gramnum+1
                    call allocham(haml(k),ndet,norb)
                    haml(k)=haml(1)
                end do
            else if (zomgflg=='n') then
                ! Read in zombie states from file
                do k=1,gramnum+1
                    if(k>1)then
                        write(gst_num,"(i2.1)")k
                        filenm="data/zom_"//trim(gst_num)//"_1001.csv"
                        inquire(file=filenm,exist=file_exists)
                        if(file_exists.eqv..false.) then
                            zstore(k,:)=zstore(1,:)
                        else 
                            call read_zombie(zstore(k,:),k)
                        end if
                    end if
                    call read_zombie(zstore(k,:),k)
                end do
                write(6,"(a)") "Zombie read in"
                ! Allocate excited state hamiltonians
                do k=1,gramnum+1
                    call allocham(haml(k),ndet,norb)
                    write(6,"(a)") "To hamiltonian gen"
                    call hamgen(haml(k),zstore(k,:),elect,ndet,an_cr,an2_cr2,1)
                    write(6,"(a)") "Hamiltonian successfully generated"
                end do
            end if
            ! Allocate dvector array
            allocate(dvecs(gramnum+1),stat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error allocating dvecs array ierr had value ",ierr
                errorflag=1
                return
            end if
            !allocate all dvectors
            do k=1 gramnum+1
                call allocdv(dvecs(k),ndet,norb)
            end do
            write(6,"(a)") "d-vector allocated"
            
            ! Allocate hamiltonian and dvector gram-schmidt specific parts
            do k=2,(gramnum+1)
                call allocgram(haml(k),(k-1),ndet)
                call allocdvgram(dvec(k),(k-1),ndet)
            end do
           ! Fill the gram-schmidt speific parts of the hamiltonian object
          
            call gram_ovrlp_fill(haml,zstore)
           
            write(6,"(a)") "Gram overlap allocated"


            

            allocate(erg(timesteps+1),stat=ierr)
            if(ierr/=0)then 
                errorflag=1
                write(0,"(a,i0)") "Error in erg allocation. ierr had value ", ierr
            end if

            call imaginary_time_prop_gs(dvecs,erg,haml,ndet,gram_ovrlp,0,0)
            
            
        else
            allocate(zstore(1,ndet),stat=ierr)
            do j=1,ndet
                call alloczf(zstore(1,j))
            end do           
            if(zomgflg=='y')then
                call genzf(zstore(1,:),ndet)
                do j=1,ndet
                    call zombiewriter(zstore(1,j),j,0)
                end do
                write(6,"(a)") "Zombie states generated"
            else if (zomgflg=='n') then


                call read_zombie(zstore(1,:))
                write(6,"(a)") "Zombie read in"
            end if
            call flush(6)
            call flush(0)
        end if 
        

        if(propflg=="y")then


            if(GDflg.eq."y")then


            end if 



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

            else 
                write(6,"(a)") "The program if here has done nothing except read in some values and then deallocate them"
                call deallocintgrl(elect)
                write(6,"(a)") "Electron integrals deallocated"
                call dealloc_oprts(an_cr)
                call dealloc_oprts(an2_cr2)
                write(6,"(a)") "creation and annihilation operators deallocated"
            end if 

        end if 
       


        write(stateno,"(i4.4)")j
                call dvec_writer(dvecs%d,ndet,j)
                call energywriter(erg,"energy_state_"//trim(stateno)//".csv",j)

    end subroutine gram_schmidt_control

    subroutine gram_ovrlp_fill(zstore)
        
        implicit none

        type(zombiest),intent(in),dimension(:,:)::zstore
        type(hamiltonian),intent(inout),dimension(:)::haml
        integer::j,k,l,plc,state

        plc=1

        do state=2,gramnum+1
            do j=1, (state-1)
                haml(state)%gs_ovrlp(j,:,:)=haml(j)%ovrlp
                haml(state)%gs_kinvh(j,:,:)=haml(j)%kinvh
                do k=1,ndet
                    do l=k,ndet
                        haml(state)%gs_ovrlp_self(j,k,l)=overlap_1(zstore(state,k)%val,zstore(j,l)%val)
                        haml(state)%gs_ovrlp_self(j,l,k)=haml(state)%gs_ovrlp_self(j,k,l)
                    end do 
                end do
            end do
        end do 

        return

    end subroutine gram_ovrlp_fill


    !Gram-Schmidt orthogonalisation 
    subroutine gs(dvecs,haml)

        implicit none
        type(dvector), intent(inout)::dvecs
        type(hamiltonian),intent(in)::haml
        real(kind=8)::numer,den
        real(kind=8),dimension(:),allocatable::dvecs_copy
        integer::states,j,k,ierr

        if (errorflag .ne. 0) return
        ierr=0
        allocate(dvecs_copy(ndet),stat=ierr)
        if(ierr/=0)then
            write(0,"(a,i0)") "Error allocating dvecs_copy array ierr had value ",ierr
            errorflag=1
            return
        end if
       
       
        do j=2,gramnum+1
            dvecs_copy=dvecs(j)%d
            do k=1, j-1
                numer=dot_product(dvecs(j)%d,matmul(gram_ovrlp(k)%ovrlps_numer((j-k),:,:),dvecs(k)%d))
                den=dot_product(dvecs(k)%d,matmul(haml(k)%ovrlp,dvecs(k)%d))
                dvecs_copy = dvecs_copy - (dvecs(k)%d*(numer/den))
            end do
            dvecs(j)%d=dvecs_copy
        end do

        return

    end subroutine gs







End Module gram_schmidt